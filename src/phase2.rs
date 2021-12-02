use crate::types::{Distance,CenterIdx,PointCount};
use crate::space::{Point,ColoredMetric};
use crate::clustering::{Clustering,Centers};

mod buckets;
use buckets::put_into_buckets;

mod flow;
use flow::{initialize_state,add_edge};

mod settle;
use settle::settle;

#[cfg(debug_assertions)]
pub(crate) mod with_sorting;

use std::collections::VecDeque;
use std::collections::HashMap;

type EdgeIdx = usize;
type BucketIdx = usize;

// An edge in the flow network. One for every center (left; in form of an gonzalez index 0,..,k-1) to every point (right). note that center appear on both sides.
// The distance is stored in d.
#[derive(Debug,Clone,Copy,PartialOrd,PartialEq)]
pub (crate) struct Edge<'a> { // Care: The ordering of attributes is important for the partial order! Ties in d are broken by left, then right
    pub d : Distance,
    left : CenterIdx, // index of the center in gonzalez
    right : &'a Point,
}

struct PendingQueues<'a> {
    queues : HashMap<(CenterIdx,BucketIdx), VecDeque<Edge<'a>>>
}

impl<'a> PendingQueues<'a> {
    fn new() -> PendingQueues<'a> {
        PendingQueues { queues: HashMap::new()}
    }

    fn push(&mut self, c: CenterIdx, b: BucketIdx, edge: Edge<'a>) {
        match self.queues.get_mut(&(c,b)) {
            Some(queue) => {queue.push_back(edge)},
            None => {
                self.queues.insert((c,b),VecDeque::from(vec![edge]));
            }
        }
    }

    fn pop(&mut self, c:CenterIdx, b: BucketIdx) -> Option<Edge<'a>> {
        match self.queues.get_mut(&(c,b)) {
            Some(queue) => queue.pop_front(),
            None => None
        }
    }

    fn is_empty(&self, c: CenterIdx, b: BucketIdx) -> bool {
        match self.queues.get(&(c,b)) {
            Some(queue) => queue.is_empty(),
            None => true
        }
    }

    fn clear(&mut self, c: CenterIdx, b: BucketIdx) {
        match self.queues.get_mut(&(c,b)) {
            Some(queue) => queue.clear(),
            None => {}
        }

    }
}

// centers_graph: nodes = centers 0,...,i (i goes from 0 to k-1) and we have an arc c_1 to c_2 iff
// there is a point x that is assigned to c_1 but could also be assigne to c_2
// centers_graph is basically the connectivity within the residual graph on the left side of the
// flow network.


/// Given a metric space and a privacy bound, make_private takes a set of ordered centers
/// 0,...,k-1, and, in a first step, determines for each prefix of centers (0,...,i) a partial clustering
/// with minimal radius that satisfy the privacy constraint, i.e., each center (0,...,i) covers exactly privacy_bound many points.
/// Afterwards it assigns the remaining points to the nearest center to provide a full clustering.
pub(crate) fn make_private<M : ColoredMetric>(space : &M, privacy_bound : PointCount, centers : &Centers) -> Vec<Clustering> { //Return value should be partialClustering

    let k = centers.m();

// create edges: care, edge.left stores the index of the gonzalez center (0,...,k-1).
    let mut edges : Vec<Edge> = Vec::with_capacity(k * space.n());
    for (j, c) in centers.get_all(space).iter().enumerate() {
        for p in space.point_iter() {
            edges.push(Edge{
                d : space.dist(c, p),
                left : j as CenterIdx,
                right : p,
            });
        }
    }

//    println!("edges: {:?}", edges);

    // step 1: Compute buckets with limit ceil(4n/k^z) (here z = 2)
    let power_of_k: u32 = 2;
    let mut buckets = put_into_buckets(edges, space.n(), k, power_of_k);

    #[cfg(debug_assertions)]
    println!("\n  - Phase 2a: Put n*k = {} edges into {} buckets, each of size at most ceil(4n/k^{}) = {}.", k*space.n(), buckets.len(), power_of_k, (4*space.n()-1)/k.pow(power_of_k)+1);

    #[cfg(debug_assertions)]
    assert!(buckets::assert_buckets_properties(&buckets, space.n(), k, power_of_k));

    #[cfg(debug_assertions)]
    for (i, bucket) in buckets.iter().enumerate() {
        let bucket_of_dist : Vec<Distance> = bucket.iter().map(|x| x.d).collect();
        println!("\t Bucket {}: {:?}", i, bucket_of_dist);
    }




    // step 2: solve flow problem
    #[cfg(debug_assertions)]
    println!("\n  - Phase 2b: Determine smallest radii and partial assignments that satisfy privacy-condition (privacy_bound = {}).", privacy_bound);

    let mut clusterings: Vec<Clustering> = Vec::with_capacity(k);

    let mut pending = PendingQueues::new();
       // pending[j][t] contains edges from buckets[t] with left = j (which is a center not considered yet)
       // note that t can also be the current bucket, but then pending only contains edges that has
       // a distance small than the settled radius of the previous center, i.e., they can be added
       // fearlessly

    let mut state = initialize_state(space.n(),k, if privacy_bound == 0 {true} else {false});


    let mut i : CenterIdx = 0; // currently processing gonzalez set containing center 0, ..., i; here, i goes from 0 to k-1
    let mut j = 0; // currently processing buckets; from 0,..., k^2-1. We have a shift by -1 compared to paper
    let mut edge_cursor : EdgeIdx = 0; // an cursor pointing to the current edge in the current bucket


//    println!("\n\n************************ Bucket {} ***********************", j);
    while i < k { // extend set of gonzalez centers one by one
        #[cfg(debug_assertions)]
        println!{"\n+++ center {} now considered!\n", i};
        assert!(j < buckets.len());
        // this is the main while-loop that deals with each center set daganzo[i] for i = 1, ..., k

        for b in 0..j {
            // here, we process edges from previous (and the current) buckets (< j) that were ignored because they
            // were adjacent to ceners not in the set that was being processed
            // We do not processed pending edges of the current bucket, as the bucket has changed,
            // and therefore pending[_][l] has been cleared

            while !pending.is_empty(i,b) {
                let e = pending.pop(i,b).unwrap();
                assert_eq!(i, e.left); // e.left should be i, (the index of the i-th center in gonzalez)
                add_edge(e, i, k, privacy_bound, &mut state);
//                println!("  Pending edge added: {:?} (pending from bucket = {}); \tmax_flow: {}", e, c,state.max_flow);
            }
        }


        while state.max_flow < (i + 1) * privacy_bound {
            assert!(j < buckets.len(), "For i = {} no flow of value (i + 1) * privacy_bound = {} found! Panic!", i, (i+1)*privacy_bound);


            if edge_cursor >= buckets[j].len() { //the current bucket has been completed
                j += 1;
//                println!("\n\n************************ Bucket {} ***********************", j);
                assert!(j < buckets.len(), "All buckets have been processed but still not all radii have been settled!");
                edge_cursor = 0;
                continue; //continue the while loop
            }
            let e = buckets[j][edge_cursor];

            let c = e.left; // t is the index of the center on the left of e
            if c > i {
                // in this case the left side (t) corresponds to a center not yet considered, so we
                // postpone the processing of the arc to the point when i == c
                pending.push(c,j,e);
            } else {
                // in this case we do add edge e.
                add_edge(e, i, k, privacy_bound, &mut state);
//                println!("  Adding: {:?} from bucket = {};\tmaxflow: {}", e,j, state.max_flow);
            }
            edge_cursor += 1;


        }


        // at this point, we have identified the bucket that settles the set S_i
        #[cfg(debug_assertions)]

        assert_eq!(state.max_flow, (i + 1) * privacy_bound, "The maximum flow value is bigger than allowed");
        // we should have equality due to the capacities of arcs (= privacy_bound) between the source and the centers in S_i

        // println!("\n+++ Center {} settles in bucket {}:\n", i, j);
        clusterings.push(settle(edge_cursor, &mut buckets[j], i, privacy_bound, &mut state, &centers, space));

        // start bucket from beginning, hence clear all pendings
        edge_cursor = 0;
        for c in i..k {
            pending.clear(c,j);
        }
        i += 1;
    }

    // assigne all unassigned points to its nearest neighbor:
    // This takes O(nk^2) time, so it is bottle-neck. TODO: Can this be improved?
    #[cfg(debug_assertions)]
    println!("\n  - For each partial assignment, assign unassigned points to their nearest center to obtain (full) assignments (= private-preserving clustering).");
    for clustering in clusterings.iter_mut() {
        clustering.fill_up(space);
    }

    clusterings
}


