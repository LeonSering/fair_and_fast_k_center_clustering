use crate::ClusteringProblem;
use crate::space::{Point,ColoredMetric};
use crate::clustering::{Clustering,Centers};

mod buckets;
use buckets::{put_into_buckets,assert_buckets_properties};

mod flow;
use flow::{initialize_state,add_edge};

mod settle;
use settle::settle;

//TEMP:
pub mod with_sorting;

use std::collections::VecDeque;

// An edge in the flow network. One for every center (left; in form of an gonzales index 0,..,k-1) to every point (right). note that center appear on both sides.
// The distance is stored in d.
#[derive(Debug,Clone,Copy,PartialOrd,PartialEq)]
pub struct Edge<'a> { // Care: The ordering of attributes is important for the partial order! Ties in d are broken by left, then right
    pub d : f32,
    left : usize, // index of the center in gonzales
    right : &'a Point,
}


// centers_graph: nodes = centers 0,...,i (i goes from 0 to k-1) and we have an arc c_1 to c_2 iff
// there is a point x that is assigned to c_1 but could also be assigne to c_2
// centers_graph is basically the connectivity within the residual graph on the left side of the
// flow network.


/// Given a metric space and a ClusteringProblem, make_private takes a set of ordered centers
/// 0,...,k-1, and determines for each prefix of centers (0,...,i) a partial clustering with minimal radius that satisfy the
/// privacy constraint, i.e., each center (0,...,i) covers exactly privacy_bound many points. 
pub fn make_private<'a>(space : &'a Box<dyn ColoredMetric>, prob : &'a ClusteringProblem, gonzales : &Centers<'a>) -> Vec<Clustering<'a>> { //Return value should be partialClustering

// create edges: care, edge.left stores the index of the gonzales center (0,...,k-1).
    let mut edges : Vec<Edge> = Vec::with_capacity(prob.k * space.n());
    for (j, c) in gonzales.iter().enumerate() {
        for p in space.point_iter() {
            edges.push(Edge{
                d : space.dist(c, p),
                left : j,
                right : p,
            });
        }
    }

//    println!("edges: {:?}", edges);

    // step 1: Compute buckets with limit ceil(4n/k^z) (here z = 1)
    let power_of_k: u32 = 1;
    let mut buckets = put_into_buckets(edges, space.n(), prob.k, power_of_k);

    println!("** Phase 2a: Put n*k = {} edges into {} buckets, each of size at most ceil(4n/k^{}) = {}.", prob.k*space.n(), buckets.len(), power_of_k, (4*space.n()-1)/prob.k+1);

    #[cfg(debug_assertions)]
    assert!(assert_buckets_properties(&buckets, space.n(), prob.k, power_of_k));

    #[cfg(debug_assertions)]
    for (i, bucket) in buckets.iter().enumerate() {
        let bucket_of_dist : Vec<f32> = bucket.iter().map(|x| x.d).collect();
        println!(" Bucket {}: {:?}", i, bucket_of_dist);
    }




    // step 2: solve flow problem
    println!("** Phase 2b: Determine smallest radii and assignment that satisfy privacy. Privacy constant privacy_bound = {}.", prob.privacy_bound);
    
    let mut clusterings: Vec<Clustering> = Vec::with_capacity(prob.k);
    

    let b = buckets.len(); // number of buckets


    // TODO: think of correct Vec_capacities
    let mut pending: Vec<Vec<VecDeque<Edge>>> = (0..prob.k).map(|_| (0..b).map(|_| VecDeque::with_capacity(prob.k*prob.k)).collect()).collect();
       // pending[j][t] contains edges from buckets[t] with left = j (which is a center not considered yet)
       // note that t can also be the current bucket, but then pending only contains edges that has
       // a distance small than the settled radius of the previous center, i.e., they can be added
       // fearlessly
    let mut state = initialize_state(space.n(),prob.k);


    let mut i = 0; // currently processing gonzales set containing center 0, ..., i; here, i goes from 0 to k-1
    let mut j = 0; // currently processing buckets; from 0,..., k^2-1. We have a shift by -1 compared to paper
    let mut edge_cursor = 0; // an cursor pointing to the current edge in the current bucket


    println!("\n\nMAKE_PRIVTE with L = {}", prob.privacy_bound);
    println!("\n\n************************ Bucket {} ***********************", j);
    while i < prob.k { // extend set of gonzales centers one by one
        println!{"\n+++ center {} now considered!\n", i};
        assert!(j < buckets.len());
        // this is the main while-loop that deals with each center set daganzo[i] for i = 1, ..., k
        for c in 0..j {
            // here, we process edges from previous (and the current) buckets (< j) that were ignored because they
            // were adjacent to ceners not in the set that was being processed
            // We do not processed pending edges of the current bucket, as the bucket has changed,
            // and therefore pending[_][l] has been cleared

            while !pending[i][c].is_empty() {
                let e = pending[i][c].pop_front().unwrap();
                assert_eq!(i, e.left); // e.left should be i, (the index of the i-th center in gonzales)
                add_edge(e, i, prob, &mut state);
                println!("  Pending edge added: {:?} (pending from bucket = {}); \tmax_flow: {}", e, c,state.max_flow);
            }
        }


        while state.max_flow < (i + 1) * prob.privacy_bound {
            assert!(j < buckets.len(), "For i = {} no flow of value (i + 1) * privacy_bound = {} found! Panic!", i, (i+1)*prob.privacy_bound);

          
            if edge_cursor >= buckets[j].len() { //the current bucket has been completed
                j += 1;
                println!("\n\n************************ Bucket {} ***********************", j);
                assert!(j < buckets.len(), "All buckets have been processed but still not all radii have been settled!");
                edge_cursor = 0; 
                continue; //continue the while loop
            }
            let e = buckets[j][edge_cursor];

            let c = e.left; // t is the index of the center on the left of e
            if c > i {
                // in this case the left side (t) corresponds to a center not yet considered, so we
                // postpone the processing of the arc to the point when i == c
                pending[c][j].push_back(e);
            } else {
                // in this case we do add edge e.
                add_edge(e, i, prob, &mut state);
                println!("  Adding: {:?} from bucket = {};\tmaxflow: {}", e,j, state.max_flow);
            }
            edge_cursor += 1;

            
        }


        // at this point, we have identified the bucket that settles the set S_i
        #[cfg(debug_assertions)]
        assert_eq!(state.max_flow, (i + 1) * prob.privacy_bound, "The maximum flow value is bigger than allowed"); // we should have equality due to the capacities of arcs (= privacy_bound) between the source and the centers in S_i
        println!("\n+++ Center {} settles in bucket {}:\n", i, j);
        //TODO Settle centers 0...i in bucket j
        clusterings.push(settle(edge_cursor, &mut buckets[j], i, prob, &mut state, &gonzales));
        
        // start bucket from beginning, hence clear all pendings
        edge_cursor = 0;
        for t in i..prob.k {
            pending[t][j].clear();
        }
        i += 1;
    }
    
    clusterings
}




