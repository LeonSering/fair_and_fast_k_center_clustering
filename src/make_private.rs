use std::collections::VecDeque;

use crate::ClusteringProblem;
use crate::clustering::Centers;
use crate::buckets::{put_into_buckets,assert_buckets_properties};

// An edge for every center (left) to every point (right). note that center appear on both sides.
// The distance is stored in d.
#[derive(Debug,Clone,Copy,PartialOrd,PartialEq)]
pub struct Edge { // Care: The ordering of attributes is important for the partial order! Ties in d are broken by left, then right
    pub d : f32,
    left : usize,
    right : usize,
}



pub struct State {
    pub center: Vec<Option<usize>>, // center[x] is the assigned center of point x
    pub reassign: Vec<Vec<VecDeque<usize>>>, // reassing[c][t] contains all points that are assigned to c but could also be assigned to t
    pub unassigned: Vec<VecDeque<Option<usize>>>, // unassigned[c] contains all points that are unassigned but could be assigned to c
    pub aux: Vec<Vec<Option<usize>>>,
    pub number_of_covered_points: Vec<usize>, // the number of points covered by center c; if this equals privacy_bound this center is "private"; if it is smaller than it is "non-private"
    pub path_in_aux_graph: Vec<Vec<bool>>, // indicates whether there is a path in aux_graph from c_1 to c_2; The arcs in aux_graph are excactly the one for which reassign is non-empty; should always be transitiv
    pub path_in_aux_graph_to_non_private: Vec<bool>, // indicates whether a center is non-private or has a path in aux_graph to a non-private center
    pub max_flow: usize, // the current value of a maximum flow

    // method "take_snapshot" and "restore_snapshot" should probably be implemented here.
    // So whenever take_snapshot is called, all later changes are recorded.
    // Calling restore_snapshot rollback all theses changes.
    // I wonder if we need to include "Pending" to the struct.
}


pub fn make_private(prob : ClusteringProblem, gonzales : Centers) { //TODO: Return value should be partialClustering

// create edges: care, edge.left stores the index of the gonzales center (0,...,k).
    let mut edges : Vec<Edge> = Vec::with_capacity(prob.k * prob.space.n());
    for (j, c) in gonzales.iter().enumerate() {
        for i in 0..prob.space.n() {
            edges.push(Edge{
                left : j,
                right : i,
                d : prob.space.dist(*c, i)});
        }
    }

//    println!("edges: {:?}", edges);

    // step 1: Compute buckets
    let buckets = put_into_buckets(edges, (4*prob.space.n())/prob.k);

    println!("** Phase 2a: Put n*k = {} edges into {} buckets, each of size at most 4n/k = {}.", prob.k*prob.space.n(), buckets.len(), (4*prob.space.n())/prob.k);

    #[cfg(debug_assertions)]
    assert!(assert_buckets_properties(&buckets, prob.space.n(), prob.k));
    #[cfg(debug_assertions)]
    for (i, bucket) in buckets.iter().enumerate() {
        let bucket_of_dist : Vec<f32> = bucket.iter().map(|x| x.d).collect();
        println!(" Bucket {}: {:?}", i, bucket_of_dist);
    }



    // step 2: solve flow prolbem
    println!("** Phase 2b: Determine smallest radii that satisfy privacy. Privacy constant privacy_bound = {}.", prob.privacy_bound);
    let b = buckets.len(); // number of buckets


    // TODO: think of correct capacities
    let mut pending: Vec<Vec<VecDeque<Edge>>> = (0..prob.k).map(|_| (0..b).map(|_| VecDeque::with_capacity(prob.k*prob.k)).collect()).collect();
       // pending[j][t] contains only edges from gonzales.centers[j] contained in buckets[t]

    let mut state = State {
        center: vec!(None; prob.space.n()), // gives for each point the index (in gonzales array) of the center it is assigned to
        reassign: (0..prob.k).map(|_| (0..prob.k).map(|_| VecDeque::with_capacity(prob.k*prob.k)).collect()).collect(),
        unassigned: (0..prob.k).map(|_| VecDeque::with_capacity(prob.k*prob.k)).collect(),
        aux: (0..prob.space.n()).map(|_| (0..prob.k).map(|_| None ).collect()).collect(), //TODO: need  to be initialised differently
        path_in_aux_graph: (0..prob.k).map(|_| (0..prob.k).map(|_| false).collect()).collect(),
        path_in_aux_graph_to_non_private: vec![true; prob.k],
        number_of_covered_points: vec![0; prob.k],
        max_flow: 0,
    };

    let mut i = 0; // currently processing gonzales set containing center 0, ..., i; here, i goes from 0 to k-1
    let mut j = 0; // currently processing buckets; from 0,..., k^2-1. We have a shift by -1 compared to paper

    while i < prob.k {
        assert!(j < buckets.len());
        // this is the main while-loop that deals with each center set daganzo[i] for i = 1, ..., k
        for l in 0..j {
            // here, we process edges from previous buckets (< j) that were ignored because they
            // were adjacent to ceners not in the set that was being processed

            while !pending[i][l].is_empty() {
                let e = pending[i][l].pop_front().unwrap();
                println!("Process pending edge TODO: {:?}", e);
                assert_eq!(i, e.left); // e.left should be i, (the index of the i-th center in gonzales)
                process_edge(e, i, prob.privacy_bound, &mut state);
            }
        }
        if state.max_flow >= (i + 1) * prob.privacy_bound {
            #[cfg(debug_assertions)]
            assert_eq!(state.max_flow, (i + 1) * prob.privacy_bound, "The maximum flow value is bigger than allowed"); // we should have equality due to the capacities of arcs (= privacy_bound) between the source and the centers in S_i
            println!("Center {} done.", i);
            i += 1;
        } else {
            while state.max_flow < (i + 1) * prob.privacy_bound {
                assert!(j < buckets.len(), "For i = {} no flow of value (i + 1) * privacy_bound = {} found! Panic!", i, (i+1)*prob.privacy_bound);
                for t in 0..i {
                    // here, we release edges of the current bucket that were "pending" because
                    // they were adjecent to centers not in the gonzales set at the point in time
                    // LEON: not sure if this is needed at all. But I guess having empty
                    // pending-queus is nice.
                    pending[t][j].clear();
                }
                println!("TODO: Take snapshot"); // that means basically to copy everthing (save the flow state) but really copying is too costly.
                for e in buckets[j].iter() {
                    // here, we (try to) process the current bucket with index j
                    let t = e.left;
                    if t > i {
                        // in this case the left side (t) is a center not yet considered, so we
                        // postpone the processing of the arc to the point when i >= t
                        pending[t][j].push_back(*e);
                    } else {
                        // in this case we do the processing.
                        println!("Process edge : {:?}", *e);
                        process_edge(*e, i, prob.privacy_bound, &mut state);
                    }

                    // TODO: I would expect that we need to check the max_flow value after each
                    // edge.
                }
                println!("Bucket {} done.", j);
                j += 1;
            }
            // at this point, we have identified the bucket that settles the set S_i
            j -= 1;
            println!("TODO settleSet({},{})", i, j); // settleSet mean to do binary research on the current bucket
            i += 1;
        }
    }
}

// input: the edge e for processing, the current gonzales set i, the current bucket j, and the
//        index t of the center of e.left (within the gonzales set)
// task: add edge e to the current flow network and look for an augmenting path to increase the
// flow by 1; then execute this augmentation
#[allow(non_snake_case)]
fn process_edge(e: Edge, i: usize, privacy_bound: usize, state: &mut State){
    // the pending part is moved to the main procedure
    let t = e.left;// maybe cloning needed
    let x = e.right; // maybe cloning needed
    match state.center[x] {
        None => {// x is not assigned yet
            // a new node correspdonding to x is added to the tail of the queue unassigned:
            state.unassigned[t].push_back(Some(x));

            // the pointer Aux[x][t] points to the index representing x at the queue unassigned[t]
            state.aux[x][t] = Some(state.unassigned[t].len()-1);
        },
        Some(center_of_x) => {
            // in this case, x is already assigned to center[x], so we have to add the new edges in
            // the auxiliary graph aux_graph
            state.reassign[center_of_x][t].push_back(x);
            state.aux[x][t] = Some(state.reassign[center_of_x][t].len()-1);
            state.path_in_aux_graph[center_of_x][t] = true;

            // we now update reachability status in graph aux_graph
            for q in 0..(i+1) {
                state.path_in_aux_graph[q][t] = state.path_in_aux_graph[q][t] || state.path_in_aux_graph[q][center_of_x];
                if !state.path_in_aux_graph_to_non_private[q] && state.path_in_aux_graph[q][center_of_x] && state.path_in_aux_graph_to_non_private[t] {
                    state.path_in_aux_graph_to_non_private[q] = true;
                }
            }
        }
    }

    // looking for the first center v (in 0,...,i) that has a path and unassigned nodes that could be
    // assigned to v
    let mut v = 0;
    while v <=i {
        if !state.unassigned[v].is_empty() && state.path_in_aux_graph_to_non_private[v] {
            break;
        } else {
            v += 1;
        }
    }

    if v == i+1{
        println!("no agumenting path found");
        return;
    }
    println!("augmenting path exists");
    state.max_flow += 1;

    // y is a point that is unassigned, i.e. there is a free arc from y to the sink. This is
    // our first arc of the augmenting path (from sink to source)
    let mut y: usize = state.unassigned[v].pop_front().expect("unassigned is empty").expect("unassigned element was None");


    // some update in updates to aux and reassign are needed (not sure yet)
    for z in 0..(i+1) {
        if z != v {
            match state.aux[y][z] {
                Some(_) => {
                    state.aux[y][z] = None;
                    state.reassign[v][z].push_back(y);
                }
                None => {}
            }
        }
    }
    // TODO: Someting with AUX

    state.center[y] = Some(v); // assign y to v
    state.number_of_covered_points[v] += 1; // v covers now one points more

    // Now it could be the case that v is private already (covers privacy_bound points already), so we need to find a new center

    while state.number_of_covered_points[v] > privacy_bound { // while v is overfull
//            #[cfg(debug_assertions)]
//            assert_eq!(state.number_of_covered_points[v], privacy_bound, "Center {} covers too many points, namely {}", v, state.number_of_covered_points[v]);

        // need to find new center w, such that arc (v,w) is in aux_graph
        // TODO: maybe we need to be more carful, s.t. we don't run in circles
        let mut w = 0;
        while w <= i {
            if !state.reassign[v][w].is_empty() {
                break;
            } else {
                w += 1;
            }
        }
        if w == i+1 {
            panic!("There is no center w such that (v, w) is in aux_graph. This contradicts the fact that B[v] == True; v = {}", v);
        }
        // w is our next center in the augmenting path. We reassign y, which means adding
        // forward arc (w,y) and backwards arc (y, v) in front of our augmenting path
        y = state.reassign[v][w].pop_front().expect("unassigned is empty");
        state.center[y] = Some(w);
        state.number_of_covered_points[w] += 1;
        // TODO: now we need to fix the reassign queues.
        // and withit the complete aux_graph structure; probably aux is needed.
        v = w;

    }

}
