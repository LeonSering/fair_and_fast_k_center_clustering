

struct ClusterProblem {
    pub k : usize, // maximal number of centers
    pub L : usize, // lower bound of representation
    //pub gamma : u16, // number of colors
    //pub a : Vec<usize>, // lower bounds for color classes
    //pub b : Vec<usize>, // upper bounds for color classes
    
    // method: test if valid: color classes is metric space are 0, ..., gamma-1; sum of a_j <= k.
    //
    // Solving algorithm as method? Maybe not.
}

//use std::collections::HashMap;

mod space;
use crate::space::{ColoredMetric,new_space_by_2dpoints_file};

mod clustering;
//use crate::clustering::{Centers,Clustering};

mod gonzales;
use crate::gonzales::gonzales_heuristic;

mod buckets;
use crate::buckets::{put_into_buckets, assert_buckets_properties};

// An edge for every center (left) to every point (right). note that center appear on both sides.
// The distance is stored in d.
#[derive(Debug,Clone,Copy,PartialOrd,PartialEq)]
pub struct Edge { // Care: The ordering of attributes is important for the partial order! Ties in d are broken by left, then right
    d : f32,
    left : usize,
    right : usize,
}

use std::collections::VecDeque;

#[allow(non_snake_case)]
pub struct State {
    pub center: Vec<Option<usize>>,
    pub reassign: Vec<Vec<VecDeque<usize>>>,
    pub unassigned: Vec<VecDeque<Option<usize>>>,
    pub aux: Vec<Vec<Option<usize>>>,
    pub number_of_covered_points: Vec<usize>,
    pub path_in_H: Vec<Vec<bool>>,
    pub path_in_H_available: Vec<bool>,
    pub max_flow: usize,
}

/////////////////////////////////////////////////////////////////////////
///////////////////////// main() ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

fn main() {


    // load test Space2D from file:
    let space : Box<dyn ColoredMetric> = new_space_by_2dpoints_file("test.2dspace",15);
    let prob = ClusterProblem{  
        k : 5, // number of center;
        L : 7, // number of points to represent;
        //gamma : space.gamma();
    };
    

    /////////////////////////////////
    // phase 1: Gonzalez heuristic //
    /////////////////////////////////

    assert!(prob.k >= 1); // we want to allow at least 1 center
    assert!(space.n() >= prob.k); // the number of points should not be less than the number of centers 

    let (gonzales, gonzales_index_by_center) = gonzales_heuristic(&space, prob.k);

    println!("** Phase 1: Determined k = {} centers by the Gonzales Heuristic: {:?}", prob.k, gonzales.iter().collect::<Vec<&usize>>());
//    println!("index of center {}: {}", gonzales.centers[2], gonzales_index_by_center.get(&gonzales.centers[2]).expect(""));
    gonzales.save_to_file("test.centers");
    
    ///////////////////////////////////////
    // phase 2: determine privacy radius //
    ///////////////////////////////////////
    let mut edges : Vec<Edge> = Vec::with_capacity(prob.k * space.n());
    for c in gonzales.iter() {
        for i in 0..space.n() {
            edges.push(Edge{
                left : *c,
                right : i,
                d : space.dist(*c, i)});
        }
    }

//    println!("edges: {:?}", edges);
    
    // step 1: Compute buckets
    let mut buckets = put_into_buckets(edges, (4*space.n())/prob.k);

    println!("** Phase 2a: Put n*k = {} edges into {} buckets, each of size at most 4n/k = {}.", prob.k*space.n(), buckets.len(), (4*space.n())/prob.k);
   
    #[cfg(debug_assertions)]
    assert!(assert_buckets_properties(&buckets, space.n(), prob.k));
    #[cfg(debug_assertions)]
    for (i, bucket) in buckets.iter().enumerate() {
        let bucket_of_dist : Vec<f32> = bucket.iter().map(|x| x.d).collect(); 
        println!(" {}. bucket: {:?}", i, bucket_of_dist);
    }



    // step 2: solve flow prolbem
    let b = buckets.len(); // number of buckets


    // TODO: think of correct capacities
    let mut pending: Vec<Vec<VecDeque<Edge>>> = (0..prob.k).map(|_| (0..b).map(|_| VecDeque::with_capacity(prob.k*prob.k)).collect()).collect();
       // pending[j][t] contains only edges from gonzales.centers[j] contained in buckets[t]

    let mut state = State {
        center: vec!(None; space.n()), // gives for each point the index (in gonzales array) of the center it is assigned to
        reassign: (0..prob.k).map(|_| (0..prob.k).map(|_| VecDeque::with_capacity(prob.k*prob.k)).collect()).collect(),
        unassigned: (0..prob.k).map(|_| VecDeque::with_capacity(prob.k*prob.k)).collect(),
        aux: (0..space.n()).map(|_| (0..prob.k).map(|_| None ).collect()).collect(), //TODO: need  to be initialised differently
        path_in_H: (0..prob.k).map(|_| (0..prob.k).map(|_| false).collect()).collect(),
        path_in_H_available: vec![false; prob.k],
        number_of_covered_points: vec![0; prob.k],
        max_flow: 0,
    };

    let mut i = 0; // currently processing gonzales set containing center 0, ..., i; here, i goes from 0 to k-1
    let mut j = 0; // currently processing buckets; from 0,..., k^2-1. We have a shift by -1 compared to paper

    while i < prob.k {
        assert!(j < buckets.len());
        // this is the main while-loop that deals with each center set daganzo[i] for i = 1, ..., k 
        for l in 0..j {
            // here, we process edges from previous buckets that were ignored because they
            // were adjacent to ceners not in the set that was being processed
            
            while !pending[i][l].is_empty() {
                let e = pending[i][l].pop_front().unwrap();
                println!("Process pending edge TODO: {:?}", e);
                assert_eq!(i, *gonzales_index_by_center.get(&e.left).expect("Cannot find gonzales index")); // e.left should be the center with index i in gonzales
                process_edge(e, i, j, i, prob.k, prob.L, &mut state);
            }
        }
        if state.max_flow == (i + 1) * prob.L {
            i += 1;
        } else {
            while state.max_flow < (i + 1) * prob.L {
                assert!(j < buckets.len(), "For i = {} no flow of value (i + 1) * L = {} found! Panic!", i, (i+1)*prob.L);
                for t in 0..i {
                    // here, we add back edges of the current bucket that were removed because
                    // they were adjecent to centers not in the set that was being processed
                     
                    while !pending[t][j].is_empty() {
                        let e = pending[t][j].pop_front().unwrap();
                        println!("Add edge to bucket {}: {:?}",j, e);
                        buckets[j].push(e);
                    }
                }
                println!("TODO: Take snapshot");
                for e in buckets[j].iter() {
                    // here, we process the current bucket with index j
                    let t = *gonzales_index_by_center.get(&e.left).expect("Cannot find gonzales index");
//                    println!("Trying to process edge TODO: {:?}, i = {}, t = {}", *e,i, t);
                    if t > i {
                        pending[t][j].push_back(*e);
                    } else {
                        println!("Process edge TODO: {:?}", *e);
                        process_edge(*e, i, j, t, prob.k, prob.L, &mut state);
                    }
                }
                j += 1;
            }
            // at this point, we have identified the bucket that settles the set S_i
            j -= 1;
            println!("TODO settleSet({},{})", i, j);
            i += 1;
        }
    }
    // input: the edge e for processing, the current gonzales set i, the current bucket j, and the
    //        index t of the center of e.left (within the gonzales set)
    fn process_edge(e: Edge, i: usize, j: usize, t: usize, k: usize, L: usize, state: &mut State){
        // the pending part is moved to the main procedure 
        let x = e.right.clone();
        match state.center[x] { 

            None => {// x is not assigned yet

                // a new node correspdonding to x is added to the tail of the queue unassigned:
                state.unassigned[t].push_back(Some(x));

                // the pointer Aux[x][t] points to the index representing x at the queue unassigned[t]
                state.aux[x][t] = Some(state.unassigned[t].len()-1);
            },
            Some(center_of_x) => {
                // in this case, x is already assigned to center[x], so we have to add the new edges in
                // the auxiliary graph H
                state.reassign[center_of_x][t].push_back(x);
                state.aux[x][t] = Some(state.reassign[center_of_x][t].len()-1);
                state.path_in_H[center_of_x][t] = true;

                // we now update reachability status in graph H
                for q in 0..k {
                    state.path_in_H[q][t] = state.path_in_H[q][t] || state.path_in_H[q][center_of_x];
                    if !state.path_in_H_available[q] && state.path_in_H[q][center_of_x] && state.path_in_H_available[t] {
                        state.path_in_H_available[q] = true;
                    }
                }
            }
        }

        // looking for the first center v that has a path and unassigned nodes that could be
        // assigned to v
        let mut v = 0;
        while v < k {
            if !state.unassigned[v].is_empty() && state.path_in_H_available[v] {
                break;
            } else {
                v += 1;
            }
        }
    
        if v == k{
            println!("no agumenting path found");
            return;
        }
        println!("Augmenting exists");
        state.max_flow += 1;
        let y: usize = state.unassigned[v].pop_front().expect("unassigned is empty").expect("unassigned element was None");

        for z in 0..k {
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

        // Now it could be the case that v has too many points, so we need to augment
        
        while state.number_of_covered_points[v] > L {
            
        }

    }

    //////////////////////////////////////////////////////////////
    // phase 3: modify assignment, s.t. sizes are multiple of L //
    //////////////////////////////////////////////////////////////
    

    //////////////////////////////////////////////////////////////////////
    // phase 4: open new centers and  determine actual set of centers C //
    //////////////////////////////////////////////////////////////////////


}
