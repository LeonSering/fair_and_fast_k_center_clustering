struct ClusterProblem{
    // metric space
    // number of colors (gamma)
    // number of points n
    // max number of clusters k
    // for each color j: upper and lower bounds for centers a_j, b_j
    // representative loweder bound L
    //
    // method: test if valid: color classes is metric space are 0, ..., gamma-1; sum of a_j <= k.
    //
    // Solving algorithm as method? Maybe not.
}

use std::collections::HashMap;
fn gonzales_heuristic(space : &Box<dyn ColoredMetric>, k : usize) -> (Centers, HashMap<usize,usize>) {

    let mut gonzales : Centers = Centers{centers : Vec::with_capacity(k)};
    let mut gonzales_index_by_center: HashMap<usize, usize> = HashMap::with_capacity(k);

    // we can add any point as first center, so lets take 0
    gonzales.centers.push(0);

    let mut dist_x_center : Vec<f32> = Vec::with_capacity(space.n()); // current distance of point x to the set of already determined centers

    for i in 0..space.n() {
        dist_x_center.push(space.dist(0,i)); // initilized to the distance of first center: 0.
        gonzales_index_by_center.insert(0,0);
    }

    for i in 1..k {
        let mut current_distance = std::f32::MIN; // maximal distance to set of centers
        let mut current_point : Option<usize> = None; // corresponing point with this max distance.
        for j in 0..space.n(){
            let dist_to_newest_center = space.dist(j, gonzales.centers[i-1]); // as distance of j to gonzales 0..i-2 is known, we only need to measure distance to newest center i-1.
            
            // update dist_x_S to now include the newest center:
            if dist_to_newest_center < dist_x_center[j] {
                dist_x_center[j] = dist_to_newest_center;
            }

            // check whether dist(j, gonzales 0..i-2) is bigger than current biggest dist.
            if dist_x_center[j] > current_distance {
                current_distance = dist_x_center[j];
                current_point = Some(j);
            }
        }

        // create new center vector including the current farthest point (= current_point):
        gonzales.centers.push(current_point.expect("No new center could be found, i.e., current_point = None"));
        gonzales_index_by_center.insert(current_point.expect("No new center"),i);
    }
    (gonzales, gonzales_index_by_center)
}


// An edge for every center (left) to every point (right). note that center appear on both sides.
// The distance is stored in d.
#[derive(Debug,Clone,Copy,PartialOrd,PartialEq)]
pub struct Edge { // Care: The ordering of attributes is important for the partial order! Ties in d are broken by left, then right
    d : f32,
    left : usize,
    right : usize,
}

use std::collections::VecDeque;
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

mod space;
pub use crate::space::*;

mod clustering;
pub use crate::clustering::*;

mod buckets;
pub use crate::buckets::*;

/////////////////////////////////////////////////////////////////////////
///////////////////////// main() ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

fn main() {
    use space::*;


    // load test Space2D from file:
    let space : Box<dyn ColoredMetric> = new_space_by_2dpoints_file("test.2dspace",15);
   
    let k : usize = 5; // number of center;
    let L : usize = 7; // number of points to represent;
    

    /////////////////////////////////
    // phase 1: Gonzalez heuristic //
    /////////////////////////////////

    assert!(k >= 1); // we want to allow at least 1 center
    assert!(space.n() >= k); // the number of points should not be less than the number of centers 

    let (gonzales, gonzales_index_by_center) = gonzales_heuristic(&space, k);

    println!("** Phase 1: Determined k = {} centers by the Gonzales Heuristic: {:?}", k, gonzales.centers);
//    println!("index of center {}: {}", gonzales.centers[2], gonzales_index_by_center.get(&gonzales.centers[2]).expect(""));
//    gonzales.save_to_file("test.centers");
    
    ///////////////////////////////////////
    // phase 2: determine privacy radius //
    ///////////////////////////////////////
    let mut edges : Vec<Edge> = Vec::with_capacity(k * space.n());
    for c in gonzales.centers.iter() {
        for i in 0..space.n() {
            edges.push(Edge{
                left : *c,
                right : i,
                d : space.dist(*c, i)});
        }
    }

//    println!("edges: {:?}", edges);
    
    // step 1: Compute buckets
    let mut buckets = put_into_buckets(edges, (4*space.n())/k);

    println!("** Phase 2a: Put n*k = {} edges into {} buckets, each of size at most 4n/k = {}.", k*space.n(), buckets.len(), (4*space.n())/k);
   
    #[cfg(debug_assertions)]
    assert!(assert_buckets_properties(&buckets, space.n(), k));
    #[cfg(debug_assertions)]
    for (i, bucket) in buckets.iter().enumerate() {
        let bucket_of_dist : Vec<f32> = bucket.iter().map(|x| x.d).collect(); 
        println!(" {}. bucket: {:?}", i, bucket_of_dist);
    }



    // step 2: solve flow prolbem
    let b = buckets.len(); // number of buckets


    // TODO: think of correct capacities
    let mut pending: Vec<Vec<VecDeque<Edge>>> = (0..k).map(|_| (0..b).map(|_| VecDeque::with_capacity(k*k)).collect()).collect();
       // pending[j][t] contains only edges from gonzales.centers[j] contained in buckets[t]


    let mut state = State {
        center: vec!(None; space.n()), // gives for each point the index (in gonzales array) of the center it is assigned to
        reassign: (0..k).map(|_| (0..k).map(|_| VecDeque::with_capacity(k*k)).collect()).collect(),
        unassigned: (0..k).map(|_| VecDeque::with_capacity(k*k)).collect(),
        aux: (0..space.n()).map(|_| (0..k).map(|_| None ).collect()).collect(), //TODO: need  to be initialised differently
        path_in_H: (0..k).map(|_| (0..k).map(|_| false).collect()).collect(),
        path_in_H_available: vec![false; k],
        number_of_covered_points: vec![0; k],
        max_flow: 0,
    };

    let mut i = 0; // currently processing gonzales set containing center 0, ..., i; here, i goes from 0 to k-1
    let mut j = 0; // currently processing buckets; from 0,..., k^2-1. We have a shift by -1 compared to paper

    while i < k {
        assert!(j < buckets.len());
        // this is the main while-loop that deals with each center set daganzo[i] for i = 1, ..., k 
        for l in 0..j {
            // here, we process edges from previous buckets that were ignored because they
            // were adjacent to ceners not in the set that was being processed
            
            while !pending[i][l].is_empty() {
                let e = pending[i][l].pop_front().unwrap();
                println!("Process pending edge TODO: {:?}", e);
                assert_eq!(i, *gonzales_index_by_center.get(&e.left).expect("Cannot find gonzales index")); // e.left should be the center with index i in gonzales
                process_edge(e, i, j, i, k, L, &mut state);
            }
        }
        if state.max_flow == (i + 1) * L {
            i += 1;
        } else {
            while state.max_flow < (i + 1) * L {
                assert!(j < buckets.len(), "For i = {} no flow of value (i + 1) * L = {} found! Panic!", i, (i+1)*L);
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
                        process_edge(*e, i, j, t, k, L, &mut state);
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
