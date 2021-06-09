use crate::ClusteringProblem;
use crate::space::{Point,ColoredMetric};
use crate::clustering::{Clustering,Centers,new_centers};

mod buckets;
use buckets::{put_into_buckets,split_at_median,assert_buckets_properties};

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

// describes the state of the flow network and of the auxiliray graph.
#[derive(Debug)]
struct State<'a>{
    pub center_of: Vec<Option<usize>>, // center[x] is the assigned center (in form of gonzales index) of point x
    reassign: Vec<Vec<VecDeque<&'a Point>>>, // reassing[c][t] contains all points that are assigned to c but could also be assigned to t
    unassigned: Vec<VecDeque<&'a Point>>, // unassigned[c] contains all points that are unassigned but could be assigned to c
    // reassign and unassigned might not be correctly updated at all times. E.g., a point in
    // reassign might not be assigned to c anymore, or not be assignable to t; 
    // or a point in unassigned might be assigned to some point or not assignable to c. 
    // So whenever an element is poped (or we check for is_empty()) it needs to be checked if it still makes sense.
    pub number_of_covered_points: Vec<usize>, // the number of points covered by center c; if this equals privacy_bound this center is "private"; if it is smaller than it is "non-private".; if the has more than that it is "overfull"
    pub path_in_centers_graph: Vec<Vec<bool>>, // indicates whether there is a path in centers_graph from c_1 to c_2; 
    pub path_in_centers_graph_to_non_private: Vec<bool>, // indicates whether a center is non-private or has a path in centers_graph to a non-private center
    pub max_flow: usize, // the current value of a maximum flow
    pub edge_present: Vec<Vec<bool>> // edge_present[c][p] denotes whether the arc from center c to point p is currently present in the flow network 
}
impl<'a> State<'a> {
    pub fn reassign_pop(&mut self, c1 : usize, c2: usize) -> Option<&'a Point> {
        self.clean_reassign(c1,c2);
        self.reassign[c1][c2].pop_front()
    }

    pub fn reassign_push(&mut self, c1 : usize, c2: usize, p: &'a Point) {
        self.reassign[c1][c2].push_back(p);
    }

//    pub fn reassing_peak(&mut self, c1: usize, c2: usize) -> Option<&&'a Point> {
//        self.clean_reassign(c1,c2);
//        self.reassign[c1][c2].front()
//    }

    pub fn reassign_is_empty(&mut self, c1: usize, c2: usize) -> bool {
        self.clean_reassign(c1,c2);
        self.reassign[c1][c2].is_empty()
    }

    fn clean_reassign(&mut self, c1: usize, c2: usize) { // removes invalid entries at the front
        while !self.reassign[c1][c2].is_empty() {
            let p = self.reassign[c1][c2].front().unwrap().idx();
            if self.center_of[p].is_none() || self.center_of[p].unwrap() != c1 || !self.edge_present[c2][p] {
                // if p is not assigned to c1 or p could not be assigned to c2..
                self.reassign[c1][c2].pop_front(); //.. discard entry
            } else {
                break;
            }
        }
    }


    pub fn unassigned_pop(&mut self, c : usize) -> Option<&'a Point> {
        self.clean_unassigned(c);
        self.unassigned[c].pop_front()
    }

    pub fn unassigned_push(&mut self, c : usize, p: &'a Point) {
        self.unassigned[c].push_back(p);
    }

//    pub fn unassigned_peak(&mut self, c: usize) -> Option<&&'a Point> {
//        self.clean_unassigned(c);
//        self.unassigned[c].front()
//    }

    pub fn unassigned_is_empty(&mut self, c: usize) -> bool {
        self.clean_unassigned(c);
        self.unassigned[c].is_empty()
    }
    
    fn clean_unassigned(&mut self, c: usize) { // removes invalid entries at the front
        while !self.unassigned[c].is_empty() {
            let p = self.unassigned[c].front().unwrap().idx();
            if self.center_of[p].is_some() || !self.edge_present[c][p] {
                // if p is assigned to some center or p could not be assigned to c..
                self.unassigned[c].pop_front(); //.. discard entry
            } else {
                break;
            }
        }
    }
}


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

    // step 1: Compute buckets with limit ceil(4n/k)
    let mut buckets = put_into_buckets(edges, (4*space.n()-1)/prob.k + 1);

    println!("** Phase 2a: Put n*k = {} edges into {} buckets, each of size at most ceil(4n/k) = {}.", prob.k*space.n(), buckets.len(), (4*space.n()-1)/prob.k+1);

    #[cfg(debug_assertions)]
    assert!(assert_buckets_properties(&buckets, space.n(), prob.k));
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
    let mut state = State {
        center_of: vec!(None; space.n()), // gives for each point the index (in gonzales array) of the center it is assigned to; at the beginning all are unassigned (= None)
        reassign: (0..prob.k).map(|_| (0..prob.k).map(|_| VecDeque::with_capacity(prob.k*prob.k)).collect()).collect(),
        unassigned: (0..prob.k).map(|_| VecDeque::with_capacity(prob.k*prob.k)).collect(),
        path_in_centers_graph: (0..prob.k).map(|i| (0..prob.k).map(|j| if i == j {true} else {false}).collect()).collect(), // in the beginning there are no arcs in centers_graph
        path_in_centers_graph_to_non_private: vec![true; prob.k], // in the beginning all centers are non_private
        number_of_covered_points: vec![0; prob.k],
        max_flow: 0,
        edge_present: vec!(vec!(false; space.n()); prob.k),
    };


    let mut i = 0; // currently processing gonzales set containing center 0, ..., i; here, i goes from 0 to k-1
    let mut j = 0; // currently processing buckets; from 0,..., k^2-1. We have a shift by -1 compared to paper
    let mut edge_cursor = 0; // an cursor pointing to the current edge in the current bucket


    println!("\n\nMAKE_PRIVTE with L = {}", prob.privacy_bound);
    println!("\n\n************************ Bucket {} ***********************", j);
    while i < prob.k { // extend set of gonzales centers one by one
        println!{"\n+++ center {} now considered!\n", i};
        assert!(j < buckets.len());
        // this is the main while-loop that deals with each center set daganzo[i] for i = 1, ..., k
        for l in 0..j {
            // here, we process edges from previous (and the current) buckets (< j) that were ignored because they
            // were adjacent to ceners not in the set that was being processed
            // We do not processed pending edges of the current bucket, as the bucket has changed,
            // and therefore pending[_][l] has been cleared

            while !pending[i][l].is_empty() {
                let e = pending[i][l].pop_front().unwrap();
                assert_eq!(i, e.left); // e.left should be i, (the index of the i-th center in gonzales)
                add_edge(e, i, prob, &mut state);
                println!("  Pending edge added: {:?} (pending from bucket = {}); \tmax_flow: {}", e, l,state.max_flow);
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

            let t = e.left; // t is the index of the center on the left of e
            if t > i {
                // in this case the left side (t) corresponds to a center not yet considered, so we
                // postpone the processing of the arc to the point when i == t
                pending[t][j].push_back(e);
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



// input: the edge e for processing, the current gonzales set (0,...,i), the current bucket j, and the
//        index t of the center of e.left (within the gonzales set)
// task: add edge e to the current flow network and look for an augmenting path to increase the
// flow by 1; then execute this augmentation
#[allow(non_snake_case)]
fn add_edge<'a>(e: Edge<'a>, i: usize, prob: &ClusteringProblem, state: &mut State<'a>){
    let t = e.left;
    let x = e.right; // maybe cloning needed
    
    if t > i { // this edge is not considered yet
//        println!("\tnot yet! Max flow: {}\n", state.max_flow);
        return;
    }

    // as the edge is now added to the flow networks we mark it that way:
    assert_eq!(state.edge_present[t][x.idx()], false, "Edge was present before");
    state.edge_present[t][x.idx()] = true;

    match state.center_of[x.idx()] {
        None => {// x is not assigned yet
            // a new node correspdonding to x is added to the tail of the queue unassigned:
            state.unassigned_push(t,x);
//            println!("\tcurrent point is unassigned");
        },
        Some(center_of_x) => {
            // in this case, x is already assigned to center[x], so we have to add the new edges in
            // the centers graph 
            state.reassign_push(center_of_x,t,x);
//            println!("\tcurrent point is assigned to {}", center_of_x);

            // There is now a path from center_of_x to t:
            state.path_in_centers_graph[center_of_x][t] = true;

            // update reachability status in centers graph
            for q in 0..(i+1) {
                state.path_in_centers_graph[q][t] = state.path_in_centers_graph[q][t] || state.path_in_centers_graph[q][center_of_x];
                state.path_in_centers_graph_to_non_private[q] = state.path_in_centers_graph_to_non_private[q] || (state.path_in_centers_graph[q][center_of_x] && state.path_in_centers_graph_to_non_private[t]);
            }
        }
    }
    augment_flow(prob,i,state);
//    println!("\tmax_flow: {}\n", state.max_flow);
}

fn remove_edge<'a>(e: Edge<'a>, i: usize, prob: &ClusteringProblem, state: &mut State<'a>) {
    let t = e.left;
    let x = e.right;
    
    if t > i { // this edge has not been considered yet, and hence is not in the network
//        println!("\n\tnot yet! Max flow: {}", state.max_flow);
        return;
    }
    
    // as the edge is now removed from the flow networks we mark it that way:
    assert_eq!(state.edge_present[t][x.idx()], true, "Edge was present before");
    state.edge_present[t][x.idx()] = false;

//    println!("x: {}, center of x: {:?}", x.idx(), state.center_of[x.idx()]);
    match state.center_of[x.idx()] {
        None => {
//            println!{"\tarc was not flow carrying (point not assigned)!"};

        }
        Some(center_of_x) => {
            if center_of_x == t { // arc e is flow carrying:
                //println!{"\tarc is flow carrying!"};
                state.max_flow -= 1;
                state.number_of_covered_points[t] -= 1;
                state.center_of[x.idx()] = None;
                
                // this unassignment indirectly also updates reassign, because whenever something is
                // poped from reassign[t][_] we check whether that x is assigned to t 
                // also x cannot be part of unassigned[t] at this point, because the moment it was
                // assigned to t it was poped from unassigned[t] (and it can only enter this queue
                // again, when edge e is added to the network)
                build_centers_graph(prob, state); // TODO: This is too expensive
                augment_flow(prob, i, state); // see if we can add a flow unit along other paths
            } else {
//                println!{"\tarc was not flow carrying!"};
            }
        }
    }
}

fn augment_flow<'a>(prob: &ClusteringProblem, i: usize, state: &mut State<'a>) {
    // We need to check if there is an augmenting path to increase the max_flow:

    // looking for the first center v (in 0,...,i) that has a path to a non_private center 
    // and has unassigned nodes that could be assigned to v
    let mut v = 0;
    while v <=i {
        if !state.unassigned_is_empty(v) && state.path_in_centers_graph_to_non_private[v] {
            break;
        } else {
            v += 1;
        }
    }

    if v == i+1{
        // in this case there is no augmenting path
//        println!("\tno agumenting path found");
        return;
    }

    // augmenting path exists
//    println!("\taugmenting path exists");

    // p is a point that is unassigned, i.e. there is a free arc from p to the sink. This is
    // our first arc of the augmenting path (from sink to source)
    let mut p = state.unassigned_pop(v).expect("unassigned should not be empty");
    
    assert_eq!(state.center_of[p.idx()], None, "WAIT: p is not unassigned");

    state.center_of[p.idx()] = Some(v); // assign p to v
    state.number_of_covered_points[v] += 1; // v covers now one points more
    state.max_flow += 1;

    // Now it could be the case that v is private already (covers privacy_bound many points already), so we need to find a new center
    // for this we do a BFS on the centers_graph (but only on the centers that has a path to a non_private center)
  
    // for the BFS we save the hops needed to reach a non_private center:
    // This implementation takes O(k^3) time
    let mut hops_to_non_private: Vec<Option<usize>> = (0..i+1).map(|c| if state.number_of_covered_points[c] < prob.privacy_bound {Some(0)} else {None} ).collect();
    for _ in 1..i+1 { // we have a maximum of i+1 hops
        for c1 in 0..i+1 {
            for c2 in 0..i+1 {

                if !state.reassign_is_empty(c1,c2) {
                    hops_to_non_private[c1] = match hops_to_non_private[c2] {
                        None => hops_to_non_private[c1],
                        Some(hops_c2) => match hops_to_non_private[c1] {
                            None => Some(hops_c2 + 1),
                            Some(hops_c1) => if hops_c2 + 1 < hops_c1 {Some(hops_c2 + 1)} else { Some(hops_c1)},
                        }
                    }
                }
            }
        }
    }

//    println!("hops to non_private: {:?}", hops_to_non_private);
//    println!("number_of_covered_points {:?}", state.number_of_covered_points);
//    
//    println!("reassign: {:?}", state.reassign);
//    println!("path_in_centers_graph: {:?}", state.path_in_centers_graph);
             
    while state.number_of_covered_points[v] > prob.privacy_bound { // while v is overfull

        // need to find new center w, such that arc (v,w) is in centers_graph, at w is closer to
        // non-private center
        let mut w = 0;
        while w <= i {
            if hops_to_non_private[w] == None || hops_to_non_private[w].unwrap() > hops_to_non_private[v].unwrap() - 1 { // w is not closer to a non-private center
                w += 1;
                continue;
            }


            if !state.reassign_is_empty(v,w) { // if there are still a elements left, there is an arc v w in centers_graph, and w leads to a non_private center
                break;
            } else {
                w += 1;
            }
        }
        assert!(w <= i,"There is no center w such that (v, w) is in centers_graph and w leads to a non_private center. This contradicts the fact that B[v] == True; v = {}", v);

        // w is our next center in the augmenting path. We reassign y, which means adding
        // forward arc (w,y) and backwards arc (y, v) in front of our augmenting path
        p = state.reassign_pop(v,w).expect("unassigned is empty");
        assert_eq!(state.center_of[p.idx()], Some(v), "WAIT: p was not assigned to v!");
        state.number_of_covered_points[v] -= 1;
        state.center_of[p.idx()] = Some(w); // reassign y to w
        state.reassign_push(w,v,p);
        state.number_of_covered_points[w] += 1;
        
        v = w;

    }
    
    // Since the connectivity in the center_graph has changed we need to update to corresponding
    // data_structures:
    build_centers_graph(prob,state); // takes O(k^3) time
    // this is done way too often. Is it possible to update these data_structures instead of
    // rebuilding? Maybe save the arcs that are added/removed.
}


// rebuilds the data structures path_in_centers_graph and path_in_centers_graph_to_non_private
// takes O(k^3) time
fn build_centers_graph(prob: &ClusteringProblem, state: &mut State) {
    // we need to rebuild path_in_centers_graph and path_in_centers_graph_to_non_private
    // TODO: This is not done properly, it takes k^3
    state.path_in_centers_graph_to_non_private = (0..prob.k).map(|c| state.number_of_covered_points[c] < prob.privacy_bound).collect();
    state.path_in_centers_graph = (0..prob.k).map(|c1| (0..prob.k).map(|c2| if c1 == c2 {true} else {false}).collect()).collect();
    for c1 in 0..prob.k {
        for c2 in 0..prob.k {
            if !state.reassign_is_empty(c1,c2) { // there is an arc from c1 to c2:
                for c3 in 0..prob.k {
                    state.path_in_centers_graph[c1][c3] = state.path_in_centers_graph[c1][c3] || state.path_in_centers_graph[c2][c3];
                    state.path_in_centers_graph_to_non_private[c1] = state.path_in_centers_graph_to_non_private[c1] || (state.path_in_centers_graph[c1][c3] && state.path_in_centers_graph_to_non_private[c3]);
                }
            }
        }
    }

}

// note that edge_cursor points at the edge that has not been added yet 
// (the edge at edge_curser-1 has been added already)
fn settle<'a, 'b>(edge_cursor: usize, bucket: &mut Vec<Edge<'a>>, i: usize, prob: &ClusteringProblem, state: &mut State<'a>, gonzales: &Centers<'b>) -> Clustering<'b>{
    let mut cursor = edge_cursor;

    println!("  Edge_cursor: {}", edge_cursor);
    println!("  Bucket to settle: {:?}\n", bucket.iter().map(|x| x.d).collect::<Vec<f32>>());
    
    let edges_present : bool;
    // first clear or fill bucket:
    if edge_cursor >= bucket.len()/2 { // in this case fill the bucket
        while cursor < bucket.len() {
            add_edge(bucket[cursor],i,prob,state);
            println!("\tInitial adding: {:?} in binary search for center: {};\tmax_flow: {}", bucket[cursor], i, state.max_flow);
            cursor += 1;
        }
        edges_present = true;
    } else { // in this case empty bucket
        while cursor > 0 {
            cursor -= 1;
            remove_edge(bucket[cursor],i,prob,state);
            println!("\tInitial removing: {:?} in binary search for center: {};\tmax_flow: {}", bucket[cursor], i, state.max_flow);
        }
        edges_present = false;
    }

    let radius = search_for_radius(edges_present, bucket, &mut cursor, i, prob, state);

    println!("\n==> Radius for center {} found: {};\tmax_flow: {}\n", i, radius, state.max_flow);
//    println!("center_of: {:?}", state.center_of);
//    println!("number_of_points_covered: {:?}", state.number_of_covered_points);

    let mut centers = new_centers(i+1);
    for c in gonzales.iter().take(i+1) {
        centers.push(c);
    }

    let clustering = Clustering {
        centers,  // centers are the gonzales centers 0,...,i
        radius, // TODO
        center_of : state.center_of.iter().map(|c| match c {
                        Some(idx) => Some(gonzales.get(*idx)),
                        None => None, 
                    }).collect(),
    };

    println!("  Bucket after settling: {:?}\n", bucket.iter().map(|x| x.d).collect::<Vec<f32>>());

    // empty bucket for the final time
    while cursor > 0 {
        cursor -= 1;
        remove_edge(bucket[cursor],i,prob,state);
        println!("\tClean up removing: {:?} in binary search for center: {};\tmaxflow: {}", bucket[cursor], i,state.max_flow);
    }

    clustering
}

// search recursevely (binary research) for the correct radius.
// Input: the edge list, current center list 0,...,i, the cluster prob, the flow state
// If edges_present == TRUE: expects that ALL edges in the list are part of the flow network
// If edges_present == FALSE: expects that NONE of the edges in the list are part of the flow network
// and expects that state.max_flow < max_flow_target
//
// Output: None if max_flow_target cannot be reached or Some(radius) if found
//   list becomes more sorted
fn search_for_radius<'a>(edges_present: bool, list: &mut Vec<Edge<'a>>, cursor : &mut usize, i : usize, prob: &ClusteringProblem, state: &mut State<'a>) -> f32 {
    println!("\n  List to settle: {:?} edges_present: {}", list.iter().map(|x| x.d).collect::<Vec<f32>>(), edges_present);
    let list_len = list.len();
    assert!(list_len > 0, "Empty list in binary search");
    if list_len == 1 {
        if !edges_present {
            add_edge(list[0],i,prob,state);
            *cursor += 1;
            println!("\n\tAdd final edge: {:?} in binary search for center: {};\tmax_flow: {}", list[0], i, state.max_flow);
        }
        assert_eq!(state.max_flow, (i+1) * prob.privacy_bound, "Something went wrong in the binary search.");
        return list[0].d;
    }

    let (mut smaller, mut bigger) = split_at_median(list);
    println!("     smaller: {:?} bigger: {:?}\n", smaller.iter().map(|x| x.d).collect::<Vec<f32>>(), bigger.iter().map(|x| x.d).collect::<Vec<f32>>());



    // take care that smaller edges are added and bigger edges are not present in the flow network
    if edges_present {
        for e in bigger.iter().rev() {
            remove_edge(*e, i, prob, state);
            println!("\tTry removing: {:?} in binary search for center: {};\tmax_flow: {}", e, i,state.max_flow);
            *cursor -= 1;
        }
    } else {
        for e in smaller.iter() {
            add_edge(*e, i, prob, state);
            println!("\tTry adding: {:?} in binary search for center: {};\tmax_flow: {}", e, i,state.max_flow);
            *cursor += 1;
        }
    }


    let radius: f32;
//    let mut left : Vec<Edge> = Vec::with_capacity(list_len);
//    let mut right : Vec<Edge> = Vec::with_capacity(list_len/2 + 1);
    if state.max_flow >= (i+1) * prob.privacy_bound { // we need to settle in smaller
        radius = search_for_radius(true, &mut smaller, cursor, i, prob, state);
    } else { // we need to settle in bigger
        radius = search_for_radius(false, &mut bigger, cursor, i, prob, state);
    }
    
    // concatenate smaller and bigger;
    list.clear();
    list.append(&mut smaller);
    list.append(&mut bigger);
    radius
}


