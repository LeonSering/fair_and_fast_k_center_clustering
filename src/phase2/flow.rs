use super::Edge;
use crate::{ClusteringProblem,space::Point};
use std::collections::VecDeque;

// describes the state of the flow network and of the centers graph.
#[derive(Debug)]
pub struct State<'a>{
    pub center_of: Vec<Option<usize>>, // center[x] is the assigned center (in form of gonzales index) of point x
    reassign: Vec<Vec<VecDeque<&'a Point>>>, // reassing[c][t] contains all points that are assigned to c but could also be assigned to t
    unassigned: Vec<VecDeque<&'a Point>>, // unassigned[c] contains all points that are unassigned but could be assigned to c
    // reassign and unassigned might not be correctly updated at all times. E.g., a point in
    // reassign might not be assigned to c anymore, or not be assignable to t; 
    // or a point in unassigned might be assigned to some point or not assignable to c. 
    // So whenever an element is poped (or we check for is_empty()) it needs to be checked if it still makes sense.
    pub number_of_covered_points: Vec<usize>, // the number of points covered by center c; if this equals privacy_bound this center is "private"; if it is smaller than it is "non-private".; if the has more than that it is "overfull"
    pub can_reach: Vec<Vec<bool>>, // indicates whether there is a path in centers_graph from c_1 to c_2; 
    pub path_to_non_private: Vec<bool>, // indicates whether a center is non-private or has a path in centers_graph to a non-private center
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

pub fn initialize_state<'a>(n: usize, k: usize) -> State<'a> {
    State {
        center_of: vec!(None; n), // gives for each point the index (in gonzales array) of the center it is assigned to; at the beginning all are unassigned (= None)
        reassign: (0..k).map(|_| (0..k).map(|_| VecDeque::with_capacity(k*k)).collect()).collect(),
        unassigned: (0..k).map(|_| VecDeque::with_capacity(k*k)).collect(),
        can_reach: (0..k).map(|i| (0..k).map(|j| if i == j {true} else {false}).collect()).collect(), // in the beginning there are no arcs in centers_graph
        path_to_non_private: vec![true; k], // in the beginning all centers are non_private
        number_of_covered_points: vec![0; k],
        max_flow: 0,
        edge_present: vec!(vec!(false; n); k),
    }

}

// input: the edge e for processing, the current gonzales set (0,...,i), the current bucket j, and the
//        index t of the center of e.left (within the gonzales set)
// task: add edge e to the current flow network and look for an augmenting path to increase the
// flow by 1; then execute this augmentation

pub fn add_edge<'a>(e: Edge<'a>, i: usize, prob: &ClusteringProblem, state: &mut State<'a>){
    let c = e.left;
    let x = e.right; // maybe cloning needed

    if c > i { // this edge is not considered yet
//        println!("\tnot yet! Max flow: {}\n", state.max_flow);
        return;
    }

    // as the edge is now added to the flow networks we mark it that way:
    assert_eq!(state.edge_present[c][x.idx()], false, "Edge was present before");
    state.edge_present[c][x.idx()] = true;

    match state.center_of[x.idx()] {
        None => {// x is not assigned yet
            // a new node correspdonding to x is added to the tail of the queue unassigned:
            state.unassigned_push(c,x);
//            println!("\tcurrent point is unassigned");
        },
        Some(center_of_x) => {
            // in this case, x is already assigned to center[x], so we have to add the new edges in
            // the centers graph
            state.reassign_push(center_of_x,c,x);
//            println!("\tcurrent point is assigned to {}", center_of_x);

            // There is now a path from center_of_x to c:
            state.can_reach[center_of_x][c] = true;

            // update reachability status in centers graph
            for d in 0..(i+1) {
                state.can_reach[d][c] = state.can_reach[d][c] || state.can_reach[d][center_of_x];
                state.path_to_non_private[d] = state.path_to_non_private[d] || (state.can_reach[d][center_of_x] && state.path_to_non_private[c]);
            }
        }
    }
    augment_flow(prob,i,state);
//    println!("\tmax_flow: {}\n", state.max_flow);
}

pub fn remove_edge<'a>(e: Edge<'a>, i: usize, prob: &ClusteringProblem, state: &mut State<'a>) {
    let c = e.left;
    let x = e.right;

    if c > i { // this edge has not been considered yet, and hence is not in the network
//        println!("\n\tnot yet! Max flow: {}", state.max_flow);
        return;
    }

    // as the edge is now removed from the flow networks we mark it that way:
    assert_eq!(state.edge_present[c][x.idx()], true, "Edge was present before");
    state.edge_present[c][x.idx()] = false;

//    println!("x: {}, center of x: {:?}", x.idx(), state.center_of[x.idx()]);
    match state.center_of[x.idx()] {
        None => {
//            println!{"\tarc was not flow carrying (point not assigned)!"};

        }
        Some(center_of_x) => {
            if center_of_x == c { // arc e is flow carrying:
                //println!{"\tarc is flow carrying!"};
                state.max_flow -= 1;
                state.number_of_covered_points[c] -= 1;
                state.center_of[x.idx()] = None;

                // this unassignment indirectly also updates reassign, because whenever something is
                // poped from reassign[t][_] we check whether that x is assigned to t
                // also x cannot be part of unassigned[t] at this point, because the moment it was
                // assigned to t it was poped from unassigned[t] (and it can only enter this queue
                // again, when edge e is added to the network)
                rebuild_reachability(prob, state); // TODO: This is too expensive
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
        if !state.unassigned_is_empty(v) && state.path_to_non_private[v] {
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
    let mut x = state.unassigned_pop(v).unwrap();

    assert_eq!(state.center_of[x.idx()], None, "WAIT: p is not unassigned");

    state.center_of[x.idx()] = Some(v); // assign p to v
    state.number_of_covered_points[v] += 1; // v covers now one points more
    state.max_flow += 1;

    // Now it could be the case that v is private already (covers privacy_bound many points already), so we need to find a new center
    // for this we do a BFS on the centers_graph (but only on the centers that has a path to a non_private center)

//    // This implementation takes O(k^3) time
//    let mut hops_to_non_private: Vec<Option<usize>> = (0..i+1).map(|c| if state.number_of_covered_points[c] < prob.privacy_bound {Some(0)} else {None} ).collect();
//    for _ in 1..i+1 { // we have a maximum of i+1 hops
//        for c1 in 0..i+1 {
//            for c2 in 0..i+1 {
//
//                if !state.reassign_is_empty(c1,c2) {
//                    hops_to_non_private[c1] = match hops_to_non_private[c2] {
//                        None => hops_to_non_private[c1],
//                        Some(hops_c2) => match hops_to_non_private[c1] {
//                            None => Some(hops_c2 + 1),
//                            Some(hops_c1) => if hops_c2 + 1 < hops_c1 {Some(hops_c2 + 1)} else { Some(hops_c1)},
//                        }
//                    }
//                }
//            }
//        }
//    }


    // compute some possible amount of hops to non private (not the minimal number of hops) doing a
    // BFS. Hence it takes O(i^2) time; note that each center goes into the queue at most once.
    let mut hops_to_non_private: Vec<Option<usize>> = (0..i+1).map(|c| if state.number_of_covered_points[c] < prob.privacy_bound {Some(0)} else {None} ).collect();
    let mut center_queue: VecDeque<usize> = VecDeque::with_capacity(i);

    for c in 0..i+1 {
        if hops_to_non_private[c].is_some() {
            center_queue.push_back(c);
        }
    }
    while !center_queue.is_empty() {
        let c1 = center_queue.pop_front().unwrap();
        for c2 in 0..i+1 {
            if hops_to_non_private[c2].is_none() && !state.reassign_is_empty(c2,c1) {
                hops_to_non_private[c2] = Some(hops_to_non_private[c1].unwrap() + 1);
                center_queue.push_back(c2);
            }
        }
    }



    for c1 in 0..i+1 {
        for c2 in 0..i+1 {

            if !state.reassign_is_empty(c1,c2) {
            }
        }
    }

//    println!("hops to non_private: {:?}", hops_to_non_private);
//    println!("number_of_covered_points {:?}", state.number_of_covered_points);
//
//    println!("reassign: {:?}", state.reassign);
//    println!("can_reach: {:?}", state.can_reach);

    while state.number_of_covered_points[v] > prob.privacy_bound { // while v is overfull

        // need to find new center w, such that arc (v,w) is in centers_graph, at w is closer to
        // non-private center
        let mut w = 0;
        while w <= i {
            if hops_to_non_private[w] != None && hops_to_non_private[w].unwrap() <= hops_to_non_private[v].unwrap() - 1 && !state.reassign_is_empty(v,w) { 
                // w is closer to a non-private center and the arc vw exists
                break;
            }
            w += 1;
        }
        assert!(w <= i,"There is no center w such that (v, w) is in centers_graph and w leads to a non_private center. This contradicts the fact that B[v] == True; v = {}", v);

        // w is our next center in the augmenting path. We reassign y, which means adding
        // forward arc (w,y) and backwards arc (y, v) in front of our augmenting path
        x = state.reassign_pop(v,w).unwrap();
        assert_eq!(state.center_of[x.idx()], Some(v), "WAIT: p was not assigned to v!");
        state.number_of_covered_points[v] -= 1;
        state.center_of[x.idx()] = Some(w); // reassign y to w
        state.reassign_push(w,v,x);
        state.number_of_covered_points[w] += 1;

        v = w;

    }

    // Since the connectivity in the center_graph has changed we need to update to corresponding
    // data_structures:
    rebuild_reachability(prob,state); // takes O(k^3) time
    // this is done way too often. Is it possible to update these data_structures instead of
    // rebuilding? Maybe save the arcs that are added/removed.
}


// rebuilds the data structures can_reach and path_to_non_private
// takes O(k^3) time
fn rebuild_reachability(prob: &ClusteringProblem, state: &mut State) {
    // we need to rebuild can_reach and path_to_non_private
    // TODO: This is not done properly, it takes k^3
    state.path_to_non_private = (0..prob.k).map(|c| state.number_of_covered_points[c] < prob.privacy_bound).collect();
    state.can_reach = (0..prob.k).map(|c1| (0..prob.k).map(|c2| if c1 == c2 {true} else {false}).collect()).collect();
    for c1 in 0..prob.k {
        for c2 in 0..prob.k {
            if !state.reassign_is_empty(c1,c2) { // there is an arc from c1 to c2:
                for c3 in 0..prob.k {
                    state.can_reach[c1][c3] = state.can_reach[c1][c3] || state.can_reach[c2][c3];
                    state.path_to_non_private[c1] = state.path_to_non_private[c1] || (state.can_reach[c1][c3] && state.path_to_non_private[c3]);
                }
            }
        }
    }

}
