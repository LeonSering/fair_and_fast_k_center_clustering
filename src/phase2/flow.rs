use super::Edge;
use crate::types::CenterIdx;
use crate::{PointCount,space::Point};
use std::collections::VecDeque;

// describes the state of the flow network and of the centers graph.
#[derive(Debug)]
pub(super) struct State<'a>{
    pub center_of: Vec<Option<CenterIdx>>, // center[x] is the assigned center (in form of gonzales index) of point x
    pub max_flow: PointCount, // the current value of a maximum flow
    reassign: Vec<Vec<VecDeque<&'a Point>>>, // reassing[c][t] contains all points that are assigned to c but could also be assigned to t
    unassigned: Vec<VecDeque<&'a Point>>, // unassigned[c] contains all points that are unassigned but could be assigned to c
    // reassign and unassigned might not be correctly updated at all times. E.g., a point in
    // reassign might not be assigned to c anymore, or not be assignable to t;
    // or a point in unassigned might be assigned to some point or not assignable to c.
    // So whenever an element is poped (or we check for is_empty()) it needs to be checked if it still makes sense.
    pub next_to_non_private: Vec<Option<CenterIdx>>, // indicates whether a center is non-private or has a path in centers_graph to a non-private center
    pub number_of_covered_points: Vec<PointCount>, // the number of points covered by center c; if this equals privacy_bound this center is "private"; if it is smaller than it is "non-private".; if the has more than that it is "overfull"
    pub edge_present: Vec<Vec<bool>>, // edge_present[c][p] denotes whether the arc from center c to point p is currently present in the flow network
//    pub target_edge: Edge, // TODO determines which edges are in the flow networks and which are not. All edges e with e < target_edge are present; all others are not.
}
impl<'a> State<'a> {
    pub fn reassign_pop(&mut self, c1 : CenterIdx, c2: CenterIdx) -> Option<&'a Point> {
        self.clean_reassign(c1,c2);
        self.reassign[c1][c2].pop_front()
    }

    pub fn reassign_push(&mut self, c1 : CenterIdx, c2: CenterIdx, p: &'a Point) {
        self.reassign[c1][c2].push_back(p);
    }

    pub fn reassign_is_empty(&mut self, c1: CenterIdx, c2: CenterIdx) -> bool {
        self.clean_reassign(c1,c2);
        self.reassign[c1][c2].is_empty()
    }

    fn clean_reassign(&mut self, c1: CenterIdx, c2: CenterIdx) { // removes invalid entries at the front
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


    pub fn unassigned_pop(&mut self, c : CenterIdx) -> Option<&'a Point> {
        self.clean_unassigned(c);
        self.unassigned[c].pop_front()
    }

    pub fn unassigned_push(&mut self, c : CenterIdx, p: &'a Point) {
        self.unassigned[c].push_back(p);
    }

    pub fn unassigned_is_empty(&mut self, c: CenterIdx) -> bool {
        self.clean_unassigned(c);
        self.unassigned[c].is_empty()
    }

    fn clean_unassigned(&mut self, c: CenterIdx) { // removes invalid entries at the front
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

pub(super) fn initialize_state<'a>(n: PointCount, k: PointCount) -> State<'a> {
    State {
        center_of: vec!(None; n), // gives for each point the index (in gonzales array) of the center it is assigned to; at the beginning all are unassigned (= None)
        max_flow: 0,
        reassign: (0..k).map(|_| (0..k).map(|_| VecDeque::new()).collect()).collect(),
        unassigned: (0..k).map(|_| VecDeque::new()).collect(),
        next_to_non_private: (0..k).map(|l| Some(l)).collect(), // in the beginning all centers are non_private
        number_of_covered_points: vec![0; k],
        edge_present: vec!(vec!(false; n); k),
    }

}

// input: the edge e for processing, the current gonzales set (0,...,i), the current bucket j, and the
//        index t of the center of e.left (within the gonzales set)
// task: add edge e to the current flow network and look for an augmenting path to increase the
// flow by 1; then execute this augmentation

pub(super) fn add_edge<'a>(e: Edge<'a>, i: CenterIdx, k: PointCount, privacy_bound:PointCount, state: &mut State<'a>){
    let c = e.left;
    let x = e.right; // maybe cloning needed

    if c > i { // this edge is not considered yet
//        println!("\tnot yet! Max flow: {}\n", state.max_flow);
        return;
    }
    // as the edge is now added to the flow networks we mark it that way:
    assert_eq!(state.edge_present[c][x.idx()], false, "Edge {:?} was present before", e);
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

            if state.next_to_non_private[center_of_x].is_none() && state.next_to_non_private[c].is_some() {
                // in this case there is also a path from center_of_x to non_private
                state.next_to_non_private[center_of_x] = Some(c);
                update_reachability(vec!(center_of_x), k, state) // update path_to_non_center in backwards bfs starting from center_of_x
            }
            // otherwise the path_to_non_private will not change at all
        }
    }
    augment_flow(k,privacy_bound,i,state); // trying to increase flow. If no augmenting path: O(k), if increase possible: O(k^2)
}

pub(super) fn remove_edge<'a>(e: Edge<'a>, i: CenterIdx, k: PointCount, privacy_bound: PointCount, state: &mut State<'a>) {
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
            rebuild_reachability(k, privacy_bound, state);
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
                rebuild_reachability(k, privacy_bound, state);
                augment_flow(k, privacy_bound, i, state); // see if we can add a flow unit along other paths
            } else {
                rebuild_reachability(k, privacy_bound, state);
//                println!{"\tarc was not flow carrying!"};
            }
        }
    }
}

fn augment_flow<'a>(k: PointCount, privacy_bound: PointCount, i: CenterIdx, state: &mut State<'a>) {
    // We need to check if there is an augmenting path to increase the max_flow:

    // looking for the first center v (in 0,...,i) that has a path to a non_private center
    // and has unassigned nodes that could be assigned to v
    let mut v = 0;
    while v <=i {
        if !state.unassigned_is_empty(v) && state.next_to_non_private[v].is_some() {
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


    while state.number_of_covered_points[v] > privacy_bound { // while v is overfull

        // need to find new center w, such that arc (v,w) is in centers_graph, at w is closer to
        // non-private center
        let w = state.next_to_non_private[v].unwrap();

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

    // Since the connectivity in the center_graph has changed we need to update to
    // path_to_non_private:
    rebuild_reachability(k,privacy_bound,state); // takes O(k^2) time
}


// rebuilds the data structures can_reach and path_to_non_private
// takes up to O(k^2) time
fn rebuild_reachability(k: PointCount, privacy_bound: PointCount, state: &mut State) {
    // we need to rebuild path_to_non_private
    // first set only non_private centers to true
    state.next_to_non_private = (0..k).map(|c| if state.number_of_covered_points[c] < privacy_bound {Some(c)} else {None}).collect();
    // then update with a backwards BFS all other centers: takes up to O(k^2) time
    update_reachability((0..k).filter(|c| state.next_to_non_private[*c].is_some()).collect(), k, state);
}

fn update_reachability(start_centers: Vec<CenterIdx>, k: PointCount, state: &mut State) {
    // do a backwards bfs
    let mut queue : VecDeque<CenterIdx> = VecDeque::with_capacity(k);
    for c in start_centers {
        queue.push_back(c);
    }
    while !queue.is_empty() {
        let w = queue.pop_front().unwrap();
        for v in 0..k {
            if !state.reassign_is_empty(v, w) && state.next_to_non_private[v].is_none() {
                state.next_to_non_private[v] = Some(w);
                queue.push_back(v);
            }
        }
    }
}


