use crate::types::{ColorCount,PointCount,CenterIdx,PointIdx,ColorIdx, Distance};
use std::collections::VecDeque;

use super::{ColorEdge,Network,PNodeIdx,CNodeIdx,ShiftedCenters};

type FlowValue = usize;

/// State represents a flow in the network. Throughout the program we keep it to be a max-flow
#[derive(Debug)]
struct State {
    // network:
    current_largest_edge: ColorEdge, // is used to indicate which edges are present in the network

    // flow:
    flow_value: FlowValue,
    flow_source_center: Vec<bool>, // true <=> the arc from source to center is saturated (flow = 1)
    flow_source_z: FlowValue, // has capacity sum_of_a
    flow_z_center: Vec<FlowValue>, // one entry for each center j (0 .. i); these arcs have capacity eta_j - 1
    point_covered_by: Vec<Option<CenterIdx>>, // edge_flow_carrying[p] = Some(j) <=> point with point-node index p is covered by center j <=> edges from j to p and arc from p to l (color of p) is saturated (flow = 1) 
    flow_color_t: Vec<FlowValue>, // one entry for each color class l (color-node index); each has a capacity of a_l - b_l
    flow_t_z: FlowValue, // has capacity k- sum_of_a
    flow_color_sink: Vec<FlowValue>, // one entry for each color class l; each has a capacity of a_l
    flow_z_sink: FlowValue, // has capacity i

}

impl State {
    /// crates a new flow of maximal flow value in the network without any color_edges
    fn new(net : &Network) -> State{
        // without any edges present the maximal flow has only flow an the path (source, z, sink)
        // of value min{sum_of_a, i + 1} (note that we have i+1 centers)
        let initial_flow_value = if net.sum_of_a < net.i + 1 {net.sum_of_a} else {net.i + 1}; // = min {sum_of_a, i}
        State{
            current_largest_edge: ColorEdge{ d: <Distance>::MIN, center: 0, point: 0, color: 0},
            flow_value: initial_flow_value, // max flow value
            flow_source_center: vec!(false; net.i+1),
            flow_source_z: initial_flow_value, 
            flow_z_center: vec!(0; net.i+1),
            point_covered_by: vec!(None;net.number_of_points),
            flow_color_t: vec!(0;net.number_of_colors),
            flow_t_z: 0,
            flow_color_sink: vec!(0;net.number_of_colors),
            flow_z_sink: initial_flow_value,
        }



    }



}


/// For each node (except for the source) the partent node is stored, i.e. the nodes that discovered it during the BFS
/// starting from the source;
/// The value is None if it has not been discovered yet.
/// In other words: Some means that the node can be reached from source in the residual network
/// None means there is no path from the source to this node in the residual network.
struct BFSTree {
    z_node : Option<Node>,
    center_nodes : Vec<Option<Node>>,
    point_nodes: Vec<Option<Node>>,
    color_nodes: Vec<Option<Node>>,
    t_node: Option<Node>,
    sink_node: Option<Node>,
}

impl BFSTree {
    /// creates an empty BFSTree, i.e., no node is discovered yet, so everything is None.
    fn new(i : CenterIdx, number_of_points : PointCount, number_of_colors : ColorCount) -> BFSTree {
        BFSTree {
            z_node : None,
            center_nodes : vec!(None; i + 1),
            point_nodes : vec!(None; number_of_points),
            color_nodes : vec!(None; number_of_colors),
            t_node : None, 
            sink_node : None,
        }

    }
}


#[derive(Debug,Clone)]
enum Node {
    Source, 
    Z,
    Center(CenterIdx),
    Point(PointIdx),
    Color(ColorIdx),
    T,
    Sink,
}

/// Given a phase4-network it add ColorEdges one by one and computes the max flow,
/// until all source-leaving arcs are saturated, which corresponds to a feasible flow in the
/// flow-network with lower flow bounds.
///
/// Returns the partial assignment of points to gonzales center, such that the number of assigned points
/// of color class l are within [a_l,b_l]; each cluster j has between 1 and eta_j assigned points
/// and the maximal distance between an assigned point and its cluster center (= gonzales center)
/// is minimized.
pub(super) fn compute_assignment_by_flow<'a>(network : &Network) -> ShiftedCenters {
    // create a new flow of maximal flow value of min{sum_of_a, i+1} in the network without
    // any color_edges
    let mut state = State::new(network);
    
    // start with an empty bfs_tree
    let mut bfs_tree = BFSTree::new(network.i,network.number_of_points, network.number_of_colors);
    // then start a bfs from the source
    bfs(&state, &network, &mut bfs_tree, Node::Source);

    // now we add one edge after the other in increasing order:
    for &edge in network.edges.iter() {

        // add edge (happens implicitly now)
        state.current_largest_edge = *edge;
        
        let p = network.point_to_node[&edge.point]; // node idx of the point of the edge
        let j = edge.center; // center idx of edge

        if bfs_tree.center_nodes[j].is_some() && bfs_tree.point_nodes[p].is_none() {
            // if the center is can be reached from the source but not the point-node, then we
            // can extend the BFS along the new edge:
            bfs_tree.point_nodes[p] = Some(Node::Center(j));
            bfs(&state, &network, &mut bfs_tree, Node::Point(p));
        }


        if bfs_tree.sink_node.is_some() {
            // There is an augmenting path!
            augment_flow(&mut state, &bfs_tree);

            // max flow has increased by 1

            // check whether max_flow saturates all source-leaving arcs (i.e., flow value
            // is sum over all a PLUS i + 1:
            if state.flow_value == network.sum_of_a + network.i + 1 {
                break;
            }

            // Now the bfs tree is outdated so it needs to be resetted
            bfs_tree = BFSTree::new(network.i,network.number_of_points, network.number_of_colors);
            bfs(&state, &network, &mut bfs_tree, Node::Source);
        }

        
//                println!("\n\tState: {:?}", state);
    }

    println!("\n  eta: {:?}, a: {:?}, b: {:?}", network.opening.eta, network.a, network.b);
    println!("  Final State: {:?}", state);
//    println!("  node index of point of last largest edge: {}", network.point_to_node[&state.current_largest_edge.point]); 


    assert_eq!(state.flow_value, network.sum_of_a + network.i + 1, "All edges added but still not all source leaving arcs saturated.");
    let mut new_centers: Vec<PNodeIdx> = Vec::with_capacity(state.flow_value - network.i + 1 + state.flow_t_z);
    let mut origins: Vec<CenterIdx> = Vec::with_capacity(state.flow_value - network.i + 1 + state.flow_t_z);
    for (p, origin) in state.point_covered_by.iter().enumerate().filter(|(_,o)| o.is_some()) {
        new_centers.push(p);
        origins.push(origin.unwrap());
    }
//    println!("\npoint_covered_by: {:?}", state.point_covered_by);
//    println!("new_centers: {:?}", new_centers);
//    println!("origins: {:?}\n", origins);

    ShiftedCenters{
        shift_radius : state.current_largest_edge.d,
        new_centers,
        origins,
    }



}

/// Continues a BFS from the starting node. The BFSTree (tree) is only expanded.
/// For starting a new BFS start with Node::Source and an new (empty) BFSTree.
fn bfs(state: &State, net: &Network, tree: &mut BFSTree, starting_node: Node) {

    // search for augmenting path:

    let mut queue: VecDeque<Node> = VecDeque::new();

    queue.push_back(starting_node);

    while !queue.is_empty() {
        let current_node = queue.pop_front().unwrap();

        match current_node {
            Node::Source => {
                if state.flow_source_z < net.sum_of_a {
                    tree.z_node = Some(Node::Source);
                    queue.push_back(Node::Z);
                }

                for j in 0..net.i+1 {
                    if !state.flow_source_center[j] {
                        tree.center_nodes[j] = Some(Node::Source);
                        queue.push_back(Node::Center(j));
                    }
                }


            }
            Node::Z => {
                if state.flow_z_sink < net.i+1 && tree.sink_node.is_none() {
                    tree.sink_node = Some(Node::Z);
                    // Augmenting path found!
                    break;
                }

                if state.flow_t_z > 0 && tree.t_node.is_none() {
                    tree.t_node = Some(Node::Z);
                    queue.push_back(Node::T);
                }

                // To be continued with edges z to center
                for j in 0..net.i+1 {
                    if state.flow_z_center[j] < net.opening.eta[j] - 1 && tree.center_nodes[j].is_none() {
                        tree.center_nodes[j] = Some(Node::Z);
                        queue.push_back(Node::Center(j));

                    }
                }

            }
            Node::Center(idx) => {
                if state.flow_z_center[idx] > 0 && tree.z_node.is_none() {
                    tree.z_node = Some(Node::Center(idx));
                    queue.push_back(Node::Z);
                }
                for &edge in net.edges_of_cluster[idx].iter() {
                    // edges are increasing, so only go until the current largest edge
                    if edge > state.current_largest_edge {
                        break;
                    }

                    let p:PNodeIdx = *net.point_to_node.get(&edge.point).unwrap();
                    if tree.point_nodes[p].is_none() && (state.point_covered_by[p].is_none() || state.point_covered_by[p].unwrap() != idx) {  
                        // point_node has not been visited yet and edge present in res network
                        tree.point_nodes[p] = Some(Node::Center(idx));
                        queue.push_back(Node::Point(p));

                    }

                }


            }

            Node::Point(idx) => {
                let l: CNodeIdx = net.point_idx_to_color_idx[idx];
                if state.point_covered_by[idx].is_none() && tree.color_nodes[l].is_none() {
                    tree.color_nodes[l] = Some(Node::Point(idx));
                    queue.push_back(Node::Color(l));
                }
                for &edge in net.edges_by_point_node[idx].iter() {
                    if edge > &state.current_largest_edge {
                        break; // edges_by_point_node are sorted; we go from small to big
                    }
                    let j = edge.center;
                    if state.point_covered_by[idx] == Some(j) && tree.center_nodes[j].is_none() {
                        tree.center_nodes[j] = Some(Node::Point(idx));
                        queue.push_back(Node::Center(j));
                    }
                     

                }

            }

            Node::Color(idx) => {
                if state.flow_color_sink[idx] < net.a[idx] && tree.sink_node.is_none() {
                    // Augmenting path found!
                    tree.sink_node = Some(Node::Color(idx));
                    break;
                }

                if state.flow_color_t[idx] < net.b[idx] - net.a[idx] && tree.t_node.is_none() {
                    tree.t_node = Some(Node::Color(idx));
                    queue.push_back(Node::T);
                }

                for &edge in net.edges_by_color_node[idx].iter() {
                    if edge > &state.current_largest_edge {
                        break;
                    }
                    let p = *net.point_to_node.get(&edge.point).unwrap();

                    if state.point_covered_by[idx] == Some(edge.center) && tree.point_nodes[p].is_none() {
                        tree.point_nodes[p] = Some(Node::Color(idx));
                        queue.push_back(Node::Point(p));
                    }
                }
            }

            Node::T => {
                if state.flow_t_z < net.k - net.sum_of_a && tree.z_node.is_none() {
                    tree.z_node = Some(Node::T);
                    queue.push_back(Node::Z);
                }

                for l in 0 .. net.number_of_colors {
                    if state.flow_color_t[l] > 0 && tree.color_nodes[l].is_none() {
                        tree.color_nodes[l] = Some(Node::T);
                        queue.push_back(Node::Color(l));
                    }
                }

            }
            Node::Sink => {
                panic!("The sink should never enter the queue.");
            }

        }

    }

}

/// Augments along a path in the residual network.
/// Panics if there is no augmenting path from Source to Sink, i.e.,
/// if tree.sink_node == None
/// Care: If backwards path starting from the sink leads to a cycles, this function
/// end in a endless loop.
fn augment_flow(state: &mut State, tree: &BFSTree) { 
    assert!(tree.sink_node.is_some(), "Cannot augment flow if there is no augmenting path!");
    
    // Augmenting path found!
    state.flow_value += 1;

    // Now augment flow along this path:
    let mut current_node = &Node::Sink; 
    let mut last_node = &Node::Sink;
    loop {
//        println!("{:?}", current_node);
        match *current_node {
            Node::Source => {
                match *last_node {
                    Node::Z => {
                        state.flow_source_z += 1;
                    }
                    Node::Center(j) => {
                        state.flow_source_center[j] = true;
                    }
                    _ => {
                        panic!("This should never happen");
                    }
                }
                break;
            }
            Node::Z => {
                match *last_node {
                    Node::Sink => {
                        state.flow_z_sink += 1;
                    }
                    Node::T => {
                        state.flow_t_z -= 1;
                    }
                    Node::Center(j) => {
                        state.flow_z_center[j] += 1;
                    }
                    _ => {
                        panic!("This should never happen");
                    }
                }
                last_node = current_node;
                current_node = tree.z_node.as_ref().unwrap();
            }
            Node::Center(idx) => {
                match *last_node {
                    Node::Z => {
                        state.flow_z_center[idx] -= 1;
                    }
                    Node::Point(p) => {
                        state.point_covered_by[p] = Some(idx);
                    }
                    _ => {
                        panic!("This should never happen");
                    }
                }
                last_node = current_node;
                current_node = tree.center_nodes[idx].as_ref().unwrap();
            }
            Node::Point(idx) => {
                match *last_node {
                    Node::Center(_) => {
                        state.point_covered_by[idx] = None;
                    }
                    Node::Color(_) => {}
                    _ => {
                        panic!("This should never happen");
                    }
                }
                last_node = current_node;
                current_node = tree.point_nodes[idx].as_ref().unwrap();
            }
            Node::Color(idx) => {
                match *last_node {
                    Node::T => {
                        state.flow_color_t[idx] += 1;
                    }
                    Node::Sink => {
                        state.flow_color_sink[idx] += 1;
                    }
                    Node::Point(_) => {}
                    _ => {
                        panic!("This should never happen");
                    }
                }
                last_node = current_node;
                current_node = tree.color_nodes[idx].as_ref().unwrap();
            }
            Node::T => {
                match *last_node {
                    Node::Color(l) => {
                        state.flow_color_t[l] -= 1;
                    }
                    Node::Z => {
                        state.flow_t_z += 1;
                    }
                    _ => {
                        panic!("This should never happen");
                    }
                }
                last_node = current_node;
                current_node = tree.t_node.as_ref().unwrap();
            }
            Node::Sink => {
                last_node = current_node;
                current_node = tree.sink_node.as_ref().unwrap();
            }
        }
    }
}
