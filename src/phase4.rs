use crate::{ClusteringProblem,Centers,Clustering,ColoredMetric,PointCount,space::{PointIdx,ColorIdx,Distance}};
use crate::clustering::CenterIdx;
use super::OpeningList;
use crate::utilities;
use std::collections::HashMap;

/// an edge between a gonzales center and a color class;
/// lablled with the point and the distance between center and point
#[derive(Debug,Clone,Copy,PartialOrd,PartialEq)]
struct ColorEdge {
    d: Distance, // distance between the gonzales center and the point
    center: CenterIdx, // index to the gonzales center
    point: PointIdx, // a point that might be a final center
    color: ColorIdx, // the color of the point
}

/// Given a metric space and a clustering problem,
/// finalize takes a vector of clusterings in which each center (except for one) covers a multple of L
/// points, and returns a vector of clusterings (which potentially more cluster as before) that also satisfy the representative constaints. 
///
pub(crate) fn finalize<'a, M : ColoredMetric>(space : &'a M, prob : &ClusteringProblem, opening_lists : Vec<Vec<OpeningList>>,  gonzales : &Centers) -> Vec<Clustering<'a>> {



    //TEMP: For test reasons we also look at flow problems with ALL edges present:
    let mut all_edges_of_cluster: Vec<Vec<ColorEdge>> = vec!(Vec::with_capacity(space.n());prob.k);
    for i in 0..prob.k {
        for p in space.point_iter() {
            let color = space.color(p);
            all_edges_of_cluster[i].push(ColorEdge{
                d: space.dist(gonzales.get(i), p),
                center: i,
                point: p.idx(),
                color});
            
        }
        all_edges_of_cluster[i].sort_by(|a,b| a.partial_cmp(b).unwrap());
    }
    // TEMP END
    
    /////////////////////////////////////////////////////////////////    
    //////// define the neighborhood of each gonzales center: ///////
    /////////////////////////////////////////////////////////////////    
    let mut edges_of_cluster: Vec<Vec<ColorEdge>> = vec!(Vec::with_capacity(prob.k);prob.k);


    for i in 0..prob.k {

        // first we make edge lists for each color class
        let mut edges_by_color: Vec<Vec<ColorEdge>> = vec!(Vec::new();space.gamma());

        // for this we create all edges from gonzales center i to all points
        for p in space.point_iter() {
            let color = space.color(p);
            edges_by_color[space.color(p)].push(ColorEdge{
                d: space.dist(gonzales.get(i), p),
                center: i,
                point: p.idx(),
                color});
            
        }


        let restricted_colors = if space.gamma() < prob.rep_interval.len() {space.gamma()} else {prob.rep_interval.len()}; // min{gamma, rep_interval.len()}
        let mut remaining_edges: Vec<ColorEdge> = Vec::new();
        let mut num_edges_to_fill = prob.k;


        for c in 0..restricted_colors {
            // take only the b smallest in each class; all others cannot play any role.
            utilities::truncate_to_smallest(&mut edges_by_color[c], prob.rep_interval[c].1);
            // the a smallest have to be present for sure, so each center can satisfy this
            // condition by itself;
            // the remaining b-a edges are collected in remaining_edges
            remaining_edges.append(&mut utilities::split_off_at(&mut edges_by_color[c], prob.rep_interval[c].0));
            num_edges_to_fill -= prob.rep_interval[c].0;
        }

        // for all colors without restriction fill them into the remaining edges:
        for c in restricted_colors..space.gamma() {
            remaining_edges.append(&mut edges_by_color[c]);
        }
        
        // so far sum(a) edges where chosen, so we fill up with the nearest k - sum(a) edges,
        // independent of the color 
        utilities::truncate_to_smallest(&mut remaining_edges,num_edges_to_fill);



        // put all edges into one list:
        for c in 0..restricted_colors {
            edges_of_cluster[i].append(&mut edges_by_color[c]);
        }
        edges_of_cluster[i].append(&mut remaining_edges);

        // there are now max k edges in edges_of_cluster[i]; so we can sort them:
        edges_of_cluster[i].sort_by(|a,b| a.partial_cmp(b).unwrap());
        println!("edges of center i: {:?}", edges_of_cluster[i]);

    }






    /////////////////////////////////////////////////////////////////    
    ///////////////////////// solve flow problem: ///////////////////
    /////////////////////////////////////////////////////////////////    


    let sum_of_a: PointCount = prob.rep_interval.iter().map(|interval| interval.0).sum();
    
    for i in 0..prob.k {

        println!("\n\n************** i = {} ******************", i);

        // consider all edges starting from centers 0 to i in increasing order:
    
        let mut edges: Vec<&ColorEdge> = Vec::with_capacity((i+1) * prob.k); // a reference to all edges with centers in S_i
        for j in 0..i+1 {
            edges.extend(edges_of_cluster[j].iter());
        }

        edges.sort_by(|a,b| a.partial_cmp(b).unwrap());

        


        // not that we only have k^2 edges so at most k^2 point and color classes can occur
        // We reindex the relevant points and colors. The new indices are called point-node and color-node index and they are
        // used in the flow network
       
        let mut point_to_node: HashMap<PointIdx, usize> = HashMap::with_capacity((i + 1) * prob.k); // maps a point index to the point-node index used in the flow network
        let mut point_counter = 0; // number of relevant points (= number of point-nodes)
        
        let mut a: Vec<PointCount> = Vec::with_capacity((i + 1) * prob.k); // lower bound given by the color-node index
        let mut b: Vec<PointCount> = Vec::with_capacity((i +1) * prob.k); // upper bound given by the color-node index
        let mut color_to_node: HashMap<ColorIdx, usize> = HashMap::with_capacity((i + 1) * prob.k); // maps a color index to the color-node index used in the flow network
        let mut color_counter = 0; // number of relevant colors (= number of color-nodes)


        let mut point_idx_to_color_idx: Vec<ColorIdx> = Vec::with_capacity((i+1) * prob.k);
        for e in edges.iter() {
            if !color_to_node.contains_key(&e.color) {
                color_to_node.insert(e.color, color_counter);
                if e.color >= prob.rep_interval.len() {
                    a.push(0);
                    b.push(<PointCount>::MAX);
                } else {
                    a.push(prob.rep_interval[e.color].0);
                    b.push(prob.rep_interval[e.color].1);
                }
                color_counter += 1;
            }
            if !point_to_node.contains_key(&e.point) {
                point_to_node.insert(e.point, point_counter);
                point_idx_to_color_idx.push(*color_to_node.get(&e.color).unwrap());
                point_counter += 1;
            }
        }



        // edges should also be referrable from their points and their color classes:
        
        let mut edges_by_color_node: Vec<Vec<&ColorEdge>> = vec!(Vec::new(); color_counter);
        let mut edges_by_point_node: Vec<Vec<&ColorEdge>> = vec!(Vec::new(); point_counter);

        for e in edges.iter() {
            edges_by_color_node[color_to_node[&e.color]].push(e);
            edges_by_point_node[point_to_node[&e.point]].push(e);
        }
        edges_by_color_node.sort_by(|a,b| a.partial_cmp(b).unwrap());
        edges_by_point_node.sort_by(|a,b| a.partial_cmp(b).unwrap());

       


        // for each eta vector we solve the flow problem individually:

        for opening in opening_lists[i].iter() {

            let network = Network {
                k : prob.k,
                opening,
                i,
                number_of_points: point_counter,
                number_of_colors: color_counter,
                sum_of_a,
                a : &a,
                b : &b,
                edges_of_cluster : &edges_of_cluster,
                point_to_node : &point_to_node,
//                color_to_node : &color_to_node,
                point_idx_to_color_idx : &point_idx_to_color_idx,
                edges_by_color_node : &edges_by_color_node,
                edges_by_point_node : &edges_by_point_node,
            };

            // without any edges present the maximal flow has only flow an the path (source, z, sink)
            // of value min{sum_of_a, i + 1} (note that we have i+1 centers)
            let initial_flow_value = if sum_of_a < i + 1 {sum_of_a} else {i + 1}; // = min {sum_of_a, i}
            let mut state = State{
                current_largest_edge: ColorEdge{ d: <Distance>::MIN, center: 0, point: 0, color: 0},
                flow_value: initial_flow_value, // max flow value
                flow_source_center: vec!(false; i+1),
                flow_source_z: initial_flow_value, 
                flow_z_center: vec!(0; i+1),
                point_covered_by: vec!(None;point_counter),
                flow_color_t: vec!(0;color_counter),
                flow_t_z: 0,
                flow_color_sink: vec!(0;color_counter),
                flow_z_sink: initial_flow_value,
            };
            
            // start with an empty bfs_tree
            let mut bfs_tree = BFSTree::new(i,network.number_of_points, network.number_of_colors);
            // then start a bfs from the source
            bfs(&state, &network, &mut bfs_tree, Node::Source);

            // now we add one edge after the other in increasing order:
            for &edge in edges.iter() {

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
                    augmenting_flow(&mut state, &bfs_tree);

                    // max flow has increased by 1

                    // check whether max_flow saturates all source-leaving arcs (i.e., flow value
                    // is sum over all a PLUS i + 1:
                    if state.flow_value == sum_of_a + i + 1 {
                        break;
                    }

                    // Now the bfs tree is outdated so it needs to be resetted
                    bfs_tree = BFSTree::new(i,network.number_of_points, network.number_of_colors);
                    bfs(&state, &network, &mut bfs_tree, Node::Source);
                }

                
//                println!("\n\tState: {:?}", state);
            }

            println!("\n  eta: {:?}, a: {:?}, b: {:?}", network.opening.eta, network.a, network.b);
            println!("  Final State: {:?}", state);


            assert_eq!(state.flow_value, sum_of_a + i + 1, "All edges added but still not all source leaving arcs saturated.");


            // TODO: Create centers from max flow
        }

    }





    let mut clusterings: Vec<Clustering> = Vec::with_capacity(prob.k);
    let mut centers = Centers::with_capacity(prob.k);
    centers.push(space.point_iter().next().unwrap());
    let center_of: Vec<Option<CenterIdx>> = space.point_iter().map(|_| None).collect();
    clusterings.push(Clustering::new(centers,center_of, space));
    clusterings
}

use std::collections::VecDeque;

struct Network<'a> {
    k : CenterIdx,
    opening : &'a OpeningList,
    i : CenterIdx,
    number_of_points : PointCount,
    number_of_colors : usize,
    sum_of_a : usize,
    a : &'a Vec<usize>,
    b : &'a Vec<usize>,
    edges_of_cluster : &'a Vec<Vec<ColorEdge>>,
    point_to_node : &'a HashMap<PointIdx, usize>,
    point_idx_to_color_idx : &'a Vec<usize>,
    edges_by_color_node : &'a Vec<Vec<&'a ColorEdge>>,
    edges_by_point_node : &'a Vec<Vec<&'a ColorEdge>>,
}

#[derive(Debug)]
struct State {
    // network:
    current_largest_edge: ColorEdge, // is used to indicate which edges are present in the network

    // flow:
    flow_value: usize,
    flow_source_center: Vec<bool>, // true <=> the arc from source to center is saturated (flow = 1)
    flow_source_z: usize, // has capacity sum_of_a
    flow_z_center: Vec<usize>, // one entry for each center j (0 .. i); these arcs have capacity eta_j - 1
    point_covered_by: Vec<Option<CenterIdx>>, // edge_flow_carrying[p] = Some(j) <=> point with point-node index p is covered by center j <=> edges from j to p and arc from p to l (color of p) is saturated (flow = 1) 
    flow_color_t: Vec<usize>, // one entry for each color class l (color-node index); each has a capacity of a_l - b_l
    flow_t_z: usize, // has capacity k- sum_of_a
    flow_color_sink: Vec<usize>, // one entry for each color class l; each has a capacity of a_l
    flow_z_sink: usize, // has capacity i

}


/// For each node the partent node is stored, i.e. the nodes that discovered it during the BFS
/// starting from the source;
/// the value is None if it has not been discovered yet.
struct BFSTree {
    z_node : Option<Node>,
    center_nodes : Vec<Option<Node>>,
    point_nodes: Vec<Option<Node>>,
    color_nodes: Vec<Option<Node>>,
    t_node: Option<Node>,
    sink_node: Option<Node>,
}

impl BFSTree {
    fn new(i : CenterIdx, number_of_points : usize, number_of_colors : usize) -> BFSTree {
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

                    let p:usize = *net.point_to_node.get(&edge.point).unwrap();
                    if tree.point_nodes[p].is_none() && (state.point_covered_by[p].is_none() || state.point_covered_by[p].unwrap() != idx) {  
                        // point_node has not been visited yet and edge present in res network
                        tree.point_nodes[p] = Some(Node::Center(idx));
                        queue.push_back(Node::Point(p));

                    }

                }


            }

            Node::Point(idx) => {
                let l: usize = net.point_idx_to_color_idx[idx];
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

fn augmenting_flow(state: &mut State, tree: &BFSTree) { 
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


