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
    //            direct_res_path_from_center_to_sink: vec!(VecDeque::new(); i),
    //            direct_res_path_from_center_to_z: vec!(VecDeque::new(); i),
    //            direct_res_path_from_center_to_center: vec!(vec!(VecDeque::new(); i); i),
    //            path_to_sink: vec!(false; i),
    //            direct_res_path_from_color_to_center: vec!(Vec::new(); i),
            };


            // now we add one edge after the other in increasing order:
            for &edge in edges.iter() {

                // add edge (happens implicitly now)
                state.current_largest_edge = *edge;

                // augment flow by one if possible
                augmenting_flow(&mut state, &network, prob);

                
                // check whether max_flow saturates all source leaving arcs:
                if state.flow_value == sum_of_a + i + 1 {
                    break;
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
    opening : &'a OpeningList,
    i : CenterIdx,
    number_of_points : PointCount,
    number_of_colors : usize,
    sum_of_a : usize,
    a : &'a Vec<usize>,
    b : &'a Vec<usize>,
    edges_of_cluster : &'a Vec<Vec<ColorEdge>>,
    point_to_node : &'a HashMap<PointIdx, usize>,
//    color_to_node : &'a HashMap<ColorIdx, usize>,
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

    // utilities:
//    direct_res_path_from_center_to_sink: Vec<VecDeque<PointIdx>>, // one queue (containing points) for each center j; non_empty <=> there exists a point p of color l with unsaturated edge (j,p) and (p,l) and unsaturated arc from l to sink.
//    direct_res_path_from_center_to_z: Vec<VecDeque<PointIdx>>, // one queue (containing points) for each center j; non_empty <=> there exists a color class point p of color l with unsaturated edge (j,p) and (p,l) and unsaturated path from l to t to z.
//    direct_res_path_from_center_to_center: Vec<Vec<VecDeque<CenterPath>>>, // non_empty <=> there is a res path from j1 to j2 that only goes over point-nodes and color-nodes and t.
//    path_to_sink: Vec<bool>, // true <=> there is some (possibly complicated) residual path from center j to the sink; is updated over a BFS over the centers-graph given by direct_res_path_from_center_to_center
//    direct_res_path_from_color_to_center: Vec<Vec<CenterIdx>>, // entry[l] contains a list of centers such that center j covers a point of color class l <=> there is a direct res part from the color-node l to center j
}

//#[derive(Clone)]
//enum CenterPath {
//    OverPoint(PointIdx),
//    OverColor(PointIdx, PointIdx),
//    OverT(PointIdx, PointIdx),
//}
//

struct Node {
    discovered_by: Option<NodeKind>,
    kind : NodeKind,
}

#[derive(Debug)]
enum NodeKind {
    Source, 
    Z,
    Center(CenterIdx),
    Point(PointIdx),
    Color(ColorIdx),
    T,
    Sink,
}

fn augmenting_flow(state: &mut State, net: &Network, prob: &ClusteringProblem) { 

    // create nodes:
    let mut z_node = Node{discovered_by:None, kind:NodeKind::Z};
    
    let mut center_nodes : Vec<Node> = Vec::with_capacity(net.i + 1);
    for j in 0..net.i+1 {
        center_nodes.push(Node{discovered_by:None, kind: NodeKind::Center(j)});
    }

    let mut point_nodes : Vec<Node> = Vec::with_capacity(net.number_of_points);
    for p in 0..net.number_of_points {
        point_nodes.push(Node{discovered_by:None, kind:NodeKind::Point(p)});
    }

    let mut color_nodes : Vec<Node> = Vec::with_capacity(net.number_of_colors);
    for l in 0..net.number_of_colors {
        color_nodes.push(Node{discovered_by:None, kind:NodeKind::Color(l)});
    }
    let mut t_node = Node{discovered_by:None, kind: NodeKind::T};
    let mut sink_discovered_by : Option<NodeKind> = None;



    // search for augmenting path:

    let mut queue: VecDeque<NodeKind> = VecDeque::new();

    queue.push_back(NodeKind::Source);

    while !queue.is_empty() {
        let current_node = queue.pop_front().unwrap();

        match current_node {
            NodeKind::Source => {
                if state.flow_source_z < net.sum_of_a {
                    z_node.discovered_by = Some(NodeKind::Source);
                    queue.push_back(NodeKind::Z);
                }

                for j in 0..net.i+1 {
                    if !state.flow_source_center[j] {
                        center_nodes[j].discovered_by = Some(NodeKind::Source);
                        queue.push_back(NodeKind::Center(j));
                    }
                }


            }
            NodeKind::Z => {
                if state.flow_z_sink < net.i+1 && sink_discovered_by.is_none() {
                    sink_discovered_by = Some(NodeKind::Z);
                    // Augmenting path found!
                    break;
                }

                if state.flow_t_z > 0 && t_node.discovered_by.is_none() {
                    t_node.discovered_by = Some(NodeKind::Z);
                    queue.push_back(NodeKind::T);
                }

                // To be continued with edges z to center
                for j in 0..net.i+1 {
                    if state.flow_z_center[j] < net.opening.eta[j] - 1 && center_nodes[j].discovered_by.is_none() {
                        center_nodes[j].discovered_by = Some(NodeKind::Z);
                        queue.push_back(NodeKind::Center(j));

                    }
                }

            }
            NodeKind::Center(idx) => {
                if state.flow_z_center[idx] > 0 && z_node.discovered_by.is_none() {
                    z_node.discovered_by = Some(NodeKind::Center(idx));
                    queue.push_back(NodeKind::Z);
                }
                for &edge in net.edges_of_cluster[idx].iter() {
                    let p:usize = *net.point_to_node.get(&edge.point).unwrap();
                    if edge <= state.current_largest_edge // edge present in network
                    && (state.point_covered_by[p].is_none() || state.point_covered_by[p].unwrap() != idx)  // edge present in res network
                    && point_nodes[p].discovered_by.is_none() { // point_node not visited yet

                        point_nodes[p].discovered_by = Some(NodeKind::Center(idx));
                        queue.push_back(NodeKind::Point(p));

                    }

                }


            }

            NodeKind::Point(idx) => {
                let l: usize = net.point_idx_to_color_idx[idx];
                if state.point_covered_by[idx].is_none() && color_nodes[l].discovered_by.is_none() {
                    color_nodes[l].discovered_by = Some(NodeKind::Point(idx));
                    queue.push_back(NodeKind::Color(l));
                }
                for &edge in net.edges_by_point_node[idx].iter() {
                    if edge > &state.current_largest_edge {
                        break; // edges_by_point_node are sorted; we go from small to big
                    }
                    let j = edge.center;
                    if state.point_covered_by[idx] == Some(j) && center_nodes[j].discovered_by.is_none() {
                        center_nodes[j].discovered_by = Some(NodeKind::Point(idx));
                        queue.push_back(NodeKind::Center(j));
                    }
                     

                }

            }

            NodeKind::Color(idx) => {
                if state.flow_color_sink[idx] < net.a[idx] && sink_discovered_by.is_none() {
                    // Augmenting path found!
                    sink_discovered_by = Some(NodeKind::Color(idx));
                    break;
                }

                if state.flow_color_t[idx] < net.b[idx] - net.a[idx] && t_node.discovered_by.is_none() {
                    t_node.discovered_by = Some(NodeKind::Color(idx));
                    queue.push_back(NodeKind::T);
                }

                for &edge in net.edges_by_color_node[idx].iter() {
                    if edge > &state.current_largest_edge {
                        break;
                    }
                    let p = *net.point_to_node.get(&edge.point).unwrap();

                    if state.point_covered_by[idx] == Some(edge.center) && point_nodes[p].discovered_by.is_none() {
                        point_nodes[p].discovered_by = Some(NodeKind::Color(idx));
                        queue.push_back(NodeKind::Point(p));
                    }
                }
            }

            NodeKind::T => {
                if state.flow_t_z < prob.k - net.sum_of_a && z_node.discovered_by.is_none() {
                    z_node.discovered_by = Some(NodeKind::T);
                    queue.push_back(NodeKind::Z);
                }

                for l in 0 .. net.number_of_colors {
                    if state.flow_color_t[l] > 0 && color_nodes[l].discovered_by.is_none() {
                        color_nodes[l].discovered_by = Some(NodeKind::T);
                        queue.push_back(NodeKind::Color(l));
                    }
                }

            }
            NodeKind::Sink => {
                panic!("The sink should never enter the queue.");
            }

        }

    }

    if sink_discovered_by.is_none() {
        // No augmenting path found!
        return;
    }

    // Augmenting path found!
    state.flow_value += 1;

    // Now augment flow along this path:
    let mut current_node = &NodeKind::Sink; 
    let mut last_node = &NodeKind::Sink;
    loop {
//        println!("{:?}", current_node);
        match *current_node {
            NodeKind::Source => {
                match *last_node {
                    NodeKind::Z => {
                        state.flow_source_z += 1;
                    }
                    NodeKind::Center(j) => {
                        state.flow_source_center[j] = true;
                    }
                    _ => {
                        panic!("This should never happen");
                    }
                }
                break;
            }
            NodeKind::Z => {
                match *last_node {
                    NodeKind::Sink => {
                        state.flow_z_sink += 1;
                    }
                    NodeKind::T => {
                        state.flow_t_z -= 1;
                    }
                    NodeKind::Center(j) => {
                        state.flow_z_center[j] += 1;
                    }
                    _ => {
                        panic!("This should never happen");
                    }
                }
                last_node = current_node;
                current_node = z_node.discovered_by.as_ref().unwrap();
            }
            NodeKind::Center(idx) => {
                match *last_node {
                    NodeKind::Z => {
                        state.flow_z_center[idx] -= 1;
                    }
                    NodeKind::Point(p) => {
                        state.point_covered_by[p] = Some(idx);
                    }
                    _ => {
                        panic!("This should never happen");
                    }
                }
                last_node = current_node;
                current_node = center_nodes[idx].discovered_by.as_ref().unwrap();
            }
            NodeKind::Point(idx) => {
                match *last_node {
                    NodeKind::Center(_) => {
                        state.point_covered_by[idx] = None;
                    }
                    NodeKind::Color(_) => {}
                    _ => {
                        panic!("This should never happen");
                    }
                }
                last_node = current_node;
                current_node = point_nodes[idx].discovered_by.as_ref().unwrap();
            }
            NodeKind::Color(idx) => {
                match *last_node {
                    NodeKind::T => {
                        state.flow_color_t[idx] += 1;
                    }
                    NodeKind::Sink => {
                        state.flow_color_sink[idx] += 1;
                    }
                    NodeKind::Point(_) => {}
                    _ => {
                        panic!("This should never happen");
                    }
                }
                last_node = current_node;
                current_node = color_nodes[idx].discovered_by.as_ref().unwrap();
            }
            NodeKind::T => {
                match *last_node {
                    NodeKind::Color(l) => {
                        state.flow_color_t[l] -= 1;
                    }
                    NodeKind::Z => {
                        state.flow_t_z += 1;
                    }
                    _ => {
                        panic!("This should never happen");
                    }
                }
                last_node = current_node;
                current_node = t_node.discovered_by.as_ref().unwrap();
            }
            NodeKind::Sink => {
                last_node = current_node;
                current_node = sink_discovered_by.as_ref().unwrap();
            }
        }
    }
}


// ideas for quicker flow augmentation:
//fn add_edge j, p, l
// 0) current_largest_edge = edge
// 1) if point_coverd_by[p] == Some(c) => directed_res_path_from_center_to_center[j][c] = true
// 2) Else (point_covered_by[p] == None:
// 3) for all center c with direct_res_path_from_color_to_center[l][c] == true: add arc from j to c 
//
// fn add edge (s_j,l):
// Check whether there is a res path from source to s_j (backwards):
// 1) check whether flow_source_center[j] = false => Path found
// 2) if not, check whether flow_z_center[j] is saturated yet => No Path
// 3) if not saturated, check whether flow_source_z is not saturated yet => Path found
// 4) if saturated, check whether flow_t_z is saturated yet => No Path
// 5) if not saturated, check whether there is a non saturated flow_color_t with a color with
//    direct_res_path_from_color_to_soure = true; => Path Found
// 6) If no such color class exists: No Path
//
// Check whether there is a res path from l to sink
// 1) check flow_color_to_sink[l] is not saturated => Path Found
// 2) if saturated, check whether flow_color_to_t[l] is saturated => No Path
// 3) if not saturated, check whether flow_t_z is saturated => No Path
// 4) if not saturated, check whether flow_z_sink is not saturated => Path found
// 5) if saturated, check whether there is some center with a non saturated flow_z_center and
//    direct_res_path_from_center_to_sink = true => Path Found
// 6) If no such center exists: No path
//
// If both paths exists: augment flow 
