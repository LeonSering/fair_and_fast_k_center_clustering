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
    
    
    // define the neighborhood of each gonzales center:
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

//        println!("restricted_colors: {}", restricted_colors);

        for c in 0..restricted_colors {
            // take only the b smallest in each class; all others cannot play any role.
            utilities::truncate_to_smallest(&mut edges_by_color[c], prob.rep_interval[c].1);
            // the a smallest have to be present for sure, so each center can satisfy this
            // condition by itself;
            // the remaining b-a edges are collected in remaining_edges
//            println!("c: {}, edges_by_color {:?}, a: {}", c, edges_by_color[c], prob.rep_interval[c].0);
            remaining_edges.append(&mut utilities::split_off_at(&mut edges_by_color[c], prob.rep_interval[c].0));
            num_edges_to_fill -= prob.rep_interval[c].0;
        }

        // for all colors without restriction fill them into the remaining edges:
        for c in restricted_colors..space.gamma() {
            remaining_edges.append(&mut edges_by_color[c]);
        }
        
//        println!("\n\n************ i = {}; center = {:?} ****************", i, gonzales.get(i));
//        for (c, cclass) in edges_by_color.iter().enumerate() {
//            println!("\nColorclass {}: {:?}", c, cclass);
//        }
//        println!("\nRemaining: {:?}", remaining_edges);

        
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

    }



    let sum_of_a: PointCount = prob.rep_interval.iter().map(|interval| interval.0).sum();


    println!("edges of center 0: {:?}", edges_of_cluster[0]);

    for i in 0..prob.k {
    
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
        for &e in edges.iter() {
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
        

//        println!("\n i = {}", i);
//        println!("number of edges: {:?}", edges.len());
//        println!("color_to_node: {:?}", color_to_node);
//        println!("a: {:?}", a);
//        println!("b: {:?}", b);




//        println!("\nedges_of_cluser[{}]: {:?}", i, edges_of_cluster[i]);
//        println!("** i = {}; openings = {:?}", i, openings);
        


        // without any edges present the maximal flow has only flow an the path (source, z, sink)
        // of value min{sum_of_a, i + 1} (note that we have i+1 centers)
        let initial_flow_value = if sum_of_a < i + 1 {sum_of_a} else {i + 1}; // = min {sum_of_a, i}
        let mut state = State{
            current_largest_edge: ColorEdge{ d: <Distance>::MIN, center: 0, point: 0, color: 0},
            flow_value: initial_flow_value,
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

        for edge in edges {
            // add edge:
            state.current_largest_edge = *edge;

            // augment flow by 1 if possible
            find_augmenting_path(&mut state, point_counter, color_counter, i, sum_of_a, &prob);
            
            // check whether max_flow saturates all source leaving arcs:
            if state.flow_value == sum_of_a + i + 1 {
                break;
            }
        }

        assert_eq!(state.flow_value, sum_of_a + i + 1, "All edges added but still not all source leaving arcs saturated.");


        // TODO: Create centers from max flow

    }





    let mut clusterings: Vec<Clustering> = Vec::with_capacity(prob.k);
    let mut centers = Centers::with_capacity(prob.k);
    centers.push(space.point_iter().next().unwrap());
    let center_of: Vec<Option<CenterIdx>> = space.point_iter().map(|_| None).collect();
    clusterings.push(Clustering::new(centers,center_of, space));
    clusterings
}

use std::collections::VecDeque;

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
    visited: bool,
    discovered_by: Option<NodeKind>,
    kind : NodeKind,
}

enum NodeKind {
    Source, 
    Z,
    Center(CenterIdx),
    Point(PointIdx),
    Color(ColorIdx),
    T,
    Sink,
}

fn find_augmenting_path(state: &mut State, point_counter: usize, color_counter: usize, i : CenterIdx, sum_of_a: usize, prob: &ClusteringProblem) { 

    // create nodes:
//    let source_node = Node{visited:true, discovered_by:None, kind:NodeKind::Source};
    let mut z_node = Node{visited:false, discovered_by:None, kind:NodeKind::Z};
    
    let mut center_nodes : Vec<Node> = Vec::with_capacity(i + 1);
    for j in 0..i+1 {
        center_nodes.push(Node{visited: false, discovered_by:None, kind: NodeKind::Center(j)});
    }

    let mut point_nodes : Vec<Node> = Vec::with_capacity(point_counter);
    for p in 0..point_counter {
        point_nodes.push(Node{visited: false, discovered_by:None, kind:NodeKind::Point(p)});
    }

    let mut color_nodes : Vec<Node> = Vec::with_capacity(color_counter);
    for l in 0..color_counter {
        color_nodes.push(Node{visited: false, discovered_by:None, kind:NodeKind::Color(l)});
    }
    let mut t_node = Node{visited: false, discovered_by:None, kind: NodeKind::T};
    let mut sink_node = Node{visited:false, discovered_by:None, kind:NodeKind::Sink};

    // find augmenting path:

    let mut queue: VecDeque<NodeKind> = VecDeque::new();

    queue.push_back(NodeKind::Source);

    while !queue.is_empty() {
        let current_node = queue.pop_front().unwrap();

        match current_node {
            NodeKind::Source => {
                if state.flow_source_z < sum_of_a {
                    z_node.visited = true;
                    z_node.discovered_by = Some(NodeKind::Source);
                    queue.push_back(NodeKind::Z);
                }

                for j in 0..i+1 {
                    if !state.flow_source_center[j] {
                        center_nodes[j].visited = true;
                        center_nodes[j].discovered_by = Some(NodeKind::Source);
                        queue.push_back(NodeKind::Center(j));
                    }
                }


            }
            NodeKind::Z => {
                if state.flow_z_sink < i+1 && !sink_node.visited {
                    sink_node.visited = true;
                    sink_node.discovered_by = Some(NodeKind::Z);
                    break;
                }

                if state.flow_t_z > 0 && !t_node.visited {
                    t_node.visited = true;
                    t_node.discovered_by = Some(NodeKind::Z);
                    queue.push_back(NodeKind::T);
                }

                // To be continued with edges z to center



            }
            NodeKind::Center(idx) => {
            }

            NodeKind::Point(idx) => {

            }

            NodeKind::Color(idx) => {

            }

            NodeKind::T => {

            }

            NodeKind::Sink => {

            }

        }

    }




}

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
