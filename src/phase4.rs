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

        


        // not that we only have k^2 edges so at most k^2 color classes can occur
        // We reindex the relevant colors. The new indices are called color-node index and they are
        // used in the flow network

        let mut a: Vec<PointCount> = Vec::with_capacity((i + 1) * prob.k); // lower bound given by the color-node index
        let mut b: Vec<PointCount> = Vec::with_capacity((i +1) * prob.k); // upper bound given by the color-node index
        let mut color_to_node: HashMap<ColorIdx, usize> = HashMap::with_capacity((i + 1) * prob.k); // maps a color index to the color-node index used in the flow network
        let mut color_counter = 0; // number of relevant colors (= number of color-nodes)
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
        }

//        println!("\n i = {}", i);
//        println!("number of edges: {:?}", edges.len());
//        println!("color_to_node: {:?}", color_to_node);
//        println!("a: {:?}", a);
//        println!("b: {:?}", b);




//        println!("\nedges_of_cluser[{}]: {:?}", i, edges_of_cluster[i]);
//        println!("** i = {}; openings = {:?}", i, openings);

        // it makes more sense to implements a general max flow algorithm
        // now we initialize the flow:


        // without any edges present the maximal flow has only flow an the path (source, z, sink)
        // of value min{sum_of_a, i}
        let initial_flow_value = if sum_of_a < i {sum_of_a} else {i}; // = min {sum_of_a, i}
        let mut state = State{
            current_largest_edge: ColorEdge{ d: <Distance>::MIN, center: 0, point: 0, color: 0},
            flow_source_center: vec!(false; i),
            flow_source_z: initial_flow_value, 
            flow_z_center: vec!(0; i),
            edge_flow_carrying: vec!(vec!(false;prob.k);i),
            flow_color_t: vec!(0;color_counter),
            flow_t_z: 0,
            flow_color_sink: vec!(0;color_counter),
            flow_z_sink: initial_flow_value,
            direct_res_path_from_color_to_source: vec!(false;color_counter),
            direct_res_path_from_center_to_sink: vec!(false; i)};


    }




    let mut clusterings: Vec<Clustering> = Vec::with_capacity(prob.k);
    let mut centers = Centers::with_capacity(prob.k);
    centers.push(space.point_iter().next().unwrap());
    let center_of: Vec<Option<CenterIdx>> = space.point_iter().map(|_| None).collect();
    clusterings.push(Clustering::new(centers,center_of, space));
    clusterings
}

struct State {
    // network:
    current_largest_edge: ColorEdge, // is used to indicate which edges are present in the network

    // flow:
    flow_source_center: Vec<bool>, // true <=> the arc from source to center is saturated (flow = 1)
    flow_source_z: usize, // has capacity sum_of_a
    flow_z_center: Vec<usize>, // one entry for each center j (0 .. i); these arcs have capacity eta_j - 1
    edge_flow_carrying: Vec<Vec<bool>>, // edge_flow_carrying[j][l] = true <=> edges_of_cluster[j][l] is saturated (flow = 1) 
    flow_color_t: Vec<usize>, // one entry for each color class l; each has a capacity of a_l - b_l
    flow_t_z: usize, // has capacity k- sum_of_a
    flow_color_sink: Vec<usize>, // one entry for each color class l; each has a capacity of a_l
    flow_z_sink: usize, // has capacity i

    // utility:
    direct_res_path_from_color_to_source: Vec<bool>, // one entry for each color class l; true <=> there exists a center s_j with a unsaturated edge to color class l, and flow_source_center[j] = false.
    direct_res_path_from_center_to_sink: Vec<bool>, // one entry for each center j; true <=> there exists a color class l with unsaturated edge (j,l) and unsaturated arc from l to sink.
}

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
