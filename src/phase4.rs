use crate::{ClusteringProblem,Centers,Clustering,ColoredMetric};
use crate::types::{PointCount,PointIdx,ColorCount,ColorIdx,Distance,CenterIdx};
use super::OpeningList;
use std::collections::HashMap;

mod neighborhood;
use neighborhood::determine_neighborhood;

mod color_flow;
use color_flow::compute_assignment_by_flow;

/// an edge between a gonzales center and a color class;
/// lablled with the point and the distance between center and point
#[derive(Debug,Clone,Copy,PartialOrd,PartialEq)]
struct ColorEdge {
    d: Distance, // distance between the gonzales center and the point
    center: CenterIdx, // index to the gonzales center
    point: PointIdx, // a point that might be a final center
    color: ColorIdx, // the color of the point
}

type PNodeIdx = usize; // type for the index of a point node (nodes in the flow network are reindexed)
type CNodeIdx = usize; // type for the index of a color node

/// Given a metric space and a clustering problem,
/// finalize takes a vector of clusterings in which each center (except for one) covers a multple of L
/// points, and returns a vector of clusterings (which potentially more cluster as before) that also satisfy the representative constaints. 
///
pub(crate) fn finalize<'a, M : ColoredMetric>(space : &'a M, prob : &ClusteringProblem, opening_lists : Vec<Vec<OpeningList>>,  gonzales : &Centers) -> Vec<Clustering<'a>> {

    let sum_of_a: PointCount = prob.rep_interval.iter().map(|interval| interval.0).sum();

    let mut shifted_centers: Vec<Vec<ShiftedCenters>> = (0..prob.k).map(|i| Vec::with_capacity(i+1)).collect();


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
    


    // define the neighborhood of each gonzales center
    // each center has to be connected to the clostest a_l ponits of each color class l
    // fill up the neighborhood to k points with the closest remaining points (but never more than b_l
    // points of color l)
    let edges_of_cluster: Vec<Vec<ColorEdge>> = determine_neighborhood(space, prob, gonzales);

    // next, we solve k^2 flow problem. One for each gonzales set and each opening-vector
    for i in 0..prob.k {
        println!("\n\n************** i = {} ******************", i);

        // except for the opening vector the network can be defined now:

        // first, collect all edges in the neighborhood of centers 0 to i.
        let mut edges: Vec<&ColorEdge> = Vec::with_capacity((i+1) * prob.k); // a reference to all edges with centers in S_i
        for j in 0..i+1 {
            edges.extend(edges_of_cluster[j].iter());
        }
        // we want to consider them in increasing order
        edges.sort_by(|a,b| a.partial_cmp(b).unwrap());

        


        // Note that we only have k^2 edges so at most k^2 point and color classes can occur.
        // We reindex the relevant points and colors. The new indices are called point-node and color-node index and they are
        // used in the flow network
       
        let mut point_to_node: HashMap<PointIdx, PNodeIdx> = HashMap::with_capacity((i + 1) * prob.k); // maps a point index to the point-node index used in the flow network
        let mut point_counter = 0; // number of relevant points (= number of point-nodes)
        
        let mut a: Vec<PointCount> = Vec::with_capacity((i + 1) * prob.k); // lower bound given by the color-node index
        let mut b: Vec<PointCount> = Vec::with_capacity((i +1) * prob.k); // upper bound given by the color-node index
        let mut color_to_node: HashMap<ColorIdx, CNodeIdx> = HashMap::with_capacity((i + 1) * prob.k); // maps a color index to the color-node index used in the flow network
        let mut color_counter = 0; // number of relevant colors (= number of color-nodes)

        let mut point_idx_to_color_idx: Vec<ColorIdx> = Vec::with_capacity((i+1) * prob.k);

        for e in edges.iter() { // there are k^2 edges
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

      
        // The network is almost done, only the opening-vector is missing.


        for opening in opening_lists[i].iter() {

            let network = Network {
                edges : &edges,
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

            // for each eta vector we solve the flow problem individually:
            shifted_centers[i].push(compute_assignment_by_flow(&network));
            
        }

        println!("\nradii of shifts with i = {}: {:?}", i, shifted_centers[i].iter().map(|c| c.shift_radius).collect::<Vec<_>>());

    }





    let mut clusterings: Vec<Clustering> = Vec::with_capacity(prob.k);
    let mut centers = Centers::with_capacity(prob.k);
    centers.push(space.point_iter().next().unwrap());
    let center_of: Vec<Option<CenterIdx>> = space.point_iter().map(|_| None).collect();
    clusterings.push(Clustering::new(centers,center_of, space));
    clusterings
}


struct Network<'a> {
    edges : &'a Vec<&'a ColorEdge>,
    k : CenterIdx,
    opening : &'a OpeningList,
    i : CenterIdx,
    number_of_points : PointCount,
    number_of_colors : ColorCount,
    sum_of_a : PointCount,
    a : &'a Vec<PointCount>,
    b : &'a Vec<PointCount>,
    edges_of_cluster : &'a Vec<Vec<ColorEdge>>,
    point_to_node : &'a HashMap<PointIdx, PNodeIdx>,
    point_idx_to_color_idx : &'a Vec<CNodeIdx>,
    edges_by_color_node : &'a Vec<Vec<&'a ColorEdge>>,
    edges_by_point_node : &'a Vec<Vec<&'a ColorEdge>>,
}

// return value of a flow problem.
// Contains the shift_radius and new_centers (with their node index)
// as well as the origin, i.e., the cluster, of each new center.
struct ShiftedCenters {
    shift_radius : Distance, // the shift radius
    new_centers : Vec<PNodeIdx>,
    origins : Vec<CenterIdx> // the original gonzales center for each new center
}

