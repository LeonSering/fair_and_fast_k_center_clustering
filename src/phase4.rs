use crate::{ClusteringProblem,Centers,ColoredMetric};
use crate::types::{PointCount,PointIdx,ColorCount,ColorIdx,Distance,CenterIdx};
use super::datastructures::{OpeningList,NewCenters};
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

// return value of a flow problem.
// Contains the assignment_radius and new_centers (with their node index)
// as well as the origin, i.e., the cluster, of each new center.
#[derive(Debug)]
struct ShiftedCenters {
    forrest_radius : Distance, // the underlying forrest_radius from phase 3
    assignment_radius : Distance, // the shift radius

    new_centers : Vec<PNodeIdx>,
    origins : Vec<CenterIdx> // the original gonzales center for each new center
}


type PNodeIdx = usize; // type for the index of a point node (nodes in the flow network are reindexed)
type CNodeIdx = usize; // type for the index of a color node

/// Given a metric space and a clustering problem,
/// phase4 takes a vector of clusterings in which each center (except for one) covers a multple of L
/// points, and returns a single list of new centers that also satisfy the representative constaints and has minimum shifting radius
///
pub(crate) fn phase4<M : ColoredMetric>(space : &M, prob : &ClusteringProblem, opening_lists : Vec<Vec<OpeningList>>,  gonzales : &Centers) -> Vec<NewCenters> {

    let sum_of_a: PointCount = prob.rep_intervals.iter().map(|interval| interval.0).sum();

    let mut shifted_centers: Vec<Option<ShiftedCenters>> = (0..prob.k).map(|_| None).collect();
    let mut node_to_point_list: Vec<Vec<PointIdx>> = Vec::with_capacity(prob.k);


    ////TEMP: For test reasons we also look at flow problems with ALL edges present:
    //let mut all_edges_of_cluster: Vec<Vec<ColorEdge>> = vec!(Vec::with_capacity(space.n());prob.k);
    //for i in 0..prob.k {
    //    for p in space.point_iter() {
    //        let color = space.color(p);
    //        all_edges_of_cluster[i].push(ColorEdge{
    //            d: space.dist(gonzales.get(i), p),
    //            center: i,
    //            point: p.idx(),
    //            color});

    //    }
    //    all_edges_of_cluster[i].sort_by(|a,b| a.partial_cmp(b).unwrap());
    //}
    //// TEMP END



    // define the neighborhood of each gonzales center
    // each center has to be connected to the clostest a_l ponits of each color class l
    // fill up the neighborhood to k points with the closest remaining points (but never more than b_l
    // points of color l)
    let edges_of_cluster: Vec<Vec<ColorEdge>> = determine_neighborhood(space, prob, gonzales);

    // next, we solve k^2 flow problem. One for each gonzales set and each opening-vector
    for i in 0..prob.k {
        // println!("\n\n************** i = {} ******************", i);

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
        let mut node_to_point: Vec<PointIdx> = Vec::with_capacity((i+1) * prob.k); // maps a point-node index to the original point index
        let mut point_counter = 0; // number of relevant points (= number of point-nodes)

        let mut a: Vec<PointCount> = Vec::with_capacity((i + 1) * prob.k); // lower bound given by the color-node index
        let mut b: Vec<PointCount> = Vec::with_capacity((i +1) * prob.k); // upper bound given by the color-node index
        let mut color_to_node: HashMap<ColorIdx, CNodeIdx> = HashMap::with_capacity((i + 1) * prob.k); // maps a color index to the color-node index used in the flow network
        let mut color_counter = 0; // number of relevant colors (= number of color-nodes)

        let mut point_idx_to_color_idx: Vec<ColorIdx> = Vec::with_capacity((i+1) * prob.k);

        for e in edges.iter() { // there are k^2 edges
            if !color_to_node.contains_key(&e.color) {
                color_to_node.insert(e.color, color_counter);
                if e.color >= prob.rep_intervals.len() {
                    a.push(0);
                    b.push(<PointCount>::MAX);
                } else {
                    a.push(prob.rep_intervals[e.color].0);
                    b.push(prob.rep_intervals[e.color].1);
                }
                color_counter += 1;
            }
            if !point_to_node.contains_key(&e.point) {
                point_to_node.insert(e.point, point_counter);
                node_to_point.push(e.point);
                point_idx_to_color_idx.push(*color_to_node.get(&e.color).unwrap());
                point_counter += 1;
            }
        }

        node_to_point_list.push(node_to_point); // save node_to_point for recreate the point index at the very end of the phase

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
        // println!("point_to_node: {:?}", point_to_node);

        let mut current_best_radius = <Distance>::MAX;
        for (j, opening) in opening_lists[i].iter().enumerate() {

            // first test if eta-vector makes sense at all:
            if sum_of_a > opening.eta.iter().sum() {
                println!("  - Cannot open enough new centers as sum_of_a = {} > eta = {}: opening_list: {}; i = {}", sum_of_a, opening.eta.iter().sum::<PointCount>(), opening, i);
                continue;
            }

            // we can skip the flow computation if the forrest_radius alone is bigger than the
            // previously best sum of forrest_radius + assignment_radius
            if current_best_radius <= opening.forrest_radius {
                #[cfg(debug_assertions)]
                println!("  - Flow-network for i = {} and open_list = {} is obsolet as the forrest_radius is bigger than forrest_radius + assignment_radius of earlier open list (= {}).", i, opening, current_best_radius);
                continue;

            }

            // now test if the exact (or less restrictive) flow problem has been solved already:
            let mut has_been_solved = false;
            for old in 0..j {
                let old_eta = &opening_lists[i][old].eta;
                let mut old_equal_or_bigger = true;
                for l in 0..old_eta.len() {
                    if opening.eta[l] > old_eta[l] {
                        old_equal_or_bigger = false;
                    }
                }
                if old_equal_or_bigger {
                    has_been_solved = true;
                }
            }
            if has_been_solved {
                #[cfg(debug_assertions)]
                println!("  - Flow-network for i = {} and opening_list = {} has been solved before (or less restricted version). Skipping this flow computation.", i, opening);
                continue;
            }

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
                point_idx_to_color_idx : &point_idx_to_color_idx,
                edges_by_color_node : &edges_by_color_node,
                edges_by_point_node : &edges_by_point_node,
            };

            // for each eta vector we solve the flow problem individually:
            let new = compute_assignment_by_flow(&network);

            if current_best_radius > new.forrest_radius + new.assignment_radius {
                current_best_radius = new.forrest_radius + new.assignment_radius;
                shifted_centers[i] = Some(new);
            }

        }

        // println!("\nradii of shifts with i = {}: {:?}", i, shifted_centers[i].iter().map(|c| c.assignment_radius).collect::<Vec<_>>());


    }


    let mut new_centers: Vec<NewCenters> = Vec::with_capacity(prob.k+1); // return vector

    // determine centers with minimal shift radius + tree radius:
    for (i,c) in shifted_centers.into_iter().enumerate() {
        let best_centers = c.unwrap();
        let mut centers : Centers = Centers::with_capacity(best_centers.new_centers.len());
        let mut new_center_idx_of_cluster = vec!(Vec::new(); i+1);
        for (l,&idx) in best_centers.new_centers.iter().enumerate() {
            let point = space.get_point(node_to_point_list[i][idx]).unwrap();
            centers.push(point);
            new_center_idx_of_cluster[best_centers.origins[l]].push(centers.m()-1);
        }

        // println!("Best centers: {:?} and in idx form grouped by cluster: {:?}", new_centers, new_center_idx_of_cluster);

        new_centers.push(NewCenters{as_points: centers, forrest_radius: best_centers.forrest_radius, assignment_radius: best_centers.assignment_radius, new_centers_of_cluster: new_center_idx_of_cluster});

    }


    new_centers
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

