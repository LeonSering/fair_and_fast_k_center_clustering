use super::datastructures::{NewCenters, OpeningList};
use crate::types::{CenterIdx, ColorCount, ColorIdx, Distance, PointCount, PointIdx};
use crate::{Centers, ClusteringProblem, ColoredMetric};
use std::collections::{hash_map, HashMap, VecDeque};

mod neighborhood;
use neighborhood::determine_neighborhood;

mod color_flow;
use color_flow::compute_assignment_by_flow;

// for parallel execution:
use rayon::ThreadPoolBuilder;
use std::sync::mpsc;

/// an edge between a gonzalez center and a color class;
/// lablled with the point and the distance between center and point
#[derive(Debug, Clone, Copy, PartialOrd, PartialEq)]
struct ColorEdge {
    d: Distance,       // distance between the gonzalez center and the point
    center: CenterIdx, // index to the gonzalez center
    point: PointIdx,   // a point that might be a final center
    color: ColorIdx,   // the color of the point
}

// return value of a flow problem.
// Contains the assignment_radius and new_centers (with their node index)
// as well as the origin, i.e., the cluster, of each new center.
#[derive(Debug)]
struct ShiftedCenters {
    forrest_radius: Distance,    // the underlying forrest_radius from phase 3
    assignment_radius: Distance, // the shift radius

    new_centers: Vec<PNodeIdx>,
    origins: Vec<CenterIdx>, // the original gonzalez center for each new center
}

type PNodeIdx = usize; // type for the index of a point node (nodes in the flow network are reindexed)
type CNodeIdx = usize; // type for the index of a color node

/// Given a metric space and a clustering problem,
/// phase4 takes a vector of clusterings in which each center (except for one) covers a multple of L
/// points, and returns a single list of new centers that also satisfy the representative constaints and has minimum shifting radius
pub(crate) fn phase4<M: ColoredMetric>(
    space: &M,
    prob: &ClusteringProblem,
    mut opening_lists: Vec<Vec<OpeningList>>,
    gonzalez: &Centers,
    thread_count: usize,
) -> (Vec<NewCenters>, Vec<usize>) {
    let sum_of_a: PointCount = prob.rep_intervals.iter().map(|interval| interval.0).sum();

    let mut shifted_centers: Vec<Option<ShiftedCenters>> = (0..prob.k).map(|_| None).collect();
    let mut counts = vec![0; prob.k]; // counter for how many flow problems are solved.
    let mut node_to_point_list: Vec<Vec<PointIdx>> = Vec::with_capacity(prob.k);

    // define the neighborhood of each gonzalez center
    // each center has to be connected to the clostest a_l ponits of each color class l
    // fill up the neighborhood to k points with the closest remaining points (but never more than b_l
    // points of color l)

    let edges_of_cluster: Vec<Vec<ColorEdge>> = determine_neighborhood(space, prob, gonzalez);

    let thread_pool = ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .unwrap();
    let mut receivers = VecDeque::with_capacity(prob.k);

    thread_pool.scope(|pool| {
        // define this scope as refs to edges_of_cluster should be moved into threads but it cannot be 'static

        // next, we solve k^2 flow problem. One for each gonzalez set and each opening-vector
        for i in (0..prob.k).rev() {
            let opening_list = opening_lists.pop().unwrap();
            let (tx, rx) = mpsc::channel();
            receivers.push_front(rx);

            let edges_of_cluster_ref = &edges_of_cluster; // take a ref in order to move it the the thread

            let problem = &prob;
            pool.spawn(move |_| {
                let (centers, node_to_point, count) =
                    shift_centers(problem, sum_of_a, i, edges_of_cluster_ref, opening_list);
                tx.send((centers, node_to_point, count)).unwrap();
            });
        }
        for i in 0..prob.k {
            let (centers, node_to_point, count) = receivers[i].recv().unwrap();
            shifted_centers[i] = centers;
            node_to_point_list.push(node_to_point); // save node_to_point for recreate the point index at the very end of the phase
            counts[i] = count;
        }
    });

    let mut new_centers: Vec<NewCenters> = Vec::with_capacity(prob.k + 1); // return vector

    // determine centers with minimal shift radius + tree radius:
    for (i, c) in shifted_centers.into_iter().enumerate() {
        let best_centers = c.unwrap();
        let mut centers: Centers = Centers::with_capacity(best_centers.new_centers.len());
        let mut new_center_idx_of_cluster = vec![Vec::new(); i + 1];
        for (l, &idx) in best_centers.new_centers.iter().enumerate() {
            let point = space.get_point(node_to_point_list[i][idx]).unwrap();
            centers.push(point);
            new_center_idx_of_cluster[best_centers.origins[l]].push(centers.m() - 1);
        }

        new_centers.push(NewCenters {
            as_points: centers,
            forrest_radius: best_centers.forrest_radius,
            assignment_radius: best_centers.assignment_radius,
            new_centers_of_cluster: new_center_idx_of_cluster,
        });
    }

    (new_centers, counts)
}

/// The following function is executed once for each gonzalez prefix C_i.
/// First the flow network network is prepared.
/// For each of the i opening_list (one for each forest) it checks wether this flow problem has
/// been solved already. Then it executes the flow algorithm to obtain new shifted centers.
/// It only return the best shifted centers for this Gonzalez prefix.
///
/// Output: shifted_centers (None if not feasible),
/// node_to_point : maps node index to point index,
/// count : number of flow problems that have been solved,
fn shift_centers(
    prob: &ClusteringProblem,
    sum_of_a: PointCount,
    i: CenterIdx,
    edges_of_cluster: &Vec<Vec<ColorEdge>>,
    opening_lists: Vec<OpeningList>,
) -> (Option<ShiftedCenters>, Vec<PointIdx>, usize) {
    // except for the opening vector the network can be defined now:

    // first, collect all edges in the neighborhood of centers 0 to i.
    let mut edges: Vec<&ColorEdge> = Vec::with_capacity((i + 1) * prob.k); // a reference to all edges with centers in S_i
    (0..i + 1).for_each(|j| {
        edges.extend(edges_of_cluster[j].iter());
    });
    // we want to consider them in increasing order
    edges.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Note that we only have k^2 edges so at most k^2 point and color classes can occur.
    // We reindex the relevant points and colors. The new indices are called point-node and color-node index and they are
    // used in the flow network

    let mut point_to_node: HashMap<PointIdx, PNodeIdx> = HashMap::with_capacity((i + 1) * prob.k); // maps a point index to the point-node index used in the flow network
    let mut node_to_point: Vec<PointIdx> = Vec::with_capacity((i + 1) * prob.k); // maps a point-node index to the original point index
    let mut point_counter = 0; // number of relevant points (= number of point-nodes)

    let mut a: Vec<PointCount> = Vec::with_capacity((i + 1) * prob.k); // lower bound given by the color-node index
    let mut b: Vec<PointCount> = Vec::with_capacity((i + 1) * prob.k); // upper bound given by the color-node index
    let mut color_to_node: HashMap<ColorIdx, CNodeIdx> = HashMap::with_capacity((i + 1) * prob.k); // maps a color index to the color-node index used in the flow network
    let mut color_counter = 0; // number of relevant colors (= number of color-nodes)

    let mut point_idx_to_color_idx: Vec<ColorIdx> = Vec::with_capacity((i + 1) * prob.k);

    for e in edges.iter() {
        // there are k^2 edges
        if let hash_map::Entry::Vacant(entry) = color_to_node.entry(e.color) {
            entry.insert(color_counter);
            if e.color >= prob.rep_intervals.len() {
                a.push(0);
                b.push(<PointCount>::MAX);
            } else {
                a.push(prob.rep_intervals[e.color].0);
                b.push(prob.rep_intervals[e.color].1);
            }
            color_counter += 1;
        }
        if let hash_map::Entry::Vacant(entry) = point_to_node.entry(e.point) {
            entry.insert(point_counter);
            node_to_point.push(e.point);
            point_idx_to_color_idx.push(*color_to_node.get(&e.color).unwrap());
            point_counter += 1;
        }
    }

    // edges should also be referrable from their points and their color classes:

    let mut edges_by_color_node: Vec<Vec<&ColorEdge>> = vec![Vec::new(); color_counter];
    let mut edges_by_point_node: Vec<Vec<&ColorEdge>> = vec![Vec::new(); point_counter];

    let mut count = 0;

    for e in edges.iter() {
        edges_by_color_node[color_to_node[&e.color]].push(e);
        edges_by_point_node[point_to_node[&e.point]].push(e);
    }
    edges_by_color_node.sort_by(|a, b| a.partial_cmp(b).unwrap());
    edges_by_point_node.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // The network is almost done, only the opening-vector is missing.

    let mut current_best_radius = <Distance>::MAX;

    let mut centers = None;

    for (j, opening) in opening_lists.iter().enumerate() {
        // we can skip the flow computation if the forrest_radius alone is bigger than the
        // previously best sum of forrest_radius + assignment_radius
        if current_best_radius <= opening.forrest_radius {
            continue;
        }

        // first test if eta-vector makes sense at all:
        if sum_of_a > opening.eta.iter().sum() {
            continue;
        }

        // now test if the exact (or less restrictive) flow problem has been solved already:
        let mut has_been_solved = false;
        for open_list in opening_lists.iter().take(j) {
            let old_eta = &open_list.eta;
            let mut old_equal_or_bigger = true;
            for (l, eta) in old_eta.iter().enumerate() {
                if opening.eta[l] > *eta {
                    old_equal_or_bigger = false;
                    break;
                }
            }
            if old_equal_or_bigger {
                has_been_solved = true;
                break;
            }
        }
        if has_been_solved {
            continue;
        }

        let network = Network {
            edges: &edges,
            k: prob.k,
            opening,
            i,
            number_of_points: point_counter,
            number_of_colors: color_counter,
            sum_of_a,
            a: &a,
            b: &b,
            edges_of_cluster,
            point_to_node: &point_to_node,
            point_idx_to_color_idx: &point_idx_to_color_idx,
            edges_by_color_node: &edges_by_color_node,
            edges_by_point_node: &edges_by_point_node,
        };

        count += 1;

        // for each eta vector we solve the flow problem individually:
        let new = compute_assignment_by_flow(&network);

        if current_best_radius > new.forrest_radius + new.assignment_radius {
            current_best_radius = new.forrest_radius + new.assignment_radius;
            centers = Some(new);
        }
    }

    (centers, node_to_point, count)
}

struct Network<'a> {
    edges: &'a Vec<&'a ColorEdge>,
    k: CenterIdx,
    opening: &'a OpeningList,
    i: CenterIdx,
    number_of_points: PointCount,
    number_of_colors: ColorCount,
    sum_of_a: PointCount,
    a: &'a Vec<PointCount>,
    b: &'a Vec<PointCount>,
    edges_of_cluster: &'a Vec<Vec<ColorEdge>>,
    point_to_node: &'a HashMap<PointIdx, PNodeIdx>,
    point_idx_to_color_idx: &'a Vec<CNodeIdx>,
    edges_by_color_node: &'a Vec<Vec<&'a ColorEdge>>,
    edges_by_point_node: &'a Vec<Vec<&'a ColorEdge>>,
}
