use super::ColorEdge;
use crate::utilities;
use crate::{Centers, ClusteringProblem, ColoredMetric};

/// define the neighborhood of each gonzalez center
/// each center has to be connected to the clostest a_l ponits of each color class l
/// fill up the neighborhood to k points with the closest remaining points (but never more than b_l
/// points of color l)
pub(super) fn determine_neighborhood<M: ColoredMetric>(
    space: &M,
    prob: &ClusteringProblem,
    gonzalez: &Centers,
) -> Vec<Vec<ColorEdge>> {
    let mut edges_of_cluster: Vec<Vec<ColorEdge>> = vec![Vec::with_capacity(prob.k); prob.k];

    (0..prob.k).for_each(|i| {
        // first we make edge lists for each color class
        let mut edges_by_color: Vec<Vec<ColorEdge>> = vec![Vec::new(); space.gamma()];

        // for this we create all edges from gonzalez center i to all points
        for p in space.point_iter() {
            let color = space.color(p);
            edges_by_color[space.color(p)].push(ColorEdge {
                d: space.dist(gonzalez.get(i, space), p),
                center: i,
                point: p.idx(),
                color,
            });
        }

        let restricted_colors = if space.gamma() < prob.rep_intervals.len() {
            space.gamma()
        } else {
            prob.rep_intervals.len()
        }; // min{gamma, rep_intervals.len()}
        let mut remaining_edges: Vec<ColorEdge> = Vec::new();
        let mut num_edges_to_fill = prob.k;

        (0..restricted_colors).for_each(|c| {
            // take only the b smallest in each class; all others cannot play any role.
            utilities::truncate_to_smallest(&mut edges_by_color[c], prob.rep_intervals[c].1);

            // the a smallest have to be present for sure, so each center can satisfy this
            // condition by itself;
            // Hence, the smallest a elements should for now remain in edges_by_color[c] and the
            // others b-a edges are collected in remaining_edges:
            let a = prob.rep_intervals[c].0;
            utilities::split_in_two_at(&mut edges_by_color[c], a);
            let bigger = edges_by_color[c].split_off(a);
            remaining_edges.extend(bigger);
            num_edges_to_fill -= a;
        });

        // for all colors without restriction fill them into the remaining edges:
        (restricted_colors..space.gamma()).for_each(|c| {
            remaining_edges.append(&mut edges_by_color[c]);
        });

        // so far sum(a) edges where chosen, so we fill up with the nearest k - sum(a) edges,
        // independent of the color
        utilities::truncate_to_smallest(&mut remaining_edges, num_edges_to_fill);

        // put all edges into one list:
        (0..restricted_colors).for_each(|c| {
            edges_of_cluster[i].append(&mut edges_by_color[c]);
        });
        edges_of_cluster[i].extend(remaining_edges);

        // there are now max k edges in edges_of_cluster[i]; so we can sort them:
        edges_of_cluster[i].sort_by(|a, b| a.partial_cmp(b).unwrap());
    });
    edges_of_cluster
}
