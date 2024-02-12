use super::{
    flow::{add_edge, initialize_state},
    Edge,
};
use crate::clustering::{Centers, Clustering};
use crate::space::ColoredMetric;
use crate::types::{Distance, PointCount};
use rayon::prelude::*;
use std::collections::VecDeque;

use std::time;

pub(crate) fn make_private_with_sorting<M: ColoredMetric>(
    space: &M,
    privacy_bound: PointCount,
    centers: &Centers,
) -> Vec<Clustering> {
    //Return value should be partialClustering

    #[cfg(debug_assertions)]
    println!("  - Phase 2 with sorting:");

    #[cfg(debug_assertions)]
    let time_start_make_private_with_sorting = time::Instant::now();
    let k = centers.m();
    // create edges: care, edge.left stores the index of the gonzalez center (0,...,k-1).
    let mut edges: Vec<Edge> = Vec::with_capacity(k * space.n());
    for (j, c) in centers.get_all(space).iter().enumerate() {
        for p in space.point_iter() {
            edges.push(Edge {
                d: space.dist(c, p),
                left: j,
                right: p,
            });
        }
    }

    #[cfg(debug_assertions)]
    let time_end_edge_creation = time::Instant::now();
    #[cfg(debug_assertions)]
    println!(
        "    - creation of {} edges takes: {:?}.",
        space.n() * k,
        time_end_edge_creation.duration_since(time_start_make_private_with_sorting)
    );

    // edges.sort_by(|a, b| a.partial_cmp(b).unwrap());
    edges.par_sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

    let mut edge_iter = edges.iter();

    #[cfg(debug_assertions)]
    let time_after_sorting = time::Instant::now();
    #[cfg(debug_assertions)]
    println!(
        "    - sorting {} edges with rayon::par_sort takes: {:?}.",
        space.n() * k,
        time_after_sorting.duration_since(time_end_edge_creation)
    );

    // step 2: solve flow problem

    let mut clusterings: Vec<Clustering> = Vec::with_capacity(k);

    let mut state = initialize_state(space.n(), k, privacy_bound == 0);

    let mut pending: Vec<VecDeque<Edge>> = (0..k).map(|_| VecDeque::new()).collect();

    let mut i = 0; // currently processing gonzalez set containing center 0, ..., i; here, i goes from 0 to k-1

    let mut _current_d = <Distance>::MIN;

    while i < k {
        // extend set of gonzalez centers one by one
        // this is the main while-loop that deals with each center set daganzo[i] for i = 1, ..., k

        while !pending[i].is_empty() {
            let e = pending[i].pop_front().unwrap();
            debug_assert_eq!(i, e.left); // e.left should be i, (the index of the i-th center in gonzalez)
            add_edge(e, i, k, privacy_bound, &mut state);
        }

        while state.max_flow < (i + 1) * privacy_bound {
            let e = *edge_iter.next().expect("no edges remaining");

            let t = e.left; // t is the index of the center on the left of e
            if t > i {
                // in this case the left side (t) corresponds to a center not yet considered, so we
                // postpone the processing of the arc to the point when i == t
                pending[t].push_back(e);
            } else {
                // in this case we do add edge e.
                add_edge(e, i, k, privacy_bound, &mut state);
            }
        }

        // create new clustering:
        let mut center_prefix = Centers::with_capacity(i + 1);
        for c in centers.get_all(space).iter().take(i + 1) {
            center_prefix.push(c);
        }

        let clustering = Clustering::new(center_prefix, state.center_of.clone(), space);

        clusterings.push(clustering);

        i += 1;
    }

    #[cfg(debug_assertions)]
    let time_after_flow_with_sorting = time::Instant::now();
    #[cfg(debug_assertions)]
    println!(
        "    - solving flow problems and settle all {} centers takes (sorting): {:?}.",
        k,
        time_after_flow_with_sorting.duration_since(time_after_sorting)
    );

    for clustering in clusterings.iter_mut() {
        clustering.fill_up(space);
    }

    #[cfg(debug_assertions)]
    let time_after_filling_up_with_sorting = time::Instant::now();

    #[cfg(debug_assertions)]
    println!(
        "    - filling up all remaining points (sorting) takes: {:?}.",
        time_after_filling_up_with_sorting.duration_since(time_after_flow_with_sorting)
    );

    clusterings
}
