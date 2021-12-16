use crate::types::{Distance,PointCount};
use crate::space::ColoredMetric;
use crate::clustering::{Clustering,Centers};
use std::collections::VecDeque;
use super::{Edge,flow::{initialize_state,add_edge}};
use rayon::prelude::*;

use std::time;

pub(crate) fn make_private_with_sorting<M : ColoredMetric>(space : &M, privacy_bound: PointCount, centers : &Centers) -> Vec<Clustering> { //Return value should be partialClustering

    let time_start_make_private_with_sorting = time::Instant::now();
    let k = centers.m();
// create edges: care, edge.left stores the index of the gonzalez center (0,...,k-1).
    let mut edges : Vec<Edge> = Vec::with_capacity(k * space.n());
    for (j, c) in centers.get_all(space).iter().enumerate() {
        for p in space.point_iter() {
            edges.push(Edge{
                d : space.dist(c, p),
                left : j,
                right : p,
            });
        }
    }

    // println!("\n  ** make_private_with_sort with privacy_bound = {}\n", privacy_bound);

   edges.sort_by(|a, b| a.partial_cmp(b).unwrap());
    // edges.par_sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
//    println!("Edges: {:?}", edges.iter().map(|x| x.d).collect::<Vec<_>>());

    let mut edge_iter = edges.iter();

    let time_after_sorting = time::Instant::now();
    println!("  - sorting {} edges takes: {:?}.", space.n() * k, time_after_sorting.duration_since(time_start_make_private_with_sorting));
    // step 2: solve flow problem

    let mut clusterings: Vec<Clustering> = Vec::with_capacity(k);


    let mut state = initialize_state(space.n(), k, if privacy_bound == 0 {true} else {false});

    let mut pending: Vec<VecDeque<Edge>> = (0..k).map(|_| VecDeque::new()).collect();

    let mut i = 0; // currently processing gonzalez set containing center 0, ..., i; here, i goes from 0 to k-1

    let mut current_d = <Distance>::MIN;
    while i < k { // extend set of gonzalez centers one by one
//        println!{"\n+++ center {} now considered!\n", i};
        // this is the main while-loop that deals with each center set daganzo[i] for i = 1, ..., k

        while !pending[i].is_empty() {
            let e = pending[i].pop_front().unwrap();
            assert_eq!(i, e.left); // e.left should be i, (the index of the i-th center in gonzalez)
            add_edge(e, i, k, privacy_bound, &mut state);
//            println!("  Pending edge added: {:?}; \tmax_flow: {}", e, state.max_flow);
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
//                println!("  Adding: {:?};\tmaxflow: {}", e, state.max_flow);
            }
            current_d = e.d;

        }


        #[cfg(debug_assertions)]
        assert_eq!(state.max_flow, (i + 1) * privacy_bound, "The maximum flow value is bigger than allowed"); // we should have equality due to the capacities of arcs (= privacy_bound) between the source and the centers in S_i

        #[cfg(debug_assertions)]
        println!("\tRadius for center {} found: {}", i, current_d);

        // create new clustering:
        let mut center_prefix = Centers::with_capacity(i+1);
        for c in centers.get_all(space).iter().take(i+1) {
            center_prefix.push(c);
        }

        let clustering = Clustering::new(center_prefix,state.center_of.clone(),space);
        #[cfg(debug_assertions)]
        assert_eq!(current_d,clustering.get_radius(), "Determined radius differs from cluster radius! This should never happen!");
        clusterings.push(clustering);

        i += 1;
    }

    let time_after_flow_with_sorting = time::Instant::now();
    println!("  - solving flow problems and settle all {} centers takes (sorting): {:?}.", k, time_after_flow_with_sorting.duration_since(time_after_sorting));

    for clustering in clusterings.iter_mut() {
        clustering.fill_up(space);
    }

    let time_after_filling_up_with_sorting = time::Instant::now();
    println!("  - filling up all remaining points (sorting) takes: {:?}.", time_after_filling_up_with_sorting.duration_since(time_after_flow_with_sorting));

    clusterings
}

