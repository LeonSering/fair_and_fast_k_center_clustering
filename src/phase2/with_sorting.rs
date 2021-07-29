use crate::ClusteringProblem;
use crate::types::Distance;
use crate::space::ColoredMetric;
use crate::clustering::{Clustering,Centers};
use std::collections::VecDeque;
use super::{Edge,flow::{initialize_state,add_edge}};

pub(crate) fn make_private_with_sorting<'a, M : ColoredMetric>(space : &'a M, prob : &'a ClusteringProblem, gonzales : &Centers<'a>) -> Vec<Clustering<'a>> { //Return value should be partialClustering

// create edges: care, edge.left stores the index of the gonzales center (0,...,k-1).
    let mut edges : Vec<Edge> = Vec::with_capacity(prob.k * space.n());
    for (j, c) in gonzales.iter().enumerate() {
        for p in space.point_iter() {
            edges.push(Edge{
                d : space.dist(c, p),
                left : j,
                right : p,
            });
        }
    }

    println!("\n  ** make_private_with_sort with privacy_bound = {}\n", prob.privacy_bound);

    edges.sort_by(|a, b| a.partial_cmp(b).unwrap());
//    println!("Edges: {:?}", edges.iter().map(|x| x.d).collect::<Vec<_>>());

    let mut edge_iter = edges.iter();

    // step 2: solve flow problem
    
    let mut clusterings: Vec<Clustering> = Vec::with_capacity(prob.k);
    

    let mut state = initialize_state(space.n(), prob.k);
    
    let mut pending: Vec<VecDeque<Edge>> = (0..prob.k).map(|_| VecDeque::with_capacity(prob.k*prob.k)).collect();

    let mut i = 0; // currently processing gonzales set containing center 0, ..., i; here, i goes from 0 to k-1

    let mut current_d = <Distance>::MIN;
    while i < prob.k { // extend set of gonzales centers one by one
//        println!{"\n+++ center {} now considered!\n", i};
        // this is the main while-loop that deals with each center set daganzo[i] for i = 1, ..., k
        
        while !pending[i].is_empty() {
            let e = pending[i].pop_front().unwrap();
            assert_eq!(i, e.left); // e.left should be i, (the index of the i-th center in gonzales)
            add_edge(e, i, prob, &mut state);
//            println!("  Pending edge added: {:?}; \tmax_flow: {}", e, state.max_flow);
        }

        while state.max_flow < (i + 1) * prob.privacy_bound {

          
            let e = *edge_iter.next().expect("no edges remaining");

            let t = e.left; // t is the index of the center on the left of e
            if t > i {
                // in this case the left side (t) corresponds to a center not yet considered, so we
                // postpone the processing of the arc to the point when i == t
                pending[t].push_back(e);
            } else {
                // in this case we do add edge e.
                add_edge(e, i, prob, &mut state);
//                println!("  Adding: {:?};\tmaxflow: {}", e, state.max_flow);
            }
            current_d = e.d;
            
        }


        // at this point, we have identified the bucket that settles the set S_i
        #[cfg(debug_assertions)]
        assert_eq!(state.max_flow, (i + 1) * prob.privacy_bound, "The maximum flow value is bigger than allowed"); // we should have equality due to the capacities of arcs (= privacy_bound) between the source and the centers in S_i

        println!("\tRadius for center {} found: {}", i, current_d);

        // create new clustering:
        let mut centers = Centers::with_capacity(i+1);
        for c in gonzales.iter().take(i+1) {
            centers.push(c);
        }

        let clustering = Clustering::new(centers,state.center_of.clone(),space);
        #[cfg(debug_assertions)]
        assert_eq!(current_d,clustering.get_radius(), "Determined radius differs from cluster radius! This should never happen!");
        clusterings.push(clustering);
        
        i += 1;
    }
    for clustering in clusterings.iter_mut() {
        clustering.fill_up(space);
    }
    clusterings
}

