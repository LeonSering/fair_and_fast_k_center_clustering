use crate::{ClusteringProblem,Clustering,clustering::{Centers,new_centers}};
use super::{Edge,buckets::split_at_median,flow::{State,add_edge,remove_edge}};

// note that edge_cursor points at the edge that has not been added yet
// (the edge at edge_curser-1 has been added already)
pub fn settle<'a, 'b>(edge_cursor: usize, bucket: &mut Vec<Edge<'a>>, i: usize, prob: &ClusteringProblem, state: &mut State<'a>, gonzales: &Centers<'b>) -> Clustering<'b>{
    let mut cursor = edge_cursor;

    println!("  Edge_cursor: {}", edge_cursor);
    println!("  Bucket to settle: {:?}\n", bucket.iter().map(|x| x.d).collect::<Vec<f32>>());

    let edges_present : bool;
    // first clear or fill bucket:
    if edge_cursor >= bucket.len()/2 { // in this case fill the bucket
        while cursor < bucket.len() {
            add_edge(bucket[cursor],i,prob,state);
            println!("\tInitial adding: {:?} in binary search for center: {};\tmax_flow: {}", bucket[cursor], i, state.max_flow);
            cursor += 1;
        }
        edges_present = true;
    } else { // in this case empty bucket
        while cursor > 0 {
            cursor -= 1;
            remove_edge(bucket[cursor],i,prob,state);
            println!("\tInitial removing: {:?} in binary search for center: {};\tmax_flow: {}", bucket[cursor], i, state.max_flow);
        }
        edges_present = false;
    }

    let radius = search_for_radius(edges_present, bucket, &mut cursor, i, prob, state);

    println!("\n==> Radius for center {} found: {};\tmax_flow: {}\n", i, radius, state.max_flow);
//    println!("center_of: {:?}", state.center_of);
//    println!("number_of_points_covered: {:?}", state.number_of_covered_points);

    let mut centers = new_centers(i+1);
    for c in gonzales.iter().take(i+1) {
        centers.push(c);
    }

    let clustering = Clustering {
        centers,  // centers are the gonzales centers 0,...,i
        radius, // TODO
        center_of : state.center_of.iter().map(|c| match c {
                        Some(idx) => Some(gonzales.get(*idx)),
                        None => None,
                    }).collect(),
    };

    println!("  Bucket after settling: {:?}\n", bucket.iter().map(|x| x.d).collect::<Vec<f32>>());

    // empty bucket for the final time
    while cursor > 0 {
        cursor -= 1;
        remove_edge(bucket[cursor],i,prob,state);
        println!("\tClean up removing: {:?} in binary search for center: {};\tmaxflow: {}", bucket[cursor], i,state.max_flow);
    }

    clustering
}

// search recursevely (binary research) for the correct radius.
// Input: the edge list, current center list 0,...,i, the cluster prob, the flow state
// If edges_present == TRUE: expects that ALL edges in the list are part of the flow network
// If edges_present == FALSE: expects that NONE of the edges in the list are part of the flow network
// and expects that state.max_flow < max_flow_target
//
// Output: None if max_flow_target cannot be reached or Some(radius) if found
//   list becomes more sorted
fn search_for_radius<'a>(edges_present: bool, list: &mut Vec<Edge<'a>>, cursor : &mut usize, i : usize, prob: &ClusteringProblem, state: &mut State<'a>) -> f32 {
    println!("\n  List to settle: {:?} edges_present: {}", list.iter().map(|x| x.d).collect::<Vec<f32>>(), edges_present);
    let list_len = list.len();
    assert!(list_len > 0, "Empty list in binary search");
    if list_len == 1 {
        if !edges_present {
            add_edge(list[0],i,prob,state);
            *cursor += 1;
            println!("\n\tAdd final edge: {:?} in binary search for center: {};\tmax_flow: {}", list[0], i, state.max_flow);
        }
        assert_eq!(state.max_flow, (i+1) * prob.privacy_bound, "Something went wrong in the binary search.");
        return list[0].d;
    }

    let (mut smaller, mut bigger) = split_at_median(list);
    println!("     smaller: {:?} bigger: {:?}\n", smaller.iter().map(|x| x.d).collect::<Vec<f32>>(), bigger.iter().map(|x| x.d).collect::<Vec<f32>>());



    // take care that smaller edges are added and bigger edges are not present in the flow network
    if edges_present {
        for e in bigger.iter().rev() {
            remove_edge(*e, i, prob, state);
            println!("\tTry removing: {:?} in binary search for center: {};\tmax_flow: {}", e, i,state.max_flow);
            *cursor -= 1;
        }
    } else {
        for e in smaller.iter() {
            add_edge(*e, i, prob, state);
            println!("\tTry adding: {:?} in binary search for center: {};\tmax_flow: {}", e, i,state.max_flow);
            *cursor += 1;
        }
    }


    let radius: f32;
//    let mut left : Vec<Edge> = Vec::with_capacity(list_len);
//    let mut right : Vec<Edge> = Vec::with_capacity(list_len/2 + 1);
    if state.max_flow >= (i+1) * prob.privacy_bound { // we need to settle in smaller
        radius = search_for_radius(true, &mut smaller, cursor, i, prob, state);
    } else { // we need to settle in bigger
        radius = search_for_radius(false, &mut bigger, cursor, i, prob, state);
    }

    // concatenate smaller and bigger;
    list.clear();
    list.append(&mut smaller);
    list.append(&mut bigger);
    radius
}


