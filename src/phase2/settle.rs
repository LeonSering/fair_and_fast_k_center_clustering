use crate::{Clustering,Centers,ColoredMetric,utilities};
use crate::types::{CenterIdx,Distance,PointCount};
use super::{Edge,EdgeIdx,flow,flow::State};

// note that edge_cursor points at the edge that has not been added yet
// (the edge at edge_curser-1 has been added already)
pub(super) fn settle<'a, M: ColoredMetric>(edge_cursor: EdgeIdx, bucket: &mut Vec<Edge<'a>>, i: CenterIdx, privacy_bound:PointCount, state: &mut State<'a>, centers: &Centers, space: &M) -> Clustering{

    let k = centers.m();
    let mut cursor = edge_cursor;

//    println!("  Edge_cursor: {}", edge_cursor);
//    println!("  Bucket to settle: {:?}\n", bucket.iter().map(|x| x.d).collect::<Vec<Distance>>());

    let edges_present : bool;
    // first clear or fill bucket:
    if edge_cursor >= bucket.len()/2 { // in this case fill the bucket
        while cursor < bucket.len() {
            flow::add_edge(bucket[cursor],i,k,privacy_bound,state);
//            println!("\tInitial adding: {:?} in binary search for center: {};\tmax_flow: {}", bucket[cursor], i, state.max_flow);
            cursor += 1;
        }
        edges_present = true;
    } else { // in this case empty bucket
        while cursor > 0 {
            cursor -= 1;
            flow::remove_edge(bucket[cursor],i,k,privacy_bound,state);
//            println!("\tInitial removing: {:?} in binary search for center: {};\tmax_flow: {}", bucket[cursor], i, state.max_flow);
        }
        edges_present = false;
    }

    let _radius = search_for_radius(edges_present, bucket, &mut cursor, i, k, privacy_bound, state);

    #[cfg(debug_assertions)]
    println!("\tRadius of partial assignment for C_{} found:\t{};\t(max_flow: {})", i, _radius, state.max_flow);
//    println!("center_of: {:?}", state.center_of);
//    println!("number_of_points_covered: {:?}", state.number_of_covered_points);

    let mut center_prefix = Centers::with_capacity(i+1);
    for c in centers.get_all(space).iter().take(i+1) {
        center_prefix.push(c);
    }

    let clustering = Clustering::new(center_prefix,state.center_of.clone(),space);
    #[cfg(debug_assertions)]
    assert_eq!(_radius,clustering.get_radius(), "Determined radius differs from cluster radius! This should never happen!");

//    println!("  Bucket after settling: {:?}\n", bucket.iter().map(|x| x.d).collect::<Vec<Distance>>());

    // empty bucket for the final time
    while cursor > 0 {
        cursor -= 1;
        flow::remove_edge(bucket[cursor],i,k,privacy_bound,state);
//        println!("\tClean up removing: {:?} in binary search for center: {};\tmaxflow: {}", bucket[cursor], i,state.max_flow);
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
fn search_for_radius<'a>(edges_present: bool, list: &mut Vec<Edge<'a>>, cursor : &mut EdgeIdx, i : CenterIdx, k : PointCount, privacy_bound: PointCount, state: &mut State<'a>) -> Distance {
//    println!("\n  List to settle: {:?} edges_present: {}", list.iter().map(|x| x.d).collect::<Vec<Distance>>(), edges_present);
    let list_len = list.len();
    assert!(list_len > 0, "Empty list in binary search");
    if list_len == 1 {
        if !edges_present {
            flow::add_edge(list[0],i,k,privacy_bound,state);
            *cursor += 1;
//            println!("\n\tAdd final edge: {:?} in binary search for center: {};\tmax_flow: {}", list[0], i, state.max_flow);
        }
        assert_eq!(state.max_flow, (i+1) * privacy_bound, "Something went wrong in the binary search.");
        return list[0].d;
    }

    let (mut smaller, mut bigger) = utilities::split_at_median(list);
//    println!("     smaller: {:?} bigger: {:?}\n", smaller.iter().map(|x| x.d).collect::<Vec<Distance>>(), bigger.iter().map(|x| x.d).collect::<Vec<Distance>>());



    // take care that smaller edges are added and bigger edges are not present in the flow network
    if edges_present {
        for e in bigger.iter().rev() {
            flow::remove_edge(*e, i, k, privacy_bound, state);
//            println!("\tTry removing: {:?} in binary search for center: {};\tmax_flow: {}", e, i,state.max_flow);
            *cursor -= 1;
        }
    } else {
        for e in smaller.iter() {
            flow::add_edge(*e, i, k, privacy_bound, state);
//            println!("\tTry adding: {:?} in binary search for center: {};\tmax_flow: {}", e, i,state.max_flow);
            *cursor += 1;
        }
    }


    let radius: Distance;
//    let mut left : Vec<Edge> = Vec::with_capacity(list_len);
//    let mut right : Vec<Edge> = Vec::with_capacity(list_len/2 + 1);
    if state.max_flow >= (i+1) * privacy_bound { // we need to settle in smaller
        radius = search_for_radius(true, &mut smaller, cursor, i, k, privacy_bound, state);
    } else { // we need to settle in bigger
        radius = search_for_radius(false, &mut bigger, cursor, i, k, privacy_bound, state);
    }

    // concatenate smaller and bigger;
    list.clear();
    list.append(&mut smaller);
    list.append(&mut bigger);
    radius
}


