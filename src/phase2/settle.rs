use super::{flow, flow::State, Edge, EdgeIdx};
use crate::types::{CenterIdx, Distance, PointCount};
use crate::{utilities, Centers, Clustering, ColoredMetric};

// note that edge_cursor points at the edge that has not been added yet
// (the edge at edge_curser-1 has been added already)
pub(super) fn settle<'a, M: ColoredMetric>(
    edge_cursor: EdgeIdx,
    bucket: &mut [Edge<'a>],
    i: CenterIdx,
    privacy_bound: PointCount,
    state: &mut State<'a>,
    centers: &Centers,
    space: &M,
) -> Clustering {
    let k = centers.m();
    let mut cursor = edge_cursor;

    let edges_present =
    // first clear or fill bucket:
    if edge_cursor >= bucket.len()/2 { // in this case fill the bucket
        while cursor < bucket.len() {
            flow::add_edge(bucket[cursor],i,k,privacy_bound,state);
            cursor += 1;
        }
        true
    } else { // in this case empty bucket
        while cursor > 0 {
            cursor -= 1;
            flow::remove_edge(bucket[cursor],i,k,privacy_bound,state);
        }
        false
    };

    let _radius = search_for_radius(
        edges_present,
        bucket,
        &mut cursor,
        i,
        k,
        privacy_bound,
        state,
    );

    let mut center_prefix = Centers::with_capacity(i + 1);
    for c in centers.get_all(space).iter().take(i + 1) {
        center_prefix.push(c);
    }

    let clustering = Clustering::new(center_prefix, state.center_of.clone(), space);

    // empty bucket for the final time
    while cursor > 0 {
        cursor -= 1;
        flow::remove_edge(bucket[cursor], i, k, privacy_bound, state);
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
fn search_for_radius<'a>(
    edges_present: bool,
    list: &mut [Edge<'a>],
    cursor: &mut EdgeIdx,
    i: CenterIdx,
    k: PointCount,
    privacy_bound: PointCount,
    state: &mut State<'a>,
) -> Distance {
    let list_len = list.len();
    debug_assert!(list_len > 0, "Empty list in binary search");
    if list_len == 1 {
        if !edges_present {
            flow::add_edge(list[0], i, k, privacy_bound, state);
            *cursor += 1;
        }
        debug_assert_eq!(
            state.max_flow,
            (i + 1) * privacy_bound,
            "Something went wrong in the binary search."
        );
        return list[0].d;
    }

    let (smaller, bigger) = utilities::split_in_half(list);

    // take care that smaller edges are added and bigger edges are not present in the flow network
    if edges_present {
        for e in bigger.iter().rev() {
            flow::remove_edge(*e, i, k, privacy_bound, state);
            *cursor -= 1;
        }
    } else {
        for e in smaller.iter() {
            flow::add_edge(*e, i, k, privacy_bound, state);
            *cursor += 1;
        }
    }

    if state.max_flow >= (i + 1) * privacy_bound {
        // we need to settle in smaller
        search_for_radius(true, smaller, cursor, i, k, privacy_bound, state)
    } else {
        // we need to settle in bigger
        search_for_radius(false, bigger, cursor, i, k, privacy_bound, state)
    }
}
