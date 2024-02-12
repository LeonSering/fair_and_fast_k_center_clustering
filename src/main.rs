extern crate ff_k_center;
use ff_k_center::*;

fn main() {
    let n = 5_000;
    let k = 8;
    let space = SpaceND::new_random(n);
    let prob = ClusteringProblem {
        k,                                                 // number of center;
        privacy_bound: n / k,                              // number of points to represent;
        rep_intervals: vec![(0, 500), (0, 500), (0, 500)], // representation interval [a,b] for each color class; for color classes without interval we subsitute [0. inf]
    };

    let (clustering, total_time) =
        compute_privacy_preserving_representative_k_center(&space, &prob, None);
    println!(
        "Radius: {}; Running time: {}s.",
        clustering.get_radius(),
        total_time
    );
}
