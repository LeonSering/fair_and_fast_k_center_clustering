extern crate ff_k_center_lib;

use ff_k_center_lib::{ClusteringProblem,compute_privacy_preserving_representative_k_center};
use ff_k_center_lib::{ColoredMetric,new_space_by_2dpoints_file};

fn main() {
    // load test Space2D from file:
    let space : Box<dyn ColoredMetric> = new_space_by_2dpoints_file("mono_space.2dspace",15);
    let prob = ClusteringProblem{
        k : 6, // number of center;
        privacy_bound : 6, // number of points to represent;
        //gamma : space.gamma();
    };

    let clustering = compute_privacy_preserving_representative_k_center(&space, &prob);
    clustering.save_to_file("output.clustering");

}
