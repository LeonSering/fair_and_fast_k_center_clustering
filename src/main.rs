extern crate ff_k_center_lib;

use ff_k_center_lib::{ClusteringProblem,compute_privacy_preserving_representative_k_center};
use ff_k_center_lib::{ColoredMetric,new_space_by_2dpoints_file};

fn main() {


    // load test Space2D from file:
    let space : Box<dyn ColoredMetric> = new_space_by_2dpoints_file("test.2dspace",15);
    let prob = ClusteringProblem{
        k : 5, // number of center;
        privacy_bound : 2, // number of points to represent;
        space : space,
        //gamma : space.gamma();
    };

    compute_privacy_preserving_representative_k_center(prob);

}
