extern crate ff_k_center_lib;
use std::fs::create_dir_all;

use ff_k_center_lib::{ClusteringProblem,compute_privacy_preserving_representative_k_center};
use ff_k_center_lib::{ColoredMetric,new_space_by_2dpoints_file};

fn main() {
    let space_name = "mono_space";


    match create_dir_all("output") {
        Err(error) => println!("Cannot create output-directory. Reason: {}", error),
        _ => {},
    }

    // load test Space2D from file:
    let space : Box<dyn ColoredMetric> = new_space_by_2dpoints_file(format!("{}.2dspace",space_name).as_str(),15);
    let prob = ClusteringProblem{
        k : 6, // number of center;
        privacy_bound : 5, // number of points to represent;
        //gamma : space.gamma();
    };

    let clustering = compute_privacy_preserving_representative_k_center(&space, &prob);
    clustering.save_to_file(format!("output/{}.clustering", space_name).as_str());


}
