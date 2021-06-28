extern crate ff_k_center_lib;
use std::fs::create_dir_all;

use ff_k_center_lib::{ClusteringProblem,compute_privacy_preserving_representative_k_center};
use ff_k_center_lib::Space2D;

fn main() {
    let space_name = "example";


    match create_dir_all("output") {
        Err(error) => println!("Cannot create output-directory. Reason: {}", error),
        _ => {},
    }

    // load test Space2D from file:
    let space = Space2D::by_file(format!("{}.2dspace",space_name).as_str(),15);
    let prob = ClusteringProblem{
        k : 6, // number of center;
        privacy_bound : 5, // number of points to represent;
        rep_interval : vec!((1,3); 4), // representation interval [a,b] for each color class; for color classes without interval we subsitute [0. inf]
    };

    let clustering = compute_privacy_preserving_representative_k_center(&space, &prob);
    clustering.save_to_file(format!("output/{}.clustering", space_name).as_str());


}
