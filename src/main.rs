extern crate ff_k_center;
use std::fs::create_dir_all;

use ff_k_center::{ClusteringProblem,compute_privacy_preserving_representative_k_center};
use ff_k_center::Space2D;

fn main() {


    match create_dir_all("output") {
        Err(error) => println!("Cannot create output-directory. Reason: {}", error),
        _ => {},
    }
    
    match create_dir_all("output/temp") {
        Err(error) => println!("Cannot create output-directory. Reason: {}", error),
        _ => {},
    }

    // load test Space2D from file:
    // let space_name = "example";
    // let space = Space2D::by_file(format!("{}.2dspace",space_name).as_str(),15);
    // let prob = ClusteringProblem{
        // k : 9, // number of center;
        // privacy_bound : 3, // number of points to represent;
        // rep_intervals : vec!((0,2),(1,2),(2,6),(3,3),(0,40),(0,500)), // representation interval [a,b] for each color class; for color classes without interval we subsitute [0. inf]
    // };
    
    let space_name = "shifting3";
    let space = Space2D::by_file(format!("{}.2dspace",space_name).as_str(),15);
    let prob = ClusteringProblem{
        k : 4, // number of center;
        privacy_bound : 7, // number of points to represent;
        rep_intervals : vec!((1,4),(0,0),(3,4)), // representation interval [a,b] for each color class; for color classes without interval we subsitute [0. inf]
    };
    
    let clustering = compute_privacy_preserving_representative_k_center(&space, &prob);
    clustering.save_to_file(format!("output/{}.clustering", space_name).as_str());


}
