extern crate ff_k_center;
use std::fs::create_dir_all;

use ff_k_center::{ClusteringProblem,OptionalParameters,compute_privacy_preserving_representative_k_center};
use ff_k_center::SpaceND;

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

    let verbose = 1; // 0: silent, 1:brief, 2:verbose

    let space_name = "classIris_crop";
    // let space = SpaceND::by_file(format!("{}.2dspace",space_name).as_str(),15,verbose);
    let n = 5_000_000;
    let k = 8;
    let space = SpaceND::new_random(n);
    let prob = ClusteringProblem{
        k, // number of center;
        privacy_bound : n/k, // number of points to represent;
        rep_intervals : vec!((0,500),(0,500),(0,500)), // representation interval [a,b] for each color class; for color classes without interval we subsitute [0. inf]
    };

    let opt_paras = OptionalParameters{
        verbose : Some(verbose),
        thread_count : None,
        phase_2_rerun : Some(false),
        phase_5_gonzalez : Some(true)
    };

    let (clustering, total_time) = compute_privacy_preserving_representative_k_center(&space, &prob, Some(opt_paras));
    clustering.save_to_file(format!("output/{}.clustering", space_name).as_str());

    println!("Running time: {}s.", total_time);

}
