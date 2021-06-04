pub struct ClusteringProblem {
    pub k : usize, // maximal number of centers
    pub privacy_bound : usize, // lower bound of representation L
    //pub gamma : u16, // number of colors
    //pub a : Vec<usize>, // lower bounds for color classes
    //pub b : Vec<usize>, // upper bounds for color classes

    // method: test if valid: color classes is metric space are 0, ..., gamma-1; sum of a_j <= k.
}

mod space;
pub use space::{Point,ColoredMetric,new_space_by_matrix,new_space_by_2dpoints_file,new_space_by_2dpoints,new_space_with_random_2dpoints}; // make it also public for external user (without the need to write ::space::)

mod clustering;
use clustering::Clustering;

mod phase1;
use phase1::gonzales_heuristic;
//mod buckets;
//use buckets::{put_into_buckets, assert_buckets_properties};

mod phase2;
use phase2::make_private;

mod phase3;
use phase3::redistribute;

pub fn compute_privacy_preserving_representative_k_center<'a>(space : &'a Box<dyn ColoredMetric>, prob : ClusteringProblem) -> Clustering::<'a> {


    assert!(prob.k >= 1); // we want to allow at least 1 center
    assert!(space.n() >= prob.k); // the number of points should not be less than the number of centers

    let gonzales = gonzales_heuristic(space, prob.k);

    println!("** Phase 1: Determined k = {} centers by the Gonzales Heuristic: {:?}", prob.k, gonzales.iter().as_slice());
//    println!("index of center {}: {}", gonzales.centers[2], gonzales_index_by_center.get(&gonzales.centers[2]).expect(""));
    gonzales.save_to_file("test.centers");


    ///////////////////////////////////////
    // phase 2: determine privacy radius //
    ///////////////////////////////////////

    let mut clustering = make_private(space, &prob, gonzales);
    // clustering is now a partial clustering

    ////////////////////////////////////////////////////////////////////
    // phase 3: redistribute assignment, s.t. sizes are multiple of L //
    ////////////////////////////////////////////////////////////////////
    clustering = redistribute(space, &prob, clustering);

    //////////////////////////////////////////////////////////////////////
    // phase 4: open new centers and  determine actual set of centers C //
    //////////////////////////////////////////////////////////////////////
    
    // TEMP: return empty clustering
    //use clustering::new_centers;
    


    clustering
}
