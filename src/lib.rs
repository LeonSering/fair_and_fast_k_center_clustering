pub struct ClusteringProblem {
    pub k : usize, // maximal number of centers
    pub privacy_bound : usize, // lower bound of representation L
    pub space : Box<dyn ColoredMetric>,
    //pub gamma : u16, // number of colors
    //pub a : Vec<usize>, // lower bounds for color classes
    //pub b : Vec<usize>, // upper bounds for color classes

    // method: test if valid: color classes is metric space are 0, ..., gamma-1; sum of a_j <= k.
}



mod space;
pub use space::{ColoredMetric,new_space_by_matrix,new_space_by_2dpoints_file,new_space_by_2dpoints}; // make it also public for external user (without the need to write ::space::)

mod clustering;
use clustering::{Centers,Clustering};

mod gonzales;
use gonzales::gonzales_heuristic;

mod buckets;
use buckets::{put_into_buckets, assert_buckets_properties};

mod make_private;
use make_private::make_private;

pub fn compute_privacy_preserving_representative_k_center(prob : ClusteringProblem) {//-> Clustering {


    assert!(prob.k >= 1); // we want to allow at least 1 center
    assert!(prob.space.n() >= prob.k); // the number of points should not be less than the number of centers

    let gonzales = gonzales_heuristic(&prob.space, prob.k);

    println!("** Phase 1: Determined k = {} centers by the Gonzales Heuristic: {:?}", prob.k, gonzales.iter().collect::<Vec<&usize>>());
//    println!("index of center {}: {}", gonzales.centers[2], gonzales_index_by_center.get(&gonzales.centers[2]).expect(""));
    gonzales.save_to_file("test.centers");


    ///////////////////////////////////////
    // phase 2: determine privacy radius //
    ///////////////////////////////////////

    make_private(prob, gonzales);

    //////////////////////////////////////////////////////////////
    // phase 3: modify assignment, s.t. sizes are multiple of L //
    //////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////
    // phase 4: open new centers and  determine actual set of centers C //
    //////////////////////////////////////////////////////////////////////


}
