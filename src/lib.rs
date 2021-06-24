//! Algorithm and utilities for computing a privacy preserving representative k-center in O(nk<sup>2</sup> + k<sup>4</sup>).
//!
//! # Quick Start
//!
//! First create a metric space with color classes by either loading them through a file or by
//! manually entering the distance matrix and colors.
//! The metric space must implement the [ColoredMetric]-trait.
//! For details on how to create or load metric spaces see: [SpaceMatrix] and [Space2D].
//!
//! Then create a [ClusteringProblem] by specifying the maximal number of clusters k, the lower bound
//! for privacy condition privacy_bound, and a list intervals (one for each color class) for the representaive conditon.
//!
//! With [compute_privacy_preserving_representative_k_center] a [Clustering] is created, which
//! contains a list of centers (see [Centers]) and an assignment of each point to a center.
//! This clustering is a 13-approximation on the optimal "privacy preserving representative k-center"-clustering.
//!


type PointCount = usize;
type ColorCount = u16;
//type Interval = (usize,usize);

/// ClusteringProblem defines for a given colored metric space a problem instance of privacy preserving representative k-clustering.
/// k is the maximal number of centers that can be opened;
/// privacy_bound specifies a lower bound on the number of clients that needs to be assigned to each center;
/// rep_interval contains an interval [a,b] for each color class, the number of centers of that
/// color must be within the interval.
pub struct ClusteringProblem {
    pub k : PointCount, // maximal number of centers
    pub privacy_bound : PointCount, // lower bound of representation L
    //pub rep_interval : Vec<Interval>, // one integer interval for each color class

    // method: test if valid: color classes is metric space are 0, ..., gamma-1; sum of a_j <= k.
}

mod space;
pub use space::{Space2D,SpaceMatrix,ColoredMetric,Point};

mod clustering;
pub use clustering::{Clustering,Centers};

mod phase1;
use phase1::gonzales_heuristic;

mod phase2;
use phase2::make_private;
use phase2::with_sorting::make_private_with_sorting; // TEMP

mod phase3;
use phase3::redistribute;

mod phase4;
use phase4::finalize;

struct OpeningList {
    eta : Vec<PointCount>
}

/// Computes a privacy preserving representative k-clustering.
/// Input: A metric space implementing the [ColoredMetric] trait and a [ClusteringProblem].
/// Output: A clustering of type [Clustering], that contains up to k centers (see [Centers]) and an assignment of
/// each client to a center.
/// The radius is a 13-approximation and the running time is O(nk<sup>2</sup> + k<sup>4</sup>). 
/// TODO: Update this doc-comment
pub fn compute_privacy_preserving_representative_k_center<'a, M : ColoredMetric>(space : &'a M, prob : &'a ClusteringProblem) -> Clustering<'a> {


    assert!(prob.k >= 1, "We have k = {}! There should be at least one center.", prob.k); // we want to allow at least 1 center
    assert!(space.n() >= prob.k, "we have n < k ({} < {})! We need more points than centers", space.n(), prob.k); // the number of points should not be less than the number of centers
    assert!(space.n() >= prob.k * prob.privacy_bound, "We have n < k * L ({} < {} * {})! We need enough points so that k center can satisfy the privacy condition.", space.n(), prob.k, prob.privacy_bound);
    
    /////////////////////////////////////////////////////////////////////////////////////
    // phase 1: use gonzales heuristic to obtain an ordered set of preliminary centers //
    /////////////////////////////////////////////////////////////////////////////////////

    let gonzales = gonzales_heuristic(space, prob.k);

    println!("** Phase 1: Determined k = {} centers by the Gonzales Heuristic: {:?}", prob.k, gonzales.iter().as_slice());
//    println!("index of center {}: {}", gonzales.centers[2], gonzales_index_by_center.get(&gonzales.centers[2]).expect(""));
    gonzales.save_to_file("test.centers");


    ///////////////////////////////////////
    // phase 2: determine privacy radius //
    ///////////////////////////////////////
    
    // TEMP:
    // clusterings is now a vector of partial clustering
    let clusterings_with_sorting : Vec<Clustering<'a>> = make_private_with_sorting(space, prob, &gonzales);

    let mut clusterings : Vec<Clustering<'a>> = make_private(space, prob, &gonzales);
    
    // TEMP:
    for (c, clustering) in clusterings.iter().enumerate() {
        let mut dist = 0.0f32;
        for p in space.point_iter() {
            if clustering.get_assignment()[p.idx()].is_some() {
                let current_d = space.dist(p,clustering.get_center(clustering.get_assignment()[p.idx()].unwrap())); 
                if current_d > dist {
                    dist = current_d;
                }
            }
        }
        println!("Clustering {}: measured radius: {}. written radius: {}, radius_with_sorting: {}", c, dist, clustering.get_radius(), clusterings_with_sorting[c].get_radius());
    }

    println!("** Phase 2: Determined k = {} radii: {:?}", prob.k, clusterings.iter().map(|clustering| clustering.get_radius()).collect::<Vec<f32>>());
    
    for i in 0..prob.k {
        let save_path = format!("output/after_phase_2_with_i_{}.clustering", i);
        clusterings[i].save_to_file(save_path.as_str());
    }

    ////////////////////////////////////////////////////////////////////
    // phase 3: redistribute assignment, s.t. sizes are multiple of L //
    ////////////////////////////////////////////////////////////////////
    
    let opening_lists = redistribute(space, prob, clusterings);



    //////////////////////////////////////////////////////////////////////
    // phase 4: open new centers and  determine actual set of centers C //
    //////////////////////////////////////////////////////////////////////

    let mut final_clusterings = finalize(space, &prob, opening_lists); 
    



    // TODO somehwere here we should pick the best clustering
    final_clusterings.pop().unwrap()
}
