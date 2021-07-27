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



/// ClusteringProblem defines for a given colored metric space a problem instance of privacy preserving representative k-clustering.
/// k is the maximal number of centers that can be opened;
/// privacy_bound specifies a lower bound on the number of clients that needs to be assigned to each center;
/// rep_interval contains an interval [a,b] for each color class, the number of centers of that
/// color must be within the interval.
pub struct ClusteringProblem {
    pub k : PointCount, // maximal number of centers
    pub privacy_bound : PointCount, // lower bound of representation L
    pub rep_interval : Vec<Interval>, // one integer interval for each color class

    // method: test if valid: color classes is metric space are 0, ..., gamma-1; sum of a_j <= k.
}
use std::fmt;
impl fmt::Display for ClusteringProblem {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Clustering-Problem with k = {}, privacy_bound = {} and representative_intervals: ",self.k,self.privacy_bound)?;
        let mut iter = self.rep_interval.iter();
        if let Some(first_interval) = iter.next() {
            write!(f, "({},{})", first_interval.0, first_interval.1)?;
            for interval in iter {
                write!(f, ", ({},{})", interval.0, interval.1)?;
            }
        }
        write!(f, ".")
    }


}

mod types;
pub use types::{PointCount,ColorCount,Interval};
use types::Distance;

mod space;
pub use space::{Space2D,SpaceMatrix,ColoredMetric,Point};

mod clustering;
pub use clustering::{Clustering,Centers};

mod utilities;

mod phase1;
use phase1::gonzales_heuristic;

mod phase2;
use phase2::make_private;
//use phase2::with_sorting::make_private_with_sorting; // TEMP

mod phase3;
use phase3::redistribute;

mod phase4;
use phase4::finalize;

mod phase5;
use phase5::phase5;

#[derive(Debug)]
struct OpeningList {
    eta : Vec<PointCount>,
    forrest_radius : Distance
}

/// Computes a privacy preserving representative k-clustering.
/// Input: A metric space implementing the [ColoredMetric] trait and a [ClusteringProblem].
/// Output: A clustering of type [Clustering], that contains up to k centers (see [Centers]) and an assignment of
/// each client to a center.
/// The radius is a 13-approximation and the running time is O(nk<sup>2</sup> + k<sup>4</sup>). 
/// TODO: Update this doc-comment
pub fn compute_privacy_preserving_representative_k_center<'a, M : ColoredMetric>(space : &'a M, prob : &'a ClusteringProblem) -> Clustering<'a> {

    println!("\n**** Solving: {}", prob);

    /////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// Assertions on the ClusteringProblem //////////////////////
    /////////////////////////////////////////////////////////////////////////////////////

    assert!(prob.k >= 1, "We have k = {}! There should be at least one center.", prob.k); // we want to allow at least 1 center
    assert!(space.n() >= prob.k, "we have n < k ({} < {})! We need more points than centers", space.n(), prob.k); // the number of points should not be less than the number of centers
    assert!(space.n() >= prob.k * prob.privacy_bound, "We have n < k * L ({} < {} * {})! We need enough points so that k center can satisfy the privacy condition.", space.n(), prob.k, prob.privacy_bound);

    for (c, interval) in prob.rep_interval.iter().enumerate() {
        assert!(interval.0 <= interval.1, "The interval of color {} is ({}, {}). Lower bound cannot be bigger than the upper bound.", c, interval.0, interval.1);
    }

    let mut sum_of_a = 0;

    // check whether the sum of the lower bounds is not bigger than k
    for (a,_) in prob.rep_interval.iter() {
        sum_of_a += a;
    }
    assert!(sum_of_a <= prob.k, "The sum of the lower bounds of the representative intervals is {}, which is larger than k = {}.", sum_of_a, prob.k);

    let mut number_of_points_of_color: Vec<PointCount> = vec!(0; space.gamma());
    for p in space.point_iter() {
        number_of_points_of_color[space.color(p)] += 1;
    }

    let restricted_colors = if space.gamma() < prob.rep_interval.len() {space.gamma()} else {prob.rep_interval.len()}; // min{gamma, rep_interval.len()}

    // check if lower bound can be satisfied
    for c in 0..restricted_colors {
        assert!(number_of_points_of_color[c] >= prob.rep_interval[c].0, "There are {} points of color {}, but we require a = {} of the centers to be of this color.", number_of_points_of_color[c], c, prob.rep_interval[c].0);
    }
    for c in restricted_colors..prob.rep_interval.len() {
        assert_eq!(prob.rep_interval[c].0, 0, "We want {} centers of color {}, but there is not a single point of that color. We have gamma = {}.", prob.rep_interval[c].0, c, space.gamma());
    }

    // check sum of the upper bounds (min{ left side of interval, number of points }) is not
    // smaller than k 
    let mut sum_of_b = 0;
    for c in 0..restricted_colors {
        let upper_bound = if prob.rep_interval[c].1 < number_of_points_of_color[c] {prob.rep_interval[c].1} else {number_of_points_of_color[c]};
        sum_of_b += upper_bound;
    }
    for c in restricted_colors..space.gamma() {
        sum_of_b += number_of_points_of_color[c];
    }
    assert!(sum_of_b >= prob.k, "The sum of the upper bounds of the representative intervals (or the number of points of that color) is {}, which is smaller than k = {}.", sum_of_b, prob.k);


    
    /////////////////////////////////////////////////////////////////////////////////////
    // phase 1: use gonzales heuristic to obtain an ordered set of preliminary centers //
    /////////////////////////////////////////////////////////////////////////////////////

    let gonzales = gonzales_heuristic(space, prob.k);

    print!("\n**** Phase 1 done: Determined k = {} centers by the Gonzales heuristic: {}.\n", prob.k, gonzales);
//    println!("index of center {}: {}", gonzales.centers[2], gonzales_index_by_center.get(&gonzales.centers[2]).expect(""));
    gonzales.save_to_file("test.centers");


    ///////////////////////////////////////
    // phase 2: determine privacy radius //
    ///////////////////////////////////////
    
    let mut clusterings : Vec<Clustering<'a>> = make_private(space, prob, &gonzales);
    
    // TEMP:
    // clusterings is now a vector of partial clustering
    #[cfg(debug_assertions)]
    {
        let clusterings_with_sorting : Vec<Clustering<'a>> = phase2::with_sorting::make_private_with_sorting(space, prob, &gonzales);
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
            println!("\tClustering {}: measured radius: {}. written radius: {}, radius_with_sorting: {}", c, dist, clustering.get_radius(), clusterings_with_sorting[c].get_radius());
        }
    }

    println!("\n**** Phase 2 done: Determined k = {} radii: {:?}", prob.k, clusterings.iter().map(|clustering| clustering.get_radius()).collect::<Vec<f32>>());
    
    for i in 0..prob.k {
        let save_path = format!("output/after_phase_2_with_i_{}.clustering", i);
        clusterings[i].save_to_file(save_path.as_str());
    }

    ////////////////////////////////////////////////////////////////////
    // phase 3: redistribute assignment, s.t. sizes are multiple of L //
    ////////////////////////////////////////////////////////////////////
    
    let (spanning_trees, opening_lists) = redistribute(space, prob, &mut clusterings);

    
    println!("\n**** Phase 3 done: Determined the following opening lists:");
    #[cfg(debug_assertions)]
    for (i,openings) in opening_lists.iter().enumerate() {
        println!("  C_{}:\t{:?}", i, openings.iter().map(|open| &open.eta).collect::<Vec<_>>());
    }


    //////////////////////////////////////////////////////////////////////
    // phase 4: open new centers and  determine actual set of centers C //
    //////////////////////////////////////////////////////////////////////

    let (forrest_radii, final_centers) = finalize(space, &prob, opening_lists, &gonzales); 

    ////////////////////////////////////////////////////////////////////////////////////////////
    // phase 5: assign point to the final point of centers and determine the final clustering //
    ////////////////////////////////////////////////////////////////////////////////////////////
    
    let final_clustering = phase5(space, prob, final_centers, spanning_trees, forrest_radii);
    // TODO: Phase 5: Assign points to the new center such that the privacy containt is satisfied.
    // We do this cluster wise.
    //
    // for each clustering (after the shifting from phase 3) do the following:
    // for each point p of the cluster store dist(p, new_centers_of_cluster), i.e. the minimum distance of p
    // to one of the new centers of this cluster).
    // Heuristics: 
    // 1) Consider a batch of points of size max{L, n/k}, that has the maximal distance and assign
    // them to the nearest new center of the cluster.
    // 2) Repeat this, but also consider the number of spare node:
    // number_of_spare_nodes = number_of_point_in_cluster - number_of_new_centers_in_cluster * L
    // 3) Whenever we want to assign a point to a cluster that is already private, then decrease
    //    the number of spare_points by one.
    // 4) If the number of spare points is 0, then never assign to private new centers anymore but
    //    instead assign to the nearest non-private center.
    // 5) As soon as all new centers are private dont do buckets anymore but instead directly
    //    assign all points to their nearest new center

    final_clustering
    




}
