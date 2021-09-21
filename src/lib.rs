//! Algorithm and utilities for computing a privacy preserving representative k-center in O(nk<sup>2</sup> + k<sup>4</sup>).
//!
//! # Quick Start
//!
//! First create a metric space with color classes by either loading them through a file or by
//! manually entering the distance matrix and colors.
//! The metric space must implement the [ColoredMetric]-trait.
//! For details on how to create or load metric spaces see: [SpaceMatrix] and [SpaceND].
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

use std::fmt;
use std::time::Instant;

extern crate num_cpus; // read the number of cpus and create that meany threads

#[derive(Clone)]
pub struct ClusteringProblem {
    pub k : PointCount, // maximal number of centers
    pub privacy_bound : PointCount, // lower bound of representation L
    pub rep_intervals : Vec<Interval>, // one integer interval for each color class

}
impl fmt::Display for ClusteringProblem {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Clustering-Problem with k = {}, privacy_bound = {} and representative_intervals: ",self.k,self.privacy_bound)?;
        let mut iter = self.rep_intervals.iter();
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

mod space;
pub use space::{SpaceND,SpaceMatrix,ColoredMetric,Point};

mod clustering;
pub use clustering::{Clustering,Centers};

mod utilities;

mod datastructures;

mod assertions;
pub use assertions::assert_clustering_problem;

mod phase1;
use phase1::gonzales_heuristic;

mod phase2;
use phase2::make_private;

mod phase3;
use phase3::redistribute;

mod phase4;
use phase4::phase4;

mod phase5;
use phase5::phase5;


mod python_interface;
use python_interface::FFKCenter;


use pyo3::prelude::*;
// use pyo3::prelude::{Python,PyModule,PyResult};
// use pyo3::proc_macro::pymodule;
// use pyo3::wrap_pyfunction;


#[pymodule]
#[pyo3(name = "ff_k_center")]
fn python_interface(_py: Python, m: &PyModule) -> PyResult<()> {
    // m.add_class::<FFKCenterModel>()?;
    m.add_class::<FFKCenter>()?;

    Ok(())
}


/// Computes a privacy preserving representative k-clustering.
/// The radius is a 13-approximation and the running time is O(nk<sup>2</sup> + k<sup>4</sup>).
///
/// # Inputs
/// * a metric space implementing the [ColoredMetric] trait;
/// * a [ClusteringProblem];
/// * a u8 indicating the length of the output. 0: silent, 1: brief, 2: verbose;
///
/// # Output
/// * a clustering of type [Clustering], that contains up to k centers (see [Centers]) and an assignment of
/// each client to a center.
/// TODO: Update this doc-comment
pub fn compute_privacy_preserving_representative_k_center<M : ColoredMetric>(space : &M, prob : &ClusteringProblem, verbose : u8) -> Clustering {
    let time_start = Instant::now();

    if verbose >= 2 {
        println!("\n**** Solving: {}", prob);
    }

    /////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////// Assertions on the ClusteringProblem //////////////////////
    /////////////////////////////////////////////////////////////////////////////////////

    match assert_clustering_problem(space, prob) {
        Err(msg) => {
            println!("\n\x1b[0;31mInvalidProblemFormulation:\x1b[0m {}", msg);
            panic!("Invalid problem. Abort!")}
        _ => {}
    }

    let time_after_assertions = Instant::now();
    if verbose >= 2 {
        println!("\n  - Assertions done (time: {:?}): ClusteringProblem seems well stated.", time_after_assertions.duration_since(time_start));
    }


    /////////////////////////////////////////////////////////////////////////////////////
    // phase 1: use gonzales heuristic to obtain an ordered set of preliminary centers //
    /////////////////////////////////////////////////////////////////////////////////////

    #[cfg(debug_assertions)]
    println!("\n**** Phase 1 ****");

    let gonzales: Centers = gonzales_heuristic(space, prob.k);

    let time_after_phase1 = Instant::now();
    if verbose >= 2{
        println!("\n  - Phase 1 done (time: {:?}): Determined k = {} centers by the Gonzales heuristic: ({}).", time_after_phase1.duration_since(time_after_assertions), prob.k, gonzales);
    } else if verbose == 1 {
        println!("  - Phase 1 done (time: {:?}): Determined k = {} centers by the Gonzales heuristic.", time_after_phase1.duration_since(time_after_assertions), prob.k);
    }

    ///////////////////////////////////////
    // phase 2: determine privacy radius //
    ///////////////////////////////////////

    #[cfg(debug_assertions)]
    println!("\n**** Phase 2 ****");
    let mut clusterings : Vec<Clustering> = make_private(space, prob, &gonzales);


    let time_after_phase2 = Instant::now();
    if verbose >= 2 {
        println!("\n  - Phase 2 done (time: {:?}): Determined k = {} clusterings:", time_after_phase2.duration_since(time_after_phase1), prob.k);
        for (i,clustering) in clusterings.iter().enumerate() {
            print!("\tC_{}:\tradius = {};\n", i, clustering.get_radius());
        }
    } else if verbose == 1 {
        println!("  - Phase 2 done (time: {:?}): Determined k = {} clusterings.", time_after_phase2.duration_since(time_after_phase1), prob.k);
    }

    #[cfg(debug_assertions)]
    {
        let clusterings_with_sorting : Vec<Clustering> = phase2::with_sorting::make_private_with_sorting(space, prob, &gonzales);

        let time_after_phase2_sorting = Instant::now();
        println!("\n  - Phase 2 (with sorting) done (time: {:?}): Determined k = {} clusterings:", time_after_phase2_sorting.duration_since(time_after_phase2), prob.k);
        for (i,clustering) in clusterings_with_sorting.iter().enumerate() {
            print!("\tC_{}:\tradius = {};\n", i, clustering.get_radius());
        }
    }

    ////////////////////////////////////////////////////////////////////
    // phase 3: redistribute assignment, s.t. sizes are multiple of L //
    ////////////////////////////////////////////////////////////////////
    let time_before_phase3 = Instant::now();

    #[cfg(debug_assertions)]
    println!("\n**** Phase 3 ****");

    let (spanning_trees, opening_lists) = redistribute(space, prob, &clusterings);


    let time_after_phase3 = Instant::now();
    if verbose >= 2 {
        println!("\n  - Phase 3 done (time: {:?}): Determined the following opening lists:", time_after_phase3.duration_since(time_before_phase3));
        for (i,openings) in opening_lists.iter().enumerate() {
            println!("\tC_{}:", i);
            for o in openings.iter() {
                println!("\t\t{};",o);
            }
        }
    } else if verbose == 1 {
        println!("  - Phase 3 done (time: {:?}): Determined opening lists.", time_after_phase3.duration_since(time_before_phase3));
    }


    //////////////////////////////////////////////////////////////////////
    // phase 4: open new centers and  determine actual set of centers C //
    //////////////////////////////////////////////////////////////////////

    #[cfg(debug_assertions)]
    println!("\n**** Phase 4 ****");


    let thread_count = num_cpus::get();
    let (new_centers, counts) = phase4(space, &prob, opening_lists, &gonzales,thread_count);

    let time_after_phase4 = Instant::now();
    let count_sum: usize = counts.iter().sum();
    if verbose >= 2 {
        println!("\n  - Phase 4 done (time: {:?} with {} threads on {} cores): Number of flow problems solved: {}. Determined the following new centers:", time_after_phase4.duration_since(time_after_phase3), thread_count, num_cpus::get(), count_sum);
        for (i,centers) in new_centers.iter().enumerate() {
            println!("\tC_{}:\t({}) \tassignment_radius: {};\tforrest_radius: {}; number of flow problems solved: {}", i, centers.as_points, centers.assignment_radius, centers.forrest_radius, counts[i]);
        }
    } else if verbose == 1 {
        println!("  - Phase 4 done (time: {:?} with {} threads on {} cores): Determined new centers. Number of flow problem solved: {}.", time_after_phase4.duration_since(time_after_phase3), thread_count, num_cpus::get(), count_sum);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    // phase 5: assign point to the final point of centers and determine the final clustering //
    ////////////////////////////////////////////////////////////////////////////////////////////

    #[cfg(debug_assertions)]
    println!("\n**** Phase 5 ****\n");

    let (best_i, final_centers, phase5_radius) = phase5(space, prob, new_centers, &mut clusterings, &spanning_trees);

    let time_after_phase5 = Instant::now();
    if verbose >= 2 {
        println!("\n  - Phase 5 done (time: {:?}): Created assignments and chose the final clustering (based on C_{}) with centers: ({}) and radius: {}.", time_after_phase5.duration_since(time_after_phase4), best_i, final_centers, phase5_radius);
    } else if verbose == 1 {
        println!("  - Phase 5 done (time: {:?}): Created assignments and chose the final clustering (based on C_{}) with radius: {}.", time_after_phase5.duration_since(time_after_phase4), best_i, phase5_radius);
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    //////////// rerun phase 2 on the final centers to obtain optimal clustering //////////
    ///////////////////////////////////////////////////////////////////////////////////////
    #[cfg(debug_assertions)]
    println!("\n**** Phase 2 Rerun ****\n");

    let final_clustering = make_private(space, &prob, &final_centers).pop().unwrap();


    let time_after_phase2rerun = Instant::now();
    if verbose >= 2 {
        println!("\n  - Rerun of phase 2 done (time: {:?}): final clustering determined with radius: {} (improvement to phase 5: {}).", time_after_phase2rerun.duration_since(time_after_phase5), final_clustering.get_radius(), final_clustering.get_radius() - phase5_radius);
    } else if verbose == 1 {
        println!("  - Rerun of phase 2 done (time: {:?}): final clustering determined with radius: {} (improvement to phase 5: {}).", time_after_phase2rerun.duration_since(time_after_phase5), final_clustering.get_radius(), final_clustering.get_radius() - phase5_radius);
    }

    ////////////////////////////////////////////////////////////////
    ///////////////////////// final print //////////////////////////
    ////////////////////////////////////////////////////////////////

    if verbose >= 2 {
        println!("\n**** Algorithm done (total time: {:?}): Privacy-preserving representative k-center computed. Number of centers: {}; final radius: {}.", time_after_phase2rerun.duration_since(time_start),final_clustering.get_centers().m(),final_clustering.get_radius());
        println!("\tTimes: phase 1: {:?}; phase 2: {:?}; phase 3: {:?}; phase 4: {:?}; phase 5: {:?}; phase 2 rerun: {:?};",
             time_after_phase1.duration_since(time_start),
             time_after_phase2.duration_since(time_after_phase1),
             time_after_phase3.duration_since(time_after_phase2),
             time_after_phase4.duration_since(time_after_phase3),
             time_after_phase5.duration_since(time_after_phase4),
             time_after_phase2rerun.duration_since(time_after_phase5));
    } else if verbose == 1 {
        println!("**** Algorithm done (total time: {:?}): Privacy-preserving representative k-center computed. Number of centers: {}; final radius: {}.", time_after_phase2rerun.duration_since(time_start),final_clustering.get_centers().m(),final_clustering.get_radius());
    }

    final_clustering

}
