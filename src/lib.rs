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
//! for privacy condition privacy_bound, and a list intervals (one for each color class) for the representative condition.
//!
//! With [compute_privacy_preserving_representative_k_center] a [Clustering] is created, which
//! contains a list of centers (see [Centers]) and an assignment of each point to a center.
//! This clustering is a 13-approximation on the optimal "privacy preserving representative k-center"-clustering.



use std::fmt;
use std::time;

extern crate num_cpus; // read the number of cpus and create that meany threads

mod types;
pub use types::{PointCount,ColorCount,Interval,DurationInSec};

mod space;
pub use space::{SpaceND,SpaceMatrix,ColoredMetric,Point};

mod clustering;
pub use clustering::{Clustering,Centers};

mod utilities;

mod datastructures;

mod assertions;
pub use assertions::assert_clustering_problem;
use assertions::assert_optional_parameters;

mod phase1;
use phase1::gonzalez_heuristic;

mod phase2;
use phase2::make_private;

mod phase3;
use phase3::algebraic_pushing;

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

/// ClusteringProblem defines for a given colored metric space a problem instance of privacy preserving representative k-clustering.
/// k is the maximal number of centers that can be opened;
/// privacy_bound specifies a lower bound on the number of clients that needs to be assigned to each center;
/// rep_interval contains an interval [a,b] for each color class, the number of centers of that
/// color must be within the interval.
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

/// A collection of optional parameters. Each can be set to None to use the default values.
/// # Optional Parameters:
/// * verbose: Option<u8> specifying the verbosity of the command line output. 0: sieltn, 1: brief,
/// 2: verbose (default: 1);
/// * thread_count: Option<u8> specifying the number of threads that are used in phase 4 and 5.
/// (default: #cores);
/// * phase_2_rerun: Option<bool> determining whether an additonal phase 2 is run at the very end
/// using the final centers in order to obtain the best privacy-preserving assignment. (default: true)
/// * phase_5_gonzalez: Option<bool> determining whether a colored-based gonzalez is run during
/// phase 5. This heuristics causes the final centers to be more spread within the same cluster.
pub struct OptionalParameters {
    pub verbose: Option<u8>,
    pub thread_count: Option<usize>,
    pub phase_2_rerun: Option<bool>,
    pub phase_5_gonzalez: Option<bool>
}

const DEFAULT_VERBOSE: u8 = 1;
const DEFAULT_PHASE2_RERUN: bool = true;
const DEFAULT_PHASE5_GONZALEZ: bool = true;


/// Computes a privacy preserving representative k-clustering.
/// The radius is a 13-approximation and the running time is O(nk<sup>2</sup> + k<sup>4</sup>).
///
/// # Inputs
/// * a metric space implementing the [ColoredMetric] trait;
/// * a [ClusteringProblem];
/// * a optional struct of optional parameters, containing verbose: Option<u8>, thread_count: Option<usize,
/// phase_2_rerun: Option<bool>, phase_5_gonzalez: Option<bool>. For None default values are filled in (see [OptionalParamters]).
/// * a u8 indicating the verbosity of the command line output. 0: silent, 1: brief, 2: verbose;
///
/// # Output
/// * a clustering of type [Clustering], that contains up to k centers (see [Centers]) and an assignment of
/// each client to a center.
/// TODO: Update this doc-comment
pub fn compute_privacy_preserving_representative_k_center<M : ColoredMetric + std::marker::Sync> (space : &M, prob : &ClusteringProblem, options: Option<OptionalParameters>) -> (Clustering, DurationInSec) {
    let time_start = time::Instant::now();


    let cpu_count = num_cpus::get();

    let default_thread_count = cpu_count;

    let optional_parameters = options.unwrap_or(OptionalParameters{
        verbose: Some(DEFAULT_VERBOSE),
        thread_count: Some(default_thread_count),
        phase_2_rerun: Some(DEFAULT_PHASE2_RERUN),
        phase_5_gonzalez: Some(DEFAULT_PHASE5_GONZALEZ)});

    match assert_optional_parameters(&optional_parameters) {
        Err(msg) => {
            println!("\n\x1b[0;31mInvalidOptionalParameters:\x1b[0m {}", msg);
            panic!("Invalid parameters. Abort!")}
        _ => {}
    }

    let verbose = optional_parameters.verbose.unwrap_or(DEFAULT_VERBOSE);
    let thread_count = optional_parameters.thread_count.unwrap_or(default_thread_count);
    let phase_2_rerun = optional_parameters.phase_2_rerun.unwrap_or(DEFAULT_PHASE2_RERUN);
    let phase_5_gonzalez = optional_parameters.phase_5_gonzalez.unwrap_or(DEFAULT_PHASE5_GONZALEZ);

    if verbose >= 1 {
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

    let time_after_assertions = time::Instant::now();
    if verbose >= 2 {
        println!("\n  - Assertions done (time: {:?}): ClusteringProblem seems well stated.", time_after_assertions.duration_since(time_start));
    }


    /////////////////////////////////////////////////////////////////////////////////////
    // phase 1: use gonzalez heuristic to obtain an ordered set of preliminary centers //
    /////////////////////////////////////////////////////////////////////////////////////

    #[cfg(debug_assertions)]
    println!("\n**** Phase 1 ****");

    let gonzalez: Centers = gonzalez_heuristic(space, prob.k);

    let time_after_phase1 = time::Instant::now();
    if verbose >= 2{
        println!("\n  - Phase 1 done (time: {:?}): Determined k = {} centers by the Gonzalez heuristic: ({}).", time_after_phase1.duration_since(time_after_assertions), prob.k, gonzalez);
    } else if verbose == 1 {
        println!("  - Phase 1 done (time: {:?}): Determined k = {} centers by the Gonzalez heuristic.\n", time_after_phase1.duration_since(time_after_assertions), prob.k);
    }

    ///////////////////////////////////////
    // phase 2: determine privacy radius //
    ///////////////////////////////////////

    #[cfg(debug_assertions)]
    println!("\n**** Phase 2 ****");
    let clusterings : Vec<Clustering> = make_private(space, prob.privacy_bound, &gonzalez);


    let time_after_phase2 = time::Instant::now();
    if verbose >= 2 {
        println!("\n  - Phase 2 done (time: {:?}): Determined k = {} clusterings:", time_after_phase2.duration_since(time_after_phase1), prob.k);
        for (i,clustering) in clusterings.iter().enumerate() {
            print!("\tC_{}:\tradius = {};\n", i, clustering.get_radius());
        }
    } else if verbose == 1 {
        println!("  - Phase 2 done (time: {:?}): Determined k = {} clusterings.\n", time_after_phase2.duration_since(time_after_phase1), prob.k);
    }
    
    let clusterings_with_sorting : Vec<Clustering> = phase2::with_sorting::make_private_with_sorting(space, prob.privacy_bound, &gonzalez);

    let time_after_phase2_sorting = time::Instant::now();
    println!("  - Phase 2 (with sorting) done (time: {:?})\n", time_after_phase2_sorting.duration_since(time_after_phase2));

    ///////////////////////////////////////////////////////////////////////////
    // phase 3: determine spanning trees and eta-vector by algebraic_pushing //
    ///////////////////////////////////////////////////////////////////////////
    let time_before_phase3 = time::Instant::now();

    #[cfg(debug_assertions)]
    println!("\n**** Phase 3 ****");

    let (spanning_trees, opening_lists) = algebraic_pushing(space, prob, &clusterings);


    let time_after_phase3 = time::Instant::now();
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
    // phase 4: open new centers and determine actual set of centers C //
    //////////////////////////////////////////////////////////////////////

    #[cfg(debug_assertions)]
    println!("\n**** Phase 4 ****");


    let (new_centers, counts) = phase4(space, prob, opening_lists, &gonzalez,thread_count);

    let time_after_phase4 = time::Instant::now();
    let count_sum: usize = counts.iter().sum();
    if verbose >= 2 {
        println!("\n  - Phase 4 done (time: {:?} with {} threads on {} cores): Number of flow problems solved: {}. Determined the following new centers:", time_after_phase4.duration_since(time_after_phase3), thread_count, cpu_count, count_sum);
        for (i,centers) in new_centers.iter().enumerate() {
            println!("\tC_{}:\t({}) \tassignment_radius: {};\tforrest_radius: {}; number of flow problems solved: {}", i, centers.as_points, centers.assignment_radius, centers.forrest_radius, counts[i]);
        }
    } else if verbose == 1 {
        println!("  - Phase 4 done (time: {:?} with {} threads on {} cores): Determined new centers. Number of flow problem solved: {}.", time_after_phase4.duration_since(time_after_phase3), thread_count, cpu_count, count_sum);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    // phase 5: assign point to the final point of centers and determine the final clustering //
    ////////////////////////////////////////////////////////////////////////////////////////////

    #[cfg(debug_assertions)]
    println!("\n**** Phase 5 ****\n");

    let (best_i, mut final_clustering, phase5_radius) = phase5(space, prob, new_centers, clusterings, spanning_trees, thread_count, phase_5_gonzalez);

    let time_after_phase5 = time::Instant::now();
    if verbose >= 2 {
        println!("\n  - Phase 5 done (time: {:?} with {} threads on {} cores): Created assignments and chose the final clustering (based on C_{}) with centers: ({}) and radius: {}.", time_after_phase5.duration_since(time_after_phase4), thread_count, cpu_count, best_i, final_clustering.get_centers(), phase5_radius);
    } else if verbose == 1 {
        println!("  - Phase 5 done (time: {:?} with {} threads on {} cores): Created assignments and chose the final clustering (based on C_{}) with radius: {}.", time_after_phase5.duration_since(time_after_phase4), thread_count, cpu_count, best_i, phase5_radius);
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    //////////// rerun phase 2 on the final centers to obtain optimal clustering //////////
    ///////////////////////////////////////////////////////////////////////////////////////
    // let mut final_clustering = clusterings.into_iter().nth(best_i).unwrap();
    let mut phase2_rerun_time_str = String::from("-");
    if phase_2_rerun {

        #[cfg(debug_assertions)]
        println!("\n**** Phase 2 Rerun ****\n");

        final_clustering = make_private(space, prob.privacy_bound, &final_clustering.get_centers()).pop().unwrap();


        let time_after_phase2rerun = time::Instant::now();
        if verbose >= 2 {
            println!("\n  - Rerun of phase 2 done (time: {:?}): final clustering determined with radius: {} (improvement to phase 5: {}).", time_after_phase2rerun.duration_since(time_after_phase5), final_clustering.get_radius(), final_clustering.get_radius() - phase5_radius);
        } else if verbose == 1 {
            println!("  - Rerun of phase 2 done (time: {:?}): final clustering determined with radius: {} (improvement to phase 5: {}).", time_after_phase2rerun.duration_since(time_after_phase5), final_clustering.get_radius(), phase5_radius - final_clustering.get_radius());
        }

        phase2_rerun_time_str = format!("{:?}", time_after_phase2rerun.duration_since(time_after_phase5));

    }

    ////////////////////////////////////////////////////////////////
    ///////////////////////// final print //////////////////////////
    ////////////////////////////////////////////////////////////////

    let finishing_time = time::Instant::now();
    let total_time = finishing_time.duration_since(time_start);

    if verbose >= 2 {
        println!("\n**** Algorithm done (total time: {:?}): Privacy-preserving representative k-center computed. Number of centers: {}; final radius: {}.", total_time,final_clustering.get_centers().m(),final_clustering.get_radius());
        println!("\tTimes: phase 1: {:?}; phase 2: {:?}; phase 3: {:?}; phase 4: {:?}; phase 5: {:?}; phase 2 rerun: {};",
             time_after_phase1.duration_since(time_start),
             time_after_phase2.duration_since(time_after_phase1),
             time_after_phase3.duration_since(time_after_phase2),
             time_after_phase4.duration_since(time_after_phase3),
             time_after_phase5.duration_since(time_after_phase4),
             phase2_rerun_time_str);
    } else if verbose == 1 {
        println!("**** Algorithm done (total time: {:?}): Privacy-preserving representative k-center computed. Number of centers: {}; final radius: {}.", total_time,final_clustering.get_centers().m(),final_clustering.get_radius());
    }

    (final_clustering, total_time.as_secs_f64())

}
