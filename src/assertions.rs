use crate::{ColoredMetric,ClusteringProblem,OptionalParameters};
use crate::types::PointCount;

/// Asserts a clustering problem.
/// Checks whether the clustering problem is feasible. If this assertion runs through, the main
/// algorithm will return a feasible clustering.
///
/// # Inputs
/// * a metric space implementing the [ColoredMetric] trait;
/// * a [ClusteringProblem];
///
/// # Panics
/// Panics if one of the following asserts fails:
/// * k must be bigger than 1.
/// * n must be bigger than k.
/// * n must be bigger or equal to k * privacy_bound.
/// * for each representative interval the lower limit must be smaller or equal to the upper limit.
/// * the sum of the lower limits must be smaller or equal to k.
/// * the number of points of a specific color must be greater or equal to the lower limit of that
/// color.
/// * k must be smaller or equal to the sum of the upper limits of the representative intervals.
pub fn assert_clustering_problem<M : ColoredMetric>(space: &M, prob: &ClusteringProblem) -> Result<(),String> {

    assert_problem_parameters(prob)?; // first check, whether prob-parameter are consistent within itself (without the space)

    if !(space.n() >= prob.k) {return Err(format!("We have n < k ({} < {})! We need more points than centers.", space.n(), prob.k))}; // the number of points should not be less than the number of centers

    if !(space.n() >= prob.k * prob.privacy_bound) {return Err(format!("We have n < k * L ({} < {} * {})! We need enough points so that k center can satisfy the privacy condition.", space.n(), prob.k, prob.privacy_bound))};



    let mut number_of_points_of_color: Vec<PointCount> = vec!(0; space.gamma());
    for p in space.point_iter() {
        number_of_points_of_color[space.color(p)] += 1;
    }

    let restricted_colors = if space.gamma() < prob.rep_intervals.len() {space.gamma()} else {prob.rep_intervals.len()}; // min{gamma, rep_intervals.len()}

    // check if lower bound can be satisfied
    for c in 0..restricted_colors {
        if !(number_of_points_of_color[c] >= prob.rep_intervals[c].0) {return Err(format!("There are {} points of color {}, but we require a = {} of the centers to be of this color.", number_of_points_of_color[c], c, prob.rep_intervals[c].0))};
    }
    for c in restricted_colors..prob.rep_intervals.len() {
        if !(prob.rep_intervals[c].0 == 0) {return Err(format!("We want {} centers of color {}, but there is not a single point of that color. We have gamma = {}.", prob.rep_intervals[c].0, c, space.gamma()))};
    }

    // check sum of the upper bounds (min{ left side of interval, number of points }) is not
    // smaller than k
    let mut sum_of_b = 0;
    for c in 0..restricted_colors {
        let upper_bound = if prob.rep_intervals[c].1 < number_of_points_of_color[c] {prob.rep_intervals[c].1} else {number_of_points_of_color[c]};
        sum_of_b += upper_bound;
    }
    for c in restricted_colors..space.gamma() {
        sum_of_b += number_of_points_of_color[c];
    }
    if !(sum_of_b >= prob.k) {return Err(format!("The sum of the upper bounds of the representative intervals (or the number of points of that color) is {}, which is smaller than k = {}.", sum_of_b, prob.k))};

    Ok(())
}

/// Assertions that does not need the metric space.
pub(crate) fn assert_problem_parameters(prob: &ClusteringProblem) -> Result<(),String> {

    if !(prob.k >= 1) {return Err(format!("We have k = {}! There should be at least one center.", prob.k))}; // we want to allow at least 1 center

    for (c, interval) in prob.rep_intervals.iter().enumerate() {
        if !(interval.0 <= interval.1) {return Err(format!("The interval of color {} is ({}, {}). Lower bound cannot be bigger than the upper bound.", c, interval.0, interval.1))};
    }

    let mut sum_of_a = 0;

    // check whether the sum of the lower bounds is not bigger than k
    for (a,_) in prob.rep_intervals.iter() {
        sum_of_a += a;
    }
    if !(sum_of_a <= prob.k) {return Err(format!("The sum of the lower bounds of the representative intervals is {}, which is larger than k = {}.", sum_of_a, prob.k))};

    Ok(())
}

/// Assertion for the optional program parameters.
pub(crate) fn assert_optional_parameters(para: &OptionalParameters) -> Result<(),String> {
    match para.verbose{
        Some(i) => { if i > 2 {return Err(String::from("The verbose parameter must be 0 (silent), 1 (brief) or 2 (verbose)."));}}
        _ => {}
    };
    match para.thread_count{
        Some(i) => { if i <= 1 {return Err(String::from("There need to be at least two threads."));}}
        _ => {}
    };
    Ok(())
}
