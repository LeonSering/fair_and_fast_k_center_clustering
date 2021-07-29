use crate::{ColoredMetric,ClusteringProblem};
use crate::types::PointCount;

/// Asserts a clustering problem.
/// Checks whether the clustering problem is feasible. If this assertions runs through the main
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
pub fn assert_clustering_problem<M : ColoredMetric>(space: &M, prob: &ClusteringProblem) {
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
}
