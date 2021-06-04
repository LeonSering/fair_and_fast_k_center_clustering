use crate::ClusteringProblem;
use crate::space::{Point,ColoredMetric};
use crate::clustering::Clustering;

/// Given a metric space and a clusting problem,
/// redistribute takes a partial clustering that satisfies the privacy constrait
/// and returns a (full) clustering in which each center (except one) covers a multiple of L
/// points.
///
pub fn redistribute<'a>(space : &Box<dyn ColoredMetric>, prob : &ClusteringProblem, clustering : Clustering<'a>) -> Clustering<'a> {
    clustering
}

