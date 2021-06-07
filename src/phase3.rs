use crate::ClusteringProblem;
use crate::space::ColoredMetric;
use crate::clustering::Clustering;

/// Given a metric space and a clusting problem,
/// redistribute takes a vector of partial clusterings, each satisfying the privacy constrait
/// and returns a vector of (full) clustering in which each center (except one) covers a multiple of L
/// points.
///
pub fn redistribute<'a>(space : &Box<dyn ColoredMetric>, prob : &ClusteringProblem, clusterings : Vec<Clustering<'a>>) -> Vec<Clustering<'a>> {
    clusterings
}

