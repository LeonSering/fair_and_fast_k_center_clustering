use crate::ClusteringProblem;
use crate::space::{Point,ColoredMetric};
use crate::clustering::Clustering;

/// Given a metric space and a clustering problem,
/// finalize takes a clustering in which each center (except for one) covers at a multple of L
/// points, and returns a clustering (which potentially more cluster as before) that also satisfy the representative constaints. 
///
pub fn finalize<'a>(space : &Box<dyn ColoredMetric>, prob : &ClusteringProblem, clustering : Clustering<'a>) -> Clustering<'a> {
    clustering
}

