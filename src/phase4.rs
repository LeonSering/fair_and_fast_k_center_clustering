use crate::ClusteringProblem;
use crate::space::ColoredMetric;
use crate::clustering::Clustering;

/// Given a metric space and a clustering problem,
/// finalize takes a vector of clusterings in which each center (except for one) covers a multple of L
/// points, and returns a vector of clusterings (which potentially more cluster as before) that also satisfy the representative constaints. 
///
pub fn finalize<'a, M : ColoredMetric>(space : &M, prob : &ClusteringProblem, clusterings : Vec<Clustering<'a>>) -> Vec<Clustering<'a>> {
    clusterings
}

