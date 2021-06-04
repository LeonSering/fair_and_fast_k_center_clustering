use crate::ClusteringProblem;
use crate::space::{Point,ColoredMetric};
use crate::clustering::{Clustering};

pub fn redistribute<'a>(space : &Box<dyn ColoredMetric>, prob : &ClusteringProblem, clustering : Clustering<'a>) -> Clustering<'a> {
    clustering
}

