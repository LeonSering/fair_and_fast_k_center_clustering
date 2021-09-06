// use pyo3::prelude::*;
use pyo3::proc_macro::{pyfunction,pyclass,pymethods};
use pyo3::prelude::PyResult;

use crate::ClusteringProblem;
use crate::clustering::Clustering;
use crate::types::{PointCount,Interval,CenterIdx,PointIdx};
use crate::assertions::assert_problem_parameters;

type Value = f64;

#[pyclass]
pub(crate) struct FFKCenter {
    // parameters
    clustering_problem: ClusteringProblem, 
    // data
    //TODO: metric space?

    // attributes
    clustering: Clustering<'static>,
}

#[pymethods]
impl FFKCenter {
    #[new]
    fn new(k: PointCount, privacy_bound : PointCount, rep_intervals : Vec<Interval>) -> PyResult<FFKCenter> {
        let problem = ClusteringProblem{
                k,
                privacy_bound,
                rep_intervals,
            };

        assert_problem_parameters(&problem);

        //TODO: Return python error if assertion fails
        
        Ok(FFKCenter{clustering_problem: problem, clustering: Clustering::new_empty()})
    }

    /// Returns a list of point-indices. The point corresponding to centers[i] is hence,
    /// X[centers[i]]. All points assigned to this center have the label i.
    #[getter]
    fn get_centers(&self) -> Vec<PointIdx> {
        self.clustering.get_centers().iter().map(|c| c.idx()).collect()
    }


    /// Returns the center index for each point. If a point is not assignet use None instead.
    #[getter]
    fn get_labels(&self) -> Vec<Option<CenterIdx>> {
        self.clustering.get_assignment().iter().map(|x| *x).collect()
    }

    /// Returns the number of centers.
    #[getter]
    fn get_number_of_centers(&self) -> PointCount {
        self.clustering.get_centers().m()
    }

    /// Executes the algorithm.
    /// Input: 2d-Array. An array of Datapoints, which are an array of dimension-many floats
    /// (Value). TODO: Use 2d-Numpy Array here.
    ///
    /// First creates a MetricSpace. Then executes the algorithm. Output are saved into the model
    /// and can be accessed via model.centers or model.labels.
    fn fit(&self, data: Vec<Vec<Value>>) -> () {
        println!("Algorithm will be executed.");
        // First create a metric space from the data
        // Then execute the algorithm and save the output clustering into self.clustering

    }

}
