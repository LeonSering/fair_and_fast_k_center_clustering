 use pyo3::prelude::*;
use pyo3::proc_macro::{pyfunction,pyclass,pymethods};
use pyo3::prelude::PyResult;
use pyo3::exceptions::PyAttributeError;

use crate::{ClusteringProblem,compute_privacy_preserving_representative_k_center};
use crate::clustering::Clustering;
use crate::types::{PointCount,Interval,CenterIdx,PointIdx,ColorIdx,Distance};
use crate::space::SpaceND;
use crate::assertions::assert_problem_parameters;


#[pyclass]
pub(crate) struct FFKCenter {
    // parameters
    clustering_problem: ClusteringProblem,

    // data
    space: Option<SpaceND>,

    // attributes
    clustering: Option<Clustering<'static>>,
}

#[pymethods]
impl FFKCenter{
    #[new]
    fn new(k: PointCount, privacy_bound : PointCount, rep_intervals : Vec<Interval>) -> PyResult<FFKCenter> {
        let problem = ClusteringProblem{
                k,
                privacy_bound,
                rep_intervals,
            };

        assert_problem_parameters(&problem);

        //TODO: Return python error if assertion fails
        
        // Python::with_gil(|py| Py::new(py,FFKCenter{clustering_problem: problem, space: None, clustering: None}).unwrap())
        Ok(FFKCenter{clustering_problem: problem, space: None, clustering: None})
    }

    /// Returns a list of point-indices. The point corresponding to centers[i] is hence,
    /// X[centers[i]]. All points assigned to this center have the label i.
    #[getter]
    fn get_centers(&self) -> PyResult<Vec<PointIdx>> {
        match &self.clustering {
            Some(clust) => Ok(clust.get_centers().iter().map(|c| c.idx()).collect()),
            None => Err(PyAttributeError::new_err("No centers computed yet!")) //TODO: Check whether correct Error Type
        }
    }


    /// Returns the center index for each point. If a point is not assignet use None instead.
    #[getter]
    fn get_labels(&self) -> PyResult<Vec<Option<CenterIdx>>> {
        match &self.clustering {
            Some(clust) => Ok(clust.get_assignment().iter().map(|x| *x).collect()),
            None => Err(PyAttributeError::new_err("No centers computed yet!")) //TODO: Check whether correct Error Type
        }
    }

    /// Returns the number of centers.
    #[getter]
    fn get_number_of_centers(&self) -> PyResult<PointCount> {
        match &self.clustering {
            Some(clust) => Ok(clust.get_centers().m()),
            None => Err(PyAttributeError::new_err("No centers computed yet!")) //TODO: Check whether correct Error Type
        }
    }

    /// Executes the algorithm.
    /// Input: 2d-Array. An array of Datapoints, which are an array of dimension-many floats
    /// (Value). TODO: Use 2d-Numpy Array here.
    /// A 1d-Array containing a color-label for each datapoint.
    ///
    /// First creates a MetricSpace. Then executes the algorithm. Output are saved into the model
    /// and can be accessed via model.centers or model.labels.
    fn fit(&mut self, data: Vec<Vec<Distance>>, colors: Vec<ColorIdx>) -> () {
        println!("Algorithm will be executed.");
        // First create a metric space from the data
        let space = SpaceND::by_ndpoints(data,colors);

        // Then execute the algorithm and save the output clustering into self.clustering
        let clustering = Some(compute_privacy_preserving_representative_k_center(&space, &self.clustering_problem));

    }

}
