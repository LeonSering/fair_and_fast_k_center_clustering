use pyo3::proc_macro::{pyclass,pymethods};
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
    prob: Option<ClusteringProblem>,

    // data
    space: Option<SpaceND>,

    // attributes
    clustering: Option<Clustering>,
}

impl FFKCenter{
    fn get_space(&self) -> &SpaceND {
        self.space.as_ref().expect("No space defined yet.")
    }
    
    fn get_prob(&self) -> &ClusteringProblem {
        self.prob.as_ref().expect("No clustering problem defined yet.")
    }
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
        Ok(FFKCenter{prob: Some(problem), space: None, clustering: None})
    }

    #[getter]
    fn get_k(&self) -> PointCount {
        self.prob.as_ref().expect("No clustering problem defined yet.").k
    }

    #[getter]
    fn get_privacy_bound(&self) -> PointCount {
        self.prob.as_ref().expect("No clustering problem defined yet.").privacy_bound
    }

    #[getter]
    fn get_rep_intervals(&self) -> Vec<Interval> {
        self.prob.as_ref().expect("No clustering problem defined yet.").rep_intervals.clone()
    }

    fn set_problem_parameters(&mut self, k: PointCount, privacy_bound: PointCount, rep_intervals: Vec<Interval>) {
        let problem = ClusteringProblem{
                k,
                privacy_bound,
                rep_intervals,
            };
        assert_problem_parameters(&problem);
        self.clustering = None;
        self.prob = Some(problem);
    }

    #[setter]
    fn set_k(&mut self, k: PointCount) {
        self.prob.as_mut().expect("No clustering problem defined yet.").k = k;
        self.clustering = None;
    }

    #[setter]
    fn set_privacy_bound(&mut self, privacy_bound: PointCount) {
        self.prob.as_mut().expect("No clustering problem defined yet.").privacy_bound = privacy_bound;
        self.clustering = None;
    }
    
    #[setter]
    fn set_rep_intervals(&mut self, rep_intervals: Vec<Interval>) {
        self.prob.as_mut().expect("No clustering problem defined yet.").rep_intervals = rep_intervals;
        self.clustering = None;
    }


    #[getter]
    fn get_data(&self) -> PyResult<Vec<Vec<Distance>>> {
        Ok(self.space.as_ref().expect("No data found.").get_positions())
    }

    #[getter]
    fn get_colors(&self) -> PyResult<Vec<ColorIdx>> {
        Ok(self.space.as_ref().expect("No data found.").get_colors())
    }
   
    /// Input:
    /// Input: 2d-Array. An array of Datapoints, which are an array of dimension-many floats
    /// (Value). TODO: Use 2d-Numpy Array here.
    /// A 1d-Array containing a color-label for each datapoint.
    ///
    /// Creates a MetricSpace with colors.
    /// Use execute() to compute a clustering.
    fn insert(&mut self, data: Vec<Vec<Distance>>, colors: Vec<ColorIdx>) -> () {
        self.space = Some(SpaceND::by_ndpoints(data,colors));
        self.clustering = None;
    }

    /// Executes the algorithm.
    /// Input: 2d-Array. An array of Datapoints, which are an array of dimension-many floats
    /// (Value). TODO: Use 2d-Numpy Array here.
    /// A 1d-Array containing a color-label for each datapoint.
    ///
    /// First creates a MetricSpace. Then executes the algorithm. Output are saved into the model
    /// and can be accessed via model.centers or model.labels.
    fn fit(&mut self, data: Vec<Vec<Distance>>, colors: Vec<ColorIdx>) -> () {
        // First create a metric space from the data
        self.space = Some(SpaceND::by_ndpoints(data,colors));

        // Then execute the algorithm and save the output clustering into self.clustering
        self.clustering = Some(compute_privacy_preserving_representative_k_center(self.get_space(), &self.get_prob()));
    }

    /// Executes the algorithm. Space (created from datapoints and colors) must be set beforehand.
    fn compute_clustering(&mut self) {
        self.clustering = Some(compute_privacy_preserving_representative_k_center(self.get_space(), &self.get_prob()));
    }

    /// Returns a list of point-indices. The point corresponding to centers[i] is hence,
    /// X[centers[i]]. All points assigned to this center have the label i.
    #[getter]
    fn get_centers(&self) -> PyResult<Vec<PointIdx>> {
        match &self.clustering {
            Some(clust) => Ok(clust.get_centers().get_all(self.get_space()).iter().map(|c| c.idx()).collect()),
            None => Err(PyAttributeError::new_err("No centers computed yet!")) //TODO: Check whether correct Error Type
        }
    }


    /// Returns the center index for each point. If a point is not assignet use None instead.
    #[getter]
    fn get_cluster_labels(&self) -> PyResult<Vec<Option<CenterIdx>>> {
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

    /// Return the radius of the current cluster.
    #[getter]
    fn get_radius(&self) -> PyResult<Distance> {
        match &self.clustering {
            Some(clust) => Ok(clust.get_radius()),
            None => Err(PyAttributeError::new_err("No centers computed yet!")) //TODO: Check whether correct Error Type
        }
    }
}
