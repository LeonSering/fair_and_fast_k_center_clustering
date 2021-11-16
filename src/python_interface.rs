use pyo3::proc_macro::{pyclass,pymethods};
use pyo3::prelude::{Python,PyResult};
use pyo3::types::PyDict;

use crate::{ClusteringProblem,compute_privacy_preserving_representative_k_center};
use crate::clustering::Clustering;
use crate::types::{PointCount,Interval,CenterIdx,PointIdx,ColorIdx,Distance};
use crate::space::{ColoredMetric,SpaceND};
use crate::assertions::assert_problem_parameters;
use pyo3::create_exception;

create_exception!(m, InvalidClusteringProblemError, pyo3::exceptions::PyException);
// create_exception!(m, ClusteringProblemMissingError, pyo3::exceptions::PyException);
// create_exception!(m, DataMissingError, pyo3::exceptions::PyException);
create_exception!(m, ClusteringMissingError, pyo3::exceptions::PyException);
create_exception!(m, DimensionMismatchError, pyo3::exceptions::PyException);

const NOPROB: &str = "No clustering problem defined yet. Run set_problem_parameters(k,privacy_bound,rep_intervals)";
const NOCLUSTERING: &str = "No clustering computed yet. Run fit(data,colors) or compute_clustering()";
const NOSPACE: &str = "No space defined yet. Run insert(data,colors) or fit(data,colors).";


use crate::phase2;
use crate::clustering::Centers;

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
        self.space.as_ref().expect(NOSPACE)
    }

    fn get_prob(&self) -> &ClusteringProblem {
        self.prob.as_ref().expect(NOCLUSTERING)
    }
}

#[pymethods]
impl FFKCenter{
    #[new]
    #[args(k,privacy_bound = "1", rep_intervals = "vec![]")]
    fn new(k: PointCount, privacy_bound : PointCount, rep_intervals : Vec<Interval>) -> PyResult<FFKCenter> {
        let problem = ClusteringProblem{
                k,
                privacy_bound,
                rep_intervals,
            };

        match assert_problem_parameters(&problem) {
            Err(msg) => Err(InvalidClusteringProblemError::new_err(msg)),
            _ => {
                Ok(FFKCenter{prob: Some(problem), space: None, clustering: None})
            }
        }
    }

    #[getter]
    fn get_k(&self) -> PointCount {
        self.prob.as_ref().expect(NOPROB).k
    }

    #[getter]
    fn get_privacy_bound(&self) -> PointCount {
        self.prob.as_ref().expect(NOPROB).privacy_bound
    }

    #[getter]
    fn get_rep_intervals(&self) -> Vec<Interval> {
        self.prob.as_ref().expect(NOPROB).rep_intervals.clone()
    }

    fn set_problem_parameters(&mut self, k: PointCount, privacy_bound: PointCount, rep_intervals: Vec<Interval>) -> PyResult<()>{
        let problem = ClusteringProblem{
                k,
                privacy_bound,
                rep_intervals,
            };
        match assert_problem_parameters(&problem) {
            Err(msg) => Err(InvalidClusteringProblemError::new_err(msg)),
            _ => {
                self.clustering = None;
                self.prob = Some(problem);
                Ok(())
            }
        }
    }

    #[setter]
    fn set_k(&mut self, k: PointCount) {
        self.prob.as_mut().expect(NOPROB).k = k;
        self.clustering = None;
    }

    #[setter]
    fn set_privacy_bound(&mut self, privacy_bound: PointCount) {
        self.prob.as_mut().expect(NOPROB).privacy_bound = privacy_bound;
        self.clustering = None;
    }

    #[setter]
    fn set_rep_intervals(&mut self, rep_intervals: Vec<Interval>) {
        self.prob.as_mut().expect(NOPROB).rep_intervals = rep_intervals;
        self.clustering = None;
    }


    #[getter]
    fn get_data(&self) -> PyResult<Vec<Vec<Distance>>> {
        Ok(self.space.as_ref().expect(NOSPACE).get_positions())
    }

    #[getter]
    fn get_colors(&self) -> PyResult<Vec<ColorIdx>> {
        Ok(self.space.as_ref().expect(NOSPACE).get_colors())
    }

    /// Input:
    /// Input: 2d-Array. An array of Datapoints, which are an array of dimension-many floats
    /// (Value).
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
    /// (Value).
    /// A 1d-Array containing a color-label for each datapoint.
    /// Optional: keyword-argument: verbose = 1 (0: silent, 1: brief, 2: verbose)
    ///
    /// First creates a MetricSpace. Then executes the algorithm. Output are saved into the model
    /// and can be accessed via model.centers or model.labels.
    #[args(data,colors,"*",verbose="1")]
    fn fit(&mut self, data: Vec<Vec<Distance>>, colors: Vec<ColorIdx>, verbose: u8) -> () {
        // First create a metric space from the data
        self.space = Some(SpaceND::by_ndpoints(data,colors));

        // Then execute the algorithm and save the output clustering into self.clustering
        self.clustering = Some(compute_privacy_preserving_representative_k_center(self.get_space(), &self.get_prob(), verbose));
    }

    /// Executes the algorithm. Space (created from datapoints and colors) must be set beforehand.
    /// Optional: keyword-argument: verbose = 1 (0: silent, 1: brief, 2: verbose)
    #[args("*",verbose="1")]
    fn compute_clustering(&mut self, verbose: u8) {
        self.clustering = Some(compute_privacy_preserving_representative_k_center(self.get_space(), &self.get_prob(), verbose));
    }

    /// Returns a list of point-indices. The point corresponding to centers[i] is hence,
    /// X[centers[i]]. All points assigned to this center have the label i.
    #[getter]
    fn get_centers(&self) -> PyResult<Vec<PointIdx>> {
        match &self.clustering {
            Some(clust) => Ok(clust.get_centers().get_all(self.get_space()).iter().map(|c| c.idx()).collect()),
            None => Err(ClusteringMissingError::new_err(NOCLUSTERING))
        }
    }


    /// Returns the cluster index (0,1,2..,m) for each point. If a point is not assignet use None instead.
    #[getter]
    fn get_cluster_labels(&self) -> PyResult<Vec<Option<CenterIdx>>> {
        match &self.clustering {
            Some(clust) => Ok(clust.get_assignment().iter().map(|x| *x).collect()),
            None => Err(ClusteringMissingError::new_err(NOCLUSTERING))
        }
    }

    /// Returns the center (in form of a point index) for each point. If a point is not assignet use None instead.
    #[getter]
    fn get_assignment(&self) -> PyResult<Vec<Option<PointIdx>>> {
        match &self.clustering {
            Some(clust) => Ok(clust.get_assignment().iter()
                              .map(|x| if x.is_none() {None} else {Some(self.get_centers().unwrap()[x.unwrap()])}).collect()),
            None => Err(ClusteringMissingError::new_err(NOCLUSTERING))
        }
    }

    /// Returns the number of centers.
    #[getter]
    fn get_number_of_centers(&self) -> PyResult<PointCount> {
        match &self.clustering {
            Some(clust) => Ok(clust.get_centers().m()),
            None => Err(ClusteringMissingError::new_err(NOCLUSTERING))
        }
    }

    /// Return the radius of the current cluster.
    #[getter]
    fn get_radius(&self) -> PyResult<Distance> {
        match &self.clustering {
            Some(clust) => Ok(clust.get_radius()),
            None => Err(ClusteringMissingError::new_err(NOCLUSTERING))
        }
    }

    /// Saves clusterig in txt-file. One line for each cluster.
    fn save_clustering_to_file(&self, file_path: &str) -> PyResult<()>{
        match &self.clustering {
            Some(clust) => { clust.save_to_file(file_path);
                Ok(())
            }
            None => Err(ClusteringMissingError::new_err(NOCLUSTERING))
        }
    }

    /// Load points and colours from txt-file.
    /// Input: file_path to txt-file in which each line represents one point: first the position as
    /// float, last entry is the color label as int. Values are separated by ',' (no comma at the
    /// end)
    /// Optional: expected can be set to the expected number of points to speed up the process.
    /// (Default: 1000).
    /// Optional: verbose = 1 (0: silent, 1: brief, 2: verbose)
    #[args(file_path, "*", expected = "1000", verbose = "1")]
    fn load_space_from_file(&mut self, file_path: &str, expected: PointCount, verbose: u8) -> PyResult<()> {
        self.space = Some(SpaceND::by_file(file_path, expected, verbose));
        self.clustering = None;
        Ok(())
    }

    /// Plot the clustering into the 2d-plane. Matplotlib must be installed.
    /// Optional: Specify the dimension for the x- and y-axes via x_dim and y_dim.
    /// Default is x_dim = 0, y_dim = 1.
    #[args("*", x_dim = "0", y_dim = "1")]
    fn plot2d(&self, x_dim: usize, y_dim: usize) -> PyResult<()> {
        if self.clustering.is_none() {
            return Err(ClusteringMissingError::new_err(NOCLUSTERING))
        }
        let dim = self.get_data().unwrap()[0].len();
        if x_dim > dim || y_dim > dim {
            return Err(DimensionMismatchError::new_err("Dimensions do not fit the data"))
        }
        Python::with_gil(|py| {
            let locals = PyDict::new(py);

            let c: Vec<ColorIdx> = self.get_colors().unwrap();
            let colormap = vec!["k","r","g","c","m","y","b","orange","lime","pink"];
            let colors:Vec<&str> = c.iter().map(|i| colormap[i % 10]).collect();

            let x: Vec<Distance> = self.get_data().unwrap().iter().map(|p| p[x_dim]).collect();
            let y: Vec<Distance> = self.get_data().unwrap().iter().map(|p| p[y_dim]).collect();

            let x_centers: Vec<Distance> = self.get_centers().unwrap().iter().map(|&c| x[c]).collect();
            let y_centers: Vec<Distance> = self.get_centers().unwrap().iter().map(|&c| y[c]).collect();
            let colors_centers: Vec<&str> = self.get_centers().unwrap().iter().map(|&c| colors[c]).collect();

            let cluster_labels = self.get_assignment().unwrap();

            locals.set_item("x", x).unwrap();
            locals.set_item("y", y).unwrap();
            locals.set_item("colors", colors).unwrap();

            locals.set_item("x_centers", x_centers).unwrap();
            locals.set_item("y_centers", y_centers).unwrap();
            locals.set_item("colors_centers", colors_centers).unwrap();

            locals.set_item("cluster_labels", cluster_labels).unwrap();
            py.run(
r#"import matplotlib.pyplot as plt
fig = plt.figure(figsize=(12,8))
plt.scatter(x, y, c = colors)
plt.axis('equal')
for p in range(0, len(cluster_labels)):
    plt.plot([x[p],x[cluster_labels[p]]], [y[p],y[cluster_labels[p]]], c = 'lightgrey', linewidth = 0.5, zorder = -1.0)
plt.scatter(x_centers,y_centers, c = colors_centers, marker = 'x', s = 300, zorder = 1.0)"#,
                None,Some(locals)).unwrap();
            });
        Ok(())
    }

    /// Computes a private assignment for the given list of centers (by point idx).
    /// Needs data to be assigned to the model and a clustering problem.
    /// Ignores the color data, as well as, k and rep_intervals and only uses data and
    /// privacy_bound.
    /// The assignment can be accessed by model.assignment and the radius via model.radius.
    /// Note that this does not compute a ff_k_center clustering.
    fn private_assignment_by_centers(&mut self, centers: Vec<PointIdx>) -> PyResult<()> {
        let k = centers.len();
        self.set_k(k);

        let space = self.get_space();
        let n = space.n();

        let prob = self.get_prob();
        if n < k * prob.privacy_bound {
            return Err(InvalidClusteringProblemError::new_err(
                    format!("Not enough points (n = {}) to satisfy the privacy_bound of {} for k = {} centers.", n, prob.privacy_bound, k)));
        }

        let centers = Centers::new(centers);

        self.clustering = Some(phase2::make_private(space, prob.privacy_bound, &centers).pop().unwrap());
        Ok(())
    }

}
