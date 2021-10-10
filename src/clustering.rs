//////////////////////////////////////////////////////////////
//////////////////// module: clustering //////////////////////
//////////////////////////////////////////////////////////////


/// Contains two structs (both are only valid together with a space
///
/// centers: a vector of centers
/// clustering: centers, radius and assignment of points to centers
///
use crate::space::Point;
use crate::types::{PointCount,PointIdx,Distance,CenterIdx};

/// A list of centers.
#[derive(Debug,Clone)]
pub struct Centers{
    centers : Vec<PointIdx>,
}

use std::fs::File;
use std::io::prelude::*;
impl Centers{

    /// Creates a new list of centers from a vector of point indices.
    pub fn new(centers: Vec<PointIdx>) -> Centers {
        Centers{centers} 
    }

    /// Returns the number of centers m.
    pub fn m(&self) -> PointCount {
        self.centers.len()
    }

    /// Return the center of index i (from 0 to m-1)
    pub fn get<'a, M: ColoredMetric>(&self, i: CenterIdx, space: &'a M) -> &'a Point {
       space.get_point(self.centers[i]).expect("Error: Center is not a point of the current space.")
    }

    /// Adds a new center to the list.
    pub fn push(&mut self, c: &Point) {
        self.centers.push(c.idx());
    }

    /// Provides an iterator of the centers.
    pub fn get_all<'a, M: ColoredMetric>(&self, space: &'a M) -> Vec<&'a Point>{
        self.centers.iter().map(|&c| space.get_point(c).expect("Error: Some center is not a point of the current space.")).collect()
    }


    /// Provides an iterator of the centers indices.
    fn iter_idx(&self) -> std::slice::Iter<PointIdx>{
        self.centers.iter()
    }

    /// Save the centers to a file specified by file_path.
    /// The output file only contains one line containing the index of the centers separated by a
    /// comma.
    ///
    /// Example:
    ///
    /// ```txt
    /// 0,19,38,29,8,17
    /// ```
    pub fn save_to_file(&self, file_path : &str){
        let mut f = File::create(file_path).expect("Cannot open file for writing centers");
        let mut text = String::new();
        for c in self.centers.iter(){
            text = text + &format!("{},",c);
        }
        text.pop(); //delete last comma
        f.write_all(text.as_bytes()).expect("Could not write into centers-file");

    }

    /// Creates a new empty list of centers. The capacity is used to allocate enough storage on the
    /// heap.
    pub fn with_capacity(capacity : PointCount) -> Centers{
        Centers{centers : Vec::with_capacity(capacity)}
    }
}
use std::fmt;
impl<'a> fmt::Display for Centers {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut iter = self.centers.iter();
        if let Some(c) = iter.next() {
            write!(f, "{}", c)?;
            for c in iter {
                write!(f, ", {}", c)?;
            }
        }
        Ok(())
    }

}



/// A (partial) clustering specified by a list of centers and an assignment of clients to centers.
pub struct Clustering{
    /// a collection of points (of type &[Point])
    centers : Centers,
    /// for each point it contains the index of the center it is assigned to. For partial clustering this value can be None
    center_of : Vec<Option<CenterIdx>>, // returns center to which point with index x belongs, if not assigned to center its None
    /// cluster (list of points) by center idx; might be a super set of the actual points. Entries
    /// are only valid if it matches the information in center_of.
    cluster: Vec<Vec<PointIdx>>,
    /// the maximal distance between a client and its center
    /// (this is not computed but assigned; so there is some redundance)
    radius : Distance,
    /// By reassigning points the radius might get outdated (its always an upper bound). If this
    /// could be the case is indicated by the following boolean:
    radius_valid : bool,

    /// the number of points covered by each center
    cluster_sizes : Vec<PointCount>,
    // alternative ways to save: reference: list of Cluster typs (consisting of Center and
    // Vector of Points

    //space : &'a Box<dyn ColoredMetric>, // reference to metric space

    // method: cluster points of given color
    // method: print to file
}

use crate::ColoredMetric;
impl Clustering{

    /// Creates a new (partial) clustering for a given list of centers and an assignment of clients
    /// to center indices;
    /// Computes the minimum cover radius, hence, it needs a reference to the metric space;
    /// Computation takes O(n) time.
    ///
    /// #Panics
    ///
    /// Panics if the size of the assignment center_of does not match the number of points
    /// space.n()
    pub fn new<M : ColoredMetric>(centers: Centers, center_of: Vec<Option<CenterIdx>>, space: &M) -> Clustering {
        assert_eq!(center_of.len(),space.n(), "The assignment center_of must have an entry for each point");

        let mut cluster: Vec<Vec<PointIdx>> = vec!(Vec::new(); centers.m());
        let mut cluster_sizes: Vec<PointCount> = vec!(0; centers.m());
        let mut radius = <Distance>::MIN;
        for p in space.point_iter() {
            let covered_by = center_of[p.idx()];
            if let Some(c) = covered_by {
                cluster[c].push(p.idx());
                cluster_sizes[c] += 1;
                let dist = space.dist(p, centers.get(c, space));
                if dist > radius {
                    radius = dist;
                }
            }

        }

        Clustering{centers, center_of, cluster, radius, radius_valid: true, cluster_sizes}
    }

    // Creates a new empty clustering without any centers.
    // pub fn new_empty() -> Clustering {
        // Clustering {
            // centers: Centers::with_capacity(1),
            // center_of: Vec::new(),
            // cluster: Vec::new(),
            // radius: <Distance>::MAX,
            // radius_valid: false,
            // cluster_sizes: Vec::new()
        // }
    // }

    /// Returns the list of centers
    pub fn get_centers(&self) -> &Centers {
        &self.centers
    }

    /// Returns the (partial) assignment of the clients;
    /// Centers are represented by their index (0 to m) in the list of centers.
    pub fn get_assignment(&self) -> &Vec<Option<CenterIdx>> {
        &self.center_of
    }

    /// Return a vector of (cloned) point-refs containing the points of the cluster of the provided center index;
    pub fn get_cluster_of<'a, M: ColoredMetric>(&self, center_idx: CenterIdx, space: &'a M) -> Vec<&'a Point> {
        assert!(center_idx < self.cluster.len(), "center_idx must be between 0 and m-1");
        self.cluster[center_idx]
            .iter().filter(|&&p| self.center_of[p] == Some(center_idx))
            .map(|&p| space.get_point(p).expect("Some point is not part of the metric space anymore")).collect()
    }

    /// Returns the maximum distance between a client and the center it is assigned to;
    /// As this is pre-computed this takes O(1) time.
    pub fn get_radius(&self) -> Distance {
        if !self.radius_valid {
            println!("Radius is outdated and needs to be computed again! Please call update_radius(space).");
        }
        return self.radius
    }

    /// Updates an outdated radius. If the radius is not outdated nothing happens.
    pub fn update_radius<M: ColoredMetric>(&mut self, space: &M) {
        if !self.radius_valid {
            self.radius = Distance::MIN;
            for p in space.point_iter() {
                if self.center_of[p.idx()].is_some() {
                    let dist = space.dist(p, self.centers.get(self.center_of[p.idx()].unwrap(), space));
                    if dist > self.radius {
                        self.radius = dist;
                    }
                }
            }
            self.radius_valid = true;
        }
    }

    /// Returns the size of each cluster, i.e., the number of clients that are assigned to each
    /// center;
    /// As this is pre-computed this takes O(1) time.
    pub fn get_cluster_sizes(&self) -> &Vec<PointCount> {
        &self.cluster_sizes
    }

    /// Returns the number of centers/clusters m.
    pub fn m(&self) -> PointCount {
        self.centers.m()
    }

    /// Returns the center (as &Point) with index center_idx.
    pub fn get_center<'a, M: ColoredMetric>(&self, center_idx: CenterIdx, space: &'a M) -> &'a Point {
        self.centers.get(center_idx, space)
    }
    /// Assigns or reassigns a point to a new cluster given by the center_idx
    pub fn assign<M : ColoredMetric>(&mut self, p : &Point, center_idx: CenterIdx, space: &M) {
        if self.center_of[p.idx()].is_some() {
            let old_center = self.center_of[p.idx()].unwrap();
            if space.dist(p, self.centers.get(old_center,space)) >= self.radius {
                #[cfg(debug_assertions)]
                println!("WARNING: Point {} is already assigned to {} and is reassigned to {}; Radius is probably wrong (too large) now! Call update_radius()", p.idx(), old_center, center_idx);
                self.radius_valid = false;
            }
            self.cluster_sizes[old_center] -= 1;
            // Note that self.cluster[old_center] does not to be modified as it is checked whenever
            // it is read
        }
        self.center_of[p.idx()] = Some(center_idx);
        self.cluster[center_idx].push(p.idx());
        self.cluster_sizes[center_idx] += 1;
        let dist = space.dist(p, self.centers.get(center_idx,space));
        if dist >= self.radius {
            self.radius = dist;
            if !self.radius_valid {
                #[cfg(debug_assertions)]
                println!("Radius is valid again!");
                self.radius_valid = true;
            }
        }
    }

    /// Assignes non assigned clients to their nearest center;
    /// Takes O(nm) time.
    pub fn fill_up<'a, M : ColoredMetric>(&mut self, space : &'a M) {
        for p in space.point_iter() {
            if self.center_of[p.idx()].is_none() {
                // in this case p is not assigned to a center yet: Assigne it to nearest center
                let mut current_dist = <Distance>::MAX;
                let mut best_center = 0;
                for (j, &center) in self.centers.get_all(space).iter().enumerate() {
                    let dist = space.dist(center,p);
                    if dist < current_dist {
                        current_dist = dist;
                        best_center = j;
                    }
                }
                self.center_of[p.idx()] = Some(best_center);
                self.cluster[best_center].push(p.idx());
                self.cluster_sizes[best_center] += 1;
                if current_dist > self.radius {
                    self.radius = current_dist;
                }
            }

        }
}


    /// Saves the clustering to the text-file specified by file_path.
    /// The output is as follows:
    /// For each center we have one line in the output file of the form
    /// center: client_1, client_2, ...,
    ///
    /// Example:
    /// ```txt
    /// 0:0,1,2,3,4,
    /// 19:13,18,19,20,32,
    /// 38:33,34,37,38,39,
    /// 29:10,11,28,29,30,
    /// 8:6,8,24,26,27,
    /// 17:14,16,17,21,25,
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the file connot be open for writing.
    pub fn save_to_file(&self, file_path : &str){

        let mut clients_of : Vec<Vec<PointIdx>> = vec!(Vec::new(); self.centers.m());

        for (point_idx,&center_idx) in self.center_of.iter().enumerate() {
            if let Some(c) = center_idx {
                clients_of[c].push(point_idx);
            }
        }

        let mut f = File::create(file_path).expect("Cannot open file for writing clustering");
        let mut text = String::new();
        for (center_idx, center) in self.centers.iter_idx().enumerate(){
            text = text + &format!("{}:",center);
            for point_idx in clients_of.get(center_idx).unwrap().iter() {
                text = text + &format!("{},", point_idx);
            }
            text = text + "\n";
        }
        text.pop(); //delete last comma
        f.write_all(text.as_bytes()).expect("Could not write into centers-file");

    }
}


