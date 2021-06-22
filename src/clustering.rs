//////////////////////////////////////////////////////////////
//////////////////// module: clustering //////////////////////
//////////////////////////////////////////////////////////////


/// Contains two structs (both are only valid together with a space
/// 
/// centers: a vector of centers
/// clustering: centers, radius and assignment of points to centers
///
use crate::{PointCount,space::{Point,PointIdx,Distance}};
pub type CenterIdx = usize;

/// A list of centers.
pub struct Centers<'a>{
    centers : Vec<&'a Point>,
}

use std::fs::File;
use std::io::prelude::*;
impl<'a> Centers<'a>{
    /// Returns the number of centers m.
    pub fn m(&self) -> PointCount {
        self.centers.len()
    }

    /// Return the center of index i (from 0 to m-1)
    pub fn get(&self, i: CenterIdx) -> &'a Point {
       self.centers[i]  // TODO: Maybe cloning? 

    }
    
    /// Adds a new center to the list.
    pub fn push(&mut self, c: &'a Point) {
        self.centers.push(c);
    }
    
    /// Provides an iterator of the centers.
    pub fn iter(&self) -> std::slice::Iter<&'a Point>{
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
            text = text + &format!("{},",c.idx());
        }
        text.pop(); //delete last comma
        f.write_all(text.as_bytes()).expect("Could not write into centers-file");

    }
    
    /// Creates a new empty list of centers. The capacity is used to allocate enough storage on the
    /// heap.
    pub fn with_capacity(capacity : PointCount) -> Centers<'static>{
        Centers{centers : Vec::with_capacity(capacity)}
    }
}




/// A (partial) clustering specified by a list of centers and an assignment of clients to centers.
pub struct Clustering<'a>{
    /// a collection of points (of type &[Point])
    pub centers : Centers<'a>,
    /// for each point it contains the center it is assigned to. For partial clustering this value can be None
    pub center_of : Vec<Option<&'a Point>>, // returns center to which point with index x belongs, if not assigned to center its None
    /// the maximal distance between a client and its center 
    /// (this is not computed but assigned; so there is some redundance)
    pub radius : Option<Distance>,
    // alternative ways to save: reference: list of Cluster typs (consisting of Center and
    // Vector of Points

    //space : &'a Box<dyn ColoredMetric>, // reference to metric space

    // method: cluster points of given color
    // method: print to file
}

use std::collections::HashMap;
impl<'a> Clustering<'a>{

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

        let mut clients_of : HashMap<&Point, Vec<PointIdx>> = HashMap::with_capacity(self.centers.m());
        for center in self.centers.iter() {
            clients_of.insert(center,Vec::new());
        }

        for (point_idx,center) in self.center_of.iter().enumerate() {
            if let Some(c) = center {
                clients_of.get_mut(c).unwrap().push(point_idx);
            }
        }

        let mut f = File::create(file_path).expect("Cannot open file for writing clustering");
        let mut text = String::new();
        for c in self.centers.iter(){
            text = text + &format!("{}:",c.idx());
            for point_idx in clients_of.get(c).expect("Error").iter() {
                text = text + &format!("{},", point_idx);
            }
            text = text + "\n";
        }
        text.pop(); //delete last comma
        f.write_all(text.as_bytes()).expect("Could not write into centers-file");

    }
}


