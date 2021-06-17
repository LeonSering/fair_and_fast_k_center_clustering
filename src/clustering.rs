//////////////////////////////////////////////////////////////
//////////////////// module: clustering //////////////////////
//////////////////////////////////////////////////////////////


/* Contains two structs (both are only valid together with a space
 * 
 * centers: a vector of centers
 * clustering: centers, radius and assignment of points to centers
 *
 */
use crate::{PointCount,space::{Point,PointIdx,Distance}};
pub type CenterIdx = usize;

pub struct Centers<'a>{
    centers : Vec<&'a Point>,
}

use std::fs::File;
use std::io::prelude::*;
impl<'a> Centers<'a>{
    pub fn m(&self) -> PointCount {
        self.centers.len()
    }

    pub fn get(&self, i: CenterIdx) -> &'a Point {
       self.centers[i]  // TODO: Maybe cloning? 

    }

    pub fn push(&mut self, c: &'a Point) {
        self.centers.push(c);
    }

    pub fn iter(&self) -> std::slice::Iter<&'a Point>{
        self.centers.iter()
    }


    // method: save to file
    pub fn save_to_file(&self, file_path : &str){
        let mut f = File::create(file_path).expect("Cannot open file for writing centers");
        let mut text = String::new();
        for c in self.centers.iter(){
            text = text + &format!("{},",c.idx());
        }
        text.pop(); //delete last comma
        f.write_all(text.as_bytes()).expect("Could not write into centers-file");

    }

    pub fn with_capacity(capacity : PointCount) -> Centers<'static>{
        Centers{centers : Vec::with_capacity(capacity)}
    }
}





pub struct Clustering<'a>{
    pub centers : Centers<'a>,
    pub radius : Distance,
    pub center_of : Vec<Option<&'a Point>>, // returns center to which point with index x belongs, if not assigned to center its None
    // alternative ways to save: reference: list of Cluster typs (consisting of Center and
    // Vector of Points

    //space : &'a Box<dyn ColoredMetric>, // reference to metric space

    // method: cluster points of given color
    // method: print to file
}

use std::collections::HashMap;
impl<'a> Clustering<'a>{
    pub fn save_to_file(&self, file_path : &str){

        let mut clients_of : HashMap<&Point, Vec<PointIdx>> = HashMap::with_capacity(self.centers.m());
        for center in self.centers.iter() {
            clients_of.insert(center,Vec::new());
        }

        for (point_idx,center) in self.center_of.iter().enumerate() {
            match center {
                Some(c) => clients_of.get_mut(c).expect("Error").push(point_idx),
                None => (),
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


