//////////////////////////////////////////////////////////////
//////////////////// module: clustering //////////////////////
//////////////////////////////////////////////////////////////


/* Contains two structs (both are only valid together with a space
 * 
 * centers: a vector of centers
 * clustering: centers, radius and assignment of points to centers
 *
 */
pub use crate::space::ColoredMetric;

pub struct Centers{
    centers : Vec<usize>,
}

use std::fs::File;
use std::io::prelude::*;
impl Centers{
    pub fn m(&self) -> usize {
        self.centers.len()
    }

    pub fn get(&self, i: usize) -> usize {
        self.centers[i].clone()

    }

    pub fn push(&mut self, c: usize) {
        self.centers.push(c);
    }

    pub fn iter(&self) -> std::slice::Iter<usize>{
        
        self.centers.iter()
    }


    // method: save to file
    pub fn save_to_file(&self, file_path : &str){
        let mut f = File::create(file_path).expect("Cannot open file for writing centers");
        let mut text = String::new();
        for c in self.centers.iter(){
            text = text + &format!("{},",c);
        }
        text.pop(); //delete last comma
        f.write_all(text.as_bytes()).expect("Could not write into centers-file");

    }
}

pub fn new_centers(capacity : usize) -> Centers{
    Centers{centers : Vec::with_capacity(capacity)}
}





pub struct Clustering<'a>{
    pub centers : Centers,
    pub radius : f32,
    pub center_of : Vec<usize>, // returns center to which point x belongs
    // alternative ways to save: reference: list of Cluster typs (consisting of Center and
    // Vector of Points

    space : &'a Box<dyn ColoredMetric>, // reference to metric space

    // method: cluster points of given color
    // method: print to file
}


