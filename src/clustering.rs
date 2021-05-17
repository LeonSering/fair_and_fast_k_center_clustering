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
    pub center_indices : Vec<usize>,
}

use std::fs::File;
use std::io::prelude::*;
impl Centers{
    pub fn m(&self) -> usize {
        self.center_indices.len()
    }


    // method: save to file
    pub fn save_to_file(&self, file_path : &str){
        let mut f = File::create(file_path).expect("Cannot open file for writing centers");
        let mut text = String::new();
        for c in self.center_indices.iter(){
            text = text + &format!("{},",c);
        }
        text.pop(); //delete last comma
        println!("{}", text);
        f.write_all(text.as_bytes()).expect("Could not write into centers-file");

    }
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


