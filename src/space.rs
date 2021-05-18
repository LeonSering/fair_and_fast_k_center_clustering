///////////////////////////////////////////////////////////////
///////////////////// module: space ///////////////////////////
///////////////////////////////////////////////////////////////

/* Module space maintains all datastructure of metric spaces where points have colors.
 *
 * All functionalities are stored in the public trait ColoredMetric:
 * - Distances can be obtained by using dist(x1 : usize, x2 : usize) -> f32
 * - Colors can be obtained by using color(x: usize) -> u8
 * - Distance to a set can be obtained by dist_set(x: usize, point_set: &Vec<usize>) -> f32
 * - The number of points can be obtained by n() -> usize
 *
 * The most general metric can be created by new_space_by_matrix
 * Input: distances_by_NxN-Matrix: [[f32; N]; N],
 *        colors: [u8; N]
 * Output: Box<dyn ColoredMetric>
 *
 */

// Trait to obtain the distance between two points.
pub trait ColoredMetric{

    fn dist(&self, x1 : usize, x2 : usize) -> f32; // returns the distance between x1 and x2

    fn color(&self, x : usize) -> u8;

    fn dist_set(&self, x : usize, point_set: &Vec<usize>) -> f32 {
        let mut current_distance = f32::MAX;
        let mut d: f32;
        for i in 0..point_set.len() {
            d = self.dist(x,point_set[i]);
            if d < current_distance {
                current_distance = d;
            }
        }
        current_distance
    }

    fn n(&self) -> usize;

    // Care this need O(n^3) time
    fn is_metric(&self) -> bool {

        // check for symmetry, non-negativity and identity of indiscernibles
        for x in 0..self.n() {
            for y in 0..self.n() {
                let dist_xy = self.dist(x,y);
                let dist_yx = self.dist(y,x);
                if dist_xy != dist_yx {
                    println!("Symmetry is violated: x={}, y={}, dist(x,y)={}, dist(y,x)={}", x, y, dist_xy, dist_yx);
                    return false;
                }

                if x == y && dist_xy != 0.0{
                    println!("Identity of indiscernibles is violated: x={}, dist(x,x)={}", x, dist_xy);
                    return false;
                }
                if x != y && dist_xy <= 0.0 {
                    println!("Non-negativity or identity of indiscernibles is violated: x={}, y={}, dist(x,y)={}", x, y, dist_xy);
                    return false;
                }
            }
        }

        // check for triangle inequility
        for x in 0..self.n() {
            for y in 0..self.n() {
                let dist_xy = self.dist(x,y);
                for z in 0..self.n() {
                    let dist_xz = self.dist(x,z);
                    let dist_zy = self.dist(z,y);
                    if dist_xy > dist_xz + dist_zy {
                        println!("Triangle inequality is violated: x={}, y={}, z={}, dist(x,y)={}, dist(x,z)={}, dist(z,y)={}, i.e., {} > {}", x, y, z, dist_xy, dist_xz, dist_zy, dist_xy, dist_xz + dist_zy);
                        return false;
                    }
                }
            }
        }
        return true;
    }
}


//////////////////// SpaceMatrix /////////////////////////

// a metric space of N point with color classes. Distances are given by NxN-array and color class
// by N-array.
struct SpaceMatrix<const N: usize>{
    distances: [[f32; N]; N], // distance matrix (later maybe implicit as distance function to avoid n^2 space)
    colors: [u8; N], // points as array or implicit? If implicity: array of color-classes.
}

// builder:
pub fn new_space_by_matrix<const N: usize>(distances : [[f32; N]; N], colors: [u8; N]) -> Box<dyn ColoredMetric> {
    Box::new(SpaceMatrix {
        distances,
        colors,
    })
}

impl<const N: usize> ColoredMetric for SpaceMatrix<N> {
    fn dist(&self, x1: usize, x2: usize) -> f32 {
        self.distances[x1][x2]
    }

    fn color(&self, x: usize) -> u8 {
        self.colors[x]
    }

    fn n(&self) -> usize {
        N
    }
}




///////////////////////// Space2D /////////////////////////

// the euclidean metric space in 2D.

// TODO: expand this to k-dim space by using "const generic"
struct Space2D{
    points : Vec<[f32;2]>, // points are saved as vector of tuples
    colors : Vec<u8>,
}

// builder
pub fn new_space_by_2dpoints(points : Vec<[f32;2]>, colors : Vec<u8>) -> Box<dyn ColoredMetric> {
    Box::new(Space2D {
        points,
        colors,
    })
}

use std::fs::File;
use std::io::{BufReader,BufRead};
pub fn new_space_by_2dpoints_file(file_path : &str, expected_number_of_points : usize) -> Box<dyn ColoredMetric> {

    let f = File::open(file_path).expect("Cannot open file to read 2dpoints.");
    let f = BufReader::new(f);


    // create vectors with initial capacity given expected number of points.
    let mut points : Vec<[f32;2]> = Vec::with_capacity(expected_number_of_points);
    let mut colors : Vec<u8> = Vec::with_capacity(expected_number_of_points);

    for line in f.lines() {
        let content = line.unwrap();
        let content: Vec<&str> = content.split(',').collect();
        points.push([content[0].parse::<f32>().expect(format!("Cannot parse x-entry to f32  on line {}.",points.len()).as_str()), content[1].parse::<f32>().expect(format!("Cannot parse y-entry to f32 on line {}.", points.len()).as_str())]);
        colors.push(content[2].parse::<u8>().expect(format!("Cannot parse color-entry to u8 on line {}",colors.len()).as_str()));

    }
    println!("Successfully loaded {} points/colors from '{}'", points.len(), file_path);
    println!("points: {:?}", points);
    println!("colors: {:?}", colors);
    Box::new(Space2D {
        points,
        colors,
    })
}

impl ColoredMetric for Space2D {
    fn dist(&self, x1: usize, x2: usize) -> f32 { // euclidean norm
        let d_squared : f32 = (self.points[x1][0] - self.points[x2][0]) * (self.points[x1][0] - self.points[x2][0]) + (self.points[x1][1] - self.points[x2][1]) *(self.points[x1][1] - self.points[x2][1]);
        d_squared.sqrt()
    }

    fn color(&self, x : usize) -> u8 {
        self.colors[x]
    }

    fn n(&self) -> usize {
        self.points.len()
    }
}



#[cfg(test)]
mod tests {
    use crate::*;
    #[test]
    fn matrix_space_3x3() { 
    // create test SpaceMatrix space with 3 points
    const N : usize = 3;
    let space_by_matrix : Box<dyn ColoredMetric> = new_space_by_matrix::<{N}>(
                    [[0.0, 2.0, 1.5],
                     [2.0, 0.0, 0.6],
                     [1.5, 0.6, 0.0]], [0, 0, 1]);
    
    assert!(space_by_matrix.is_metric());
    assert_eq!(space_by_matrix.dist(1,2),0.6);
    assert_eq!(space_by_matrix.color(2),1);
    assert_eq!(space_by_matrix.n(),3);
    }

    #[test]
    fn point_space_3_points() {
        // create test Space2D:
        let space_by_points : Box<dyn ColoredMetric> = new_space_by_2dpoints(vec!([0.0,0.0], [1.5,1.1], [1.0,0.5]), vec!(0,0,1));

        assert!(space_by_points.is_metric());


        assert_eq!(space_by_points.dist(1,2),f32::sqrt(0.25+0.36));
        assert_eq!(space_by_points.color(2),1);
        assert_eq!(space_by_points.n(),3);
        println!("dist between 1 and 2: {}", space_by_points.dist(1,2));
    }
}
