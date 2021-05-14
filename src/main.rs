struct ClusterProblem{
    // metric space
    // number of colors (gamma)
    // number of points n
    // max number of clusters k
    // for each color j: upper and lower bounds for centers a_j, b_j
    // representative loweder bound L
    //
    // method: test if valid: color classes is metric space are 0, ..., gamma-1; sum of a_j <= k.
    //
    // Solving algorithm as method? Maybe not.
}



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
mod space { //TODO: put space module in separate file.
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
    pub fn new_space_by_2dpoint_file(file_path : &str, expected_number_of_points : usize) -> Box<dyn ColoredMetric> {

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
}

//////////////////////////////////////////////////////////////
//////////////////// module: clustering //////////////////////
//////////////////////////////////////////////////////////////

/* Contains two structs (both are only valid together with a space
 * 
 * centers: a vector of centers
 * clustering: centers, radius and assignment of points to centers
 *
 */
mod clustering{
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







    pub struct Clustering{
        centers : Centers,
        // array of centers / center indices
        // radius
        // assignment of points to centers
        // reference to space
        // method: cluster points of given color
    }

}

/////////////////////////////////////////////////////////////////////////
///////////////////////// main() ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

fn main() {
    use space::*;

    /* 
    // create test SpaceMatrix space with 3 points
    const N : usize = 3;
    let space_by_matrix : Box<dyn ColoredMetric> =new_space_by_matrix::<{N}>(
                    [[0.0, 2.0, 1.5],
                     [2.0, 0.0, 0.6],
                     [1.5, 0.6, 0.0]], [0, 0, 1]);
    

    println!("dist between 1 and 2: {}", space_by_matrix.dist(1,2));
    println!("color of 2: {}", space_by_matrix.color(2));
    println!("number of points: {}", space_by_matrix.n());
    assert!(space_by_matrix.is_metric());
    */

    /*
    // create test Space2D:
    let space_by_points : Box<dyn ColoredMetric> = new_space_by_2dpoints(vec!([0.0,0.0], [1.5,1.1], [1.0,0.5]), vec!(0,0,1));

    println!("dist between 1 and 2: {}", space_by_points.dist(1,2));
    println!("color of 2: {}", space_by_points.color(2));
    */

    // load test Space2D from file:
    let space : Box<dyn ColoredMetric> = new_space_by_2dpoint_file("test.2dspace",15);
   
    for i in 2..8 {
        println!("dist between 1 and {}: {}", i, space.dist(1,i));
    }
    let set : Vec<usize> = vec![2,3,4,5,6,7];
    
    /*println!("set: {:?}", set);
    println!("dist between 1 and [2,3,4,5,6,7]: {}", space.dist_set(1,&set));
    println!("set: {:?}", set);
    println!("color of 2: {}", space.color(2));
    println!("number of points: {}", space.n());*/





    // phase 1: Gonzalez heuristic
    
    use clustering::*;
    let k : usize = 5; // number of center;
    
    assert!(k >= 1); // we want to allow at least 1 center
    assert!(space.n() >= k); // the number of points should not be less than the number of centers 

    let mut gonzales : Vec<Centers> = Vec::with_capacity(k+1); // Gonzales[i] has i centers
    
    // gonzales[0] is the empty set
    gonzales.push(Centers {
        center_indices : vec!()
    });

    // for gonzales[1] we can add any point as first center, so lets take 0
    gonzales.push(Centers {
        center_indices : vec!(0)
    });

    let mut dist_x_center : Vec<f32> = Vec::with_capacity(space.n());

    for i in 0..space.n() {
        dist_x_center.push(space.dist(0,i));
    }

    for i in 2..k+1 {
        let mut current_distance = std::f32::MIN;
        let mut current_point : Option<usize> = None;
        println!("Iteration {}", i);
        for j in 0..space.n(){
            let dist_to_newest_center = space.dist(j, gonzales[i-1].center_indices[i-2]);
            // update dist_x_S to now include the newest center
            if dist_to_newest_center < dist_x_center[j] {
                dist_x_center[j] = dist_to_newest_center;
            }
            if dist_x_center[j] > current_distance {
                current_distance = dist_x_center[j];
                current_point = Some(j);
            }
        }
        let mut centers : Vec<usize> = gonzales[i-1].center_indices.clone();
        centers.push(current_point.expect("No new center could be found, i.e., current_point = None"));
        gonzales.push(Centers {
            center_indices : centers, 
        });
    }

    for i in 0..k+1 {
        println!("gonzales[{}]: {:?}", i, gonzales[i].center_indices);
    }

    gonzales[5].save_to_file("test.centers");
    // Output: S[0] to S[k]
    
    // phase 2: determine privacy radius
    
    // phase 3: modify assignment, s.t. sizes are multiple of L
    
    // phase 4: open new centers and  determine actual set of centers C 


}
