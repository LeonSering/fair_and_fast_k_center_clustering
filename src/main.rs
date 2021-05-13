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

/* Module space maintains all datastructure of metric spaces where points have colors.
 * 
 * All functionalities are stored in the public trait ColoredMetric:
 * - Distances can be obtained by using dist(x1 : usize, x2 : usize) -> f32
 * - Colors can be obtained by using color(x: usize) -> u8 
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

        let f = File::open(file_path).expect("Cannot open file!");
        let f = BufReader::new(f);
        
        
        // create vectors with initial capacity given expected number of points.
        let mut points : Vec<[f32;2]> = Vec::with_capacity(expected_number_of_points);
        let mut colors : Vec<u8> = Vec::with_capacity(expected_number_of_points);
        
        println!("capacity: {}", points.capacity());

        for line in f.lines() {
            let content = line.unwrap();
            let content: Vec<&str> = content.split(',').collect();
            points.push([content[0].parse::<f32>().expect(format!("Cannot parse x-entry to f32  on line {}.",points.len()).as_str()), content[1].parse::<f32>().expect(format!("Cannot parse y-entry to f32 on line {}.", points.len()).as_str())]);
            colors.push(content[2].parse::<u8>().expect(format!("Cannot parse color-entry to u8 on line {}",colors.len()).as_str()));
                
            println!("capacity: {}", points.capacity());
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
    }
}




struct Clustering{
    // array of centers / center indices
    // radius
    // assignment of points to centers
    // method: cluster points of given color
}

fn main() {
    use space::*;

    /*
    // create test SpaceMatrix space with 3 points
    const N : usize = 3;
    let space_by_matrix : Box<dyn ColoredMetric> =new_space_by_matrix::<{N}>([[0.0, 1.0, 2.0],
                     [1.0, 0.0, 3.0],
                     [2.0, 3.0, 0.0]], [0, 0, 1]);
    

    println!("dist between 1 and 2: {}", space_by_matrix.dist(1,2));
    println!("color of 2: {}", space_by_matrix.color(2));
    */

    /*
    // create test Space2D:
    let space_by_points : Box<dyn ColoredMetric> = new_space_by_2dpoints(vec!([0.0,0.0], [1.5,1.1], [1.0,0.5]), vec!(0,0,1));

    println!("dist between 1 and 2: {}", space_by_points.dist(1,2));
    println!("color of 2: {}", space_by_points.color(2));
    */

    // load test Space2D from file:
    let space_by_file : Box<dyn ColoredMetric> = new_space_by_2dpoint_file("test.2dspace",15);
    
    println!("dist between 1 and 2: {}", space_by_file.dist(1,2));
    println!("color of 2: {}", space_by_file.color(2));

    
    // phase 1: Gonzalez heuristic
    // Output: S[0] to S[k]
    
    // phase 2: determine privacy radius
    
    // phase 3: modify assignment, s.t. sizes are multiple of L
    
    // phase 4: open new centers and  determine actual set of centers C 


}
