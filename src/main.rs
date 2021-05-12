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


// Trait to obtain the distance between two points.
pub trait Metric{
    fn dist(&self, x1 : usize, x2 : usize) -> f32; // returns the distance between x1 and x2
}


// Trait to obtain the color class of a point.
pub trait Colored{
    fn color(&self, x : usize) -> u8;
}

// a metric space of N point with color classes. Distances are given by NxN-array and color class
// by N-array.
pub struct Space<const N: usize>{
    distances: [[f32; N]; N], // distance matrix (later maybe implicit as distance function to avoid n^2 space)
    colors: [u8; N], // points as array or implicit? If implicity: array of color-classes.
}

/*pub fn new_space(distances : [[f32; N]; N], colors: [u8; N]) -> Space<const N> {
    Space<const N> {
        distances,
        colors,
    }
}*/
impl<const N: usize> Metric for Space<N> {
    fn dist(&self, x1: usize, x2: usize) -> f32 {
        self.distances[x1][x2]
    }
}

impl<const N: usize> Colored for Space<N> {
    fn color(&self, x: usize) -> u8 {
        self.colors[x]
    }
}

struct Clustering{
    // array of centers / center indices
    // radius
    // assignment of points to centers
    // method: cluster points of given color
}

fn main() {

    // create test metric space with 3 points

    const N : usize = 3;
    let space1 : Space<{N}> = Space {
        distances : [[0.0, 1.0, 2.0],
                     [1.0, 0.0, 3.0],
                     [2.0, 3.0, 0.0]],
        colors : [0, 0, 1],
    };

    println!("dist between 1 and 2: {}", space1.dist(1,2));
    println!("color of 2: {}", space1.color(2));
    // phase 1: Gonzalez heuristic
    // Output: S[0] to S[k]
    
    // phase 2: determine privacy radius
    
    // phase 3: modify assignment, s.t. sizes are multiple of L
    
    // phase 4: open new centers and  determine actual set of centers C


}
