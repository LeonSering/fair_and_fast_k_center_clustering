///////////////////////////////////////////////////////////////
///////////////////// module: space ///////////////////////////
///////////////////////////////////////////////////////////////

/// Module space maintains all datastructure of metric spaces where points have colors.
///
/// All functionalities are stored in the public trait ColoredMetric:
/// - Distances can be obtained by using dist(x1 : PointIdx, x2 : PointIdx) -> Distance
/// - Colors can be obtained by using color(x: PointIdx) -> ColorIdx
/// - Distance to a set can be obtained by dist_set(x: PointIdx, point_set: &Vec<PointIdx>) -> Distance
/// - The number of points can be obtained by n() -> PointCount 
///
/// The most general metric can be created by new_space_by_matrix
/// Input: distances_by_NxN-Matrix: [[Distance; N]; N],
///        colors: [ColorCount; N]
/// Output: SpaceMatrix 
///
/// Other builder functions are new_space_by_2dpoints and new_space_by_2dpoints_file
/// These create Euclidean metrics in the plane, either by loading a file or by a array of typles [(x,y)].
///

use crate::types::{PointIdx,ColorIdx,PointCount,ColorCount,Distance};

/// A point of a metric space. Their only attribute is an index, which can be obtained by idx().
/// Points are create by metric spaces and can only be accessed via the point_inter() method of the
/// [ColoredMetric] trait.
#[derive(Debug,PartialOrd,PartialEq,Eq,Hash)]
pub struct Point{ // a point in the metric
    index : PointIdx,
}

impl Point {
    /// Return the index of the point.
    pub fn idx(&self) -> PointIdx {
        self.index
    }
}

/// Trait for a metric space with colored points.
pub trait ColoredMetric{
    
    /// Returns the distance between two points x1 and x2.
    fn dist(&self, x1 : &Point, x2 : &Point) -> Distance; // returns the distance between x1 and x2

    /// Returns the the color of point x.
    fn color(&self, x : &Point) -> ColorIdx;

    /// Returns the distance between a point x and a set of points.
    fn dist_set(&self, x : &Point, point_set: Vec<&Point>) -> Distance {
        let mut current_distance = Distance::MAX;
        for p in point_set {
            let d = self.dist(x,p);
            if d < current_distance {
                current_distance = d;
            }
        }
        current_distance
    }

    fn get_closest<'a>(&self, x : &Point, point_set: &Vec<&'a Point>) -> (Distance, &'a Point) {
        assert!(point_set.len() > 0, "there has to be at least one point in the point set.");
        let mut current_distance = Distance::MAX;
        let mut current_closest: &Point = point_set[0];
        for p in point_set {
            let d = self.dist(x,p);
            if d < current_distance {
                current_distance = d;
                current_closest = p;
            }
        }
        (current_distance, current_closest)
    }
    
    /// Return the number of points in the metric space.
    fn n(&self) -> PointCount; // return number of points

    /// Provides an iterator of all points of the metric space. This is the only way to access the
    /// points.
    fn point_iter(&self) -> std::slice::Iter<Point>;
    
    /// Return a reference (wrapped in Some) to the point with provided index.
    /// If there is not point with this index, return None.
    fn get_point(&self, idx: PointIdx) -> Option<&Point>;

    /// Returns the number of color classes that are present in the metric space. More precesely it
    /// return the highest color value + 1.
    fn gamma(&self) -> ColorCount;

    /// Checks whether the dist function satisfy the metric properties in which case true is
    /// returned.
    /// Care: this need O(n<sup>3</sup>) time.
    fn is_metric(&self) -> bool {

        // check for symmetry, non-negativity and identity of indiscernibles
        for x in self.point_iter() {
            for y in self.point_iter() {
                let dist_xy = self.dist(x,y);
                let dist_yx = self.dist(y,x);
                if dist_xy != dist_yx {
                    println!("Symmetry is violated: x={:?}, y={:?}, dist(x,y)={}, dist(y,x)={}", *x, *y, dist_xy, dist_yx);
                    return false;
                }

                if x == y && dist_xy != 0.0{
                    println!("Identity of indiscernibles is violated: x={:?}, dist(x,x)={}", *x, dist_xy);
                    return false;
                }
                if x != y && dist_xy <= 0.0 {
                    println!("Non-negativity or identity of indiscernibles is violated: x={:?}, y={:?}, dist(x,y)={}", *x, *y, dist_xy);
                    return false;
                }
            }
        }

        // check for triangle inequility
        for x in self.point_iter() {
            for y in self.point_iter() {
                let dist_xy = self.dist(x,y);
                for z in self.point_iter() {
                    let dist_xz = self.dist(x,z);
                    let dist_zy = self.dist(z,y);
                    if dist_xy > dist_xz + dist_zy {
                        println!("Triangle inequality is violated: x={:?}, y={:?}, z={:?}, dist(x,y)={}, dist(x,z)={}, dist(z,y)={}, i.e., {} > {}", *x, *y, *z, dist_xy, dist_xz, dist_zy, dist_xy, dist_xz + dist_zy);
                        return false;
                    }
                }
            }
        }
        return true;
    }
}


//////////////////// SpaceMatrix /////////////////////////

/// A general finite metric space with color classes. </br> 
/// Distances are given by a symetric distance matrx of type nxn; </br>
/// Colors are given by a list of size n; </br>
/// Implements the [ColoredMetric] trait.
pub struct SpaceMatrix{
    distances: Vec<Vec<Distance>>, // distance matrix (later maybe implicit as distance function to avoid n^2 space)
    points: Vec<Point>,
    colors: Vec<ColorIdx>, // points as array or implicit? If implicity: array of color-classes.
    gamma: ColorCount,
}


impl SpaceMatrix {
    /// Creates a new [SpaceMatrix];
    /// distances has to be a 2d-vector that represents a quadratic matrix satisfying the metric properties
    /// (symmetry, non-negativity, identity of indiscernibles and the triangle inequality);
    /// colors is a vector that stores the color of each point;
    ///
    /// To assert the metric properties, the run time is O(n^3).
    ///
    /// # Panics
    ///
    /// Panics, if distances is not a quadratic matrix;
    /// Panics, if the length of the color classes does not match the number of points given by distances;
    /// Panics, if distances does not satisfy the metric properties;
    pub fn new(distances: Vec<Vec<Distance>>, colors: Vec<ColorIdx>) -> SpaceMatrix {
        // check if distances is a qudratic matrix:
        let number_of_rows = distances.len();
        let mut points: Vec<Point> = Vec::with_capacity(number_of_rows);
        for (i, row) in distances.iter().enumerate() {
            assert_eq!(row.len(), number_of_rows, "Matrix is not qudratic. row {} has {} entries; number of rows: {}", i, row.len(), number_of_rows);
            points.push(Point{index:i});
        }
        assert_eq!(number_of_rows, colors.len(), "Number of points: {} do not match number of colors: {}", number_of_rows, colors.len());
        let gamma = colors.iter().max().expect("No maximal color found") + 1;

        let space = SpaceMatrix {
            distances,
            colors,
            points,
            gamma,
        };

        assert!(space.is_metric(), "distances do not satisfy the metric properties.");
        space
    }
    
    /// Creates a new [SpaceMatrix] as in new but the distances and color data are received by arrays.
    ///
    /// # Example
    /// ```rust
    /// use ff_k_center::{Point,SpaceMatrix,ColoredMetric};
    /// let space = SpaceMatrix::new_by_array(
    ///                 [[0.0, 2.0, 1.5],
    ///                  [2.0, 0.0, 0.6],
    ///                  [1.5, 0.6, 0.0]], [0, 0, 1]);
    /// let points : Vec<&Point> = space.point_iter().collect(); 
    /// assert!(space.is_metric());
    /// assert_eq!(space.dist(points[1],points[2]),0.6);
    /// assert_eq!(space.color(points[2]),1);
    /// assert_eq!(space.n(),3);
    /// ``` 
    pub fn new_by_array<const N: PointCount>(distances : [[Distance; N]; N], colors: [ColorIdx; N]) -> SpaceMatrix {

        // this is kind of silly, to doulbe initilize the points, but I don't see a way to make this
        // without unsafe rust.
        //let mut points = [Point{index : 0}; N];
        let mut points : Vec<Point> = Vec::with_capacity(N);
        for i in 0..N {
            points.push(Point{index:i});
        }
        let distances : Vec<Vec<Distance>> = distances.iter().map(|row| row.to_vec()).collect();
        let colors = colors.to_vec();
        let gamma = colors.iter().max().expect("No maximal color found") + 1;
        SpaceMatrix {
            distances,
            colors,
            points,
            gamma,
        }
    }
}

impl ColoredMetric for SpaceMatrix {
    fn dist(&self, x1: &Point, x2: &Point) -> Distance {
        self.distances[x1.idx()][x2.idx()]
    }

    fn color(&self, x: &Point) -> ColorIdx {
        self.colors[x.idx()]
    }

    fn n(&self) -> PointCount {
        self.points.len() 
    }

    fn get_point(&self, idx: PointIdx) -> Option<&Point> {
        if idx >= self.points.len() {
            return None;
        }
        Some(&self.points[idx])
    }

    fn point_iter(&self) -> std::slice::Iter<Point> {
        self.points.iter()
    }
    
    fn gamma(&self) -> ColorCount {
        self.gamma
    }

}




///////////////////////// Space2D /////////////////////////
type Position = Vec<f32>;
// the euclidean metric space in 2D.

// TODO: expand this to k-dim space by using "const generic"

/// A metric space in the euklidean space of some dim N. Implements the [ColoredMetric] trait.
/// Beside the color it also stores the position of type Vec<f32> of each point.
/// The distance is computed by the Eucleadean metric.
pub struct SpaceND{
    points : Vec<Point>,
    positions : Vec<Position>,
    colors : Vec<ColorIdx>,
    gamma : ColorCount,
}
use rand::Rng;
use std::fs::File;
use std::io::{BufReader,BufRead};

impl SpaceND {
    /// Creates a new metric space of type [SpaceND].
    /// Input: An vector of positions of type Vec<f32> and a vector of colors of type u16.
    ///
    /// # Panics
    ///
    /// Panics if positions is an empty vector.
    /// Panics if the two vectors have different size.
    /// Panics if there are positions of different lengths.
    ///
    /// # Example
    ///
    /// ```rust
    /// use ff_k_center::{Point,SpaceND,ColoredMetric};
    /// let space_by_points = SpaceND::by_ndpoints(vec!(vec![0.0,0.0], vec![1.5,1.1], vec![1.0,0.5]), vec!(0,0,1));
    /// 
    /// assert!(space_by_points.is_metric());
    /// 
    /// let points : Vec<&Point> = space_by_points.point_iter().collect(); 
    /// 
    /// assert_eq!(space_by_points.dist(points[1],points[2]),<f32>::sqrt(0.25+0.36));
    /// assert_eq!(space_by_points.color(points[2]),1);
    /// assert_eq!(space_by_points.n(),3);
    /// println!("dist between 1 and 2: {}", space_by_points.dist(points[1],points[2]));
    /// ```
    ///
    pub fn by_ndpoints(positions : Vec<Position>, colors : Vec<ColorIdx>) -> SpaceND {
        assert!(positions.len() > 0, "There need to be at least one position given.");

        assert_eq!(positions.len(),colors.len(),"The number of points in position must equal the number of colors!");

        let gamma = colors.iter().max().expect("No maximal color found") + 1;

        let dimension = positions[0].len();

        for i in 1..positions.len() {
            assert_eq!(positions[i].len(), dimension, "Dimension is {}, but positions[{}] has only {} entries", dimension, i, positions[i].len());
        }

        SpaceND {
            points : (0..positions.len()).map(|i| Point{index : i}).collect(),
            positions,
            colors,
            gamma,
        }
    }

    /// Crates a new metric space of type [SpaceND] of dimension 2.
    /// It containt n random points in the [-100,100]x[-100,100] box with random colors from [1..10]
    pub fn new_random(n : PointCount) -> SpaceND {
        let mut rng = rand::thread_rng();
        let positions = (0..n).map(|_| vec![rng.gen_range(-100.0f32..100.0f32), rng.gen_range(-100.0f32..100.0f32)]).collect();
        let colors = (0..n).map(|_| rng.gen_range(0..10)).collect();
        SpaceND::by_ndpoints(positions, colors)
    }

    /// Loads a new metric space of type [SpaceND] from a file.
    /// The expected_number_of_points is used to allocate enough storage.
    /// File_path must point into a text-file that stores a tuple in each line, separated by a comma.
    /// The first entries of each line must be of type f32 specifying the position;
    /// the last entry must be a non-negative integer specifying the color (of type u16).
    /// 
    /// Example:
    /// ```txt
    /// -8.19,-7.88,0
    /// -8.06,-6.58,0
    /// -7.3,-6.9,0
    /// -5.97,-8.26,0
    /// ```
    ///
    /// # Panics
    ///
    /// Panics if the file cannot be open of if it cannot parse the triplets.
    pub fn by_file(file_path : &str, expected_number_of_points : PointCount) -> SpaceND {

        let f = File::open(file_path).expect("Cannot open file to read ndpoints.");
        let f = BufReader::new(f);


        // create vectors with initial capacity given expected number of points.
        let mut positions : Vec<Position> = Vec::with_capacity(expected_number_of_points);
        let mut colors : Vec<ColorIdx> = Vec::with_capacity(expected_number_of_points);

        let mut dim = 0;

        for line in f.lines() {
            let content = line.unwrap();
            let content: Vec<&str> = content.split(',').collect();
            if dim == 0 {
                dim = content.len() - 1;
            } else {
                assert_eq!(dim, content.len() - 1, "Line {} has the wrong number of entries.", positions.len());
            }

            positions.push((0..dim).map(|i| content[i].parse::<Distance>().expect(format!("Cannot parse entry {} to f32  on line {}.", i, positions.len()).as_str())).collect());
            colors.push(content[dim].parse::<ColorIdx>().expect(format!("Cannot parse color-entry to u16 on line {}",colors.len()).as_str()));

        }
        println!("\n**** Successfully loaded {} points/colors (dimension: {}) from '{}'", positions.len(), dim, file_path);
        
        #[cfg(debug_assertions)]
        {
            print!("    positions:");
            let mut counter = 0;
            for p in positions.iter() {
                if counter % 10 == 0 {
                    print!("\n\t{}-{}:\t",counter, <usize>::min(counter+9,positions.len()));
                }
                print!("{:?} ", p);
                counter += 1;
            }
            println!("\n    colors: {:?}", colors);
        }
        let gamma = colors.iter().max().expect("No maximal color found") + 1;
        SpaceND {
            points : (0..positions.len()).map(|i| Point{index : i}).collect(),
            positions,
            colors,
            gamma,
        }
    }
}

impl ColoredMetric for SpaceND {
    fn dist(&self, x1: &Point, x2: &Point) -> Distance { // euclidean norm
        let dim = self.positions[x1.idx()].len();
        let d_squared : Distance = (0..dim).map(|i| (self.positions[x1.idx()][i] - self.positions[x2.idx()][i]) * (self.positions[x1.idx()][i] - self.positions[x2.idx()][i])).sum();
        d_squared.sqrt()
    }

    fn color(&self, x : &Point) -> ColorIdx {
        self.colors[x.idx()]
    }

    fn get_point(&self, idx: PointIdx) -> Option<&Point> {
        if idx >= self.points.len() {
            return None;
        }
        Some(&self.points[idx])
    }

    fn point_iter(&self) -> std::slice::Iter<Point> {
        self.points.iter()
    }

    fn n(&self) -> PointCount {
        self.points.len()
    }

    fn gamma(&self) -> ColorCount {
        self.gamma
    }
    
}

