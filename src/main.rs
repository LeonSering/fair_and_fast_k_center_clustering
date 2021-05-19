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


fn gonzales_heuristic(space : &Box<dyn ColoredMetric>, k : usize) -> Vec<Centers> {

    let mut gonzales : Vec<Centers> = Vec::with_capacity(k+1); // Gonzales[i] has i centers
    
    // gonzales[0] is the empty set
    gonzales.push(Centers {
        centers : vec!()
    });

    // for gonzales[1] we can add any point as first center, so lets take 0
    gonzales.push(Centers {
        centers : vec!(0)
    });

    let mut dist_x_center : Vec<f32> = Vec::with_capacity(space.n()); // current distance of point x to the set of already determined centers

    for i in 0..space.n() {
        dist_x_center.push(space.dist(0,i)); // initilized to the distance of first center: 0.
    }

    for i in 2..k+1 {
        let mut current_distance = std::f32::MIN; // maximal distance to set of centers
        let mut current_point : Option<usize> = None; // corresponing point with this max distance.
        for j in 0..space.n(){
            let dist_to_newest_center = space.dist(j, gonzales[i-1].centers[i-2]); // as distance of j to gonzales[i-2] is known, we only need to measure distance to newest center.
            
            // update dist_x_S to now include the newest center:
            if dist_to_newest_center < dist_x_center[j] {
                dist_x_center[j] = dist_to_newest_center;
            }

            // check wether dist(j, gonzales[i-1]) is bigger than current biggest dist.
            if dist_x_center[j] > current_distance {
                current_distance = dist_x_center[j];
                current_point = Some(j);
            }
        }

        // create new center vector including the current farthest point (= current_point):
        let mut centers : Vec<usize> = gonzales[i-1].centers.clone(); //TODO: Maybe cloning not needed. Just save one list of centers
        centers.push(current_point.expect("No new center could be found, i.e., current_point = None"));
        gonzales.push(Centers {
            centers : centers, 
        });
    }
    gonzales
}


// An edge for every center (left) to every point (right). note that center appear on both sides.
// The distance is stored in d.
#[derive(Debug,Clone,Copy,PartialOrd,PartialEq)]
pub struct Edge { // Care: The ordering of attributes is important for the partial order! Ties in d are broken by left, then right
    d : f32,
    left : usize,
    right : usize,
}

mod space;
pub use crate::space::*;

mod clustering;
pub use crate::clustering::*;

mod buckets;
pub use crate::buckets::*;

/////////////////////////////////////////////////////////////////////////
///////////////////////// main() ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

fn main() {
    use space::*;


    // load test Space2D from file:
    let space : Box<dyn ColoredMetric> = new_space_by_2dpoints_file("test.2dspace",15);
   
    let k : usize = 5; // number of center;
    

    /////////////////////////////////
    // phase 1: Gonzalez heuristic //
    /////////////////////////////////

    assert!(k >= 1); // we want to allow at least 1 center
    assert!(space.n() >= k); // the number of points should not be less than the number of centers 

    let gonzales = gonzales_heuristic(&space, k);

    println!("** Phase 1: Determined k = {} centers by the Gonzales Heuristic: {:?}", k, gonzales[k].centers);

//    gonzales[5].save_to_file("test.centers");
    
    ///////////////////////////////////////
    // phase 2: determine privacy radius //
    ///////////////////////////////////////
    let mut edges : Vec<Edge> = Vec::with_capacity(k * space.n());
    for c in gonzales[k].centers.iter() {
        for i in 0..space.n() {
            edges.push(Edge{
                left : *c,
                right : i,
                d : space.dist(*c, i)});
        }
    }

//    println!("edges: {:?}", edges);
    
    // step 1: Compute buckets
    let buckets = put_into_buckets(edges, (4*space.n())/k);

    println!("** Phase 2a: Put n*k = {} edges into {} buckets, each of size at most 4n/k = {}.", k*space.n(), buckets.len(), (4*space.n())/k);
   
    #[cfg(debug_assertions)]
    assert!(assert_buckets_properties(&buckets, space.n(), k));
    #[cfg(debug_assertions)]
    for (i, bucket) in buckets.iter().enumerate() {
        let bucket_of_dist : Vec<f32> = bucket.iter().map(|x| x.d).collect(); 
        println!(" {}. bucket: {:?}", i, bucket_of_dist);
    }

     


    //////////////////////////////////////////////////////////////
    // phase 3: modify assignment, s.t. sizes are multiple of L //
    //////////////////////////////////////////////////////////////
    

    //////////////////////////////////////////////////////////////////////
    // phase 4: open new centers and  determine actual set of centers C //
    //////////////////////////////////////////////////////////////////////


}
