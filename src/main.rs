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



mod space;
pub use crate::space::*;

mod clustering;
pub use crate::clustering::*;

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




    /////////////////////////////////
    // phase 1: Gonzalez heuristic //
    /////////////////////////////////

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

    let mut dist_x_center : Vec<f32> = Vec::with_capacity(space.n()); // current distance of point x to the set of already determined centers

    for i in 0..space.n() {
        dist_x_center.push(space.dist(0,i)); // initilized to the distance of first center: 0.
    }

    for i in 2..k+1 {
        let mut current_distance = std::f32::MIN; // maximal distance to set of centers
        let mut current_point : Option<usize> = None; // corresponing point with this max distance.
        for j in 0..space.n(){
            let dist_to_newest_center = space.dist(j, gonzales[i-1].center_indices[i-2]); // as distance of j to gonzales[i-2] is known, we only need to measure distance to newest center.
            
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
    
    ///////////////////////////////////////
    // phase 2: determine privacy radius //
    ///////////////////////////////////////
    
    // step 1: Compute buckets
    


    //////////////////////////////////////////////////////////////
    // phase 3: modify assignment, s.t. sizes are multiple of L //
    //////////////////////////////////////////////////////////////
    

    //////////////////////////////////////////////////////////////////////
    // phase 4: open new centers and  determine actual set of centers C //
    //////////////////////////////////////////////////////////////////////


}
