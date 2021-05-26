use crate::clustering::{Centers,new_centers};
use crate::space::ColoredMetric;

pub fn gonzales_heuristic(space : &Box<dyn ColoredMetric>, k : usize) -> Centers {

    let mut gonzales : Centers = new_centers(k);//Centers{centers : Vec::with_capacity(k)};

    // we can add any point as first center, so lets take 0
    gonzales.push(0);

    let mut dist_x_center : Vec<f32> = Vec::with_capacity(space.n()); // current distance of point x to the set of already determined centers

    for i in 0..space.n() {
        dist_x_center.push(space.dist(0,i)); // initilized to the distance of first center: 0.
    }

    for i in 1..k {
        let mut current_distance = std::f32::MIN; // maximal distance to set of centers
        let mut current_point : Option<usize> = None; // corresponing point with this max distance.
        for j in 0..space.n(){
            let dist_to_newest_center = space.dist(j, gonzales.get(i-1)); // as distance of j to gonzales 0..i-2 is known, we only need to measure distance to newest center i-1.
            
            // update dist_x_S to now include the newest center:
            if dist_to_newest_center < dist_x_center[j] {
                dist_x_center[j] = dist_to_newest_center;
            }

            // check whether dist(j, gonzales 0..i-2) is bigger than current biggest dist.
            if dist_x_center[j] > current_distance {
                current_distance = dist_x_center[j];
                current_point = Some(j);
            }
        }

        // create new center vector including the current farthest point (= current_point):
        gonzales.push(current_point.expect("No new center could be found, i.e., current_point = None"));
    }
    gonzales
}

