use crate::clustering::Centers;
use crate::space::{Point,ColoredMetric};
use crate::types::{PointCount,Distance};

pub(crate) fn gonzalez_heuristic<M : ColoredMetric>(space : &M, k : PointCount) -> Centers {

    let mut gonzalez : Centers = Centers::with_capacity(k);//Centers{centers : Vec::with_capacity(k)};

    let first_center = space.point_iter().as_slice().first().expect("No points in metric.");
    // we can add any point as first center, so lets take 0
    gonzalez.push(first_center);

    let mut dist_x_center : Vec<Distance> = Vec::with_capacity(space.n()); // current distance of point x to the set of already determined centers


    for i in space.point_iter() {
        dist_x_center.push(space.dist(first_center,i)); // initilized to the distance of first center: 0.
    }

    for i in 1..k {
        let mut current_distance = <Distance>::MIN; // maximal distance to set of centers
        let mut current_point : Option<&Point> = None; // corresponing point with this max distance.
        for (j, p) in space.point_iter().enumerate(){
            let dist_to_newest_center = space.dist(p, gonzalez.get(i-1,space)); // as distance of j to gonzalez 0..i-2 is known, we only need to measure distance to newest center i-1.

            // update dist_x_S to now include the newest center:
            if dist_to_newest_center < dist_x_center[j] {
                dist_x_center[j] = dist_to_newest_center;
            }

            // check whether dist(j, gonzalez 0..i-2) is bigger than current biggest dist.
            if dist_x_center[j] > current_distance {
                current_distance = dist_x_center[j];
                current_point = Some(p);
            }
        }

        // create new center vector including the current farthest point (= current_point):
        gonzalez.push(current_point.expect("No new center could be found, i.e., current_point = None"));
    }
    gonzalez
}

