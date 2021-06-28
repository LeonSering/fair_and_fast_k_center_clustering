use crate::{ClusteringProblem,Centers,Clustering,ColoredMetric,Point};
use crate::clustering::CenterIdx;
use super::OpeningList;

/// Given a metric space and a clustering problem,
/// finalize takes a vector of clusterings in which each center (except for one) covers a multple of L
/// points, and returns a vector of clusterings (which potentially more cluster as before) that also satisfy the representative constaints. 
///
pub(crate) fn finalize<'a, M : ColoredMetric>(space : &'a M, prob : &ClusteringProblem, opening_lists : Vec<Vec<OpeningList>>,  gonzales : &Centers) -> Vec<Clustering<'a>> {


    // first we make lists for each color class
    let mut points_of_color: Vec<Vec<&Point>> = vec!(Vec::new();space.gamma());


    for p in space.point_iter() {
        points_of_color[space.color(p)].push(p);
        
    }

//    for cclass in points_of_color.iter() {
//        println!("Colorclass: {:?}", cclass);
//    }

    // TODO: Maybe at list of points of colors with a = 0
    
    
    // then we define the colored neighborhood of each gonzales center:
    
    let mut points_around: Vec<Vec<&Point>> = vec!(Vec::new();prob.k);
    for (i, c) in gonzales.iter().enumerate() {
        println!("i = {}; center = {:?}", i, c);

        // first determine distance to the a_i nearest neighbor of the color class (median algorithm)
        // take all these points
        // second, determine the b_i - a_i nearest neighbor of all color classes
        // take the remaining points from the union of that
        //
        // Care: b_i could be larger than the number of points in the color class
        
    }
    


    for (i, openings) in opening_lists.iter().enumerate() {
        //
//        println!("** i = {}; openings = {:?}", i, openings);


        // then we define the network

    }




    let mut clusterings: Vec<Clustering> = Vec::with_capacity(prob.k);
    let mut centers = Centers::with_capacity(prob.k);
    centers.push(space.point_iter().next().unwrap());
    let center_of: Vec<Option<CenterIdx>> = space.point_iter().map(|_| None).collect();
    clusterings.push(Clustering::new(centers,center_of, space));
    clusterings
}

