use crate::{ClusteringProblem,Centers,Clustering,ColoredMetric,space::{PointIdx,ColorIdx,Distance}};
use crate::clustering::CenterIdx;
use super::OpeningList;
use crate::utilities;

/// an edge between a gonzales center and a color class;
/// lablled with the point and the distance between center and point
#[derive(Debug,Clone,Copy,PartialOrd,PartialEq)]
struct ColorEdge {
    d: Distance, // distance between the gonzales center and the point
    center: CenterIdx, // index to the gonzales center
    point: PointIdx, // a point that might be a final center
    color: ColorIdx, // the color of the point
}

/// Given a metric space and a clustering problem,
/// finalize takes a vector of clusterings in which each center (except for one) covers a multple of L
/// points, and returns a vector of clusterings (which potentially more cluster as before) that also satisfy the representative constaints. 
///
pub(crate) fn finalize<'a, M : ColoredMetric>(space : &'a M, prob : &ClusteringProblem, opening_lists : Vec<Vec<OpeningList>>,  gonzales : &Centers) -> Vec<Clustering<'a>> {



    //TEMP: For test reasons we also look at flow problems with ALL edges present:
    let mut all_edges_of_cluster: Vec<Vec<ColorEdge>> = vec!(Vec::with_capacity(space.n());prob.k);
    for i in 0..prob.k {
        for p in space.point_iter() {
            let color = space.color(p);
            all_edges_of_cluster[i].push(ColorEdge{
                d: space.dist(gonzales.get(i), p),
                center: i,
                point: p.idx(),
                color});
            
        }
        all_edges_of_cluster[i].sort_by(|a,b| a.partial_cmp(b).unwrap());
    }
    // TEMP END
    
    
    // define the neighborhood of each gonzales center:
    let mut edges_of_cluster: Vec<Vec<ColorEdge>> = vec!(Vec::with_capacity(prob.k);prob.k);

    for i in 0..prob.k {

        // first we make edge lists for each color class
        let mut edges_by_color: Vec<Vec<ColorEdge>> = vec!(Vec::new();space.gamma());

        // for this we create all edges from gonzales center i to all points
        for p in space.point_iter() {
            let color = space.color(p);
            edges_by_color[space.color(p)].push(ColorEdge{
                d: space.dist(gonzales.get(i), p),
                center: i,
                point: p.idx(),
                color});
            
        }


        let restricted_colors = if space.gamma() < prob.rep_interval.len() {space.gamma()} else {prob.rep_interval.len()}; // min{gamma, rep_interval.len()}
        let mut remaining_edges: Vec<ColorEdge> = Vec::new();
        let mut num_edges_to_fill = prob.k;

        for c in 0..restricted_colors {
            // take only the b smallest in each class; all others cannot play any role.
            utilities::truncate_to_smallest(&mut edges_by_color[c], prob.rep_interval[c].1);
            // the a smallest have to be present for sure, so each center can satisfy this
            // condition by itself;
            // the remaining b-a edges are collected in remaining_edges
            remaining_edges.append(&mut utilities::split_off_at(&mut edges_by_color[c], prob.rep_interval[c].0));
            num_edges_to_fill -= prob.rep_interval[c].0;
        }

        // for all colors without restriction fill them into the remaining edges:
        for c in restricted_colors..space.gamma() {
            remaining_edges.append(&mut edges_by_color[c]);
        }
        
//        println!("\n\n************ i = {}; center = {:?} ****************", i, gonzales.get(i));
//        for (c, cclass) in edges_by_color.iter().enumerate() {
//            println!("\nColorclass {}: {:?}", c, cclass);
//        }
//        println!("\nRemaining: {:?}", remaining_edges);

        
        // so far sum(a) edges where chosen, so we fill up with the nearest k - sum(a) edges,
        // independent of the color 
        utilities::truncate_to_smallest(&mut remaining_edges,num_edges_to_fill);



        // put all edges into one list:
        for c in 0..restricted_colors {
            edges_of_cluster[i].append(&mut edges_by_color[c]);
        }
        edges_of_cluster[i].append(&mut remaining_edges);

        // there are now max k edges in edges_of_cluster[i]; so we can sort them:
        edges_of_cluster[i].sort_by(|a,b| a.partial_cmp(b).unwrap());

//        println!("\nedges_of_cluser[{}]: {:?}", i, edges_of_cluster[i]);




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

