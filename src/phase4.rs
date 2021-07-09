use crate::{ClusteringProblem,Centers,Clustering,ColoredMetric};
use crate::types::{PointCount,PointIdx,ColorIdx,Distance,CenterIdx};
use super::OpeningList;
use crate::utilities;
use std::collections::HashMap;

mod color_flow;
use color_flow::compute_flow;

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
    
    /////////////////////////////////////////////////////////////////    
    //////// define the neighborhood of each gonzales center: ///////
    /////////////////////////////////////////////////////////////////    
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
        println!("edges of center i: {:?}", edges_of_cluster[i]);

    }






    /////////////////////////////////////////////////////////////////    
    ///////////////////////// solve flow problem: ///////////////////
    /////////////////////////////////////////////////////////////////    


    let sum_of_a: PointCount = prob.rep_interval.iter().map(|interval| interval.0).sum();
    
    for i in 0..prob.k {

        println!("\n\n************** i = {} ******************", i);

        // consider all edges starting from centers 0 to i in increasing order:
    
        let mut edges: Vec<&ColorEdge> = Vec::with_capacity((i+1) * prob.k); // a reference to all edges with centers in S_i
        for j in 0..i+1 {
            edges.extend(edges_of_cluster[j].iter());
        }

        edges.sort_by(|a,b| a.partial_cmp(b).unwrap());

        


        // not that we only have k^2 edges so at most k^2 point and color classes can occur
        // We reindex the relevant points and colors. The new indices are called point-node and color-node index and they are
        // used in the flow network
       
        let mut point_to_node: HashMap<PointIdx, usize> = HashMap::with_capacity((i + 1) * prob.k); // maps a point index to the point-node index used in the flow network
        let mut point_counter = 0; // number of relevant points (= number of point-nodes)
        
        let mut a: Vec<PointCount> = Vec::with_capacity((i + 1) * prob.k); // lower bound given by the color-node index
        let mut b: Vec<PointCount> = Vec::with_capacity((i +1) * prob.k); // upper bound given by the color-node index
        let mut color_to_node: HashMap<ColorIdx, usize> = HashMap::with_capacity((i + 1) * prob.k); // maps a color index to the color-node index used in the flow network
        let mut color_counter = 0; // number of relevant colors (= number of color-nodes)


        let mut point_idx_to_color_idx: Vec<ColorIdx> = Vec::with_capacity((i+1) * prob.k);
        for e in edges.iter() {
            if !color_to_node.contains_key(&e.color) {
                color_to_node.insert(e.color, color_counter);
                if e.color >= prob.rep_interval.len() {
                    a.push(0);
                    b.push(<PointCount>::MAX);
                } else {
                    a.push(prob.rep_interval[e.color].0);
                    b.push(prob.rep_interval[e.color].1);
                }
                color_counter += 1;
            }
            if !point_to_node.contains_key(&e.point) {
                point_to_node.insert(e.point, point_counter);
                point_idx_to_color_idx.push(*color_to_node.get(&e.color).unwrap());
                point_counter += 1;
            }
        }



        // edges should also be referrable from their points and their color classes:
        
        let mut edges_by_color_node: Vec<Vec<&ColorEdge>> = vec!(Vec::new(); color_counter);
        let mut edges_by_point_node: Vec<Vec<&ColorEdge>> = vec!(Vec::new(); point_counter);

        for e in edges.iter() {
            edges_by_color_node[color_to_node[&e.color]].push(e);
            edges_by_point_node[point_to_node[&e.point]].push(e);
        }
        edges_by_color_node.sort_by(|a,b| a.partial_cmp(b).unwrap());
        edges_by_point_node.sort_by(|a,b| a.partial_cmp(b).unwrap());

       


        // for each eta vector we solve the flow problem individually:

        for opening in opening_lists[i].iter() {

            let network = Network {
                k : prob.k,
                opening,
                i,
                number_of_points: point_counter,
                number_of_colors: color_counter,
                sum_of_a,
                a : &a,
                b : &b,
                edges_of_cluster : &edges_of_cluster,
                point_to_node : &point_to_node,
//                color_to_node : &color_to_node,
                point_idx_to_color_idx : &point_idx_to_color_idx,
                edges_by_color_node : &edges_by_color_node,
                edges_by_point_node : &edges_by_point_node,
            };

            compute_flow(opening,&network,&edges);


            // TODO: Create centers from max flow
        }

    }





    let mut clusterings: Vec<Clustering> = Vec::with_capacity(prob.k);
    let mut centers = Centers::with_capacity(prob.k);
    centers.push(space.point_iter().next().unwrap());
    let center_of: Vec<Option<CenterIdx>> = space.point_iter().map(|_| None).collect();
    clusterings.push(Clustering::new(centers,center_of, space));
    clusterings
}


struct Network<'a> {
    k : CenterIdx,
    opening : &'a OpeningList,
    i : CenterIdx,
    number_of_points : PointCount,
    number_of_colors : usize,
    sum_of_a : usize,
    a : &'a Vec<usize>,
    b : &'a Vec<usize>,
    edges_of_cluster : &'a Vec<Vec<ColorEdge>>,
    point_to_node : &'a HashMap<PointIdx, usize>,
    point_idx_to_color_idx : &'a Vec<usize>,
    edges_by_color_node : &'a Vec<Vec<&'a ColorEdge>>,
    edges_by_point_node : &'a Vec<Vec<&'a ColorEdge>>,
}



