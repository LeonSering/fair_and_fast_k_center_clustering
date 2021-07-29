use crate::ClusteringProblem;
use crate::datastructures::{NewCenters,RootedSpanningTree};
use crate::types::{Distance,CenterIdx,PointCount};
use crate::space::{Point,ColoredMetric};
use crate::clustering::{Centers,Clustering};
use crate::utilities;

use std::collections::HashMap;



/// A point-center-pair with computed distance:
#[derive(Debug,Clone,PartialOrd,PartialEq,Copy)]
struct PointCenterLink<'a> {
    dist : Distance,
    point : &'a Point,
    center_idx : CenterIdx,
}

pub(crate) fn phase5<'a, M : ColoredMetric>(space : &'a M, prob : &ClusteringProblem, final_centers: Vec<NewCenters<'a>>, clusterings: &mut Vec<Clustering<'a>>, spanning_trees: &Vec<RootedSpanningTree>) -> Clustering<'a> {


    let mut current_best_clustering: Vec<Option<CenterIdx>> = Vec::new();
    let mut current_best_centers: Centers = Centers::with_capacity(0);
    let mut current_best_radius = <Distance>::MAX;

    for (i, new_centers) in final_centers.iter().enumerate() {
        let threshold = new_centers.forrest_radius;

        println!("\n*** i = {} *** \t spanning_tree: {:?}", i, spanning_trees[i].get_edges(threshold));
        println!("  new_centers:{:?}", new_centers);

        // first realize the shifting of phase 3 but this time really shift the points:
        point_shifting(space, prob.privacy_bound, &mut clusterings[i], &spanning_trees[i], threshold);
        println!("Cluster sizes after point_shifting: {:?}", clusterings[i].get_cluster_sizes());

        // now the clusters are given, now we really open the new centers and assign the points to
        // them:
        let mut new_clustering: Vec<Option<CenterIdx>> = vec!(None; space.n());
        let mut new_radius = <Distance>::MIN;

        for cluster in 0..i+1 {
            let points = clusterings[i].get_cluster_of(cluster);
            let center_idx = new_centers.new_centers_of[cluster].clone();
            // let center_pt = center_idx.iter().map(|idx| new_centers.centers.get(*idx)).collect();
            println!("\n  ** Cluster: {} with {} center and {} points", cluster, center_idx.len(), points.len());

            let mut link_to_nearest: Vec<PointCenterLink> = Vec::with_capacity(points.len());
            for point in points.iter() {
                let mut current_closest = center_idx[0];
                let mut current_dist = <Distance>::MAX;
                for &idx in center_idx.iter() {
                    let dist = space.dist(point, new_centers.centers.get(idx));
                    if dist <= current_dist {
                        current_closest = idx;
                        current_dist = dist;
                    }
                }
                link_to_nearest.push(PointCenterLink{point,center_idx:current_closest,dist:current_dist});
            }
            // Creates a link from each point to its nearest center (reference by index):
            // let mut link_to_nearest: Vec<PointCenterLink> = points.iter().map(|p| (p, space.get_closest(p, &center_pt))).map(|(point, (dist, center_idx))| PointCenterLink{point,center_idx,dist}).collect();

            // each new center has to have at least privacy_bound many points. There might be
            // points to spare:
            let mut to_spare = points.len() - center_idx.len() * prob.privacy_bound;

            // we count how many points have been assigned to each center so far:
            let mut num_points_assigned_to: HashMap<CenterIdx, PointCount> = HashMap::with_capacity(center_idx.len());
            for &c in center_idx.iter() {
                num_points_assigned_to.insert(c, 0);
            }

            // we also track the list of centers that are not private yet:
            let mut non_private_centers = center_idx.clone();

            // goal: consider batches of size L and assign the farthest to their nearest new
            // center.
            // As L could be very small, we would end up sorting the points. For that reason the
            // batch size has to have a size of at least O(n/k). That way we look at the farthest
            // batch at most O(k) times, which would lead to a running time of O(nk). (And doing
            // this k times).

            let tradeoff_parameter = 10; // high value means better radius but longer running time.
            let safety_n_divided_k = space.n() / (prob.k * tradeoff_parameter);
            let batch_size = if prob.privacy_bound > safety_n_divided_k {prob.privacy_bound} else {safety_n_divided_k};



            while link_to_nearest.len() > 0 {
                let split_pos = if link_to_nearest.len() >= batch_size {link_to_nearest.len()-batch_size} else {0};
                let batch = utilities::split_off_at(&mut link_to_nearest, split_pos);
                // println!("\tbatch_size: {}, batch.len(): {}, batch: {:?}", batch_size, batch.len(), batch.iter().map(|x| x.dist).collect::<Vec<_>>());

                for entry in batch {
                    let mut chosen_center = center_idx[0]; // we first set it to the closest center
                    let mut dist_to_chosen_center = <Distance>::MAX;
                    if *num_points_assigned_to.get(&entry.center_idx).unwrap() < prob.privacy_bound || to_spare > 0{
                        chosen_center = entry.center_idx;
                        dist_to_chosen_center = entry.dist;
                        if *num_points_assigned_to.get(&entry.center_idx).unwrap() >= prob.privacy_bound {
                            // in this case the center was private already so we use one of our
                            // spare counter
                            to_spare -= 1;
                        }
                    } else {
                        // in this case we cannot assign the point to its closest center.
                        // instead assign it to the closest non-private center.

                        for &c in non_private_centers.iter() {
                            let dist = space.dist(entry.point, new_centers.centers.get(c));
                            if dist < dist_to_chosen_center {
                                chosen_center = c;
                                dist_to_chosen_center = dist;
                            }

                        }
                        
                    }
                    new_clustering[entry.point.idx()] = Some(chosen_center);
                    if dist_to_chosen_center > new_radius {
                        new_radius = dist_to_chosen_center;
                    }


                    let counter = num_points_assigned_to.get_mut(&chosen_center).unwrap();
                    *counter += 1;
                    if *counter == prob.privacy_bound {
                        // center is now_private
                        let idx = non_private_centers.iter().position(|c| *c == chosen_center).unwrap();
                        non_private_centers.remove(idx);
                    }
                    // println!("\n\tspare: {}, chosen center (idx): {}, non_private_centers, {:?}", to_spare, chosen_center, non_private_centers);
                    // println!("\tnum_points_assigned_to: {:?}", num_points_assigned_to);
                    


                } // end of entry

            } // end of batch

        } // and of cluster 
        println!("\n\tnew_radius: {}", new_radius);

        //TEMP:
        let clustering = Clustering::new(new_centers.centers.clone(), new_clustering.clone(), space);
        let save_path = format!("output/temp/after_phase_5_for_i_{}.clustering", i);
        clustering.save_to_file(save_path.as_str());



        // now only save the best clustering (over all k+1 gonzales sets) depending on the minimum radius
        if new_radius < current_best_radius {
            current_best_centers = new_centers.centers.clone();
            current_best_clustering = new_clustering;
            current_best_radius = new_radius;
        }


    } // end of gonzales set



    Clustering::new(current_best_centers, current_best_clustering, space)
}

fn point_shifting<'a, M : ColoredMetric>(space : &'a M, privacy_bound: PointCount, clustering : &mut Clustering<'a>, spanning_tree: &RootedSpanningTree, threshold: Distance){

    let mut cluster_sizes: Vec<PointCount> = clustering.get_cluster_sizes().clone();

    // this stacks contains all nodes in the order discovered by an BFS starting from the root;
    // we push point in the reverse ordering (therefore a stack)
    let mut stack = spanning_tree.build_bfs_stack(threshold);
    
    while !stack.is_empty() {
        let node = stack.pop().unwrap();
        let potential_up_edge = spanning_tree.get_edge(node,threshold);
        let surplus = cluster_sizes[node] % privacy_bound;
        println!("cluster_sizes:{:?}, \tsurplus:{:?}", cluster_sizes, surplus);
        if surplus > 0 {
            match potential_up_edge {
                None => {} // node is a root
                Some(up_edge) => {
                    let parent = up_edge.up;
                    // TODO: We could be more conservative and only hand_over as many points as
                    // needed.
                    hand_over(space, clustering, node, parent, surplus);
                    cluster_sizes[parent] += surplus;
                    cluster_sizes[node] -= surplus;
                }
            }
        }

    }

}


fn hand_over<'a, M: ColoredMetric>(space: &'a M, clustering: &mut Clustering<'a>, supplier: CenterIdx, consumer: CenterIdx, number: PointCount) {
    let mut point_center_pairs: Vec<PointCenterLink> = clustering.get_cluster_of(supplier).iter().cloned().map(|p| 
            PointCenterLink { dist: space.dist(p, clustering.get_centers().get(consumer)),
                              point: p,
                              center_idx: consumer }).collect();
    utilities::truncate_to_smallest(&mut point_center_pairs, number);
    println!("Point_center_pairs after truncate:{:?}", point_center_pairs);
    for PointCenterLink{dist:_, point: p, center_idx:_} in point_center_pairs {
        clustering.assign(p, consumer, space);
    }
}
