use crate::ClusteringProblem;
use crate::datastructures::{NewCenters,RootedSpanningTree};
use crate::types::{Distance,CenterIdx,PointCount};
use crate::space::{Point,ColoredMetric};
use crate::clustering::{Clustering,Centers};
use crate::utilities;

use std::collections::{VecDeque,HashMap};

// for parallel execution:
use std::sync::mpsc;
use rayon::ThreadPoolBuilder;

/// A point-center-pair with computed distance:
#[derive(Debug,Clone,PartialOrd,PartialEq,Copy)]
struct PointCenterLink<'a> {
    dist : Distance,
    point : &'a Point,
    center_idx : CenterIdx,
}

pub(crate) fn phase5<M : ColoredMetric + std::marker::Sync>(space : &M, prob : &ClusteringProblem, centers_list: Vec<NewCenters>, clustering_list: &mut Vec<Clustering>, spanning_trees: &Vec<RootedSpanningTree>, thread_count: usize) -> (CenterIdx, Centers, Distance) {


    let mut best_radius = <Distance>::MAX;
    let mut best_i = 0;

    let thread_pool = ThreadPoolBuilder::new().num_threads(thread_count).build().unwrap();
    let mut receivers = VecDeque::with_capacity(prob.k);

    // We need to define a scope as space and prob are moved to threads and we cannot make their
    // life-time static due to the nature of the python interface. Hence the scope defines the
    // range in which the parallel threads have to be finished.
    thread_pool.scope(|pool| {

        let mut clustering_iter = clustering_list.into_iter();
        let mut spanning_trees_iter = spanning_trees.into_iter();
        for (i, new_centers) in centers_list.iter().enumerate() {
            let threshold = new_centers.forrest_radius;

            let (tx, rx) = mpsc::channel();
            receivers.push_back((i,rx));

            let clustering = clustering_iter.next().unwrap();
            let spanning_tree = spanning_trees_iter.next().unwrap();

            // copy the read non-mut references from prob and space to move them into the thread
            let prob_ref = (&prob).clone();
            let space_ref = (&space).clone();
            pool.spawn(move |_| {
                // first realize the shifting of phase 3 but this time really shift the points:
                point_shifting(space_ref, prob.privacy_bound, clustering, &spanning_tree, threshold);

                // now assign the points of each cluster the new centers and compute the radius:
                let radius = assign_points_to_new_centers(space_ref, prob_ref, i, &clustering, new_centers);
                tx.send(radius).unwrap();
            });
        }
    });

    for (i, receiver) in receivers.into_iter() {
        let radius = receiver.recv().unwrap();
        #[cfg(debug_assertions)]
        println!("  - C_{}: radius: {}", i, radius);

        // now only save the best clustering (over all k+1 gonzales sets) depending on the minimum radius
        if radius < best_radius {
            best_radius = radius;
            best_i = i;
        }

    }
    (best_i, centers_list[best_i].as_points.clone(), best_radius)
}

fn point_shifting<M : ColoredMetric>(space : &M, privacy_bound: PointCount, clustering : &mut Clustering, spanning_tree: &RootedSpanningTree, threshold: Distance){

    let mut cluster_sizes: Vec<PointCount> = clustering.get_cluster_sizes().clone();

    // this stacks contains all nodes in the order discovered by an BFS starting from the root;
    // we push point in the reverse ordering (therefore a stack)
    let mut stack = spanning_tree.build_bfs_stack(threshold);

    while !stack.is_empty() {
        let node = stack.pop().unwrap();
        let potential_up_edge = spanning_tree.get_edge(node,threshold);
        let surplus = cluster_sizes[node] % privacy_bound;
        // println!("\tcluster_sizes:{:?}, \tsurplus:{:?}", cluster_sizes, surplus);
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


fn hand_over<M: ColoredMetric>(space: &M, clustering: &mut Clustering, supplier: CenterIdx, consumer: CenterIdx, number: PointCount) {
    let mut point_center_pairs: Vec<PointCenterLink> = clustering.get_cluster_of(supplier,space).iter().cloned().map(|p|
            PointCenterLink { dist: space.dist(p, clustering.get_centers().get(consumer,space)),
                              point: p,
                              center_idx: consumer }).collect();
    utilities::truncate_to_smallest(&mut point_center_pairs, number);
    // println!("Point_center_pairs after truncate:{:?}", point_center_pairs);
    for PointCenterLink{dist:_, point: p, center_idx:_} in point_center_pairs {
        clustering.assign(p, consumer, space);
    }
}


fn assign_points_to_new_centers<M : ColoredMetric>(space: &M, prob: &ClusteringProblem, i : CenterIdx, old_clustering: &Clustering, centers: &NewCenters) -> Distance {

    // now the clusters are given, now we really open the new centers and assign the points to
    // them:
    let mut new_assignment: Vec<Option<CenterIdx>> = vec!(None; space.n());
    let mut radius = <Distance>::MIN;
    for cluster in 0..i+1 {
        let points = old_clustering.get_cluster_of(cluster,space);
        let center_idx = centers.new_centers_of_cluster[cluster].clone();
        // let center_pt = center_idx.iter().map(|idx| new_centers.centers.get(*idx)).collect();
        // println!("\n  ** Cluster: {} with {} center and {} points", cluster, center_idx.len(), points.len());

        // Creates a link from each point to its nearest center (reference by index):
        let mut link_to_nearest: Vec<PointCenterLink> = Vec::with_capacity(points.len());
        for point in points.iter() {
            let mut current_closest = center_idx[0];
            let mut current_dist = <Distance>::MAX;
            for &idx in center_idx.iter() {
                let dist = space.dist(point, centers.as_points.get(idx,space));
                if dist <= current_dist {
                    current_closest = idx;
                    current_dist = dist;
                }
            }
            link_to_nearest.push(PointCenterLink{point,center_idx:current_closest,dist:current_dist});
        }

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
                        let dist = space.dist(entry.point, centers.as_points.get(c,space));
                        if dist < dist_to_chosen_center {
                            chosen_center = c;
                            dist_to_chosen_center = dist;
                        }

                    }

                }
                new_assignment[entry.point.idx()] = Some(chosen_center);
                if dist_to_chosen_center > radius {
                    radius = dist_to_chosen_center;
                }


                let counter = num_points_assigned_to.get_mut(&chosen_center).unwrap();
                *counter += 1;
                if *counter == prob.privacy_bound {
                    // center is now_private
                    let idx = non_private_centers.iter().position(|c| *c == chosen_center).unwrap();
                    non_private_centers.remove(idx);
                }
            } // end of entry

        } // end of batch

    } // and of cluster
    radius

}
