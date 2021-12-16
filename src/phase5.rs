use crate::ClusteringProblem;
use crate::datastructures::{NewCenters,RootedSpanningTree};
use crate::types::{Distance,CenterIdx,PointCount,ColorIdx};
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

/// assigns clients to new centers by first executing the pushing from phase 3 (this time for real)
/// and then assigning point within these clusters to the new centers (from phase 4) by some
/// heuristcs.
/// # Input:
/// * centers_list: the new centers computed by phase 4 (one set of centers for each k)
/// * clustering_list: the assignments created by phase 2
/// * spanning_trees: the spanning_trees in the centers class that were created by phase 3 and chosen
/// (for each k) in phase 4.
/// * phase_5_gonzalez: boolean indicated if a colored gonzalez is run within each cluster (after the real pushing but
/// before the assignment). This leads to more wide-spread centers within a cluster.
///
/// # Output:
/// * i: CenterIdx - The index i of the final centers. I.e. the final centers are based on the
/// Gonzalez prefix C_i.
/// * final_centers: Centers - The final chosen centers.
/// * r: Distance - The radius r of the final centers computed by the heuristic assignment.
pub(crate) fn phase5<M : ColoredMetric + std::marker::Sync>(space : &M, prob : &ClusteringProblem, centers_list: Vec<NewCenters>, clustering_list: Vec<Clustering>, spanning_trees: Vec<RootedSpanningTree>, thread_count: usize, phase_5_gonzalez: bool) -> (CenterIdx, Clustering, Distance) {

    let mut best_radius = <Distance>::MAX;
    let mut best_clustering = None;
    let mut best_i = 0;

    let thread_pool = ThreadPoolBuilder::new().num_threads(thread_count).build().unwrap();
    let mut receivers = VecDeque::with_capacity(prob.k);

    // We need to define a scope as space and prob are moved to threads and we cannot make their
    // life-time static due to the nature of the python interface. Hence the scope defines the
    // range in which the parallel threads have to be finished.
    thread_pool.scope(|pool| {

        let mut clustering_iter = clustering_list.into_iter();
        let mut spanning_trees_iter = spanning_trees.into_iter();
        for (i, new_centers) in centers_list.into_iter().enumerate() {
            let threshold = new_centers.forrest_radius;

            let (tx, rx) = mpsc::channel();
            receivers.push_back((i,rx));

            let mut clustering = clustering_iter.next().unwrap();
            let spanning_tree = spanning_trees_iter.next().unwrap();

            // copy the read non-mut references from prob and space to move them into the thread
            let prob_ref = (&prob).clone();
            let space_ref = (&space).clone();
            pool.spawn(move |_| {
                // first realize the shifting of phase 3 but this time really shift the points:
                point_shifting(space_ref, prob.privacy_bound, &mut clustering, &spanning_tree, threshold);

                let spread_centers = if phase_5_gonzalez {spread_centers_within_cluster(space_ref, &clustering, &new_centers)} else {new_centers};


                // now assign the points of each cluster the new centers and compute the radius:
                let (radius, new_clustering) = assign_points_to_new_centers(space_ref, prob_ref, i, clustering, spread_centers);

                tx.send((radius, new_clustering)).unwrap();
            });
        }
    });

    for (i, receiver) in receivers.into_iter() {
        let (radius, clustering) = receiver.recv().unwrap();
        #[cfg(debug_assertions)]
        println!("  - C_{}: radius: {}", i, radius);

        // now only save the best clustering (over all k+1 gonzalez sets) depending on the minimum radius
        if radius < best_radius {
            best_radius = radius;
            best_i = i;
            best_clustering = Some(clustering);
        }

    }
    (best_i, best_clustering.unwrap(), best_radius)
}


fn point_shifting<M : ColoredMetric>(space : &M, privacy_bound: PointCount, clustering : &mut Clustering, spanning_tree: &RootedSpanningTree, threshold: Distance){

    let mut cluster_sizes: Vec<PointCount> = clustering.get_cluster_sizes().clone();

    // this stacks contains all nodes in the order discovered by an BFS starting from the root;
    // we push point in the reverse ordering (therefore a stack)
    let mut stack = spanning_tree.build_bfs_stack(threshold);

    while !stack.is_empty() {
        let node = stack.pop().unwrap();
        let potential_up_edge = spanning_tree.get_edge(node,threshold);
        let surplus;
        if privacy_bound == 0 {
            surplus = 0;
        } else {
            surplus = cluster_sizes[node] % privacy_bound;
        }
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

fn spread_centers_within_cluster<M : ColoredMetric>(space: &M, pushed_clustering: &Clustering, phase4_centers: &NewCenters) -> NewCenters {

    let mut all_spread_centers = Centers::with_capacity(phase4_centers.as_points.m()); // here we collect all spread centers of all clusters
    let mut spread_centers_of_cluster = Vec::with_capacity(phase4_centers.new_centers_of_cluster.len());

    for j in 0 .. pushed_clustering.m() { // go through all clusters
        let clients = pushed_clustering.get_cluster_of(j,space);
        let old_centers : Vec<&Point> = phase4_centers.new_centers_of_cluster.get(j).unwrap().iter().map(|&idx| phase4_centers.as_points.get(idx, space)).collect();

        // we now apply colored Gonzalez, meaning that the new centers should have the exact same
        // colors as the old centers in the cluster


        // first extract the color information and determine the first center to be the
        // phase4_center that is closest to the original Gonzalez center (from phase 2):
        let mut colors_to_be_opened : HashMap<ColorIdx,PointCount> = HashMap::with_capacity(old_centers.len());
        let mut first_center = None;
        let mut min_dist_to_original = Distance::MAX;
        for &c in old_centers.iter() {
            *colors_to_be_opened.entry(space.color(c)).or_insert(0) += 1;

            let dist_to_original = space.dist(c, pushed_clustering.get_center(j,space));

            if dist_to_original < min_dist_to_original {
                min_dist_to_original = dist_to_original;
                first_center = Some(c);
            }

        }

        let mut spread_centers = Centers::with_capacity(phase4_centers.as_points.m());
        spread_centers.push(first_center.unwrap());
        *colors_to_be_opened.get_mut(&space.color(first_center.unwrap())).unwrap() -= 1;

        let mut dist_x_center : Vec<Distance> = vec!(Distance::MAX; clients.len());

        for r in 1..old_centers.len() { // one center opened already old_centers.len() - 1 to go
            let mut candidate_dist = Distance::MIN;
            let mut candidate : Option<&Point> = None;


            for (s, client) in clients.iter().enumerate() {
                let dist_to_newest_center = space.dist(client, spread_centers.get(r-1,space));

                if dist_to_newest_center < dist_x_center[s] {
                    dist_x_center[s] = dist_to_newest_center;
                }

                let color = space.color(client);

                if dist_x_center[s] > candidate_dist && colors_to_be_opened.contains_key(&color) && *colors_to_be_opened.get(&color).unwrap() > 0 {
                    candidate_dist = dist_x_center[s];
                    candidate = Some(client);
                }
            }

            // it can be the case that there were no client in the cluster of an acceptable colors,
            // i.e., at least one phase4_center is not part of the cluster.
            if candidate.is_none() {

                // look for a center within the phase4_centers:
                for old_center in old_centers.iter() {
                    let dist_to_chosen_centers = space.dist_set(old_center, spread_centers.get_all(space));

                    let color = space.color(old_center);
                    if dist_to_chosen_centers > candidate_dist && colors_to_be_opened.contains_key(&color) && *colors_to_be_opened.get(&color).unwrap() > 0 {
                        candidate_dist = dist_to_chosen_centers;
                        candidate = Some(old_center);

                    }

                }

            }

            spread_centers.push(candidate.unwrap());
            *colors_to_be_opened.get_mut(&space.color(candidate.unwrap())).unwrap() -=1;

        }

        let index_start = all_spread_centers.m();
        for c in spread_centers.get_all(space) {
            all_spread_centers.push(c);
        }
        let index_end = all_spread_centers.m();
        spread_centers_of_cluster.push((index_start..index_end).collect());
    }


    NewCenters {
        as_points : all_spread_centers,
        forrest_radius : phase4_centers.forrest_radius,
        assignment_radius: phase4_centers.assignment_radius,
        new_centers_of_cluster : spread_centers_of_cluster
    }

}


fn assign_points_to_new_centers<M : ColoredMetric>(space: &M, prob: &ClusteringProblem, i : CenterIdx, old_clustering: Clustering, centers: NewCenters) -> (Distance, Clustering) {

    // the clusters and centers are final now, time to finally assign the points to
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

        let mut remaining_links : &mut [PointCenterLink] = &mut link_to_nearest;
        while remaining_links.len() > 0 {
            let split_pos = if remaining_links.len() >= batch_size {remaining_links.len()-batch_size} else {0};
            let (remainder ,batch) = utilities::ordering_split_in_two_at(remaining_links, split_pos,true);
            remaining_links = remainder;
    
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
    (radius, Clustering::new(centers.as_points, new_assignment, space))

}
