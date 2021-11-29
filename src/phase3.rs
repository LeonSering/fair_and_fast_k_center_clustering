use crate::ClusteringProblem;
use crate::datastructures::{OpeningList,RootedSpanningTree,UpEdge};
use crate::types::{Distance,CenterIdx,PointCount};
use crate::space::ColoredMetric;
use crate::clustering::Clustering;







/// Given a metric space and a clusting problem,
/// redistribute takes a vector of (full) clusterings, each satisfying the privacy constrait
/// and returns, for each gonzales set and each of the i+1 spanning forests, the number eta of
/// centers that can be opened in the neighborhood of each center.
/// It does not update the clusterings (shifting points around).
pub(crate) fn algebraic_pushing<M : ColoredMetric>(space : &M, prob : &ClusteringProblem, clusterings : &Vec<Clustering>) -> (Vec<RootedSpanningTree>, Vec<Vec<OpeningList>>) {

    #[cfg(debug_assertions)]
    println!("\n  - Phase 3a: Computing spanning trees.");
    let spanning_trees = compute_spanning_trees(space, prob, &clusterings);


    let mut all_opening_lists: Vec<Vec<OpeningList>> = Vec::with_capacity(prob.k);

    for (i,spanning_tree) in spanning_trees.iter().enumerate() {
        #[cfg(debug_assertions)]
        println!("\tspanning_tree for C_{}: {}", i, spanning_tree);
//        #[cfg(debug_assertions)]
//        println!("\tspanning_tree[{}] sorted distances: {:?}", i, spanning_trees[i].get_sorted_dist().iter().collect::<Vec<_>>());
//        #[cfg(debug_assertions)]
//        println!("\tspanning_tree[{}] edges_to_children: {:?}", i, spanning_trees[i].edges_to_children.iter().collect::<Vec<_>>());
        let mut opening_lists: Vec<OpeningList> = Vec::with_capacity(i+1);
        for &threshold in [0.0].iter().chain(spanning_tree.get_sorted_dist().iter()) {
            // println!("* forrest: {:?}", spanning_tree.get_edges(threshold));
            // first threshold is 0.0 so we obtain the forest without edges
            opening_lists.push(algebraic_shifting(prob.privacy_bound, &clusterings[i], i, &spanning_tree, threshold));
            // println!("\t final eta: {:?}", opening_lists[opening_lists.len()-1]);
        }
        all_opening_lists.push(opening_lists);
    }
    (spanning_trees, all_opening_lists)
}


fn algebraic_shifting(privacy_bound: PointCount, clustering : &Clustering, i : CenterIdx, spanning_tree: &RootedSpanningTree, threshold: Distance) -> OpeningList {


    let mut cluster_sizes: Vec<PointCount> = clustering.get_cluster_sizes().clone();
    let mut eta: Vec<Option<PointCount>> = vec!(None;i+1); // starts with the size of each cluster and is then converted to be the muliplier of L

    // this stacks contains all nodes in the order discovered by an BFS starting from the root;
    // we push point in the reverse ordering (therefore a stack)
    let mut stack = spanning_tree.build_bfs_stack(threshold);

    while !stack.is_empty() {
        let node = stack.pop().unwrap();
        // println!("cluster_sizes:{:?}, \teta:{:?}", cluster_sizes, eta);
        let potential_up_edge = spanning_tree.get_edge(node,threshold);
        eta[node] = Some(cluster_sizes[node] / privacy_bound);
        let surplus = cluster_sizes[node] - eta[node].unwrap()*privacy_bound;
        if surplus > 0 {
            match potential_up_edge {
                None => {} // node is a root
                Some(up_edge) => {
                    let parent = up_edge.up;
                    cluster_sizes[parent] += surplus;
                    cluster_sizes[node] -= surplus;
                }
            }
        }

    }

    OpeningList{eta: eta.into_iter().map(|e| e.unwrap()).collect(), forrest_radius : threshold}
}





////// Spanning Tree Computation ////////

/// Given a set of clusterings (one for each gonzales set), return a minimum spanning tree on the
/// centers, one for each clustering
fn compute_spanning_trees<'a, M : ColoredMetric>(space : &M, prob : &ClusteringProblem, clusterings : &Vec<Clustering>) -> Vec<RootedSpanningTree> {

    let mut spanning_trees : Vec<RootedSpanningTree> = Vec::with_capacity(prob.k);

    for clustering in clusterings.iter() {
        let i = clustering.m() - 1; // the number of centers; i.e. we consider gonzales set S_i (counting from 0 to k-1)

        // first lets do the algorithm of prim to compute the minimum spanning tree.
        // takes O(i^2) time

        let mut spanning_tree: Vec<UpEdge> = Vec::with_capacity(i+1);

        let mut dist_to_tree: Vec<Option<Distance>> = Vec::with_capacity(i+1); // None means beeing part of the tree
        let mut closest: Vec<CenterIdx> = Vec::with_capacity(i+1); // denotes the cloeset center that is already in the tree


        let mut edges_to_children: Vec<Vec<UpEdge>> = vec!(Vec::new(); i+1);

        // we start with center 0 as singleton spanning tree
        closest.push(0); // closest to 0 is 0;
        dist_to_tree.push(None); // 0 is our tree at the beginning
        for c in 1..i+1 {
            // priority_queue.push(c,Reverse(Priority{d : space.dist(clustering.get_center(0), clustering.get_center(c))}));
            dist_to_tree.push(Some(space.dist(clustering.get_center(0,space), clustering.get_center(c,space))));
            closest.push(0);
        }

        for _ in 1..i+1 {
            // let c1 = argmin dist_to_tree[c], i.e., its the center that is closest to the tree
            // (but not part of it yet)
            let (c1,_) = dist_to_tree.iter().enumerate().filter(|(_,dist)| dist.is_some()).min_by(|(_,dist1),(_,dist2)| dist1.as_ref().unwrap().partial_cmp(dist2.as_ref().unwrap()).unwrap()).unwrap();

            // add c1 to the tree by adding the arc from c1 to the closest center within the tree
            let edge = UpEdge{down : c1, up : closest[c1], d: dist_to_tree[c1].unwrap()};
            edges_to_children[closest[c1]].push(edge.clone());
            spanning_tree.push(edge);

            // now the distances to the tree (of centers not yet part of it) are outdated and needs
            // to be updated
            dist_to_tree[c1] = None;
            for c2 in 1..i+1 {
                if dist_to_tree[c2].is_some() { // if it is none it is already part of the tree
                    let new_d = space.dist(clustering.get_center(c1,space), clustering.get_center(c2,space));
                    if new_d < dist_to_tree[c2].unwrap() {
                        dist_to_tree[c2] = Some(new_d);
                        closest[c2] = c1;
                    }
                }


            }
        }

        // sort edges such that spanning_tree[5] refers to the edge from node 5 to the parent of 5
        // (node 0 has not parent so the value is 0
        spanning_tree.sort_by(|a,b| a.down.cmp(&b.down));
        let mut spanning_tree_edges : Vec<Option<UpEdge>> = Vec::with_capacity(i+1);
        spanning_tree_edges.push(None); // The Root (Center 0) has no UpEdge
        spanning_tree_edges.extend(spanning_tree.into_iter().map(|e| Some(e)));



        spanning_trees.push(RootedSpanningTree::new(spanning_tree_edges, edges_to_children));


    }
    spanning_trees
}

