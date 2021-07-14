use crate::{ClusteringProblem,OpeningList};
use crate::types::{Distance,CenterIdx,PointCount};
use crate::space::ColoredMetric;
use crate::clustering::Clustering;

/// a pointer from a gonzales center to its parent in a rooted forest 
#[derive(Debug,Clone)]
struct UpEdge {
    d : Distance,
    up : CenterIdx, // the parent of down in the rooted forest
    down : CenterIdx, // the child of up in the rooted forest
}

struct RootedSpanningTree {
    edges: Vec<Option<UpEdge>>, // vector of size i; for each node s_j \in S_i it either points to the parent node or is None for the root 
    edges_to_children: Vec<Vec<UpEdge>>, // for each node we have a list of edges pointing to the children
//    sorted_distances: Option<Vec<Distance>>, // a list of all distances in edges in increasing order
}


impl RootedSpanningTree{

    /// returns a reference to the TreeEdge, if its distance is not greater than the threshold
    fn get_edge(&self, node: CenterIdx, threshold : Distance) -> Option<&UpEdge> {
        if let Some(ref edge) = self.edges[node] {
            // node has an up-going edge
            if edge.d  <= threshold {Some(edge)} else {None}
        } else {
            // node has no up-going edge
            None 
        }
    }

    fn get_edges(&self, threshold: Distance) -> Vec<&UpEdge> {
        self.edges.iter().filter_map(|e| e.as_ref()).filter(|e| e.d <= threshold).collect()
    }

    fn get_edges_to_children(&self, node: CenterIdx, threshold: Distance) -> Vec<&UpEdge> {
        self.edges_to_children[node].iter().filter(|e| e.d <= threshold).collect()
    }

    fn get_sorted_dist(&self) -> Vec<Distance> {
        let mut distances : Vec<Distance> = self.edges.iter().filter_map(|o| o.as_ref()).map(|e| e.d).collect();
        distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
        distances
    }
}

use std::collections::VecDeque;

/// Given a metric space and a clusting problem,
/// redistribute takes a vector of (full) clusterings, each satisfying the privacy constrait
/// and returns, for each gonzales set and each of the i+1 spanning forests, the number eta of
/// centers that can be opened in the neighborhood of each center. 
pub(crate) fn redistribute<'a, M : ColoredMetric>(space : &M, prob : &ClusteringProblem, clusterings : Vec<Clustering<'a>>) -> Vec<Vec<OpeningList>> {
    println!("\n  ** Phase 3a: Computing spanning trees");
    let spanning_trees = compute_spanning_trees(space, prob, &clusterings);

    let mut all_opening_lists: Vec<Vec<OpeningList>> = Vec::with_capacity(prob.k);

    for (i,spanning_tree) in spanning_trees.iter().enumerate() {
        println!("** spanning_tree[{}]: {:?}", i, spanning_trees[i].edges);
//        #[cfg(debug_assertions)]
//        println!("\tspanning_tree[{}] sorted distances: {:?}", i, spanning_trees[i].get_sorted_dist().iter().collect::<Vec<_>>());
//        #[cfg(debug_assertions)]
//        println!("\tspanning_tree[{}] edges_to_children: {:?}", i, spanning_trees[i].edges_to_children.iter().collect::<Vec<_>>());
        let clustering = &clusterings[i];
        let mut opening_lists: Vec<OpeningList> = Vec::with_capacity(i+1);


        
        for &threshold in [-1.0].iter().chain(spanning_tree.get_sorted_dist().iter()) {
            // first threshold is -1.0 so we obtain the forest without edges
//            println!("\tthreshold: {}", threshold);

//            println!("\tForest: {:?}", spanning_tree.get_edges(threshold));
//
            let mut cluster_sizes: Vec<PointCount> = clustering.get_cluster_sizes().clone();
            let mut eta: Vec<Option<PointCount>> = vec!(None;i+1); // starts with the size of each cluster and is then converted to be the muliplier of L
//            let mut pushed: Vec<bool> = vec!(false;i+1);


            let mut stack: Vec<CenterIdx> = Vec::with_capacity(i+1); // elements poped here are treated

            // look for roots:
            let mut is_root= vec!(true; i+1);
            for e in spanning_tree.get_edges(threshold).iter() {
                is_root[e.down] = false;
            }

//            println!("roots:{:?}", is_root);
            let mut queue: VecDeque<CenterIdx> = VecDeque::with_capacity(i+1); // this queue is only for building the stack
            let mut pushed = vec!(false; i+1); // true if node has been pushed to the BFS queue 
            for j in (0..i+1).filter(|l| is_root[*l]) {
                queue.push_back(j);
                pushed[j] = true;
            }

            // we now need to do a BFS starting from the roots, and push each discovered node to
            // the stack. Hence, by emptying the stack we go the BFS in reverse. 
            while !queue.is_empty() {
                let node = queue.pop_front().unwrap();

                // we push the node to the stack
                stack.push(node);


                for &edge in spanning_tree.get_edges_to_children(node, threshold).iter() {
                    let next = edge.down;
                    if pushed[next] == false { 
                        pushed[next] = true; 
                        queue.push_back(next);
                    }
                }




            }
                


            while !stack.is_empty() {
                let node = stack.pop().unwrap();
//                println!("cluster_sizes:{:?}, \teta:{:?}", cluster_sizes, eta);
                let potential_up_edge = spanning_tree.get_edge(node,threshold);
                eta[node] = Some(cluster_sizes[node] / prob.privacy_bound);
                cluster_sizes[node] -= eta[node].unwrap()*prob.privacy_bound;
                match potential_up_edge {
                    None => {} // node is a root
                    Some(up_edge) => {
                        let parent = up_edge.up;
                        cluster_sizes[parent] += cluster_sizes[node];
                        cluster_sizes[node] = 0;
                    }
                }

            }
//            println!("final cluster_sizes:{:?}, eta:{:?}", cluster_sizes, eta);

            opening_lists.push(OpeningList{eta: eta.into_iter().map(|e| e.unwrap()).collect()});
        }
        all_opening_lists.push(opening_lists);
    }
    all_opening_lists
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
            dist_to_tree.push(Some(space.dist(clustering.get_center(0), clustering.get_center(c))));
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
                    let new_d = space.dist(clustering.get_center(c1), clustering.get_center(c2));
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

        

        spanning_trees.push(RootedSpanningTree{edges : spanning_tree_edges, edges_to_children});


    }
    spanning_trees
}

