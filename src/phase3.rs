use crate::ClusteringProblem;
use crate::space::{Distance,ColoredMetric};
use crate::clustering::{CenterIdx,Clustering};
use super::{OpeningList,PointCount};

/// a pointer from a gonzales center to its parent in a rooted forest 
#[derive(Debug)]
struct UpEdge {
    d : Distance,
    up : CenterIdx, // the parent of down in the rooted forest
    down : CenterIdx, // the child of up in the rooted forest
}

struct RootedSpanningTree {
    edges: Vec<Option<UpEdge>>, // vector of size i; for each node s_j \in S_i it either points to the parent node or is None for the root 
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

    fn get_sorted_dist(&self) -> Vec<Distance> {
        let mut distances : Vec<Distance> = self.edges.iter().filter_map(|o| o.as_ref()).map(|e| e.d).collect();
        distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
        distances
    }
//    fn get_sorted_dist_iter(&mut self) -> std::slice::Iter<Distance> {
//        if self.sorted_distances.is_some() {
//            return self.sorted_distances.as_ref().unwrap().iter()
//        }
//        let mut distances : Vec<Distance> = self.edges.iter().filter_map(|o| o.as_ref()).map(|e| e.d).collect();
//        distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
//        self.sorted_distances = Some(distances);
//        self.sorted_distances.as_ref().unwrap().iter()
//    }
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
//        println!("** spanning_tree[{}]: {:?}", i, spanning_trees[i].edges);
        #[cfg(debug_assertions)]
        println!("\tspanning_tree[{}] sorted distances: {:?}", i, spanning_trees[i].get_sorted_dist().iter().collect::<Vec<_>>());
        let clustering = &clusterings[i];
        let mut opening_lists: Vec<OpeningList> = Vec::with_capacity(i+1);


        
        for &threshold in [-1.0].iter().chain(spanning_tree.get_sorted_dist().iter()) {
            // first threshold is -1.0 so we obtain the forest without edges
//            println!("\tthreshold: {}", threshold);

//            println!("\tForest: {:?}", spanning_tree.get_edges(threshold));

            let mut eta: Vec<PointCount> = clustering.get_cluster_sizes().clone(); // starts with the size of each cluster and is then converted to be the muliplier of L
            let mut pushed: Vec<bool> = vec!(false;i+1);



            // look for leafs:
            let mut is_leaf= vec!(true; i+1);
            for e in spanning_tree.get_edges(threshold).iter() {
                is_leaf[e.up] = false;
            }

//            println!("leafs:{:?}", is_leaf);
            let mut queue: VecDeque<CenterIdx> = VecDeque::with_capacity(i+1);
            for j in (0..i+1).filter(|l| is_leaf[*l]) {
                queue.push_back(j);
                pushed[j] = true;
            }
//            println!("queue with leafs:{:?}", queue);

            while !queue.is_empty() {
//                println!("eta:{:?}", eta);
                let node = queue.pop_front().unwrap();
                let potential_up_edge = spanning_tree.get_edge(node,threshold);
                match potential_up_edge {
                    None => { // node is a root
                        eta[node] = eta[node] / prob.privacy_bound;
                    }
                    Some(up_edge) => {
                        let parent = up_edge.up;
                        let remainder = eta[node] % prob.privacy_bound;
                        eta[parent] += remainder;
                        eta[node] = eta[node] / prob.privacy_bound;
                        if !pushed[parent] {
                            queue.push_back(parent);
                            pushed[parent] = true;
                        }
                    }
                }

            }
//            println!("final eta:{:?}", eta);

            opening_lists.push(OpeningList{eta});
        }
        all_opening_lists.push(opening_lists);
    }
    all_opening_lists
}


////// Spanning Tree Computation ////////

use priority_queue::PriorityQueue;
use std::cmp::{Reverse,Ordering};

#[derive(Copy,Clone,PartialEq, PartialOrd)]
struct Priority {
    d : Distance,
}

impl Eq for Priority {}

impl Ord for Priority {
    fn cmp(&self, other: &Self) -> Ordering {
        // small d means high priority
        if self.d < other.d {
            Ordering::Less
        } else if self.d > other.d{
            Ordering::Greater
        }
        else {
            Ordering::Equal
        }
    }
}
/// Given a set of clusterings (one for each gonzales set), return a minimum spanning tree on the
/// centers, one for each clustering 
fn compute_spanning_trees<M : ColoredMetric>(space : &M, prob : &ClusteringProblem, clusterings : &Vec<Clustering>) -> Vec<RootedSpanningTree> {

    let mut spanning_trees : Vec<RootedSpanningTree> = Vec::with_capacity(prob.k);

    for clustering in clusterings.iter() {
        let i = clustering.m() - 1; // the number of centers; i.e. we consider gonzales set S_i (counting from 0 to k-1)

        // first lets do the algorithm of prim to compute the minimum spanning tree.
        let mut priority_queue: PriorityQueue<CenterIdx,Reverse<Priority>> = PriorityQueue::with_capacity(i+1);

        let mut spanning_tree: Vec<UpEdge> = Vec::with_capacity(i);
        let mut closest: Vec<CenterIdx> = Vec::with_capacity(i+1); // denotes the cloest center that is already in the tree

        // we start with center 0 as singleton spanning tree
        closest.push(0); // closest to 0 is 0;
        for c in 1..i+1 {
            priority_queue.push(c,Reverse(Priority{d : space.dist(clustering.get_center(0), clustering.get_center(c))}));
            closest.push(0);
        }

        while let Some((c1,Reverse(Priority{d:dist_to_closest}))) = priority_queue.pop() {
            spanning_tree.push(UpEdge{down : c1, up : closest[c1], d: dist_to_closest});
            for (&mut c2,Reverse(prio)) in priority_queue.iter_mut() {


                let new_d = space.dist(clustering.get_center(c1), clustering.get_center(c2));
                if new_d < prio.d {
                    prio.d = new_d;
                    closest[c2] = c1;
                }
                
            }
        }
//        println!("Clustering {}: Spanning tree: {:?}", i, spanning_tree);
        spanning_tree.sort_by(|a,b| a.down.cmp(&b.down));
        let mut spanning_tree_edges : Vec<Option<UpEdge>> = Vec::with_capacity(i+1);
        spanning_tree_edges.push(None); // The Root (Center 0) has no UpEdge
        spanning_tree_edges.extend(spanning_tree.into_iter().map(|e| Some(e)));

        spanning_trees.push(RootedSpanningTree{edges : spanning_tree_edges});

        }
    spanning_trees
}

