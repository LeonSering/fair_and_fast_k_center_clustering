use crate::ClusteringProblem;
use crate::space::{Distance,ColoredMetric};
use crate::clustering::{CenterIdx,Clustering};


/// a pointer from a gonzales center to its parent in a rooted forest 
#[derive(Debug)]
struct UpEdge {
    d : Distance,
    up : CenterIdx, // the parent of down in the rooted forest
    down : CenterIdx, // the child of up in the rooted forest
}

struct RootedSpanningTree {
    edges: Vec<Option<UpEdge>>, // vector of size i; for each node s_j \in S_i it either points to the parent node or is None for the root 
    sorted_distances: Option<Vec<Distance>>, // a list of all distances in edges in increasing order
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

    fn get_sorted_dist_iter(&mut self) -> std::slice::Iter<Distance> {
        if self.sorted_distances.is_some() {
            return self.sorted_distances.as_ref().unwrap().iter()
        }
        let mut distances : Vec<Distance> = self.edges.iter().filter_map(|o| o.as_ref()).map(|e| e.d).collect();
        distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
        self.sorted_distances = Some(distances);
        self.sorted_distances.as_ref().unwrap().iter()
    }
}


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

/// Given a metric space and a clusting problem,
/// redistribute takes a vector of (full) clusterings, each satisfying the privacy constrait
/// and returns for each gonzales set a minimum spanning tree produced by the Kruskal
/// algorithm.
pub fn redistribute<'a, M : ColoredMetric>(space : &M, prob : &ClusteringProblem, clusterings : Vec<Clustering<'a>>) -> Vec<Clustering<'a>> {

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
            priority_queue.push(c,Reverse(Priority{d : space.dist(clustering.centers.get(0), clustering.centers.get(c))}));
            closest.push(0);
        }

        while let Some((c1,Reverse(Priority{d:dist_to_closest}))) = priority_queue.pop() {
            spanning_tree.push(UpEdge{down : c1, up : closest[c1], d: dist_to_closest});
            for (&mut c2,Reverse(prio)) in priority_queue.iter_mut() {


                let new_d = space.dist(clustering.centers.get(c1), clustering.centers.get(c2));
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

        spanning_trees.push(RootedSpanningTree{edges : spanning_tree_edges,  sorted_distances : None});

        println!("spanning_tree[{}]: {:?}", i, spanning_trees[i].edges);
        println!("spanning_tree[{}] sorted distances: {:?}", i, spanning_trees[i].get_sorted_dist_iter().collect::<Vec<_>>());
        }
    clusterings
}

