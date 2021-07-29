use crate::types::{PointCount,Distance,CenterIdx};
use crate::clustering::Centers;
use std::collections::VecDeque;

/// Phase 3 return a vector (one for each gonzales set) of this type (together with the computed
/// MSTs -> input of Phase 5).
/// It is the input of Phase 4.
/// It contains the eta-vector which determines how many centers can be opened in each gonzales
/// cluster.
#[derive(Debug)]
pub(crate) struct OpeningList {
    pub eta : Vec<PointCount>,
    pub forrest_radius : Distance
}

/// Phase 4 returns a vector (one for each gonzales set) of this type.
/// This together with the MSTs of Phase 3 are the input of Phase 5.
/// It contains the forrest and the assignment radius (which was minimized in the sum)
/// for each gonzales cluster it contains the new centers that needs to be open within this cluster
#[derive(Debug)]
pub(crate) struct NewCenters<'a> {
    pub centers : Centers<'a>,
    pub forrest_radius : Distance,
    pub assignment_radius : Distance,
    pub new_centers_of : Vec<Vec<CenterIdx>>,
}

/// Phase 3 return a vector (one for each gonzales set) of this type.
/// Together with the new centers of Phase 4 this is the input of Phase 5.
/// It contains a minimum spanning tree (MST) on the gonzales centers.
/// It has a designated root (always center with index 0) and each edge is pointing towards this
/// root.
#[derive(Debug)]
pub(crate) struct RootedSpanningTree {
    edges: Vec<Option<UpEdge>>, // vector of size i; for each node s_j \in S_i it either points to the parent node or is None for the root 
    edges_to_children: Vec<Vec<UpEdge>>, // for each node we have a list of edges pointing to the children
}

impl RootedSpanningTree{

    pub fn new(edges: Vec<Option<UpEdge>>, edges_to_children: Vec<Vec<UpEdge>>) -> RootedSpanningTree {
        RootedSpanningTree{edges,edges_to_children}
    }

    /// returns a reference to the TreeEdge, if its distance is not greater than the threshold
    pub fn get_edge(&self, node: CenterIdx, threshold : Distance) -> Option<&UpEdge> {
        if let Some(ref edge) = self.edges[node] {
            // node has an up-going edge
            if edge.d  <= threshold {Some(edge)} else {None}
        } else {
            // node has no up-going edge
            None 
        }
    }

    pub fn get_edges(&self, threshold: Distance) -> Vec<&UpEdge> {
        self.edges.iter().filter_map(|e| e.as_ref()).filter(|e| e.d <= threshold).collect()
    }

    pub fn get_edges_to_children(&self, node: CenterIdx, threshold: Distance) -> Vec<&UpEdge> {
        self.edges_to_children[node].iter().filter(|e| e.d <= threshold).collect()
    }

    pub fn get_sorted_dist(&self) -> Vec<Distance> {
        let mut distances : Vec<Distance> = self.edges.iter().filter_map(|o| o.as_ref()).map(|e| e.d).collect();
        distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
        distances
    }

    pub fn build_bfs_stack(&self, threshold: Distance) -> Vec<CenterIdx> {
        let i = self.edges.len()-1;
        
        let mut stack: Vec<CenterIdx> = Vec::with_capacity(i+1); // elements poped here are treated
        // look for roots:
        let mut is_root= vec!(true; i+1);
        for e in self.get_edges(threshold).iter() {
            is_root[e.down] = false;
        }

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


            for &edge in self.get_edges_to_children(node, threshold).iter() {
                let next = edge.down;
                if pushed[next] == false { 
                    pushed[next] = true; 
                    queue.push_back(next);
                }
            }
        }
        stack
    }
}


/// A pointer from a gonzales center to its parent in a rooted forest.
#[derive(Debug,Clone)]
pub(crate) struct UpEdge {
    pub d : Distance,
    pub up : CenterIdx, // the parent of down in the rooted forest
    pub down : CenterIdx, // the child of up in the rooted forest
}