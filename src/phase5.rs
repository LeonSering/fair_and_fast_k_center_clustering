use crate::ClusteringProblem;
use crate::phase3::{RootedSpanningTree,build_bfs_stack};
use crate::types::{Distance,CenterIdx,PointCount};
use crate::space::{Point,ColoredMetric};
use crate::clustering::{Centers,Clustering};
use crate::utilities;



/// A point-center-pair with computed distance:
#[derive(Debug,Clone,PartialOrd,PartialEq,Copy)]
struct PointCenterLink<'a> {
    dist : Distance,
    point : &'a Point,
    center_idx : CenterIdx,
}

pub(crate) fn phase5<'a, M : ColoredMetric>(space : &'a M, prob : &ClusteringProblem, final_centers: Vec<Centers>, spanning_trees: Vec<RootedSpanningTree>, thresholds: Vec<Distance>) -> Clustering<'a> {



    // TEMP: Create empty clustering
    let mut centers = Centers::with_capacity(prob.k);
    centers.push(space.point_iter().next().unwrap());
    let center_of: Vec<Option<usize>> = space.point_iter().map(|_| None).collect();
    Clustering::new(centers,center_of, space)
}

fn point_shifting<'a, M : ColoredMetric>(space : &'a M, privacy_bound: PointCount, clustering : &mut Clustering<'a>, i : CenterIdx, spanning_tree: &RootedSpanningTree, threshold: Distance){

    let mut cluster_sizes: Vec<PointCount> = clustering.get_cluster_sizes().clone();

    // this stacks contains all nodes in the order discovered by an BFS starting from the root;
    // we push point in the reverse ordering (therefore a stack)
    let mut stack = build_bfs_stack(i, spanning_tree, threshold);
    
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
    println!("Point_center_pairs before truncate:{:?}", point_center_pairs);
    utilities::truncate_to_smallest(&mut point_center_pairs, number);
    println!("Point_center_pairs after truncate:{:?}", point_center_pairs);
    for PointCenterLink{dist:_, point: p, center_idx:_} in point_center_pairs {
        clustering.assign(p, consumer, space);
    }
}
