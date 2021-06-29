use crate::{PointCount,space::Distance};
use super::Edge;
use crate::utilities;

// Input: unsorted list of edges, upper bound on the size of a bucket
// Output: list of buckets of size <= ceil(4n/k^4) edges; with property of Lemma 3;
pub(super) fn put_into_buckets(mut list: Vec<Edge>, n: PointCount, k: PointCount, power_of_k: u32) -> Vec<Vec<Edge>> {
    let bucket_size_limit = (4*n-1)/(k.pow(power_of_k)) + 1; // ceil(4n/k^z)

    let mut smaller_buckets : Vec<Vec<Edge>> = Vec::with_capacity(list.len()/2 + 1);
    let mut bigger_buckets : Vec<Vec<Edge>> = Vec::with_capacity(list.len()/2 + 1);


//    println!("\nList before: {:?}", list.iter().map(|e| e.d).collect::<Vec<_>>());
    let (smaller, bigger) = utilities::split_at_median(&mut list);
    
//    println!("smaller: {:?}", smaller.iter().map(|e| e.d).collect::<Vec<_>>());
//    println!("bigger: {:?}", bigger.iter().map(|e| e.d).collect::<Vec<_>>());
//    println!("List after: {:?}", list.iter().map(|e| e.d).collect::<Vec<_>>());

    if smaller.len() <= bucket_size_limit {
        smaller_buckets.push(smaller);
    } else {
        smaller_buckets.append(&mut put_into_buckets(smaller, n, k, power_of_k));
    }

    if bigger.len() <= bucket_size_limit {
        bigger_buckets.push(bigger);
    } else {
        bigger_buckets.append(&mut put_into_buckets(bigger, n, k, power_of_k));
    }

    smaller_buckets.append(&mut bigger_buckets);
    smaller_buckets
}


// asserts that there are at most k^5 buckets, each of size at most ceil(4n/k^4); the distance of
// each element of a bucket j is bigger than all elements in bucket < j and smaller than all
// elements in buckets > j.
pub(super) fn assert_buckets_properties(buckets: &Vec<Vec<Edge>>, n: PointCount, k: PointCount, power_of_k: u32) -> bool {

    assert!(buckets.len() <= k.pow(power_of_k+1));

    let size_limit = (4*n-1)/k.pow(power_of_k)+1;
    let mut d_of_last: Distance = <Distance>::MIN;
    for bucket in buckets.iter() {
        assert!(bucket.len() <= size_limit);
        let bucket_of_dist : Vec<Distance> = bucket.iter().map(|x| x.d).collect();
        let mut d_of_current = d_of_last.clone();
        for d in bucket_of_dist.iter() {
            assert!(*d >= d_of_last);
            if *d > d_of_current {
                d_of_current = *d;
            }
        }
        d_of_last = d_of_current;
    }
    true
}

#[cfg(test)]
mod tests {
    use rand::Rng;
    use super::*;
    use crate::space::{Point,Space2D,ColoredMetric};

    #[test]
    fn bucket_test() {
        let n = 100;
        let k = 20;
        // create n*k edges with dublicates:
        let mut rng = rand::thread_rng();
        let vals: Vec<Distance> = (0..(n*k/2)).map(|_| 100.0*rng.gen::<Distance>()).collect();
        let space = Space2D::new_random(n*k);
        let points : Vec<&Point> = space.point_iter().collect();
        let list: Vec<Edge> = (0..(n*k)).map(|i| Edge{
            left: i,
            right: points[i],
            d: vals[rng.gen_range(0..n*k/2)]}).collect(); // do a list with dublicates
        let power_of_k = 2;
        let buckets = put_into_buckets(list, n, k, power_of_k);
        assert!(assert_buckets_properties(&buckets, n, k, power_of_k));
    }
}


