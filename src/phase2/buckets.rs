use super::Edge;
use crate::types::PointCount;
use crate::utilities;

use rayon::ThreadPoolBuilder;

// Input: unsorted list of edges, upper bound on the size of a bucket
// Output: list of at most k^3 buckets (as mutable slices) of size <= ceil(4n/k^2) edges; with property the property
// that d(e_1) < d(e_2) for all e_1 in bucket[i] and e_2 in bucket[j] for i < j.
pub(super) fn put_into_buckets<'a, 'b>(
    list: &'b mut [Edge<'a>],
    n: PointCount,
    k: PointCount,
    power_of_k: u32,
    thread_count: usize
) -> Vec<&'b mut [Edge<'a>]> {

    let thread_pool = ThreadPoolBuilder::new().num_threads(thread_count).build().unwrap();

    let bucket_size_limit = (4 * n - 1) / (k.pow(power_of_k)) + 1; // ceil(4n/k^z)
    thread_pool.install(move || put_into_buckets_recursion(list, bucket_size_limit))
}


pub(super) fn put_into_buckets_recursion<'a, 'b>(
    list: &'b mut [Edge<'a>],
    bucket_size_limit: PointCount,
) -> Vec<&'b mut [Edge<'a>]> {

    if list.len() <= bucket_size_limit {
        return vec![list];
    }

    let (smaller, bigger) = utilities::split_in_half(list);


    let (mut smaller_buckets, bigger_buckets) = rayon::join(
        move || put_into_buckets_recursion(smaller, bucket_size_limit),
        move || put_into_buckets_recursion(bigger, bucket_size_limit));

    smaller_buckets.extend(bigger_buckets);
    smaller_buckets
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::space::{ColoredMetric, Point, SpaceND};
    use crate::types::Distance;
    use rand::Rng;

    // asserts that there are at most k^3 buckets, each of size at most ceil(4n/k^2); the distance of
    // each element of a bucket j is bigger than all elements in bucket < j and smaller than all
    // elements in buckets > j.
    pub(super) fn assert_buckets_properties(
        buckets: &Vec<&mut [Edge]>,
        n: PointCount,
        k: PointCount,
        power_of_k: u32,
    ) -> bool {
        assert!(buckets.len() <= k.pow(power_of_k + 1));

        let size_limit = (4 * n - 1) / k.pow(power_of_k) + 1;
        let mut d_of_last: Distance = <Distance>::MIN;
        for bucket in buckets.iter() {
            assert!(bucket.len() <= size_limit);
            let bucket_of_dist: Vec<Distance> = bucket.iter().map(|x| x.d).collect();
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

    #[test]
    fn bucket_test() {
        let n = 100;
        let k = 20;
        // create n*k edges with dublicates:
        let mut rng = rand::thread_rng();
        let vals: Vec<Distance> = (0..(n * k / 2))
            .map(|_| 100.0 * rng.gen::<Distance>())
            .collect();
        let space = SpaceND::new_random(n * k);
        let points: Vec<&Point> = space.point_iter().collect();
        let mut list: Vec<Edge> = (0..(n * k))
            .map(|i| Edge {
                left: i,
                right: points[i],
                d: vals[rng.gen_range(0..n * k / 2)],
            })
            .collect(); // do a list with dublicates
        let power_of_k = 2;
        let buckets = put_into_buckets(&mut list, n, k, power_of_k,8);
        assert!(assert_buckets_properties(&buckets, n, k, power_of_k));
    }
}
