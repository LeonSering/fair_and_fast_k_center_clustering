use super::Edge;

// Input: unsorted list of edges, upper bound on the size of a bucket
// Output: list of buckets of size <= ceil(4n/k^4) edges; with property of Lemma 3;
pub fn put_into_buckets(mut list: Vec<Edge>, n: usize, k: usize, power_of_k: u32) -> Vec<Vec<Edge>> {
    let bucket_size_limit = (4*n-1)/(k.pow(power_of_k)) + 1; // ceil(4n/k^z)

    let mut smaller_buckets : Vec<Vec<Edge>> = Vec::with_capacity(list.len()/2 + 1);
    let mut bigger_buckets : Vec<Vec<Edge>> = Vec::with_capacity(list.len()/2 + 1);

    let (smaller, bigger) = split_at_median(&mut list);
//    let median = median_of_medians(&list, (list.len() - 1) / 2);
//    
//    let mut smaller : Vec<Edge> = Vec::with_capacity(list.len()/2 + 1);
//    let mut bigger : Vec<Edge> = Vec::with_capacity(list.len()/2 + 1); 
//    
//    while list.len() > 0 {
//        let edge = list.pop().unwrap();
//        if edge <= median {
//            smaller.push(edge);
//        } else { // note that ties in edge.d are broken according to edge.left then edge.right 
//            bigger.push(edge);
//        }
//    }

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

// split bucket at the median into (left, right). All elements in left are smaller or equal to all
// elements on the right.
pub fn split_at_median<'a>(list: &mut Vec<Edge<'a>>) -> (Vec<Edge<'a>>, Vec<Edge<'a>>) {
    let median = median_of_medians(&list, (list.len() - 1) / 2);
    let mut smaller : Vec<Edge> = Vec::with_capacity(list.len()/2 + 1);
    let mut bigger : Vec<Edge> = Vec::with_capacity(list.len()/2 + 1); 
    
    for edge in list.iter(){
        if *edge <= median {
            smaller.push(*edge);
        } else { // note that ties in edge.d are broken according to edge.left then edge.right 
            bigger.push(*edge);
        }
    }
    (smaller, bigger)
}

// input: list of unsorted edges; an integer pos;
// output: the value of the element that would have pos as index if list was sorted.
fn median_of_medians<'a>(list: & Vec<Edge<'a>>, pos : usize) -> Edge<'a> {

    let chunks = list.chunks(5);
    let mut sublist: Vec<Edge> = Vec::with_capacity(list.len()/5);
    for chunk in chunks {
        let mut chunk = chunk.to_vec();
        chunk.sort_by(|a, b| a.partial_cmp(b).unwrap());
        sublist.push(chunk[(chunk.len()-1)/2]); // take the median: element with index = floor(length / 2)
    }
//    println!("sublist: {:?}", sublist);
    let pivot: Edge;
    if sublist.len() <= 5 {
        sublist.sort_by(|a, b| a.partial_cmp(b).unwrap());
        pivot = sublist[(sublist.len()-1) / 2];
    } else {
        pivot = median_of_medians(&sublist, (sublist.len()-1)/2); 
    }

//    println!("Pivot: {:?}", pivot);

    let mut left: Vec<Edge> = Vec::with_capacity((list.len()*7) /10);
    let mut right: Vec<Edge> = Vec::with_capacity((list.len()*7)/10);
    let mut pivot_seen: bool = false;
    for i in list {
        if *i < pivot {
            left.push(*i);
        } else if *i > pivot {
            right.push(*i);
        } else { // *i == pivot (if it is the second number with this property, push it to right)
            if pivot_seen {
                right.push(*i);
            } else {
                pivot_seen = true;
            }
        }
    }
    let k: usize = left.len();
//    println!("list.len(): {}, pos: {}, left.len(): {}, right.len(): {}", list.len(), pos, left.len(), right.len());

    if pos < k {
        return median_of_medians(&left,pos);
    } else if pos > k {
        return median_of_medians(&right,pos-k-1);
    }
//    println!("Final Pivot:{:?}", pivot);
    pivot
}

// asserts that there are at most k^5 buckets, each of size at most ceil(4n/k^4); the distance of
// each element of a bucket j is bigger than all elements in bucket < j and smaller than all
// elements in buckets > j.
pub fn assert_buckets_properties(buckets: &Vec<Vec<Edge>>, n: usize, k: usize, power_of_k: u32) -> bool {

    assert!(buckets.len() <= k.pow(power_of_k+1));

    let size_limit = (4*n-1)/k.pow(power_of_k)+1;
    let mut d_of_last: f32 = <f32>::MIN;
    for bucket in buckets.iter() {
        assert!(bucket.len() <= size_limit);
        let bucket_of_dist : Vec<f32> = bucket.iter().map(|x| x.d).collect();
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
    use crate::space::{Point,new_space_with_random_2dpoints};
    #[test]
    fn median_test() {
        let n = 100;
        let k = 20;
        // create n*k edges with dublicates:
        let mut rng = rand::thread_rng();
        let vals: Vec<f32> = (0..(n*k/2)).map(|_| 100.0*rng.gen::<f32>()).collect();
        let space = new_space_with_random_2dpoints(n*k);
        let points : Vec<&Point> = space.point_iter().collect();
        let mut list: Vec<Edge> = (0..(n*k)).map(|i| Edge{
            left: i,
            right: points[i],
            d: vals[rng.gen_range(0..100)]}).collect(); // do a list with dublicates
        let median = median_of_medians(&list, (list.len() - 1)/2);

        list.sort_by(|a, b| a.partial_cmp(b).unwrap());

        assert_eq!(median, list[(list.len()-1)/2]);
        //println!("list: {:?}", list);
        println!("our median: {}; median by sorting: {}", median.d, list[(list.len()-1)/2].d);
    }

    #[test]
    fn bucket_test() {
        let n = 100;
        let k = 20;
        // create n*k edges with dublicates:
        let mut rng = rand::thread_rng();
        let vals: Vec<f32> = (0..(n*k/2)).map(|_| 100.0*rng.gen::<f32>()).collect();
        let space = new_space_with_random_2dpoints(n*k);
        let points : Vec<&Point> = space.point_iter().collect();
        let list: Vec<Edge> = (0..(n*k)).map(|i| Edge{
            left: i,
            right: points[i],
            d: vals[rng.gen_range(0..100)]}).collect(); // do a list with dublicates
        let size_limit = (4*n)/k;
        let buckets = put_into_buckets(list, size_limit);
        assert!(assert_buckets_properties(&buckets, n, k));
    }
}


