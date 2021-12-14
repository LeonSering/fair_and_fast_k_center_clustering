// pub(super) fn split_at_median<E: Clone + PartialOrd>(list: &mut Vec<E>) -> (Vec<E>, Vec<E>) {
    // let middle = (list.len() - 1 ) /2;
    // split_at_pos(list, middle)
// }

/// finds the element that would be on position pos in a sorted list.
/// all elements in the left list are smaller or equal to all
/// elements of the right list.
/// list is empty afterwards.
// pub(super) fn split_at_pos<E: Clone + PartialOrd>(list: &mut Vec<E>, pos: usize) -> (Vec<E>, Vec<E>) {
    // if pos >= list.len() {
        // #[cfg(debug_assertions)]
        // println!("Care: list has only length {}, but is asked to be split at pos = {}. Everything will be in smaller.", list.len(), pos);
        // let mut smaller: Vec<E> = Vec::with_capacity(list.len());
        // let bigger: Vec<E> = Vec::new();
        // smaller.append(list);
        // return (smaller, bigger)
    // }
    // median_of_medians(list, pos);

    // let bigger = list.split_off(pos+1);
    // let mut smaller: Vec<E> = Vec::with_capacity(list.len());
    // smaller.append(list);
    // (smaller, bigger)
// }

/// deletes element from the list, such that only the t smallest elements remain (unsorted)
pub(super) fn truncate_to_smallest<E: Clone + PartialOrd + Send>(list: &mut Vec<E>, t: usize){
    if t == 0 {
        list.clear();
        return;
    }
    if t >= list.len() {
        #[cfg(debug_assertions)]
        {
            if t > list.len() {
                println!("Care: list has only length {}, but is asked to be truncated to a length of t = {}. Nothing will happen.", list.len(), t);
            }
        }
        return;
    }
    let (m, _) = median_of_medians(list, t-1);
    list.push(m);
}

/// splits off the half largest elements. Hence, afterwords the list consists of the half
/// smallest elements and the larger half of elements is returned.
/// For list of a odd number of elements the list.len() = return_list.len()+1;
/// i.e., the median remains in the list.
pub(super) fn split_off_at_median<E: Clone + PartialOrd + Send>(list: &mut Vec<E>) -> Vec<E>{
    let middle = (list.len() - 1)/2; //note that we take the lower median
    split_off_at(list, middle+1) //splits of all larger elements (median itself remains in list)
}

/// splits off all elements except for the t smallest. Hence, afterwords the list consists of the t
/// smallest elements and the remaining elements are returned.
pub(super) fn split_off_at<E: Clone + PartialOrd + Send>(list: &mut Vec<E>, t: usize) -> Vec<E>{
    if t == 0 {
        let mut bigger = Vec::with_capacity(list.len());
        bigger.append(list);
        return bigger;
    }
    if t >= list.len() {
        #[cfg(debug_assertions)]
        if t > list.len() {
            println!("Care: list has only length {}, but is asked to split off the t = {} smallest elements. Everything is split up.", list.len(), t);
        }
        return list.split_off(list.len())
    }
    let (m, bigger) = median_of_medians(list, t-1);
    list.push(m);
    bigger
}

use rayon::prelude::*;
/// input: list of unsorted edges; an integer pos;
/// output: 
/// the value of the element m that would have pos as index if list was sorted.
/// a list of (n-pos-1)-many elements that are larger than (or equal to) m
/// the provided list is reduced to pos-many the elements smaller than (or equal to) m
///
/// Note: The element m is neither in the provided list nore in the returned list.
fn median_of_medians<E: Clone + PartialOrd + Send>(list: &mut Vec<E>, pos : usize) -> (E, Vec<E>) {
    assert!(pos <= list.len(), "Cannot compute the element at position pos = {} if the list has only length = {}", pos, list.len());


    if list.len() < 5 {
        list.par_sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        let bigger = list.split_off(pos+1);
        let m = list.pop().unwrap();
        return (m,bigger); 
    }

    let chunks = list.chunks(5);
    let mut sublist: Vec<E> = Vec::with_capacity(list.len()/5+1);
    for chunk in chunks {
        let mut chunk = chunk.to_vec();
        chunk.par_sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        sublist.push(chunk[(chunk.len()-1)/2].clone()); // take the median: element with index = floor(length / 2)
    }
    
    let middle = (sublist.len()-1)/2;
    let (m, _) = median_of_medians(&mut sublist, middle);
    let median = m;


    // let mut left: Vec<E> = Vec::with_capacity((list.len()*7) /10);
    // let mut right: Vec<E> = Vec::with_capacity((list.len()*7)/10);
    // let mut median_seen: bool = false;

   
    // lets reorganize the list in place into two halves: Left smaller than median, and right bigger than median.
    let mut i = 0;
    let mut j = list.len()-1;

    let mut i_med = 0; // if the list is only 1 element long it has to be the median
    while i < list.len() && list[i] <= median {
        if list[i] == median {
            i_med = i;
        }
        i += 1;
    }
    while j > 0 && list[j] > median {
        j -= 1;
    }

    while i < j {
        list.swap(i, j);
        // println!("swap: i: {}, j: {}", i, j);
        while i < list.len() && list[i] <= median {
            if list[i] == median {
                i_med = i;
            }
            i += 1;
        }
        while j > 0 && list[j] > median {
            j -= 1;
        }
    }
    // println!("after: {}, i: {}, j: {}", list.len(),i, j);


    let mut right = list.split_off(i);
    list.swap(i_med, i-1);
    let mut pivot = list.pop().unwrap();
    


    let k = list.len();
    let final_right: Vec<E>;
    if pos < k {
        let (m, mut bigger) = median_of_medians(list,pos);
        pivot = m;
        bigger.push(median);
        bigger.extend(right);
        final_right = bigger;
    } else if pos > k {
        let (m, bigger) = median_of_medians(&mut right,pos-k-1);
        pivot = m;
        final_right = bigger;
        list.extend(right);
        list.push(median);
    } else {
        final_right = right;
    }
    (pivot, final_right)
}



#[cfg(test)]
mod tests {
    use rand::Rng;
    use super::*;

    fn create_random_list_with_dublicates(n: usize, range: usize) -> Vec<usize> {
        let mut rng = rand::thread_rng();
        let vals: Vec<usize> = (0..range).map(|_| rng.gen_range(0..10*range)).collect();
        (0..n).map(|_| vals[rng.gen_range(0..range)]).collect() // do a list with dublicates
    }

    #[test]
    fn median_test() {
        let mut list = create_random_list_with_dublicates(10000, 1000);
            // create n*k edges with dublicates:
        let mut clone = list.clone();
        let middle = (list.len() - 1)/2;
        let (median,_) = median_of_medians(&mut list, middle);

        clone.sort();
        assert_eq!(median, clone[(clone.len()-1)/2]);
        println!("our median: {}; median by sorting: {}", median, list[(list.len()-1)/2]);
    }

    #[test]
    fn truncate_to_shortest_test() {
        let mut list = create_random_list_with_dublicates(1000, 1000);
        let mut clone = list.clone();

        let t = 200;
        truncate_to_smallest(&mut list, t);

        clone.sort();
        clone.truncate(t);
        println!("t={} smallest element. List: {:?}; By sorting: {:?}", t,list, clone);
        list.sort();
        assert_eq!(list, clone);
    }
}
