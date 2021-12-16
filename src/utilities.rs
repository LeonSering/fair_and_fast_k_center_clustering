
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
    ordering_split_at(list, t-1, true);
    list.truncate(t);
}


pub(super) fn ordering_split_in_half<'a, E: Clone + PartialOrd + Send>(list: &'a mut [E], be_precise: bool) -> (&'a mut [E], &'a mut [E]) {
    ordering_split_in_two_at(list, (list.len() - 1) / 2 + 1, be_precise)
}

pub(super) fn ordering_split_in_two_at<'a, E: Clone + PartialOrd + Send>(list: &'a mut [E], t: usize, be_precise: bool) -> (&'a mut [E], &'a mut [E]) {
    if 0 < t && t < list.len() {
        ordering_split_at(list, t, be_precise);
    }
    #[cfg(debug_assertions)]
    if t > list.len() {
        println!("Care: list has only length {}, but is asked to split off the t = {} smallest elements. Everything is split up.", list.len(), t);
    }
    list.split_at_mut(t)
}


// pub(super) fn ordering_split_to_get_median<E: Clone + PartialOrd + Send>(list: &mut [E], be_precise: bool) -> E {
    // let middle = (list.len() - 1)/2; //note that we take the lower median
    // ordering_split_to_get(list, middle, be_precise)
// }


pub(super) fn ordering_split_to_get<E: Clone + PartialOrd + Send>(list: &mut [E], pos : usize, be_precise: bool) -> E {
    let (slices, i) = ordering_split_at(list, pos, be_precise);
    slices[i].last().unwrap().clone()
}


pub(super) fn ordering_split_at_median<'a, E: Clone + PartialOrd + Send>(list: &'a mut [E], be_precise: bool) -> (Vec<&'a mut [E]>, usize) {
    let middle = (list.len() - 1)/2; //note that we take the lower median
    ordering_split_at(list, middle, be_precise)
}


// use rayon::prelude::*;
/// input: slice of unsorted elements; 
/// an integer pos;
/// a boolean be_precise: true: the split has to be exactly at pos; false: the split can be off by
/// a bit
///
/// Suppose input slice goes from 0..n.
/// Let m be the value of the element that would have pos as index if the list was sorted.
///
/// output: 
/// * A vector of slices:
/// All elements of a slice are smaller to all elements of all slices with highter index.
/// * An index i:
/// The last element of the slice i is element m. All other elements of this slice are smaller than
/// m. That means that the size of all slices with index <= i combined is pos+1. 
pub(super) fn ordering_split_at<'a, E: Clone + PartialOrd + Send>(list: &'a mut [E], pos : usize, be_precise: bool) -> (Vec<&'a mut [E]>, usize) {
    let n = list.len();
    assert!(pos <= n, "Cannot compute the element at position pos = {} if the list has only length = {}", pos, list.len());

    if list.len() < 5 {
        list.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let (smaller, bigger) = list.split_at_mut(pos+1);
        return (vec!(smaller,bigger), 0); 
    }

    let chunks = list.chunks_mut(5);
    let mut sublist: Vec<E> = Vec::with_capacity(n/5+1);
    for chunk in chunks {
        chunk.sort_by(|a, b| a.partial_cmp(b).unwrap());
        sublist.push(chunk[(chunk.len()-1)/2].clone()); // take the median: element with index = floor(length / 2)
        // sublist.push(median_for_list_of_atmost_five(chunk));
    }
    
    // let pivot_index = (sublist.len()-1)/2;
    let pivot_index = pos/5;
    let pivot = ordering_split_to_get(&mut sublist, pivot_index, false);


    // let mut left: Vec<E> = Vec::with_capacity((list.len()*7) /10);
    // let mut right: Vec<E> = Vec::with_capacity((list.len()*7)/10);
    // let mut median_seen: bool = false;

   
    // lets reorganize the list in place into two parts: Left smaller than pivot (including pivot at last entry), and right bigger than pivot.
    let mut i = 0;
    let mut j = list.len()-1;

    let mut i_med = 0; // if the list is only 1 element long it has to be the median
    while i < list.len() && list[i] < pivot {
        i += 1;
    }
    while j > 0 && list[j] >= pivot {
        if list[j] == pivot {
            i_med = j;
        }
        j -= 1;
        
    }

    while i < j {
        list.swap(i, j);
        while i < list.len() && list[i] < pivot {
            i += 1;
        }
        while j > 0 && list[j] >= pivot {
            if list[j] == pivot {
                i_med = j;
            }
            j -= 1;
        }
    }
    list.swap(i_med, i);
    let (smaller, bigger) = list.split_at_mut(i+1);
   
    let mut slices: Vec<&mut [E]> = Vec::new();
    let slice_index;

    // println!("list.len: {}, smaller.len: {}, bigger.len: {}, pos: {}, i: {}", n, smaller.len(), bigger.len(), pos, i);

    if be_precise {
        let k = smaller.len()-1;
        if pos < k {
            let (slices_of_smaller, i) = ordering_split_at(smaller,pos, true);
            slices.extend(slices_of_smaller);
            slices.push(bigger);
            slice_index = i;
        } else if pos > k {
            let (slices_of_bigger, i) = ordering_split_at(bigger,pos-k-1, true);
            slices.push(smaller);
            slices.extend(slices_of_bigger);
            slice_index = 1 + i;
        } else {
            slices.push(smaller);
            slices.push(bigger);
            slice_index = 0;
        }
    } else { // if not precise we just keep the split as it is
        slices.push(smaller);
        slices.push(bigger);
        slice_index = 0;
    }
    (slices, slice_index)
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
        let median = ordering_split_to_get(&mut list, middle, true);

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
