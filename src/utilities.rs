pub(super) fn split_at_median<E: Clone + PartialOrd + Copy>(list: &mut Vec<E>) -> (Vec<E>, Vec<E>) {
    let middle = (list.len() - 1 ) /2;
    split_at_pos(list, middle)
}

/// finds the element that would be on position pos in a sorted list.
/// all elements in the left list are smaller or equal to all
/// elements of the right list.
/// list is empty afterwards.
pub(super) fn split_at_pos<E: Clone + PartialOrd + Copy>(list: &mut Vec<E>, pos: usize) -> (Vec<E>, Vec<E>) {
    if pos >= list.len() {
        #[cfg(debug_assertions)]
        println!("Care: list has only length {}, but is asked to be split at pos = {}. Everyting will be in smaller.", list.len(), pos);
        let mut smaller: Vec<E> = Vec::with_capacity(list.len());
        let bigger: Vec<E> = Vec::new();
        smaller.append(list);
        return (smaller, bigger)
    }
    median_of_medians(list, pos);

    let bigger = list.split_off(pos+1);
    let mut smaller: Vec<E> = Vec::with_capacity(list.len());
    smaller.append(list);
    (smaller, bigger)
}

/// deletes element from the list, such that only the t smallest elements remain (unsorted)
pub(super) fn truncate_to_smallest<E: Clone + PartialOrd + Copy>(list: &mut Vec<E>, t: usize){
    if t >= list.len() {
        #[cfg(debug_assertions)]
        println!("Care: list has only length {}, but is asked to be truncated to a length of t = {}. Nothing will happen.", list.len(), t);
        return;
    }
    median_of_medians(list, t);
    list.truncate(t);
}

/// splits of all elements except for the t smallest. Hence, afterwords the list consists of the t
/// smallest elements and the remaining elements are returned.
pub(super) fn split_off_at<E: Clone + PartialOrd + Copy>(list: &mut Vec<E>, t: usize) -> Vec<E>{
    if t >= list.len() {
        #[cfg(debug_assertions)]
        println!("Care: list has only length {}, but is asked to split of the t = {} smallest elements. Everything is split up.", list.len(), t);
//        let mut split: Vec<E> = Vec::with_capacity(list.len());
//        split.append(list);
//        return split;
        return list.split_off(list.len())
    }
    median_of_medians(list, t);
    list.split_off(t)
}

/// input: list of unsorted edges; an integer pos;
/// output: the value of the element that would have pos as index if list was sorted.
fn median_of_medians<E: Clone + PartialOrd + Copy>(list: &mut Vec<E>, pos : usize) -> E {
    assert!(pos <= list.len(), "Cannot compute the element at position pos = {} if the list has only length = {}", pos, list.len());

    let chunks = list.chunks(5);
    let mut sublist: Vec<E> = Vec::with_capacity(list.len()/5+1);
    for chunk in chunks {
        let mut chunk = chunk.to_vec();
        chunk.sort_by(|a, b| a.partial_cmp(b).unwrap());
        sublist.push(chunk[(chunk.len()-1)/2]); // take the median: element with index = floor(length / 2)
    }
    let median: E;
    if sublist.len() <= 5 {
        sublist.sort_by(|a, b| a.partial_cmp(b).unwrap());
        median = sublist[(sublist.len()-1) / 2];
    } else {
        let middle = (sublist.len()-1)/2;
        median = median_of_medians(&mut sublist, middle); 
    }


    let mut left: Vec<E> = Vec::with_capacity((list.len()*7) /10);
    let mut right: Vec<E> = Vec::with_capacity((list.len()*7)/10);
    let mut median_seen: bool = false;
    while !list.is_empty() {
        let i = list.pop().unwrap(); // empties out list
        if i < median {
            left.push(i);
        } else if i > median {
            right.push(i);
        } else { // *i == pivot (if it is the second number with this property, push it to right)
            if median_seen {
                right.push(i);
            } else {
                median_seen = true;
            }
        }
    }
    let mut pivot = median;
    let k = left.len();
    if pos < k {
        pivot = median_of_medians(&mut left,pos);
    } else if pos > k {
        pivot = median_of_medians(&mut right,pos-k-1);
    }

    list.append(&mut left);
    list.push(median);
    list.append(&mut right);
    pivot 
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
        let median = median_of_medians(&mut list, middle);

        // TODO Create two lists to compare

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

        // TODO Create two lists to compare

        clone.sort();
        clone.truncate(t);
        println!("t={} smallest element. List: {:?}; By sorting: {:?}", t,list, clone);
        list.sort();
        assert_eq!(list, clone);
    }
}
