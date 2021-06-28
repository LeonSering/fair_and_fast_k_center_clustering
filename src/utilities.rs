pub(super) fn split_at_median<E: Clone + PartialOrd + Copy>(list: &mut Vec<E>) -> (Vec<E>, Vec<E>) {
    let middle = (list.len() - 1 ) /2;
    split_at_pos(list, middle)
}

/// finds the element that would be on position pos in a sorted list.
/// all elements in the left list are smaller or equal to all
/// elements of the right list.
/// list is empty afterwards.
pub(super) fn split_at_pos<E: Clone + PartialOrd + Copy>(list: &mut Vec<E>, pos: usize) -> (Vec<E>, Vec<E>) {
    median_of_medians(list, pos);

    let bigger = list.split_off(pos+1);
    let mut smaller: Vec<E> = Vec::with_capacity(list.len());
    smaller.append(list);
    (smaller, bigger)
}

/// input: list of unsorted edges; an integer pos;
/// output: the value of the element that would have pos as index if list was sorted.
pub fn median_of_medians<E: Clone + PartialOrd + Copy>(list: &mut Vec<E>, pos : usize) -> E {

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
    #[test]
    fn median_test() {
        let n = 100;
        let k = 20;
        // create n*k edges with dublicates:
        let mut rng = rand::thread_rng();
        let vals: Vec<usize> = (0..(n*k/2)).map(|_| rng.gen_range(0..100)).collect();
        let mut list: Vec<usize> = (0..(n*k)).map(|_| vals[rng.gen_range(0..n*k/2)]).collect(); // do a list with dublicates
        let mut clone = list.clone();
        let middle = (list.len() - 1)/2;
        let median = median_of_medians(&mut list, middle);

        // TODO Create two lists to compare

        clone.sort();
        assert_eq!(median, clone[(clone.len()-1)/2]);
        println!("our median: {}; median by sorting: {}", median, list[(list.len()-1)/2]);
    }
}
