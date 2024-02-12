/// deletes element from the list, such that only the t smallest elements remain (unsorted)
pub(super) fn truncate_to_smallest<E: Clone + PartialOrd + Send>(list: &mut Vec<E>, t: usize) {
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
    list.select_nth_unstable_by(t, |a, b| a.partial_cmp(b).unwrap());
    list.truncate(t);
}

pub(super) fn split_in_half<E: Clone + PartialOrd + Send>(list: &mut [E]) -> (&mut [E], &mut [E]) {
    split_in_two_at(list, (list.len() - 1) / 2 + 1)
}

pub(super) fn split_in_two_at<E: Clone + PartialOrd + Send>(
    list: &mut [E],
    t: usize,
) -> (&mut [E], &mut [E]) {
    if 0 < t && t < list.len() {
        list.select_nth_unstable_by(t, |a, b| a.partial_cmp(b).unwrap());
        return list.split_at_mut(t);
    }
    #[cfg(debug_assertions)]
    if t > list.len() {
        println!("Care: list has only length {}, but is asked to split off the t = {} smallest elements. Everything is split up.", list.len(), t);
    }
    list.split_at_mut(t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    fn create_random_list_with_dublicates(n: usize, range: usize) -> Vec<usize> {
        let mut rng = rand::thread_rng();
        let vals: Vec<usize> = (0..range).map(|_| rng.gen_range(0..10 * range)).collect();
        (0..n).map(|_| vals[rng.gen_range(0..range)]).collect() // do a list with dublicates
    }

    #[test]
    fn split_in_two_test() {
        let n = 10000;
        let mut list = create_random_list_with_dublicates(n, 1000);
        // create n*k edges with dublicates:
        let mut sorted_list = list.clone();
        let pos = (list.len() - 1) / 2;
        let (smaller, bigger) = split_in_two_at(&mut list, pos);

        sorted_list.sort();
        assert_eq!(smaller.len(), pos);
        assert_eq!(bigger.len(), n - pos);

        let pivot = list[pos];
        assert_eq!(pivot, sorted_list[pos]);

        (0..pos).for_each(|i| {
            assert!(list[i] <= pivot);
        });
        (pos..n).for_each(|i| {
            assert!(list[i] >= pivot);
        });
    }

    #[test]
    fn truncate_to_shortest_test() {
        let mut list = create_random_list_with_dublicates(1000, 1000);
        let mut sorted_list = list.clone();

        let pos = 200;
        truncate_to_smallest(&mut list, pos);

        assert_eq!(list.len(), pos);

        sorted_list.sort();
        sorted_list.truncate(pos);
        list.sort();
        assert_eq!(list, sorted_list);
    }
}
