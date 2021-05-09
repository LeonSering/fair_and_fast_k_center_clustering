struct cluster_problem{
    // metric space
    // number of colors (gamma)
    // number of points n
    // max number of clusters k
    // for each color j: upper and lower bounds for centers a_j, b_j
    // representative loweder bound L
    //
    // method: test if valid: color classes is metric space are 0, ..., gamma-1; sum of a_j <= k.
    //
    // Solving algorithm as method? Maybe not.
}

struct space{
    // distance matrix (later maybe implicit as distance function to avoid n^2 space)
    // points as array or implicit? If implicity: array of color-classes.
}

struct clustering{
    // array of centers / center indices
    // radius
    // assignment of points to centers
    // method: cluster points of given color
}

fn main() {

    // phase 1: Gonzalez heuristic
    // Output: S[0] to S[k]
    
    // phase 2: determine privacy radius
    
    // phase 3: modify assignment, s.t. sizes are multiple of L
    
    // phase 4: open new centers and  determine actual set of centers C





}
