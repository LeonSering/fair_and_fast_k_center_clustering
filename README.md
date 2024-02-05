# Code of ICML submission "Fair and Fast k-Center Clustering for Data Summarization"

This is our Rust implementation for our Priv-Rep-kC algorithm that runs in $O(nk^2 + k^5)$ time.

The code is a Rust library that implements a python interface, leading to the following two options to run our algorithm:

- It can be imported into another Rust project.
- After installing the python wheel, the algorithm can be executed from any python script or jupyter notebook.

As the python interface is most convenient, we focus on this way to run our code.

## Requirements

- python >=3.6
- matplotlib (optional)

## Installation of the python wheel

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install the python-wheel:

```bash
pip install ff_k_center-1.2.2-cp36-abi3-manylinux_2_12_x86_64.manylinux2010_x86_64.whl
```

Moreover, we recommend to install matplotlib to get plots of a computed clustering:

```bash
pip install matplotlib
```

## Usage via python interface

The following python code presents all functions of the python interface.

```python
from ff_k_center import FFKCenter
import numpy as np

model = FFKCenter(4,1,[(1,13)]) # parameters: k, privacy_bound = 1, rep_intervals = []
pos = np.array([[1,2],[0,0],[6,7],[0,-1],[-2,-3],[0,-1.2],[1,3],[6,1]])
colors = np.array([0,0,0,1,1,0,0,1])


model.fit(pos, colors, verbose = 2, thread_count = 6, phase_2_rerun = False, phase_5_gonzalez = True) # positions of points; followed by colors; Optional: verbose (0: silent, 1: brief, 2: verbose; default: 1); thread_count (default: #cores); phase_2_rerun (default: True); phase_5_gonzalez (default: False);

print("\nr^*: ", model.radius) # the final radius of the assignment
print("C^*: ", model.centers) # Chosen centers by the point-index.
print("phi^*: ", model.assignment) # for each point the center it is assigned to.
print("running time: ", model.running_time) # running time in sec
```

You can also use the jupyter notebook ```showcase/showcase.ipynb``` for a playground with a toy example for testing.

Alternatively, you can use ```experiments/computational_study.ipynb``` for recreating our experiments. (Warning: As we created the plots in the paper by iterating over many different parameters values this can take several hours, even though a single computation is fast.) The data sets used in the manuscript are included and can be found in ```experiments/data/```.

By running ```python_scripts/create_2d_space.py``` it is possible to easy create new toy examples by left clicking into the figure. Use number keys '0' to '9' to change the color and 'Del' to go into deletion mode. The data set is saved after closing the window.

## Uninstall python wheel

Use the following command to uninstall the python wheel:

```bash
pip uninstall ff_k_center
```

## Building the python wheel by compiling the code

There is no need to recreate the python wheel as it is included in the provided code.
Nevertheless, to be able to verify that the python wheel indeed corresponds to the actual code, we provide information here on how it can be compiled from the Rust code.
This can be done via the tool [maturin](https://github.com/PyO3/maturin).

The most convenient way is to use the docker image [konstin2/maturin](https://hub.docker.com/r/konstin2/maturin) (this may take a couple of minutes):

```bash
docker run --rm -v $(pwd):/io konstin2/maturin build --manylinux=2010 --release
```

Afterwards the python wheel can be found in ```target/wheels/```.

## Usage as Rust library

Make sure that the provided Rust crate is listed as a ```dependency``` in your Cargo.toml file with properly indicated path.
You can load all functionality to a new Rust project by

```rust
extern crate ff_k_center;
use ff_k_center::*;
```

The code below shows an example that creates a data set of 5000 random points on which Priv-Rep-kC is solved with our algorithm.

```rust
extern crate ff_k_center;
use ff_k_center::*;

fn main() {
    let n = 5_000;
    let k = 8;
    let space = SpaceND::new_random(n);
    let prob = ClusteringProblem{
        k, // number of centers;
        privacy_bound : n/k, // number of points to represent;
        rep_intervals : vec!((0,500),(0,500),(0,500)), // representation interval [a,b] for each color class; for color classes without interval we subsitute [0. inf]
    };

    let (clustering, total_time) = compute_privacy_preserving_representative_k_center(&space, &prob, None);
    println!("Radius: {}; Running time: {}s.", clustering.get_radius(), total_time);
}
```

For more details on how to load a data set, save the computed clustering, set optional parameters, and more, we refer to the Rust documentation of the provided library, which can be created and opened via:

```bash
cargo doc --open
```
