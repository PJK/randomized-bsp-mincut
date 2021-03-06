# Communication-Efficient Randomized Minimum Cuts

This is a proof-of-concept implementation of the "Communication-Efficient Randomized Minimum Cuts" paper by

- [Pavel Kalvoda](https://github.com/PJK)
- [Lukas Gianinazzi](https://github.com/glukas)
- [Alessandro de Palma](https://github.com/AleDepo93)

## Overview
The main C++ MPI application is located in `src`. It implements both the sparse and the dense algorithm, as well as the sequential base cases.

We also provide our experimental setup. The inputs were (partly) generated by the scripts in `input_generators`, which use the following software

 - [NetworkX](https://networkx.github.io/) (see below)
 - [PaRMAT](https://github.com/farkhor/PaRMAT)

We have also used the generators provided by Levine along with "Experimental Study of Minimum Cut Algorithms", and the real-world inputs as described in the paper. The large dense graph is generated in memory.
 
The experiments were executed on the [CSCS Dora](http://user.cscs.ch/computing_systems/piz_dora/index.html) Cray XC40. We provide the automation scripts in `experiment_runner`. Please note that these are specific for the system and our FS layout and are only included for completeness. Before executing, be sure to adapt them to your setup.
 
## Building and running
The following dependencies are required:

 - Boost 1.61 
 - C++11 compiler
 - MPI2 libraries and runtime
 - CMake 3.4+

The specific versions of software used in our experiments are detailed in the paper.

Configure and execute the build:
```
buildir=$(mktemp -d ~/tmp_build.XXXX)
pushd $buildir
cmake -DCMAKE_BUILD_TYPE=Release <PATH-TO-SOURCES>
make -j 8
popd
```

The executables will be ready in `$buildir/src/executables`. Of particular interest are the following:

- `karger_stein` -- sequential K-S implementation used in the base case
- `square_root` -- our implementation
-  `transition` -- utility to compute allocation of work among nodes without executing the actual computation

All of the executables are documented with the input parameters and accept inputs generated by our generators or transformed using the provided utilities.

Finally, execute the code using e.g.
```
mpiexec -n 48 <PATH>/src/executables/square_root 0.95 input.in 0
```
or equivalent for your platform. This will print the value of the resulting cut, timing information, and possibly more fine-grained profiling data.

### Generators and utilities:

Interesting, albeit small graphs can be generated using the tools in `utilities`. The `cut` utility can also independently verify correctness (keep in mind that our algorithm is Monte-Carlo and all randomness is controlled by the seed).

#### Setup

Get Python 3, do `pip3 install networkx`

#### Usage

 - `cut`: Cuts a graphs from stdin
 - `generate`: Generates graphs

##### Examples

Complete graph on 100 vertices with weights uniformly distributed between 1 and 42
```
./utils/generate.py 'complete_graph(N)' 100 --weight 42 --randomize
```

Erdos-Renyi graph on 100 vertices with 0.2 edge creation probability
```
./utils/generate.py 'fast_gnp_random_graph(N, P)' 100 --prob 0.2
```

Watts–Strogatz small-world graph with `k = 6, p = 0.4` (short parameters possible)
```
./utils/generate.py 'connected_watts_strogatz_graph(N, K, P)' 100 -p 0.4 -k 6
```

Help
```
./utils/generate.py -h
```

Checking the cut value
```
./utils/generate.py 'cycle_graph(N)' 100 | ./utils/cut.py
```

Also see (the list of available generators)[https://networkx.github.io/documentation/networkx-1.10/reference/generators.html]

## License

Communication-Efficient Randomized Minimum Cuts    
Copyright (C) 2015, 2016  Pavel Kalvoda, Lukas Gianinazzi, Alessandro de Palma

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.