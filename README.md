# sw-ising
Ising model simulation using the Swendsen-Wang algorithm.

## Getting Started

The `RunSW` script demonstrates basic usage of the core `SwendsenWangIsing` function.
The `RunSWParallel` script gives an example of how the [MATLAB Parallel Computing Toolbox](https://www.mathworks.com/help/distcomp/) can be used to accelerate simulation of multiple temperatures.

## Prerequisites

This code was tested in MATLAB R2017b.
As the [`graph`](https://www.mathworks.com/help/matlab/ref/graph.html) and [`conncomp`](https://www.mathworks.com/help/matlab/ref/graph.conncomp.html) functions were first introduced in R2015b, that is likely the oldest version of MATLAB in which this code will run.
The simulation code itself - `SwendsenWangIsing` - depends only on the MATLAB standard library. The `RunSW` and `RunSWParallel` scripts depend on the utilities in the `utils` folder, and `RunSWParallel` additionally uses [`parfor`](https://www.mathworks.com/help/distcomp/parfor.html) from the [MATLAB Parallel Computing Toolbox](https://www.mathworks.com/help/distcomp/).

## Usage

`SwendsenWangIsing` takes five input arguments:
1. `N`: The total dimensionality of the system.
2. `T`: The temperature of the system.
3. `J`: The `N` by `N` interaction matrix defining the quadratic energy form. Currently, only non-negative interactions, which make it energetically favorable for spins to align, are supported. This implementation is designed for the case in which `J` is sparse, that is, that each node interacts only with a small number of other nodes.  
3. `nIter`: The number of iterations for which the simulation is run.
4. `displayIter`: An optional argument indicating the interval at which to print updates to the terminal. If `displayIter` is not specified, it defaults to zero, and no updates are displayed.

Once the simulation is complete, `SwendsenWangIsing` returns up to five output arguments:
1. `E`: An `nIter` by 1 vector of the energy at each iteration.
2. `M`: An `nIter` by 1 vector of the magnetization at each iteration.
3. `x`: An `N` by 1 vector of the final spin state.
4. `PRNGState`: The initial state of the MATLAB pseudo-random number generator, including the seed.
5. `X`: An `nIter` by `N` matrix containing the spin states at every iteration. If only four output arguments are requested, only the current spin state is stored, increasing the speed and memory efficiency of the simulation.

The `RunSW` and `RunSWParallel` scripts demonstrate the usage of `SwendsenWangIsing`, including examples of how to construct the inputs and plot the outputs.

## References

1. Swendsen, Robert H., and Jian-Sheng Wang. "Nonuniversal critical dynamics in Monte Carlo simulations." _Physical Review Letters_ 58, no. 2 (1987): 86. [(Link)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.58.86)

2. Wang, Jian-Sheng, and Robert H. Swendsen. "Cluster Monte Carlo algorithms." _Physica A: Statistical Mechanics and its Applications_ 167, no. 3 (1990): 565-579. [(Link)](https://www.sciencedirect.com/science/article/pii/037843719090275W)

3. Gilbert, John R., Cleve Moler, and Robert Schreiber. "Sparse matrices in MATLAB: Design and implementation." _SIAM Journal on Matrix Analysis and Applications_ 13, no. 1 (1992): 333-356. [(Link)](https://www.mathworks.com/help/pdf_doc/otherdocs/simax.pdf)


## License

This project is licensed under the [MIT License](LICENSE.txt).
