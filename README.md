# SOFT
SOFT (for "Synchrotron-detecting Orbit Following Toolkit") is a synthetic synchrotron diagnostic that can be used to study the synchrotron radiation emitted by runaway electrons in tokamaks. By taking the effect of electron orbits into account, SOFT can accurately reproduce synchrotron images and spectra as obtained in experiments, assuming the runaway electron distribution function is known.

![Synthetic synchrotron image](https://github.com/hoppe93/SOFT/raw/master/examples/scone.png "Synthetic synchrotron image")

## Third-party libraries
SOFT utilizes a number of third-party libraries. Some are required for the code to at all compile, while some can be optionally linked to for additional functionality. Software that **must** be installed on your system is

- A compiler supporting OpenMP (such as gcc)
- GNU Scientific Library version 2.0 or later ([https://www.gnu.org/software/gsl/])

Optional software that is recommended, but not required, to compile, is

- [HDF5](https://support.hdfgroup.org/HDF5/), for reading/writing data in this format.
- [Matlab](https://www.mathworks.com/products/matlab.html), for reading/writing data in *.mat* format.
- MPI, for MPI support. Allows running SOFT across multiple computers. See for example [MPICH](https://www.mpich.org/) or [OpenMPI](https://www.open-mpi.org/).

## Compilation
1. The first step to compiling SOFT is to configure the source code by running CMake. This is done by entering the build directory and running cmake:
```bash
$ cd build/
$ cmake ../
```
By default, CMake will look for HDF5 and Matlab, and if they are not found on the system, configuration will fail. To disable, HDF5 and Matlab support, run instead
```bash
$ cmake ../ -DUSE_HDF5=OFF -DUSE_MATLAB=OFF
```
Similarly, MPI support can be *enabled* using
```bash
$ cmake ../ -DUSE_MPI=ON
```

2. Once SOFT is configured, the code can be compiled with
```bash
$ make
```
in the build directory. The resulting binary can be found under build/src/.

3. (Optional) To install SOFT on your system, run
```bash
$ make install
```
from the build directory. This will copy the build/src/soft binary to /usr/bin.

## Documentation
Documentation of SOFT is available in two forms:
1. The official SOFT manual, located under docs/. Documenting the details of implementation in SOFT.
2. The SOFT beginner's guide, found at [https://ft.nephy.chalmers.se/~hoppe/soft/beginner.html]

## Who's behind SOFT?
SOFT is developed by the [Plasma Theory group](https://ft.nephy.chalmers.se/) at the Department of Physics, Chalmers University of Technology.

## How to cite
When refering to SOFT, please cite the SOFT paper
```
[1] M. Hoppe, O. Embréus, A. Tinguely, R. Granetz, A. Stahl and T. Fülöp, "Synthetic synchrotron diagnostics for runaway electrons in tokamaks", In progress.
```
