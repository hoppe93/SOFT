Compiling
=========
SOFT is written in C, and as such is straightforward to setup on a Linux system. While SOFT hasn't
been tested on any other system, it should be possible compile and run on for example Windows and
Mac with some additional effort.

Dependencies
------------
SOFT depends on a number other technologies, some of which are required for compilation, while
others can be compiled in optionally. Technologies that are absolutely mandatory in order to
compile SOFT are

* `CMake <https://cmake.org/>`_, for preparing necessary build files.
* A C compiler with OpenMP support (such as `gcc <https://gcc.gnu.org/>`_).
* `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_, for various mathematical
  operations. If a version of GSL older than 2.0 is used, the GSL extension
  `interp2d <https://github.com/diazona/interp2d>`_ must also be installed.

A number of libraries are also optional for compilation, and can be compiled in for additional
functionality. The optional libraries are

* `HDF5 <https://www.hdfgroup.org/HDF5/>`_ for reading/writing data in HDF5 format.
* `MATLAB <https://www.mathworks.com/products/matlab.html>`_, for reading/writing data in
  MATLAB's \*.mat format.
* An MPI library, such as `MPICH <https://www.mpich.org/>`_ or
  `OpenMPI <https://www.open-mpi.org/>`_. Compiling in support for MPI allows running SOFT across
  multiple computers, such as on a supercomputer cluster.

Obtaining the code
------------------
You may clone the latest build from the `SOFT GitHub repository <https://github.com/hoppe93/SOFT>`_
via the command line::

  $ git clone https://github.com/hoppe93/SOFT.git

or if you have your ssh keys configured with GitHub::

  $ git clone git@github.com:hoppe93/SOFT.git

Compiling
---------
Once the SOFT source code has been obtained and all required and desired dependencies have been
installed, navigate to the directory cloned from GitHub::
  
  $ cd SOFT

Next, to compile SOFT, create a ``build`` directory, navigate to it, run CMake followed by make,
using the following set of commands::

  $ mkdir build
  $ cd build
  $ cmake ../ -DUSE_HDF5=ON -DUSE_MATLAB=ON -DUSE_MPI=OFF
  $ make

If the build was successful, the SOFT binary will be found under ``build/src/soft``. The flags
starting with ``-D`` specify configuration options, and in the command above we see that in this
case SOFT would be configured with HDF5 and MATLAB support, but without MPI support. This is the
default, and would have happened even if those flags were not specified. To enable/disable
compilation for either of these libraries, simply specify ``ON``/``OFF`` as appropriate in the above.

Usage
-----
All configuration of a SOFT run is done in a separate script file, commonly referred to as a
``pi`` file (for *Particle Information*). As such, running SOFT is as simple as ::

  $ ./soft pi

assuming the ``pi`` file has been setup appropriately. There are a large number of options that can
be specified in the ``pi`` file, and for this reason the details of using SOFT are left to the
:ref:`sec:HowToRun`.
