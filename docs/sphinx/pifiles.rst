.. _sec:HowToRun:

How to run SOFT
===============
All configuration of a SOFT run is done in a separate configuration file, commonly referred to
as a ``pi`` file. In this section the basic structure of a ``pi`` file will be explained in detail.
For detailed information about which options can be set, pleaes consult the
`SOFT manual <https://github.com/hoppe93/SOFT/master/docs/manual/Manual.pdf>`_.

Examples
--------
The best way to learn how to set up run scripts for SOFT is to see examples of such run scripts.
A basic ``pi`` file can look like the following::

  # Basic SOFT pi file
  useequation=guiding-center-relativistic
  usetool=sycamera
  
  # Specify magnetic field
  magnetic_field=circular   # Use analytic magnetic field
  magnetic circular { B0=5; major_radius=0.68; minor_radius=0.22; safety_factor=1; }
  domain_has_outer_wall=no  # Remove outer device walls to prevent from blocking radiation

  # Set phase-space
  particles {
      t=0,-1
      rdyn=0.84,1000
      p=3e7,3e7,1
      pitch=0.15,0.15,1
  }

  # Specify properties for the sycamera tool
  tool sycamera {
      aperture=0.006                 # Side length (in m) of (square) detector
      cone=delta                     # Use the cone model (not full angular distribution)
      direction=0,1,0                # Normal vector of detector surface (not necessarily normalized)
      position=0,-1.069,0            # Position vector of detector, relative tokamak point of symmtetry
      product=image                  # Output a synchrotron image when done
      radiation=synchrotron_spectrum # Take spectrum of radiation into account
      spectrum=5e-7,1e-6             # Detector spectral range
      toroidal_resolution=3500       # Number of steps in toroidal integral
      vision_angle=2.0               # Size of field-of-view
  }

  # Specify properties for the 'image' sycout
  sycout image {
      pixels=1000
      name=image.dat
  }

The settings available for SOFT are many more, and for a detailed list of which settings are
available, please consult the
`SOFT manual <https://github.com/hoppe93/SOFT/master/docs/manual/Manual.pdf>`_. Further examples
of ``pi`` files for different purposes are:

* `distpi <http://ft.nephy.chalmers.se/~hoppe/soft/examples/distpi>`_ -- Illustrates how SOFT can
  be run together with a runaway distribution function.
* `hollowpi <http://ft.nephy.chalmers.se/~hoppe/soft/examples/hollowpi>`_ -- An example of
  simulating a hollow electron beam.
* `simplepi <http://ft.nephy.chalmers.se/~hoppe/soft/examples/simplepi>`_ -- The basic example
  shown above, setting just the most important options.
* `orbitpi <http://ft.nephy.chalmers.se/~hoppe/soft/examples/orbitpi>`_ -- Shows how to use the
  orbit following part of SOFT to simulate particle orbits.

Basic syntax
------------
Options in a ``pi`` file are specified by first giving the name of the option, followed by an equal
sign, followed by the value to assign to the option. White-space around the equal sign is ignored.
Typically, everything between the equal sign and the end-of-line marker is considered part of the
assigned value, except for any white-space coming either directly after the equal sign, or directly
before the end-of-line-marker. It is however possible to put several settings on the same line by
separating them with semi-colons (``;``).

Comments can be given by preceding the comment text with a hashtag symbol (``#``). Any text
following the hashtag on the same line will be ignored. Note that comments *cannot* be ended
with a semi-colon.

Some options should be assigned vectors of data, such as the ``direction`` and ``rdyn`` options
(among others) in the above example. Each component of the vector must be separated by a comma,
and any white-space surrounding commas is ignored. Note that all floating-point numbers can be
specified using either decimal form (i.e. ``1000`` or ``0.68``) or C scientific notation
(i.e. ``5e-7``).

Environments
------------
Some options in SOFT are considered *global* and are specified directly in the file, such as for
example ``useequation`` and ``usetool`` in the example above. Many options are however specific
to certain modules of SOFT, and they are instead specified inside the appropriate option
*environment*.

There are four different environments in SOFT, all of which are syntactically similar. With the
exception of the ``particles`` environment (which really just sets what could be considered
global options), they are also conceptually similar.

The ``magnetic``, ``tool`` and ``sycout`` environments specify options for a particular SOFT
module, and the name of the module must be specified in the environment *header*. The settings
are then wrapped within curly brackets (``{`` and ``}``) and given to the specified module.
Note that even if an environment for a module is present in the configuration file, it does
not mean that the module will automatically used. Other options must be set to enable modules.

The basic syntax for an environment ``environment`` configuring the module named ``module`` is::

  environment module {
      ...
  }

The ``particles`` environment does not require any module name to specified.

magnetic
^^^^^^^^
The ``magnetic`` environment specifies settings for the magnetic equilibrium to use, as well as
the surrounding walls. Currently, there are two different so called *magnetic handler* modules
that can be used. The first and simplest is the ``circular`` magnetic handler which implements
simple analytic circular magnetic field with a constant safety factor. The second magnetic
handler, named ``numeric``, allows the specification of a magnetic field numerically from for
example an HDF5 or MATLAB \*.mat file.

particles
^^^^^^^^^
The ``particles`` environment sets a number of options relating to the phase-space of the run.
Since it is necessarily tied to the *particles* module of SOFT, the module name part of the 
environment specification given above should be omitted.

In addition to specifying the bounds of and number of points in phase-space, the ``particles``
environment can also be used to specify a different mass or charge of the simulated particle
species.

.. note:: The ``orbit`` tool for tracing particle orbits only allows simulating a single point of
          phase-space at a time, and can otherwise give rise to some very anonymous errors.

tool
^^^^
The ``tool`` environment sets the options for particular tool. A tool, in SOFT, is a module which
receives information about a computed orbit and processes it. Currently, there are two tools in
SOFT, and these are the ``orbit`` and ``sycamera`` tools. The ``orbit`` tool simply traces a
particle or guiding-center orbit, keeps track of a few addiational parameters, and outputs it all
to a CSV file.

The ``sycamera`` tool is the synchrotron camera (or rather detector) tool which gives SOFT its name.
A large part of the SOFT code is dedicated to this module, and the options set by this tool include
for example the type of synchrotron radiation model to use, the number of toroidal steps to take,
various detector properties among many other things.

sycout
^^^^^^
Due to the great versatility of the ``sycamera`` tool, the types of output that could be obtained
it are numerous. Since each of the output types requires its own set of settings, a separate
environment for specifying settings to the output handler of the ``sycamera`` tool was created.

The ``sycout`` environment thus specifies settings of a ``sycamera`` output handler module. To
date there are five different *sycout* modules, namely

+--------------+----------------------------------------------------------------------------------+
| Module name  | Description                                                                      |
+==============+==================================================================================+
| green        | Generates a *Green's function* which relates the distribution of runaways to the |
|              | resulting spectrum or image. (Can) allow fast computation of image/spectrum.     |
+--------------+----------------------------------------------------------------------------------+
| image        | Generates a synthetic synchrotron image.                                         |
+--------------+----------------------------------------------------------------------------------+
| space3d      | Stores 3D information about all particles contributing to a synchrotron image    |
|              | allows visualizing the corresponding surface-of-visibility.                      |
+--------------+----------------------------------------------------------------------------------+
| spectrometer | Generates a spectrum curve.                                                      |
+--------------+----------------------------------------------------------------------------------+
| topview      | Stores information about where particles where located in the *xy*-plane when    |
|              | when they emitted towards the detector. Allows visualizing the toroidal          |
|              | distribution of particles that are visible to the detector.                      |
+--------------+----------------------------------------------------------------------------------+
