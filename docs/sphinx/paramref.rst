.. _paramref:

Parameter reference
===================
There are a number of settings that can be specified in a ``pi`` file, and each of the SOFT modules
introduces its own set of options. In this section a complete list of all the options that can be
set in a ``pi`` file are given.

.. contents:: Contents
   :local:

Global options
--------------

.. option:: debug

   | **Default value:** 0
   | **Example line:** ``debug=1``
   | **Allowed values:** 0 or 1

   If set to 1, debug output will be generated and written to ``stdout`` during the run. Default
   value is 0.

.. option:: domain_has_outer_wall

   | **Default value:** yes
   | **Example line:** ``domain_has_outer_wall=no``
   | **Allowed values:** ``yes`` or ``no``

   If set to ``no``, ignores all points of the wall/separatrix outside :math:`R = R_m`, where
   :math:`R_m` denotes the radial coordinate of the magnetic axis. This will allows the placement
   of a detector outside the device. The mid-pole will still be present to block out radiation.

.. option:: interptimestep

   | **Default value:**
   | **Example line:**
   | **Allowed values:**

   TODO

.. option:: magnetic_field

   | **Default value:** None
   | **Example line:** ``magnetic_field=numeric``
   | **Allowed values:** ``circular`` and ``numeric``

   Specifies the name of the magnetic field handler module to use. Either ``circular`` or
   ``numeric``.

.. option:: maxtimestep

   | **Default value:** None
   | **Example line:** ``maxtimestep=1e-11``
   | **Allowed values:** Any positive real value

   Sets the maximum allowed size of a timestep in the equation solver (whichever it may be).
   If the adaptive timestep becomes larger than this, it is automatically adjusted to this
   value. By default there is no limit on how long the timestep can be.

.. option:: nodrifts

   | **Default value:** no
   | **Example line:** ``nodrifts=yes``
   | **Allowed values:** ``yes`` or ``no``

   If set to ``yes``, ignores the drift terms in the first-order guiding-center equations of
   motion (effectively solving the zeroth-order guiding-center equations of motion). This option
   only influences behaviour of the code when the guiding-center equations of motion are solved.
   By default the value of this option is ``no`` so that the drift terms are kept.

.. option:: progress

   | **Default value:** 0
   | **Example line:** ``progress=10``
   | **Allowed values:** Any non-negative integer

   Specifies how many times during the run SOFT should print information about the current progress.
   Information will be printed in uniform steps as particles (defined as points in phase-space) are
   completed.

.. option:: threads

   | **Default value:** Number of threads suggested by OpenMP
   | **Example line:** ``threads=3``
   | **Allowed values:** Any positive integer (no upper limit)

   Overrides the number of threads started by each (MPI) process. By default, SOFT will start
   the number of threads indicated by the ``OMP_NUM_THREADS`` environment variable in each
   process.

.. option:: tolerance

   | **Default value:** 1e-12
   | **Example line:** ``tolerance=4e-13``
   | **Allowed values:** Any positive real number

   Specifices the tolerance in the RKF45 solver. The default tolerance is set by the tool used in
   the run. The ``orbit`` tool defaults to a tolerance of :math:`10^{-7}`, while the ``sycamera``
   defaults to a tolerance of :math:`10^{-12}`.

.. option:: useequation

   | **Default value:** None
   | **Example line:** ``useequation=guiding-center-relativistic``
   | **Allowed values:** ``guiding-center``, ``guiding-center-relativistic``, ``particle``, ``particle-relativistic``.

   Determines which set of equations of motion to solve. Note that the ``sycamera`` tool requires
   that the (relativistic) guiding-center equations of motion be solved. Possible values for this
   option are ``particle``, ``particle-relativistic``, ``guiding-center`` and
   ``guiding-center-relativistic``.

.. option:: usetool

   | **Default value:** None
   | **Example line:** ``usetool=sycamera``
   | **Allowed values:** ``orbit``, ``sycamera``

   Sets the name of the tool to use. Can either be ``orbit`` (which traces orbits), or ``sycamera``
   (which computes various synchrotron-radiation quantities for runaway electrons).

Particle settings
-----------------

.. option:: charge

   | **Default value:** One electron charge (i.e. ``-1``)
   | **Example line:** ``charge=4``
   | **Allowed values:** ``orbit``, ``sycamera``

   The charge of the particle to simulate, in units of the elementary charge (:math:`e = 9.109\times 10^{-19}\,\mathrm{C}`).
   The default value is -1, i.e. the electron charge.

.. option:: cospitch

   | **Default value:** None
   | **Example line:** ``cospitch=1,0.95,100``
   | **Allowed values:** A number :math:`\in [0,1]`; A number :math:`\in[0,1]`; any positive integer

   Specifies the range of cosines of the particle's pitch anle with which to initiate orbits. The
   first argument specifies the first value in the range to give to particles, while the second
   argument argument specifies the last value in the range. The third argument specifies the total
   number of values to simulate. Example: ``cospitch = 0.999,0.97,10``, while initiate ten particles
   with cosine of the pitch angle values between 0.97 and 0.999.

.. option:: gc_position

   | **Default value:** Yes
   | **Example line:** ``gc_position=no``
   | **Allowed values:** ``yes`` or ``no``
   
   If set to ``yes``, assumes that the position given specifies the guiding-center position when
   solving the guiding-center equations of motion. If set to ``no``, the program instead assumes
   that the particle position is specified and compensates accordingly when solving the
   guiding-center equations of motion. Has no effect when solving the full particle orbit.

.. option:: mass

   | **Default value:** One electron mass (:math:`0.000548579909\,\mathrm{u}`)
   | **Example line:** ``mass=2``
   | **Allowed values:** Any positive real number

   The particle mass in unified atomic mass units (u). The default value is 0.000548579909,
   corresponding to the electron mass.

.. option:: p

   | **Default value:** None
   | **Example line:** ``p=1e6,1.2e7,10``
   | **Allowed values:** Any real number; any real number; any positive integer

   Specifies the range of momenta with which to initiate orbits. The first argument specifies
   the first momentum value to give to particles while the second argument specifies the last
   momentum value. The third argument specifies the total number of momentum values to simulate.
   Example: ``p = 3e7,4e7,5``.

.. option:: pitch

   | **Default value:** None
   | **Example line:** ``pitch=0.05,0.15,14``
   | **Allowed values:** A number :math:`\in [0,\pi]`; a number :math:`\in [0,\pi]`; any positive integer

   Specifies the range of pitch angles with which to initiate orbits. The first argument specifies
   the first pitch angle to give to particles while the second argument specifies the last
   pitch angle. The third argument specifies the total number of pitch angles to simulate.
   Example: ``pitch = 0.03,0.25,15``.

.. option:: ppar

   | **Default value:** None
   | **Example line:** ``ppar=1e6,1.2e7,14``
   | **Allowed values:** Any real number; any real number; any positive integer

   Specifies the range of parallel momenta with which to initiate orbits. The first argument specifies
   the first parallel momentum to give to particles while the second argument specifies the last
   momentum value. The third argument specifies the total number of momentum values to simulate.
   Example: ``ppar = 3e7,4e7,5``.

.. option:: pperp

   | **Default value:** None
   | **Example line:** ``pperp=1e6,1.2e7,14``
   | **Allowed values:** Any real number; any real number; any positive integer

   Specifies the range of perpendicular momenta with which to initiate orbits. The first argument specifies
   the first perpendicular momentum to give to particles while the second argument specifies the last
   momentum value. The third argument specifies the total number of momentum values to simulate.
   Example: ``pperp = 3e6,7e6,15``.

.. option:: r

   | **Default value:** None
   | **Example line:** ``r=0.68,0.84,14``
   | **Allowed values:** Any real number inside device; any real number inside device; any positive integer

   Specifies the range of radii with which to initiate orbits. The first argument specifies
   the first radius to give to particles while the second argument specifies the last
   radius. The third argument specifies the total number of radii to simulate.
   Example: ``r = 0.68,0.84,80``.

.. option:: rdyn

   | **Default value:** None
   | **Example line:** ``rdyn=0.84,14``
   | **Allowed values:** Any real number inside device; any positive integer

   Specifies the outermost radius at which to initiate orbits, as well as the number of radii
   to drop particles on. The innermost radius is automatically set as the magnetic axis, and
   particles will only be dropped at a radius in the interval if their "effective magnetic axis"
   radial location is less than the currently simulated. The "effective magnetic axis" arises
   due to orbit drifts, and if it's presence is not properly accounted for, weird bright or
   dark spots will show up in synchrotron image (when orbit drifts are taken into account).
   Example: ``rdyn = 0.84,80``.

.. option:: t

   | **Default value:** ``0,-1``
   | **Example line:** ``t=0,1e-6``
   | **Allowed values:** Any real number; any real number

   The first argument of this parameter specifies the reference time. For most purposes this
   parameter is most conveniently set to 0. The second argument specifies the end time, at
   which point an orbit should be considered finished and no longer followed. If the second
   argument is less than the reference time (the first argument), the orbit will be followed
   for one full *poloidal* orbit, or until the simulation clock is greater than minus the
   end time.

Magnetic settings
-----------------
Two different magnetic handler modules are provided with SOFT. These are the ``circular`` module,
implementing a simple analytical magnetic field with a circular cross-section and constant
safety factor, as well as the ``numeric`` module, which loads 2D numeric magnetic fields.


Performance-wise, the ``numeric`` module is somewhat slower than the ``circular`` model, due to
that the former interpolates the 2D magnetic field with a cubic spline. The difference is however
only about a factor of two.

circular
^^^^^^^^
.. option:: B0

   | **Default value:** ``1``
   | **Example line:** ``B0=5.2``
   | **Allowed values:** Any real number

   Specifies the magnetic field strength on the magnetic axis, i.e. on the circle
   :math:`R = R_{\mathrm{m}}, Z = 0`. In units of Tesla.

.. option:: major_radius

   | **Default value:** ``1``
   | **Example line:** ``major_radius=2``
   | **Allowed values:** Any positive real number

   Specifies the major radius of the tokamak. In units of meter.

.. option:: minor_radius

   | **Default value:** ``1``
   | **Example line:** ``minor_radius=1``
   | **Allowed values:** Any real number

   Specifies the minor radius of the device. In units of meter. This parameter only influences
   the location of the walls of the tokamak, and does not affect the magnetic field.

.. option:: safety_factor

   | **Default value:** ``1``
   | **Example line:** ``B0=1``
   | **Allowed values:** Any real number

   The safety factor, or :math:`q`-factor of the tokamak magnetic field. In this analytical
   model of the magnetic field, the safety factor is a constant.

numeric
^^^^^^^
.. option:: axis

   | **Default value:** *Set in equilibrium file*
   | **Example line:** ``axis=0.68,-0.002``
   | **Allowed values:** Any positive real number; any real number

   Specifies the location of the magnetic axis in a poloidal plane. The first coordinate
   specifies the major radial location (:math:`R`) of the axis, and the second coordinate specifies
   the vertical location (:math:`Z`) of the axis. SOFT requires the magnetic equilibrium
   data file to give this value, but under some circumstances it may be desirable to
   override the value set in the equilibrium file, in which case this parameter can be used.

.. option:: file

   | **Default value:** None
   | **Example line:** ``file=/path/to/magnetic/equilibrium.mat``
   | **Allowed values:** Any real number

   Specifies the name of the file containing the magnetic equilibrium data to use. The
   format that this file must have is described under :ref:`magnetic`. The format of
   the file is determined by analyzing the file name extension. All file formats supported
   by the SOFT file interface can be used.

.. option:: format

   | **Default value:** ``auto``
   | **Example line:** ``format=mat``
   | **Allowed values:** ``auto``, ``hdf5`` or ``mat``

   Overrides the format specifier for the magnetic equilibrium data file. ``auto``
   is the default, which causes SOFT to determine the file format based on the filename
   extension. ``hdf5`` causes SOFT to interpret the data file as an HDF5 file. ``mat``
   causes SOFT to interpret the data file as a Matlab MAT file.

.. option:: wall

   | **Default value:** ``any``
   | **Example line:** ``wall=separatrix``
   | **Allowed values:** ``any``, ``separatrix``, ``wall``

   Specifies which type of wall should be used. Equilibrium data files can contain two types
   of "walls", namely the actual tokamak wall cross-section or the separatrix/last closed flux surface.
   SOFT only requires one of these two types to be present in the data file, and with ``any`` set,
   the tokamak wall will be first be considered, but if it's not present in the file, the separatrix
   will be used instead. The ``wall`` and ``separatrix`` options forces the use of either of
   the two types. *The wall is the structure beyond which particles will be considered as lost
   and no longer followed.*

sycout settings
---------------
A **sycout** (short for *SYnchrotron Camera OUTput*) is an output module that is
coupled to the ``sycamera`` tool of SOFT. Currently the following sycouts are available:

* **green** -- Generates a Green's function
* **image** -- Generates a camera image 
* **space3d** -- Stores 3D information about the contributions to an image
* **spectrometer** -- Generates a spectrum
* **topview** -- Stores X and Y coordinates of contributions to an image. Creates a top-down "map" of contributions.

green
^^^^^
The ``green`` sycout allows you to generate Green's functions for images, spectra or
any kind of function you can imagine. Green's functions are sometimes also known as
weight functions and are essentially mappings from a distribution function to a quantity
such as an image, spectrum or combination thereof.

*Instructions on how to use this sycout are available under :ref:`geomkern`.*

.. option:: format

   | **Default value:** Auto-determined from output filename extension
   | **Example line:** ``format=mat``
   | **Allowed values:** ``h5``, ``hdf5``, ``mat``, ``out``, ``sdt``

   Overrides the default setting for what file format to store the output in.
   If not set, the output file format is determined based on the filename extension
   of the output file. ``h5`` and ``hdf5`` forces HDF5 output. ``mat`` forces
   Matlab MAT output. ``out`` and ``sdt`` forces SOFT self-descriptive text (SDT)
   format output (text-based).

.. option:: function

   | **Default value:** None
   | **Example line:** ``function=r12ij``
   | **Allowed values:** Any (non-repeating) combination of the characters ``1``, ``2``, ``i``, ``j``, ``r``, ``w``

   Sets the shape and contents of the Green's function. A more detailed description
   of how this option works can be found under :ref:`geomkern`.

.. option:: output

   | **Default value:** None
   | **Example line:** ``output=outputfile.mat``
   | **Allowed values:** Any non-line-breaking string

   Sets the name of the output file. The format of the output file is determined based
   on the extension part of this setting unless the ``format`` option has also
   been specified. *By extension is meant everything that comes after the last dot (.).*

.. option:: pixels

   | **Default value:** None
   | **Example line:** ``pixels=520``
   | **Allowed values:** Any positive integer

   Sets the number of pixels of the image, i.e. the number of elements in each of the ``i``
   and ``j`` dimensions. Only required if either ``i`` or ``j`` appears in the ``function`` option.

.. option:: suboffseti
.. option:: suboffsetj

   | **Default value:** 0
   | **Example line:** ``suboffseti=20``
   | **Allowed values:** Any non-negative integer

   Green's functions for images tend to become quite large, and in many cases much of the
   Green's function is zero and provides no interesting information. In these cases, a subset
   of the image can be stored so that the correct wide-angle image distortion is still present.
   These offset parameters specify the offsets in the i and j directions respectively from
   which the image that is to be stored should start.

.. option:: subpixels

   | **Default value:** *Same as ``pixels``*
   | **Example line:** ``subpixels=30``
   | **Allowed values:** Any positive integer

   Specifies the number of pixels in each of the i and j directions of the subset image.
   Since the subset image must lie within the full image, ``suboffseti``+``subpixels`` and
   ``suboffsetj``+``subpixels`` must both be less than or equal to ``pixels``.

image
^^^^^
The ``image`` sycout generates a camera image.

.. option:: brightness

   | **Default value:** ``intensity``
   | **Example line:** ``brightness=histogram``
   | **Allowed values:** ``bw``, ``histogram``, ``intensity``

   Specifies how pixels should be colored. With ``bw`` (for black-and-white), pixels are
   simply marked if they receive a contribution. Thus, if any radiation hits the pixel
   during the run, the pixel will contain the value 1 at the end of the run and 0 otherwise.

   The ``histogram`` option specifies that each hit in a pixel should increase the value
   of the pixel by 1. The radiation intensity reaching the pixel is not considered.

   The ``intensity`` option takes the emitted radiation intensity into account, including
   spectral effects (if enabled through other options).

.. option:: name

   | **Default value:** None
   | **Example line:** ``name=output-file.mat``
   | **Allowed values:** Any string allowed by the underlying file system

   Specifies the name of the file to which the output will be written. The output
   is written through the SOFT file interface which means it will be either in
   a HDF5 file, a Matlab MAT file or a SOFT SDT (Self-Descriptive Text) format.
   The file format is determined based on the filename extension. For HDF5, use
   either *.h5* or *.hdf5*, for Matlab MAT use *.mat*, and for SDT any other extension
   (though *.sdt* is recommended).

.. option:: pixels

   | **Default value:** None
   | **Example line:** ``pixels=300``
   | **Allowed values:** Any positive integer

   Sets the number of pixels in the image. Images are always square and have the same
   number of pixels in the x (i) direction as in the y (j) direction.

space3d
^^^^^^^
The ``space3d`` can be used to store 3D data about the points of space that contribute
to an image. *Will be described in more detail in a separate section*

.. option:: output

   | **Default value:** None
   | **Example line:** ``output=name-of-outputfile.mat``
   | **Allowed values:** Any string allowed by the underlying file system

   Name of the file to which output should be written. The ``space3d`` module
   uses the SOFT file interface, meaning output can be written in either
   HDF5, Matlab MAT or SOFT SDT (Self-Descriptive Text) format. The format
   of the output file is determined based on the filename extension. For HDF5
   use *.h5* or *.hdf5*, for Matlab MAT use *.mat*, and for SDT use any other
   extension (though *.sdt* is recommended).

.. option:: pixels

   | **Default value:** None
   | **Example line:** ``pixels=300``
   | **Allowed values:** Any positive integer

   When ``type=pixels``, sets the number of pixels in each direction of the
   bounding box. A value of for example 100 means that there will be a total
   of :math:`100\times 100\times 100 = 1\,000\,000` "pixels" in the box.

.. option:: point0

   | **Default value:** None
   | **Example line:** ``point0=.40,-.75,.20``
   | **Allowed values:** Any real number; any real number; any real number

   Specifies one of the two defining edge points of the bounding box. 

.. option:: point1

   | **Default value:** None
   | **Example line:** ``point1=.63,-.15,-.20``
   | **Allowed values:** Any real number; any real number; any real number

   Specifies one of the two defining edge points of the bounding box. 

.. option:: type

   | **Default value:** None
   | **Example line:** ``type=pixels``
   | **Allowed values:** ``pixels``, ``real``

   Specifies the type of 3D object to store. ``pixels`` divides the
   bounding box into a number of smaller boxes and collects the contribution
   in each of those (the number of boxes is determined by the ``pixels`` option).
   This 3D type is fixed in size and is represented as a simple 3D array.

   The ``real`` type stores the real location of each particle that contributes
   to the image. This 3D grows in size with the number of particles that hit
   the detector, and is stored as a sparse matrix. It's usually very difficult
   to determine the final size of this 3D type, but it gives much more detailed
   data and can sometimes be the more space-efficient option.

spectrometer
^^^^^^^^^^^^

topview
^^^^^^^

