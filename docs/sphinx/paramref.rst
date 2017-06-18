.. _paramref:

Parameter reference
===================
There are a number of settings that can be specified in a ``pi`` file, and each of the SOFT modules
introduces its own set of options. In this section a complete list of all the options that can be
set in a ``pi`` file are given.

Global options
--------------

.. option:: debug

   If set to 1, debug output will be generated and written to ``stdout`` during the run. Default
   value is 0.

.. option:: domain_has_outer_wall

   If set to ``no``, ignores all points of the wall/separatrix outside :math:`R = R_m`, where
   :math:`R_m` denotes the radial coordinate of the magnetic axis. This will allows the placement
   of a detector outside the device. The mid-pole will still be present to block out radiation.

.. option:: interptimestep

   TODO

.. option:: magnetic_field

   Specifies the name of the magnetic field handler module to use. Either ``circular`` or
   ``numeric``.

.. option:: maxtimestep

   Sets the maximum allowed size of a timestep in the equation solver (whichever it may be).
   If the adaptive timestep becomes larger than this, it is automatically adjusted to this
   value. By default there is no limit on how long the timestep can be.

.. option:: nodrifts

   If set to ``yes``, ignores the drift terms in the first-order guiding-center equations of
   motion (effectively solving the zeroth-order guiding-center equations of motion). This option
   only influences behaviour of the code when the guiding-center equations of motion are solved.
   By default the value of this option is ``no`` so that the drift terms are kept.

.. option:: threads

   Overrides the number of threads started by each (MPI) process. By default, SOFT will start
   the number of threads indicated by the ``OMP_NUM_THREADS`` environment variable in each
   process.

.. option:: tolerance

   Specifices the tolerance in the RKF45 solver. The default tolerance is set by the tool used in
   the run. The ``orbit`` tool defaults to a tolerance of :math:`10^{-7}`, while the ``sycamera``
   defaults to a tolerance of :math:`10^{-12}`.

.. option:: useequation

   Determines which set of equations of motion to solve. Note that the ``sycamera`` tool requires
   that the (relativistic) guiding-center equations of motion be solved. Possible values for this
   option are ``particle``, ``particle-relativistic``, ``guiding-center`` and
   ``guiding-center-relativistic``.

.. option:: usetool

   Sets the name of the tool to use. Can either be ``orbit`` (which traces orbits), or ``sycamera``
   (which computes various synchrotron-radiation quantities for runaway electrons).
