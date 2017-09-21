.. _magnetic:

Magnetic equilibria
===================
There are currently two magnetic handler modules available for SOFT. The ``circular`` handler
implements a simple circular magnetic field with a constant safety factor, and is somewhat faster
than the alternative. The ``numeric`` allows the magnetic field to be loaded from numeric data,
which is interpolated. This handler is often the desired one as it allows complicated magnetic
geometries to be simulated.

Analytic circular
-----------------
The ``circular`` magnetic handler implements the magnetic field:

.. math::
   \boldsymbol{B}(r,\theta) = \frac{B_0}{1-(r/R_m)\cos\theta} \left(
   \frac{r}{q(r)R_m}\hat{\boldsymbol{\theta}} - \hat{\boldsymbol{\phi}} \right)

where :math:`r` is the minor radius, :math:`\theta` the poloidal angle, :math:`B_0` the magnetic
field strength on the magnetic axis (:math:`r = 0`), :math:`R_m` is the major radius, :math:`q`
is the safety factor, :math:`\hat{\boldsymbol{\theta}}` is a unit vector in the poloidal
direction and :math:`\hat{\boldsymbol{\phi}}` is a unit vector in the toroidal direction. While
this formula allows arbitrary *q*-profiles, SOFT currently only implements this magnetic field
with a linear *q*-profile.

The magnetic field shown above has three free parameters, namely the field strength :math:`B_0`,
the tokamak major radius :math:`R_m` and safety factor :math:`q(r) = q_0`. These parameters
must be specified by the user, and are set by specifying the options ``B0``, ``major_radius``
and ``safety_factor`` respectively in the ``magnetic circular`` environment. For SOFT to be able
to determine when a particle escapes confinement and hits the wall, the minor radius of the
device must also be specified. A circular cross section is assumed. All options for the
*circular* magnetic handler are set according to ::

  magnetic circular {
      B0 = 5
      major_radius = 0.68
      minor_radius = 0.22
      safety_factor = 1
  }

Numeric
-------
One of the great strengths of SOFT is that magnetic equilibrias can be specified as numeric data,
allowing complicated magnetic configurations, and in particular, experimentally measured data,
to be plugged into SOFT. Specifying a numeric equilibrium in the ``pi`` file is as simple as ::

  magnetic_handler=numeric
  magnetic numeric {
      name=/path/to/magnetic/equilibrium.mat
  }

Currently, the equilibrium data can be stored in either a HDF5, (MATLAB) MAT or SDT file. Both
HDF5 and MATLAB files can be created easily with user-friendly tools such as Python or MATLAB,
while SDT (for *Semi-Descriptive Text*) is a SOFT-specific text-based format which is likely the
best choice if the magnetic equilibrium is generated using a small C/C++ program which is
difficult to interface with HDF5 or MATLAB.

Since SOFT assumes the magnetic field to be toroidally symmetric, the magnetic field components
in a poloidal plane must be specified. SOFT uses a cylindrical coordinate system for specifying
the magnetic field, so that :math:`\boldsymbol{B} = B_r \hat{\boldsymbol{r}} + B_z\hat{\boldsymbol{z}} + B_\phi \hat{\boldsymbol{\phi}}`,
where :math:`B_r \hat{\boldsymbol{r}}` denotes the component radially out from the point of
symmetry of the tokamak, :math:`B_z\hat{\boldsymbol{z}}` denotes the component in the vertical
direction, and :math:`B_\phi\hat{\boldsymbol{\phi}}` denotes the component in the toroidal
direction, perpendicular to the poloidal plane in which the magnetic field is given.

Variables
^^^^^^^^^
Both HDF5, MATLAB and SDT files have a *variable* concept where data within the file is
named. Because of this, SOFT looks for certain variables in the datasets, loads them and
gives them meaning in the code. The following variables must be present in all SOFT
magnetic equilibrium files:

+----------------+------------------+---------------------------------------------------------------+
| Variable       | Type             | Description                                                   |
+================+==================+===============================================================+
| ``Br``         | m-by-n matrix    | Radial component of magnetic field (radius-by-z).             |
+----------------+------------------+---------------------------------------------------------------+
| ``Bphi``       | m-by-n matrix    | Toroidal component of magnetic field (radius-by-z).           |
+----------------+------------------+---------------------------------------------------------------+
| ``Bz``         | m-by-n matrix    | Vertical component of magnetic field (radius-by-z).           |
+----------------+------------------+---------------------------------------------------------------+
| ``desc``       | String           | A longer description of the equilibrium. Must be present, but |
|                |                  | may be empty.                                                 |
+----------------+------------------+---------------------------------------------------------------+
| ``maxis``      | 1-by-2 vector    | Specifies the location of the magnetic axis in the            |
|                |                  | :math:`(R, z)`-plane.                                         |
+----------------+------------------+---------------------------------------------------------------+
| ``name``       | String           | Name of the equilibrium. Must be present, but may be empty.   |
+----------------+------------------+---------------------------------------------------------------+
| ``r``          | 1-by-m vector    | List radial points in which the components of the magnetic    |
|                |                  | are given.                                                    |
+----------------+------------------+---------------------------------------------------------------+
| ``separatrix`` | 2-by-many vector | List of contour points marking the separatrix in the          |
|                |                  | :math:`(R, z)`-plane.                                         |
+----------------+------------------+---------------------------------------------------------------+
| ``wall``       | 2-by-many vector | List of contour points marking the bounds of the device in    |
|                |                  | the poloidal plane.                                           |
+----------------+------------------+---------------------------------------------------------------+
| ``z``          | 1-by-n vector    | List of vertical points in which the components of the        |
|                |                  | magnetic field are given.                                     |
+----------------+------------------+---------------------------------------------------------------+

.. note:: Only one of the ``separatrix`` and ``wall`` variables is required to be present in the
          equilibrium file. Both may be present, and in that case the domain contour to use can be
          specified as an additional option to the ``numeric`` magnetic handler. By default the
          ``wall`` contour will be used if available.
