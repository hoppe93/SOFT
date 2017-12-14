Polarization information
========================
When ``radiation=synchrotron_spectrum`` SOFT will also store information about the
polarization of the detected radiation. Using the ``polimage`` and ``polspectrometer``
sycouts, it is possible to generate output files containing the polarization
information in image or spectrum format. In this section usage and interpretation
of the data will be briefly be discussed.

What information does SOFT store?
---------------------------------
SOFT stores the four `Stokes parameters <https://en.wikipedia.org/wiki/Stokes_parameters>`_,
:math:`S`, :math:`Q`, :math:`U` and :math:`V`, averaged over the relevant parameters
(depending on which model is being used). The emitted synchrotron power per unit frequency,
per unit solid angle, can be written in terms of the two quantities :math:`A_\parallel` and
:math:`A_\perp` as

.. math::
   \frac{\mathrm{d}^2 P}{\mathrm{d}\omega\mathrm{d}\Omega} \propto \left| -\boldsymbol{\epsilon}_\parallel A_\parallel + \boldsymbol{\epsilon}_\perp A_\perp \right|^2.

where :math:`\boldsymbol{\epsilon}_\parallel` is a vector corresponding to polarization in
the gyration plane, and :math:`\boldsymbol{\epsilon}_\perp` to polarization in the plane
orthogonal to that. It can be shown that the Stokes parameters can be expressed using
:math:`A_\parallel` and :math:`A_\perp` through

.. math::
   I &\propto A_\parallel^2 + A_\perp^2,\\
   Q &\propto \left( A_\perp^2 - A_\parallel^2 \right)\cos 2\beta,\\
   U &\propto \left( A_\perp^2 - A_\parallel^2 \right)\sin2\beta,\\
   V &\propto 2A_\parallel A_\perp \cos 2\beta.

The angle :math:`\beta` is the angle between the plane of parallel polarization and
the plane in which the horizontal polarization is *measured*. The first Stokes parameter,
:math:`I`, is just the intensity of the radiation as obtained also from the SOFT
``image`` sycout.

The fourth Stokes parameter :math:`V` is often quoted as identically zero in the literature,
a result stemming from that the object :math:`A_\parallel A_\perp` is odd in the
angle :math:`\psi` between the guiding-center's emission cone and a line-of-sight. When
averaged over all emission angles, the contribution to :math:`V` therefore cancels
identically. In the angular and spectral distribution implemented in SOFT however, we
do not neglect the finite emission width, and therefore obtain a finite contribution to
the :math:`V` parameter, since it is possible for only part of the emission cone to
overlap the detector (corresponding to "cut-offs" in the integration over emission angle).

*For a derivation of the full* :math:`\frac{\mathrm{d}^2 P}{\mathrm{d}\omega\mathrm{d}\Omega}` *,
see for example Jackson's "Electrodynamics", Landau-Lifshitz "The Classical Theory
of Fields" or Mathias Hoppe's Master's thesis* (`link <http://ft.nephy.chalmers.se/publications/Hoppe_Masters_Thesis_Final.pdf>`_).

File format
-----------
The ``polimage`` sycout of SOFT outputs a variable-based file (such as SDT, HDF5 or
Matlab) containing the following variables:

+-----------------------+------------------------------------------------+
| Variable              | Description                                    |
+-----------------------+------------------------------------------------+
| ``detectorPosition``  | Vector specifying the position of the detector |
+-----------------------+------------------------------------------------+
| ``detectorDirection`` | Central viewing direction of the detector      |
+-----------------------+------------------------------------------------+
| ``detectorVisang``    | Vision angle of the detector                   |
+-----------------------+------------------------------------------------+
| ``StokesI``           | Stokes :math:`I` parameter                     |
+-----------------------+------------------------------------------------+
| ``StokesQ``           | Stokes :math:`Q` parameter                     |
+-----------------------+------------------------------------------------+
| ``StokesU``           | Stokes :math:`U` parameter                     |
+-----------------------+------------------------------------------------+
| ``StokesV``           | Stokes :math:`V` parameter                     |
+-----------------------+------------------------------------------------+
| ``wall``              | Wall data used for the simulation              |
+-----------------------+------------------------------------------------+
