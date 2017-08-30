Distribution functions
======================
In SOFT, distribution functions depend on three variables, namely the major radius :math:`\rho`
at which the guiding-center orbit was initiated, the momentum :math:`p` of the particle, as well
as the cosine of the pitch angle :math:`\xi = \cos\theta_\mathrm{p}` in the outer midplane.

File format
-----------
Distribution functions are given to SOFT as Matlab MAT-files. SOFT expects the following
variables to be present in the file:

+-----------------+----------------------------------------------------------------------------------------------------------------------+
| Name            | Description                                                                                                          |
+-----------------+----------------------------------------------------------------------------------------------------------------------+
| ``description`` | String describing the distribution function.                                                                         |
+-----------------+----------------------------------------------------------------------------------------------------------------------+
| ``f``           | Actual distribution function. An :math:`n_r`-by-:math:`n_pn_{\xi}` matrix (see below).                               |
+-----------------+----------------------------------------------------------------------------------------------------------------------+
| ``name``        | String naming the distribution function.                                                                             |
+-----------------+----------------------------------------------------------------------------------------------------------------------+
| ``p``           | Vector containing points of momentum. Size 1-by-:math:`n_p`.                                                         |
+-----------------+----------------------------------------------------------------------------------------------------------------------+
| ``punits``      | String describing the units of ``p``. Either ``ev``, ``normalized`` or ``si``.                                       |
+-----------------+----------------------------------------------------------------------------------------------------------------------+
| ``r``           | Vector containing radial points. Size 1-by-:math:`n_r`.                                                              |
+-----------------+----------------------------------------------------------------------------------------------------------------------+
| ``xi``          | Vector containing (cosine of) pitch angle points. Size 1-by-:math:`n_{\xi}`.                                         |
+-----------------+----------------------------------------------------------------------------------------------------------------------+

The most important variable in a SOFT distribution function file is ``f``, which is the actual
distribution function. The variable is stored as a matrix with each row representing a momentum-space
distribution function, i.e. with the radial coordinate changing along the row index.

Each row of ``f`` corresponds to a momentum-space distribution function, shaped as one long
:math:`n_p\times n_{\xi}` vector. The elements are ordered into :math:`n_{\xi}` groups of :math:`n_p`
elements, so that the first :math:`n_p` elements of the vector corresponds to holding :math:`\xi`
fixed and varying :math:`p`.

The ``name`` and ``description`` variables are fairly arbitrary and are only included to provide
the user with basic information about the distribution function.

The ``p``, ``r`` and ``xi`` variables are vectors consisting of :math:`n_p`, :math:`n_r` and :math:`n_{\xi}`
elements respectively. Together, the vectors specify the grid in momentum, radius and (cosine of)
pitch angle on which the distribution function is defined.

To allow users to specify momentum coordinates in the units most convenient for them, and more
importantly to prevent mix-ups of used units, the variable ``punits`` must be provided specifying
the units used for the momentum variable. Allowed values are ``ev`` (for momentum in :math:`\text{eV}/mc`),
``normalized`` (for :math:`p\equiv\gamma\beta`, where :math:`\gamma` is the electron's Lorentz factor and
:math:`\beta` is the electron's speed normalized to the speed of light) and ``si`` (for SI units, i.e.
:math:`\text{kg}\cdot\text{m/s}`).

.. note::
   Even though CODE is commonly used to generate distribution functions for SOFT, plain CODE
   distribution functions are not directly compatible with SOFT. The distribution function
   given as output by CODE consists of a set of Legendre polynomial coefficients used in evaluating the
   distribution function :math:`f(p,\xi)`. SOFT on the other hand requires the function values
   to be already evaluated.

Helper tools for CODE/NORSE
---------------------------
A nice graphical helper tool has been developed for analyzing CODE/NORSE distributions and
generating distributions readable by SOFT. The tool is called
`codeviz <https://github.com/hoppe93/codeviz>`_ and is available on GitHub.
