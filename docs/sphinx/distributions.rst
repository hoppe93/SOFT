Distribution functions
======================
In SOFT, distribution functions depend on three variables, namely the major radius :math:`\rho`
at which the guiding-center orbit was initiated, the momentum :math:`p` of the particle, as well
as the cosine of the pitch angle :math:`\xi = \cos\theta_\mathrm{p}` in the outer midplane.

File format
-----------
SOFT distribution functions are text files containing a set of matrices. Each matrix contains
the momentum-space distribution function at a specific radius, with momentum :math:`p` along
the :math:`x`-axis (i.e. varying with columns) and cosine of the pitch angle :math:`\xi`
along the :math:`y`-axis (i.e. varying with rows). Each matrix corresponds to a particular
major radius :math:`\rho`.

In the beginning of the file the phase-space grid is specified. The first three rows give the
bounds and number of points in major radius :math:`\rho`, :math:`\xi` and momentum :math:`p`
respectively. The next three rows are vectors containing all the grid points of phase-space,
allowing the distribution function to be given on a non-uniform grid. Numbers are always
separated by non-line-breaking white-space characters (i.e. tab or space).

.. note:: Radial position and momentum are given in SI units for historic reasons.

An example of a SOFT distribution function with numbers replaced by human-readable text is shown below. ::

  r0   rn    nr
  xi0  xin   nxi
  p0   pn    np
  r0   r1   r2   r3   r4 ...
  xi0  xi1  xi2  xi3  xi4 ...
  p0   p1   p2   p3   p4 ...
  f(r0,xi0,p0)   f(r0,xi0,p1)   f(r0,xi0,p2) ...
  f(r0,xi1,p0)   f(r0,xi1,p1)   f(r0,xi1,p2) ...
  f(r0,xi2,p0)   f(r0,xi2,p1)   f(r0,xi2,p2) ...
  .
  .
  .
  f(rn,xi0,p0)   f(rn,xi0,p1)   f(rn,xi0,p2) ...
  f(rn,xi1,p0)   f(rn,xi1,p1)   f(rn,xi1,p2) ...
  f(rn,xi2,p0)   f(rn,xi2,p1)   f(rn,xi2,p2) ...
  .
  .
  .

With numbers, we instead have ::

  6.8e-1     9.2e-1    200
  0          1         5
  5.36e-21   2.68e-20  5
  0.68       0.68121   0.68241 ...
  6.1232e-17 0.38268   0.70711 ...
  5.36e-21   1.072e-20 1.608e-20 ...
  0.0010962  0.0010962 0.0010962 ...
  0.0010962  0.0010962 0.0010962 ...
  0.0010962  0.0010962 0.0010962 ...
  .
  .
  .

Helper tools for CODE/NORSE
---------------------------
A nice graphical helper tool has been developed for analyzing CODE/NORSE distributions and
generating distributions readable by SOFT. The tool is called
`codeviz <https://github.com/hoppe93/codeviz>`_ and is available on GitHub.
