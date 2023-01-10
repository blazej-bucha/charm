=====================================
Some definitions to make things clear
=====================================


Notation summary
================

.. list-table:: Notation
   :header-rows: 1
   :widths: 20 80

   * - Symbol
     - Meaning

   * - :math:`\varphi \in \big[-\frac{\pi}{2}, \frac{\pi}{2} \big]`
     - Spherical latitude

   * - :math:`\theta \in \big[0, \pi \big]`
     - Spherical co-latitude

   * - :math:`\lambda \in \big[0, 2\pi \big)`
     - Spherical longitude

   * - :math:`r > 0`
     - Spherical radius

   * - :math:`n = 0, 1, 2, \dots`
     - Spherical harmonic degree

   * - :math:`m = 0, 1, 2, \dots, n`
     - Spherical harmonic order

   * - :math:`n_{\max} = \mathbb{N}_0`
     - Maximum spherical harmonic degree of the expansion

   * - :math:`P_n`
     - Un-normalized Legendre polynomial of degree :math:`n`

   * - :math:`\bar{P}_{nm}`
     - Fully-normalized associated Legendre function of the first kind of
       degree :math:`n` and order :math:`m`

   * - :math:`P_{nmj}`
     - Fourier coefficient of the fully-normalized associated Legendre function
       of degree :math:`n` and order :math:`m`. :math:`j` is a coefficient
       related to the wave-number :math:`k` as :math:`j =
       \left[\frac{k}{2}\right]_{\mathrm{floor}}`

   * - :math:`\bar{Y}_{nm}`
     - Real-valued :math:`4\pi`-fully-normalized surface spherical harmonic of
       degree :math:`n` and order :math:`m`

   * - :math:`\{ \bar{C}_{nm},\, \bar{S}_{nm} \}`
     - Real-valued :math:`4\pi`-fully-normalized surface spherical harmonic
       coefficients of degree :math:`n` and order :math:`m`

   * - :math:`\mu`
     - Scalling parameter of spherical harmonic coefficients (for instance, the
       geocentric gravitational constant)

   * - :math:`R > 0`
     - Radius of the reference sphere to which spherical harmonic coefficients
       refer to.

   * - :math:`f`
     - Real-valued function on the unit sphere

   * - :math:`\tilde{f}`
     - Mean value of a real-valued function :math:`f` over a grid cell

   * - :math:`\sigma`
     - Unit sphere

   * - :math:`\Delta \sigma`
     - Area of a grid cell on the unit sphere


.. _points_cells:

Evaluation points and evaluation cells
======================================

CHarm can work with evaluation points and evaluation cells.

* An evaluation **point** is given by the spherical latitude, :math:`\varphi
  \in [-\frac{\pi}{2}, \frac{\pi}{2}]`, the spherical longitude, :math:`\lambda
  \in [0, 2\pi)`, and the spherical radius, :math:`r > 0`.

* An evaluation **cell** is given the minimum and the maximum latitude,
  :math:`\varphi^{\mathrm{min}}` and :math:`\varphi^{\mathrm{max}}`,
  respectively, the minimum and the maximum longitude,
  :math:`\lambda^{\mathrm{min}}` and :math:`\lambda^{\mathrm{max}}`,
  respectively, and a *single* spherical radius :math:`r` that is *constant*
  over the cell.

Occasionally, spherical co-latitude, :math:`\theta = \frac{\pi}{2} - \varphi`,
:math:`\theta \in [0, \pi]`, will be used instead of the spherical latitude.

.. note::
   Angular inputs/outputs (if any) must always be provided in radians in CHarm.

Evaluation points/cells can either be organized as a grid or as scattered
points/cells.

* **Grid of points** is defined by

  * ``nlat`` latitudes over a single grid meridian,

  * ``nlon`` longitudes over a single grid latitude parallel, and

  * ``nlat`` spherical radii over a single grid meridian.

  The grid has in total ``nlat * nlon`` points.  The spherical radius ``r`` may
  vary with the latitude, but is implicitly consider as constant in the
  longitudinal direction (hence the ``nlat`` number of radii).

  CHarm does not store all ``nlat * nlon`` grid coordinates, but only the grid
  boundaries to save some memory.

* **Grid of cells** is defined by

  * ``nlat`` minimum cell latitudes over a single grid meridian,

  * ``nlat`` maximum cell latitudes over a single grid meridian,

  * ``nlon`` minimum longitudes over a single grid latitude parallel,

  * ``nlon`` maximum longitudes over a single grid latitude parallel, and

  * ``nlat`` spherical radii over a single grid meridian.

  The grid has in total ``nlat * nlon`` cells.  Given that the spherical radius 
  (1) is implicitly constant over a grid cell and (2) may vary only with 
  latitude, there is ``nlat`` spherical radii.

  Again, CHarm stores only the grid boundaries.

* **Scattered points** are defined by

  * ``nlat == nlon`` latitudes,

  * ``nlat == nlon`` longitudes, and

  * ``nlat == nlon`` spherical radii.

  The total number of scattered points is ``nlat == nlon``.

* **Scattered cells** are defined by

  * ``nlat == nlon`` minimum latitudes,

  * ``nlat == nlon`` maximum latitudes,

  * ``nlat == nlon`` minimum longitudes,

  * ``nlat == nlon`` maximum longitudes,

  * ``nlat == nlon`` spherical radii.

  The total number of scattered cells is ``nlat == nlon``.

.. note::

   Routines implementing the definitions of evaluation points/cells are
   gathered in :ref:`charm_crd` and :ref:`pyharm_crd`.


Surface spherical harmonics
===========================

Real-valued surface spherical harmonics :math:`\bar{Y}_{nm}(\varphi, \lambda)`
of degree :math:`n` and order :math:`m` are defined as (e.g., Hofmann-Wellenhof
and Moritz, 2005)

.. math::
   \bar{Y}_{nm}(\varphi, \lambda) = \bar{P}_{nm}(\sin \varphi)
   \begin{cases}
   \cos(m\, \lambda)\, {,}\\
   \sin(m\, \lambda)\, {.}
   \end{cases}

Here,

.. math::
   \bar{P}_{nm}(\sin \varphi) =
   \begin{cases}
   \sqrt{(2n + 1)} \, P_n(\sin\varphi){,} &m = 0 {,}\\
   \sqrt{2 (2n + 1) \dfrac{(n - m)!}{(n + m)!}} \, \left(1 -
   \sin^2\varphi\right)^{m / 2} \, \dfrac{\mathrm{d}^m
   P_n(\sin\varphi)}{\mathrm{d} (\sin\varphi)^m} {,} \quad &0 < m \leq n {,}
   \end{cases}

are the fully-normalized associated Legendre functions of the first kind and

.. math::
   P_n(\sin\varphi) = \dfrac{1}{2^n \, n!} \, \dfrac{\mathrm{d}^n}{\mathrm{d}
   (\sin\varphi)^n} \left(\sin^2\varphi - 1 \right)^n

are the (un-normalized) Legendre polynomials (:math:`m = 0`, so the order is
omitted from the notation).

.. note::
   Applied is the geodetic :math:`4\pi` full normalization.  Neither other
   normalizations nor complex spherical harmonics are supported (yet?).

.. note::
   The numerical evaluation of Legendre functions is performed after Fukushima
   (2012), so spherical harmonics can be safely evaluated up to high degrees
   and orders (tens of thousands and even well beyond).


Spherical harmonic analysis
===========================

Assume a harmonic function :math:`f(r, \varphi, \lambda)` given on a sphere
with the radius :math:`r`.  By surface spherical harmonic analysis,

.. math::
   \left.\begin{aligned}
   \bar{C}_{nm} \\
   \bar{S}_{nm}
   \end{aligned}\right\}
   = \dfrac{1}{4 \pi} \, \dfrac{R}{\mu} \, \left( \dfrac{r}{R} \right)^n
   \displaystyle\iint_{\sigma} f(r, \varphi, \lambda) \,
   \bar{Y}_{nm}(\varphi,\lambda) \, \mathrm{d}\sigma {,}

it is possible to compute its spherical harmonic coefficients
:math:`\{ \bar{C}_{nm},\, \bar{S}_{nm} \}`.  The coefficients are normalized by
the :math:`\mu` constant and, if :math:`r \neq R`, they are
additionally rescaled **from** the data sphere with the radius :math:`r` **to**
the reference sphere with the radius :math:`R`.

* If :math:`r = R = \mu = 1`, one arrives at the surface spherical harmonic
  analysis that is well-known from the literature.

* If :math:`r = R`, then :math:`f` does not even have to be harmonic.

* In geosciences, :math:`\mu` frequently represents the geocentric
  gravitational constant and :math:`R` stands for the equatorial radius of the
  Earth.

If :math:`f` is band-limited (that is, it has a finite spherical harmonic
expansion) and is sampled at suitable grid points, these equations can be
computed rigorously (analytically).  Examples of such quadrature, employed in
CHarm, are the Gauss--Legendre quadrature (Sneeuw, 1994) and the
Driscoll--Healy quadrature (Driscoll and Healy, 1994).

In CHarm, the coefficients can also be computed from mean values
:math:`\tilde{f}` of :math:`f` given over grid cells.  In that case, however,
the quadratures are no longer exact.

.. note::

   Routines for spherical harmonic analysis are gathered in
   :ref:`charm_sha` and :ref:`pyharm_sha`.


Spherical harmonic synthesis
============================

Assume that surface spherical harmonic coefficients :math:`\{ \bar{C}_{nm},\,
\bar{S}_{nm} \}` of a harmonic function :math:`f(r, \varphi,\lambda)` are
available up to degree :math:`n_{\mathrm{max}}` and are scaled to a constant
:math:`\mu` and a reference sphere with the radius :math:`R`.  Then, it is
possible to reconstruct a point value of :math:`f(r, \varphi,\lambda)` for any
:math:`(r > R, \varphi,\lambda)` *exactly* by solid spherical harmonic
synthesis,

.. math::

   \displaystyle f(r, \varphi,\lambda) = \frac{\mu}{R} \, \sum_{n
   = 0}^{n_{\max}} \left( \frac{R}{r} \right)^{n + 1} \, \sum_{m = 0}^{n}
   \left( \bar{C}_{nm}\, \cos(m \, \lambda) + \bar{S}_{nm} \, \sin(m \,
   \lambda) \right) \, \bar{P}_{nm}(\sin\varphi){.}

If the evaluation points :math:`(r, \varphi,\lambda)` form a grid (as defined
in :ref:`points_cells`), highly efficient FFT-based algorithms can be employed
(e.g., Colombo, 1981; Sneeuw, 1994; Jekeli et al, 2007; Rexer and Hirt, 2015).
CHarm takes advantage of these algorithms in order to achieve efficient
grid-wise numerical computations.

In addition to point values of :math:`f`, CHarm computes also mean values of
:math:`f` over cells (as defined in :ref:`points_cells`):

.. math::

   \displaystyle \tilde{f}(r, \varphi_{\mathrm{min}},\varphi_{\mathrm{max}},
   \lambda_{\mathrm{min}},\lambda_{\mathrm{max}}) = \frac{1}{\Delta \sigma}
   \int\limits_{\varphi = \varphi_{\mathrm{min}}}^{\varphi_{\mathrm{max}}}
   \int\limits_{\lambda = \lambda_{\mathrm{min}}}^{\lambda_{\mathrm{max}}} f(r,
   \varphi,\lambda) \, \mathrm{d}\lambda \, \cos\varphi \, \mathrm{d} \varphi

and

.. math::

   \displaystyle \tilde{f}(r(\varphi, \lambda), \varphi_{\mathrm{min}},
   \varphi_{\mathrm{max}}, \lambda_{\mathrm{min}},\lambda_{\mathrm{max}}) =&
   \frac{1}{\Delta \sigma} \int\limits_{\varphi
   = \varphi_{\mathrm{min}}}^{\varphi_{\mathrm{max}}} \int\limits_{\lambda
   = \lambda_{\mathrm{min}}}^{\lambda_{\mathrm{max}}} f(r(\varphi,\lambda),
   \varphi,\lambda) \\
   &\times \mathrm{d}\lambda \, \cos\varphi \, \mathrm{d} \varphi {,}

where :math:`\Delta\sigma` is the size of the cell on the unit sphere.

Note that in the latter equation, :math:`f(r(\varphi,\lambda),\varphi,\lambda)`
is defined on an irregular surface given by a spherical radius
:math:`r(\varphi,\lambda)` that varies with the latitude and longitude.  This
computation is unique to CHarm and cannot be found in any other publicly
available package or library.

.. note::

   Routines for spherical harmonic synthesis are gathered in
   :ref:`charm_shs` and :ref:`pyharm_shs`.


References
==========

* Colombo OL (1981) Numerical methods for harmonic analysis on the
  sphere. Report No. 310, Department of Geodetic Science and Surveying, The
  Ohio State University, Columbus, Ohio, 140 pp

* Driscoll, J. R., Healy, D. M. (1994) Computing Fourier transforms and
  convolutions on the 2-sphere. Advances in Applied Mathematics 15:202-250

* Fukushima T (2012) Numerical computation of spherical harmonics of arbitrary
  degree and order by extending exponent of floating point numbers. Journal of
  Geodesy 86:271--285, doi: 10.1007/s00190-011-0519-2

* Hofmann-Wellenhof B, Moritz H (2005) Physical Geodesy. Springer, Wien, New
  York, 403 pp

* Jekeli C, Lee JK, Kwon JH (2007) On the computation and approximation of
  ultra-high-degree spherical harmonic series. Journal of Geodesy 81:603--615,
  doi: 10.1007/s00190-006-0123-z

* Rexer M, Hirt C (2015) Ultra-high-degree surface spherical harmonic analysis
  using the Gauss--Legendre and the Driscoll/Healy quadrature theorem and
  application to planetary topography models of Earth, Mars and Moon. Surveys
  in Geophysics 36:803--830, doi: 10.1007/s10712-015-9345-z

* Sneeuw N (1994) Global spherical harmonic analysis by least-squares and
  numerical quadrature methods in historical perspective. Geophysical Journal
  International 118:707--716
