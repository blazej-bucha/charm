# IMPORTANT: To run this program, CHarm must be compiled with the MPFR
# support, that is, using the "--with-mpfr" installation flag.  If this is your
# case, manually delete (or comment out) all lines between the following
# marks:
#
# @@@
#
# lines to be deleted or commented out
#
# !!!


import pyharm as ph


# INPUTS
# =============================================================================
# Shape of the gravitating body
# -----------------------------------------------------------------------------
# Upper integration limit in spherical radius (topography)
# .............................................................................
# Path to spherical harmonic coefficients defining the shape of the Moon
path_shcs_shape = '../../data/input/MoonTopo2600p_to10-tbl.txt'


# Maximum harmonic degree of the lunar topography defined by coefficients in
# "path_shcs_shape"
shape_nmax = 10
# .............................................................................


# Lower integration limit in spherical radius (sphere)
# .............................................................................
# Radius of the reference sphere
rref = 1728000.0
# .............................................................................
# -----------------------------------------------------------------------------


# Density of the gravitating body
# -----------------------------------------------------------------------------
# Path to spherical harmonic coefficients defining the zero- and first-order
# polynomial density coefficients */
path_shcs_density = ['../../data/input/moon-rho0_to10.tbl',
                     '../../data/input/moon-rho1_to10.tbl']


# Maximum harmonic degrees of the polynomial density coefficients in
# "path_shcs_density".
density_nmax = [10, 10]


# Order of the density polynomial
density_order = len(path_shcs_density) - 1
# -----------------------------------------------------------------------------


# Other inputs
# -----------------------------------------------------------------------------
# Newton's gravitational constant ("kg**-1 * m**3 * s**-2")
gc = 6.67430e-11


# Mass of the Moon ("kg")
mass = 7.346e22


# Minimum and maximum topography powers.  In theory, "pmax" should be
# infinitely large, so higher "pmax" means more complete forward modelling (but
# also longer computation times).
pmin, pmax = 1, 10


# Maximum harmonic degree of the output gravitational potential.  Due to the
# non-linear relation between topography and its implied gravitational field,
# the potential series is spectrally unlimited even if "shape_nmax" is finite
# (the only exception is "shape_nmax = 0").  Therefore, this value should be as
# high as reasonably possible in order to mitigate forward modelling errors.
nmax_potential = 100
# -----------------------------------------------------------------------------
# =============================================================================


# =============================================================================
# Read the shape coefficients
shape_shcs = ph.shc.Shc.from_file('tbl', path_shcs_shape, shape_nmax)


# Read the density coefficients
density_shcs = [''] * (density_order + 1)
for i in range(density_order + 1):
    density_shcs[i] = ph.shc.Shc.from_file('tbl', path_shcs_density[i],
                                           density_nmax[i])


# Global gravity forward modelling using 3D density
# -----------------------------------------------------------------------------
potential_global_shcs = ph.gfm.global_density_3d(shape_shcs, shape_nmax, rref,
                                                 density_shcs, density_nmax,
                                                 density_order,
                                                 gc, mass,
                                                 pmin, pmax,
                                                 nmax_potential,
                                                 None, None, None)


# Now we can compute, say, the gravitational potential using
# "potential_global_shcs".  Let's do this in the Gauss--Legendre grid that
# passes above all lunar masses.
grd = ph.crd.PointGridGL(potential_global_shcs.nmax,
                         potential_global_shcs.r + 25000.0)
vgfm = ph.shs.point(grd, potential_global_shcs, potential_global_shcs.nmax)
# -----------------------------------------------------------------------------


# In CHarm, gravity forward modelling using lateral density is fairly similar
# to using 3D density, so we skip it. */


# Global gravity forward modelling using constant density
# -----------------------------------------------------------------------------
# Constant density of the lunar crust
density_const = 2550.0


potential_global_shcs = ph.gfm.global_density_const(shape_shcs, shape_nmax,
                                                    rref,
                                                    density_const,
                                                    gc, mass,
                                                    pmin, pmax,
                                                    nmax_potential)


# You could now similarly synthesize the gravitational potential, this time
# being implied by the constant density model
# -----------------------------------------------------------------------------


# @@@
# Spatially restricted gravity forward modelling of near-zone masses using 3D
# density
# -----------------------------------------------------------------------------
# Minimum and maximum order of the radial derivatives of the output quantity
kmin, kmax = 0, 3


# Radius of the sphere, on which evaluation points reside
r = rref + 25000.0


# Integration radius in radians
psi0 = 0.2


# Order of the potential derivative ("0" for potential, "1" for quantities
# related to gravitational vector elements, and "2" for quantities related
# to gravitational tensor elements)
u = 0


# Order of the potential derivative with respect to the spherical distance
v = 0


# We want to integrate all masses up to distance "psi0" from evalution
# points.  To integrate masses beyond "psi0", set "zone" to
# "ph.gfm.FAR_ZOZE". */
zone = ph.gfm.NEAR_ZONE


# Number of bits to represent significands of floating points numbers used to
# evaluate
# truncation coefficients
nbits = 512


potential_cap_shcs = ph.gfm.cap_density_3d(shape_shcs, shape_nmax, rref,
                                           density_shcs, density_nmax,
                                           density_order,
                                           gc, mass,
                                           pmin, pmax,
                                           kmin, kmax,
                                           r,
                                           psi0,
                                           u, v,
                                           zone,
                                           nbits,
                                           nmax_potential)
# -----------------------------------------------------------------------------


# Spatially restricted gravity forward modelling of near-zone masses using
# constant density
# -----------------------------------------------------------------------------
potential_cap_shcs = ph.gfm.cap_density_const(shape_shcs, shape_nmax, rref,
                                              density_const,
                                              gc, mass,
                                              pmin, pmax,
                                              kmin, kmax,
                                              r,
                                              psi0,
                                              u, v,
                                              zone,
                                              nbits,
                                              nmax_potential)
# -----------------------------------------------------------------------------


# Let's compute some Molodensky's truncation coefficients, which are used
# internally whenever calling the "charm_gfm_cap_density_*" routines.
# -----------------------------------------------------------------------------
q = ph.gfm.cap_q(rref, r, psi0, nmax_potential, pmax, kmin, kmax,
                 density_order, ph.gfm.NEAR_ZONE, ph.gfm.Q00, nbits)
# -----------------------------------------------------------------------------
# !!!


print('Great, all done!')
# =============================================================================
