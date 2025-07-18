import pyharm as ph
import numpy as np


# INPUTS
# =============================================================================
# Path to the input file with spherical harmonic coefficients in the 'gfc'
# format
shcs_file = '../../data/input/EGM96-degree10.gfc'
# =============================================================================


# =============================================================================
# Now read all coefficients in "shcs_file"
shcs = ph.shc.Shc.from_file('gfc', shcs_file)


# Compute the Gauss--Legendre point grid for the "nmax" given by the
# maximum degree stored in "shcs" and for a radius equal to the reference
# sphere, to which the spherical harmonic coefficients are scaled to.  We
# intentionally use here the Gauss--Legendre grid instead of the
# Driscoll--Healy grids in order to avoid the inaccuracies due to the
# singularities at the poles (see the documentation).
grd = ph.crd.PointGridGL(shcs.nmax, shcs.r)


# Compute the full first-order gradient (vector)
fx, fy, fz = ph.shs.point_grad1(grd, shcs, shcs.nmax)


# Compute the magnitude of the gravitational acceleration (no contribution due
# to the centrifugal force is considered here)
g = np.sqrt(fx**2 + fy**2 + fz**2)


# Now compute the full second-order gradient (tensor)
fxx, fxy, fxz, fyy, fyz, fzz = ph.shs.point_grad2(grd, shcs, shcs.nmax)


# Check whether "fxx + fyy + fzz" is zero within numerical errors
print(f'The largest error of the gravitational tensor trace is '
      f'{np.abs(fxx + fyy + fzz).max()} s**-2')


print('Great, all done!')
# =============================================================================
