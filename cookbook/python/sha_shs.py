import numpy as np
import pyharm as ph


# INPUTS
# =============================================================================
# Define the path to an input text file with spherical harmonic coefficients.
# For details on the structure of the text file, see the documentation.
shcs_in_file = '../../data/input/EGM96-degree10-mtx.txt'


# Maximum harmonic degree of coefficients to read.  The same degree is used
# later in the harmonic synthesis and harmonic analysis.
nmax = 10
# =============================================================================


# =============================================================================
print('Closed-loop experiment -- Spherical harmonic synthesis and analysis')
print('===========================')


# Read spherical harmonic coefficients from the input text file.
shcs = ph.shc.Shc.from_file('mtx', shcs_in_file, nmax)


# Now let's say we do not want to use the zero-degree term.  This can be
# achieved simply by setting the respective "C00" coefficient to zero.
shcs.set_coeffs(n=0, m=0, c=0.0)


# Compute the Gauss--Legendre grid for a given "nmax" on a sphere that is
# "1000.0" metres above the reference sphere of spherical harmonic
# coefficients.  PyHarm represents the Gauss--Legendre grid via the
# "PointGridGL" class.  For the Driscoll--Healy grids, PyHarm has the
# "PointGridDH1" and "PointGridDH2" classes (see the documentation).
grd_pnt = ph.crd.PointGridGL(nmax, shcs.r + 1000.0)


# Perform the synthesis.  Since the Gauss--Legendre grid resides "1000.0"
# metres above the reference sphere of the coefficients, performed is
# solid spherical harmonic synthesis.
f = ph.shs.point(grd_pnt, shcs, nmax)


# Note on parallelization
# -----------------------
#
# If you compiled the library with OpenMP parallelization enabled, you can
# control the number of threads being used either by an environmental variable
# or in your Python code, for instance, like this:
#
# from threadpoolctl import threadpool_limits
# with threadpool_limits(limits=2):
#     f = ph.shs.point(grd_pnt, shcs, nmax)


# Print some synthesized values
print('Print some synthesized values of the signal...')
i = 0
j = 0
print(f'f({i}, {j}) = {f[i, j]}')


i = grd_pnt.lat.size // 2
j = grd_pnt.lat.size // 2
print(f'f({i}, {j}) = {f[i, j]}')
print('')


# Now use the synthesized signal and compute back its spherical harmonic
# coefficients by harmonic analysis.  The output coefficients in "shcs2" should
# be the same as the input ones in "shcs".  Note that the signal "f" resides on
# a sphere that is "1000.0" metres above the reference sphere "shcs.r".  The
# last input parameter therefore allows user to define the reference radius of
# the output coefficients.  In this case, this value is set to "shcs.r" to
# ensure that "shcs2.r" will be the same as "shcs.r".
shcs2 = ph.sha.point(grd_pnt, f, nmax, shcs.mu, shcs.r)


# Now check whether "shcs" and "shcs2" are the same by computing their
# difference degree amplitudes
dda = shcs.dda(shcs2)


# Print some difference degree amplitudes
n = 0
print('Now print the difference degree amplitudes.  These should be'
      'very small, say, at the order of 1e-18 or less...')
print(f'Difference degree amplitude for harmonic degree {n} = {dda[n]}')
n = 4
print(f'Difference degree amplitude for harmonic degree {n} = {dda[n]}')
n = 10
print(f'Difference degree amplitude for harmonic degree {n} = {dda[n]}')


# We will not need some of the variables anymore, so let's free the memory by
# deleting them
del grd_pnt, shcs2, f, dda


print('===========================')
# =============================================================================


# =============================================================================
print('')
print('Solid spherical harmonic synthesis at scattered points')
print('===========================')


# Create a "PointSctr" class instance to store three scattered points.  The
# instance will be created from three numpy arrays that hold the coordinates of
# the scattered points.
lat      = np.array([0.1, 0.436231, -0.9651], dtype=np.float64)
lon      = np.array([0.0, 1.53434, 4.2316], dtype=np.float64)
r        = np.array([1000.0, 2000.0, 3000.0], dtype=np.float64) + shcs.r
sctr_pnt = ph.crd.PointSctr.from_arrays(lat, lon, r)


# Do the synthesis
f = ph.shs.point(sctr_pnt, shcs, nmax)


# It's really this easy!


# Print the synthesized values
print(f'{f}')


del sctr_pnt, f
print('===========================')
# =============================================================================


# =============================================================================
print('')
print('Solid spherical harmonic synthesis at a custom grid of points')
print('===========================')

# Define some grid points as numpy arrays.  Then create a "PointGrid" class
# instance that is designed to store custom point grids.
nlat    = 5
nlon    = 10
lat     = np.linspace(np.pi / 2.0, -np.pi / 2.0, nlat)
lon     = np.linspace(0.0, 2.0 * np.pi, nlon)
r       = np.zeros((nlat,), dtype=np.float64) + shcs.r + np.arange(nlat)
grd_pnt = ph.crd.PointGrid.from_arrays(lat, lon, r)


# Do the synthesis
f = ph.shs.point(grd_pnt, shcs, nmax)


# Print the synthesized values
print(f'{f}')


del grd_pnt, f
print('===========================')
# =============================================================================


# =============================================================================
print('')
print('Solid spherical harmonic synthesis at scattered cells')
print('===========================')


# Now we define evaluation cells to synthesize area mean values.  To this end,
# we specify numpy arrays of minimum and maximum latitudes and longitudes and
# an array of spherical radii.  Using these arrays, we create a "CellSctr"
# class instance.
latmax    = np.array([0.323413435, -0.90234320952, 0.0])
latmin    = latmax - np.array([0.234323, 0.4456, np.pi / 2.0])
lonmin    = np.array([0.123456789, 4.3445234, 0.0])
lonmax    = lonmin + np.array([1.3235, 0.01, 2.0 * np.pi])
r         = np.zeros((latmin.size,), dtype=np.float64) + shcs.r + \
            np.arange(latmax.size) + 1000.0
cell_sctr = ph.crd.CellSctr.from_arrays(latmin, latmax, lonmin, lonmax, r)


# Do the synthesis
f = ph.shs.cell(cell_sctr, shcs, nmax)


# Print the synthesized values
print(f'{f}')


del cell_sctr, f
print('===========================')
# =============================================================================


# =============================================================================
print('')
print('Solid spherical harmonic synthesis at a grid of cells')
print('===========================')


# Define the grid cells and create a "CellGrid" class instance from numpy
# arrays.
nlat     = 15
nlon     = 30
latmax   = np.linspace(np.pi / 2.0, -np.pi / 2.0, nlat, False)
latmin   = np.linspace(latmax[1], -np.pi / 2.0, nlat)
lonmin   = np.linspace(0.0, 2.0 * np.pi, nlon, False)
lonmax   = np.linspace(lonmin[1], 2.0 * np.pi, nlon)
r        = np.zeros((nlat,), dtype=np.float64) + shcs.r
grd_cell = ph.crd.CellGrid.from_arrays(latmin, latmax, lonmin, lonmax, r)


# Do the synthesis
f = ph.shs.cell(grd_cell, shcs, nmax)


# Print some of the synthesized values
# -----------------------------------------------------------------------------
i = 0
j = 10
print(f'f({i}, {j}) = {f[i, j]}')


i = 4
j = 3
print(f'f({i}, {j}) = {f[i, j]}')
# -----------------------------------------------------------------------------

print('===========================')
# =============================================================================


# =============================================================================
print('')
print('Spherical harmonic analysis of block-mean values in cells')
print('===========================')


# In this example, we use the "grd_cell" instance and the signal "f" from the
# previous example


# Do the analysis using the approximate quadrature
shcs2 = ph.sha.cell(grd_cell, f, nmax, ph.sha.CELL_AQ, shcs.mu, shcs.r)


# Print some of the computed coefficients.  Note that the harmonic analysis
# with block-mean values in cells is *not* exact, hence the coefficients will
# not be equal to the input ones.
n = [2, 9, 9, 9]  # List of harmonic degrees
m = [0, 0, 4, 9]  # List of harmonic orders
cnm, snm = shcs2.get_coeffs(n, m)
for i, (n1, m1) in enumerate(zip(n, m)):
    print(f'C({n1}, {m1}) = {cnm[i]}')
    print(f'S({n1}, {m1}) = {snm[i]}')


del grd_cell, shcs2, f, shcs


print('===========================')
# =============================================================================
