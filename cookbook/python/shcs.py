import pyharm as ph


# INPUTS
# =============================================================================
# Define the path to an input "gfc" text file with spherical harmonic
# coefficients.  For details on the structure of the "gfc" file, see the
# documentation.
shcs_in_file = '../../data/input/EGM96-degree10.gfc'


# Define the path to an output binary file with spherical harmonic
# coefficients.  For details on the structure of the output binary file, see
# the documentation of the "shc" module of "CHarm".
shcs_out_file = '../../data/output/EGM96-degree10.shcs'
# =============================================================================


# =============================================================================
# Read the spherical harmonic coefficients using the factory method "from_file"
# of the "Shc" class from the "ph.shc" module.
# -----------------------------------------------------------------------------
shcs = ph.shc.Shc.from_file('gfc', shcs_in_file)
# -----------------------------------------------------------------------------


# Now print some more or less randomly chosen spherical harmonic coefficients
# -----------------------------------------------------------------------------
# Let's start with zonal coefficients
n = 9
m = 0
cnm, snm = shcs.get_coeffs(n, m)
print(f'C({n}, {m}) = {cnm}')
print(f'S({n}, {m}) = {snm}')


# Now some tesseral coefficients
n = 9
m = 4
cnm, snm = shcs.get_coeffs(n, m)
print(f'C({n}, {m}) = {cnm}')
print(f'S({n}, {m}) = {snm}')


# And finally some sectorial coefficients
n = 9
m = 9
cnm, snm = shcs.get_coeffs(n, m)
print(f'C({n}, {m}) = {cnm}')
print(f'S({n}, {m}) = {snm}')
# -----------------------------------------------------------------------------


# Now let's save the coefficients to a binary file and then read them back to
# a new "Shc" class instance.  Just for the fun...
# -----------------------------------------------------------------------------
shcs.to_file('bin', shcs_out_file)


# Read back the coefficients to a new class instance called "shcs2"
shcs2 = ph.shc.Shc.from_file('bin', shcs_out_file)
# -----------------------------------------------------------------------------


# Compute degree variances from "shcs"
# -----------------------------------------------------------------------------
dv = shcs.dv()


n = 0
print('\n')
print(f'Degree variance for harmonic degree {n} = {dv[n]}')
n = 4
print(f'Degree variance for harmonic degree {n} = {dv[n]}')
n = 10
print(f'Degree variance for harmonic degree {n} = {dv[n]}')
# -----------------------------------------------------------------------------


# Now check whether "shcs" and "shcs2" contain the same coefficients by
# computing difference degree variances
# -----------------------------------------------------------------------------
ddv = shcs.ddv(shcs2)


n = 0
print('\n')
print(f'Difference degree variance for harmonic degree {n} = {ddv[n]}')
n = 4
print(f'Difference degree variance for harmonic degree {n} = {ddv[n]}')
n = 10
print(f'Difference degree variance for harmonic degree {n} = {ddv[n]}')
# -----------------------------------------------------------------------------


print('Great, all done!')
# =============================================================================
