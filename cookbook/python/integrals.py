import pyharm as ph


# Type of the first spherical harmonic function, its harmonic degree and order,
# respectively
i1 = 0
n1 = 287
m1 = 122


# Type of the second spherical harmonic function, its harmonic degree and
# order, respectively
i2 = 1
n2 = 34
m2 = 9


# Minimum and maximum co-latitudes of the integration domain
cltmin = 0.1
cltmax = 0.9


# Minimum and maximum longitudes of the intration domain
lonmin = 0.5
lonmax = 0.6


# Compute the Fourier coefficients
nmax = max(n1, n2)
pnmj = ph.leg.fourier_coeffs(nmax)


# Compute the integral of a product of two spherical harmonics
iy = ph.integ.yi1n1m1yi2n2m2(cltmin, cltmax, lonmin, lonmax,
                             i1, n1, m1, i2, n2, m2, pnmj)


# Print the value of the integral
print(f'The integral of the product of two spherical harmonics '
      f'i1 = {i1}, n1 = {n1}, m1 = {m1}, i2 = {i2}, n2 = {n2}, m2 = {m2} '
      f'is {iy}')


# Now compute the integral of a product of two Legendre functions over
# a restricted domain
ip = ph.integ.pn1m1pn2m2(cltmin, cltmax, n1, m1, n2, m2, pnmj)


print(f'The integral of the product of two Legendre functions '
      f'n1 = {n1}, m1 = {m1}, n2 = {n2}, m2 = {m2} is {ip}')
