import pyharm as ph


# Maximum harmonic degree to compute the Fourier coefficients of Legendre 
# functions
nmax = 500


# Create a "ph.leg.Pnmj" class instance using the factory method called 
# "from_zeros"
pnmj = ph.leg.Pnmj.from_zeros(nmax, ph.leg.PNMJ_ORDER_MNJ)


# Compute the Fourier coefficients
pnmj.coeffs()


# Print some Fourier coefficients
# Harmonic degree
n = 123
# Harmonic order
m = 23
# Wave-number-related variable
j = 12
k = pnmj.j2k(n, j)
c = pnmj.get_coeff(n, m, j)
print(f'Fourier coefficients for degree {n}, order {m} and wave-number {k} = '
      f'{c}')


n = 360
m = 358
j = 101
k = pnmj.j2k(j, k)
print(f'Fourier coefficients for degree {n}, order {m} and wave-number {k} = '
      f'{c}')


print('Great, all done!')
