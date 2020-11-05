from numpy import *

kz = 0.0

nk = 30

v1 = array([1,-1])
v2 = array([1, 1])

for n1 in arange(nk):
  for n2 in arange(nk):
    k = n1*v1/nk + n2*v2/nk

    print (f'{k[0]:10.6f} {k[1]:10.6f} {kz:10.6f}')
