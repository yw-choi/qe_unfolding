from numpy import *
import matplotlib.pyplot as plt
from multiprocessing import Pool
from functools import partial
from scipy import interpolate
from scipy.spatial import distance_matrix
import time
from concurrent.futures import ThreadPoolExecutor

debug = True

def main():

  print ('Calculating intensity map from grid data')
  print ('Author: ywchoi, 2020')

  print ('[Setup]')
  # crystal structure
  at_tetra = array([[8.80625,0.0,0.0],
    [0.0,8.80625,0.0],
    [0.0,0.0,12.7127]]).T # Ang
  bg_tetra = 2*pi*linalg.inv(at_tetra).T # Ang^{-1}
  SC = array([[1,1,0],[-1,1,0],[0,0,2]], dtype=float64)

  at_cubic = matmul(at_tetra, linalg.inv(SC))
  bg_cubic = 2*pi*linalg.inv(at_cubic).T # Ang^{-1}

  # control parameters
  read_intensity = False
  nprocs = 16
  qz_frac = 0.0

  # broadenings
  sigk = 0.03  # Ang^{-1}
  sigE = 0.03  # eV

  print (f' read_intensity = {read_intensity}')
  print (f' nprocs = {nprocs}')
  print (f' qz_frac = {qz_frac}')
  print (f' sigk = {sigk}')
  print (f' sigE = {sigE}')

  print (f'[Energy and q-grid]')
  # q-grid and energy grid
  dimq = 100
  nq = dimq*dimq
  # qxlim = [-0.5*bg_cubic[0,0], 1.5*bg_cubic[0,0]]
  # qylim = [-0.5*bg_cubic[1,0], 1.5*bg_cubic[1,0]]
  qxlim = [-0.5*bg_tetra[0,0], 1.5*bg_tetra[0,0]]
  qylim = [-0.5*bg_tetra[1,1], 1.5*bg_tetra[1,1]]
  qx = linspace(qxlim[0], qxlim[1], dimq)
  qy = linspace(qylim[0], qylim[1], dimq)
  qz = qz_frac*bg_tetra[2,2]

  qpoints = zeros((nq,3))
  iq = 0
  for ix in range(dimq):
    for iy in range(dimq):
      qpoints[iq,0] = qx[ix]
      qpoints[iq,1] = qx[iy]
      qpoints[iq,2] = qz
      iq += 1

  print (f' dimq = {dimq}')

  energies = [-.4]
  ne = len(energies)

  print (f' energies = {energies}')

  ef = 2.848630612
  enk = 13.605*loadtxt('./enk.dat') - ef
  nbnd, nk = shape(enk)
  wnk = loadtxt('wnk.dat')
  # wnk = ones((nbnd,nk))

  klist = loadtxt('./kmesh')
  # k-points in cartesian coordinate (Ang^{-1})
  kpoints = einsum('xi,ki->kx', bg_tetra, klist[:,:3])

  print ('[Expanding data]')
  print (f' before: nbnd, nk = {nbnd} {nk}')
  kpoints, enk, wnk = expand_data(bg_tetra, kpoints, enk, wnk)
  nbnd, nk = shape(enk)
  print (f' after: nbnd, nk = {nbnd} {nk}')

  dimk = int32(sqrt(nk))

  print ('[Grid data to intensity]')
  intensity = grid2intensity(qpoints, energies, kpoints, enk, wnk, sigk, sigE, nprocs)

  print ('[Plotting]')
  ie = 0
  energy = energies[ie]

  qx = qpoints[:,0].reshape((dimq,dimq))
  qy = qpoints[:,1].reshape((dimq,dimq))
  kx = kpoints[:,0].reshape((dimk,dimk))
  ky = kpoints[:,1].reshape((dimk,dimk))

  z = intensity[:,ie].reshape((dimq,dimq))

  plt.pcolormesh(qx, qy, z, cmap='afmhot')

  for ibnd in range(nbnd):
    if enk[ibnd,:].max()<energy:
      continue
    if enk[ibnd,:].min()>energy:
      break
    ek = enk[ibnd,:].reshape((dimk, dimk))
    plt.contour(kx, ky, ek, levels=[energy], colors='white', linestyles='solid')

  plt.xlim(qxlim[0], qxlim[1])
  plt.ylim(qylim[0], qylim[1])
  plt.savefig('intensity.png')
  plt.show()

  return

def grid2intensity(qpoints, energies, kpoints, enk, wnk, sigk, sigE, nprocs):

  if debug: start = time.time()

  nq = len(qpoints)
  ne = len(energies)
  nbnd, nk = enk.shape

  intensity = zeros((nq, ne))

  if debug: start1 = time.time()
  dk = distance_matrix(qpoints, kpoints)
  if debug: end1 = time.time(); print (f' dk = {end1-start1} s')
  #
  if debug: start1 = time.time()
  I_qk = 1/(sigk*sqrt(2*pi))*exp(-0.5*(dk*dk)/(sigk*sigk))
  if debug: end1 = time.time(); print (f' I_qk = {end1-start1} s')
  #
  if debug: start1 = time.time()
  de = einsum('e,nk->enk', energies, ones((nbnd,nk))) \
      - einsum('e,nk->enk', ones(ne), enk)
  if debug: end1 = time.time(); print (f' de = {end1-start1} s')
  #
  if debug: start1 = time.time()
  I_enk = 1/(sigE*sqrt(2*pi))*exp(-0.5*(de*de)/(sigE*sigE))
  I_enk = einsum('nk,enk->enk', wnk, I_enk)
  if debug: end1 = time.time(); print (f' I_enk = {end1-start1} s')
  #
  if debug: start1 = time.time()
  with ThreadPoolExecutor(max_workers=nprocs) as executor:
    jobs = executor.map(lambda I_k: einsum('k,enk->e', I_k, I_enk), I_qk)
    for iq, intensity_q in enumerate(jobs):
      intensity[iq,:] = intensity_q
  if debug: end1 = time.time(); print (f' intensity = {end1-start1} s')
  #
  if debug: end = time.time(); print (f' grid2intensity: {end-start} s')

  return intensity

def expand_data(bg, kpoints, enk, wnk):

  nbnd, nk = shape(enk)

  nmin, nmax = -1, 1

  lsc = nmax-nmin+1
  nsc = lsc*lsc
  nktot = nk*nsc

  kpoints_all = zeros((nktot,3))
  enk_all = zeros((nbnd,nktot))
  wnk_all = zeros((nbnd,nktot))

  dimk = int(sqrt(nk))
  ik = 0
  for i1 in arange(nmin,nmax+1):
    for ik1 in range(dimk):
      for i2 in arange(nmin,nmax+1):
        for ik2 in range(dimk):

          ik0 = ik1*dimk+ik2

          k = i1*bg[:,0]+i2*bg[:,1]+kpoints[ik0,:]
          kpoints_all[ik,:] = k
          enk_all[:,ik] = enk[:,ik0]
          wnk_all[:,ik] = wnk[:,ik0]

          ik += 1
  return kpoints_all, enk_all, wnk_all

if __name__ == '__main__':
  main()
