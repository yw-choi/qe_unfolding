
## unfold.x

Unfolding supercell band structure obtained from Quantum ESPRESSO.

Young Woo Choi (ywchoi@yonsei.ac.kr)

## Compilation

Modify make.inc to fit your system environment. 
You may use the same make.inc for compiling Quantum ESPRESSO distribution.
But, you must specify TOPDIR, which is the top directory of your QE(>6.4).

## Usage

1. Do supercell scf calculation to converge the charge density.
2. Do nscf or band calculations for k-points along the primitive Brillouin zone (See K_POINTS section in test/pw.bands.in)
3. Fill up the input file for unfold.x: Here the SC variable refers to the supercell matrix such that [t1 t2 t3]^T = SC [a1 a2 a3]^T.

```
&inputpp
  outdir = './'
  prefix = 'graphene'
  first_k = 0 
  last_k  = 0
  first_band = 1
  last_band = 18
  SC(1,:) = 2, 0, 0
  SC(2,:) = 0, 2, 0
  SC(3,:) = 0, 0, 1
/
```

4. mpirun -np X unfold.x -i unfold.in 
5. You get enk.dat (KS eigenvalues from supercell calculation) and wnk.dat (Unfolding weights)
6. Plot unfolded band structure (refer test/plot_unfold.py)
