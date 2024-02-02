
# Change Log
All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).
 
## [2.0] - 2024-02-01
 
Release of [`ECOMAN2.0-geodynamics`](https://github.com/mfaccenda/ECOMAN2.0-geodynamics.git)
 
### Major changes relative to ECOMAN1.0

- `D_REX_S`: 
   - postPerovskite fabric 
   - new MATLAB file for plotting fabrics (M-/J-index, fraction of anisotropy classes)

- `D_REX_M`: 
  - postPerovskite fabric 
  - MPI parallelization 
  - definition of max\_strain above which no more aggregate deformation and advection
  - elastic tensors for isotropic aggregates computed without averaging random orientatons (acs0) but with average. For isotropic rocktypes (3 for D-REX; 2,3,4,5 for SBFTEX), directly taken from databases
  
- `EXEV`: 
  - new DEMelastic: use the DEM theory to compute elastic tensor of a two-phase composite made by ellipsoidal inclusions in a background matrix
  - new MATLAB files readDEMelastic.m and readDEMviscous.m to plot elastic and viscous tensor properties in a Flinn diagram
  
- `VIZTOMO`:
  - possibility to plot fraction of anisotropy classes (triclinic, monoclinic, orthorhombic, tetragonal, hexagonal)
  
- `SKS-SPLIT`: now included in ECOMAN2.0-seismology.SKS-SPLIT
