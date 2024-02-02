# ECOMAN - Exploring the COnsequences of Mechanical ANisotropy

ECOMAN is a software package that (i) simulates the development of strain-induced mantle fabrics (LPO+SPO) and (ii) tests the effects of the mechanical (elastic and viscous) anisotropy associated with these fabrics on seismic imaging and mantle convection. It includes programs that:

1) estimate strain/stress-induced rock fabrics (LPO and SPO) and their elastic and viscous anisotropic mechanical properties (ECOMAN-geodynamics: D-REX_S, D-REX_M, EXEV),

2) post-process the simulated rock fabrics for visualisation of their isotropic/anisotropic mechanical properties and deformational history (ECOMAN-geodynamics: VIZTOMO, VIZVISC), and format the elastic tensors generating input files for seismological synthetics (ECOMAN-geodynamics: VIZTOMO), and

3) test the elastic response of anisotropic media by performing seismological forward/inverse modelling and, in particular, isotropic and anisotropic seismic tomographies on synthetic and real seismic datasets (ECOMAN-seismology: SKS-SPLIT, PSI).

ECOMAN is supported by the ERC StG 758199 NEWTON

# General information and Installation

ECOMAN-geodynamics includes:

1) D-REX_S: a program designed for modelling the evolution of strain-induced LPO fabrics and related elastic properties of a single, two-mineral phases mantle aggregate, as a function of the imposed flow field, amount of strain, crystal plasticity, P-T conditions and additional effects related to SPO fabrics.

2) D-REX_M: a program that computes the evolution of the LPO and related elastic properties of multiple, two-mineral phases mantle aggregates, as a function of the single crystal plastic and elastic properties, and of the flow field, deformation mechanisms and P-T conditions resulting from 2D-3D geodynamic simulations.

3) EXEV: includes routines to compute SPO fabrics using Effective Medium Theory modelling for a multi-component layered system (Smoothed Transversely Isotropic Long-Wavelength Equivalent, STILWE; Backus, 1962) or a two-component system with similar ellipsoidal inclusions in a uniform background matrix (Differential Effective Medium, DEM; e.g., Mainprice, 1997).

4) VIZTOMO: processes the D-REX_M output for the visualisation of the aggregates' elastic and deformational history properties.

5) VIZVISC: processes the D-REX_M output for the visualisation of the aggregates' viscous and deformational history properties.


ECOMAN-geodynamics programs are mostly written in Fortran, and where most of the routines are parallelized with shared memory architecture (OpenMP), providing good scalability with increasing number of cores. In addition, D-REX_M is parallelized with a hybrid MPI-OpenMP architecture. 
As a result, ECOMAN-geodynamics requires installation of the Intel Fortran compilers and HDF5 libraries.

Software compilation: from the directory of each software, execute ./bash_compile

More detailed information and instructions are provided in the user [manual](https://newtonproject.geoscienze.unipd.it/wp-content/uploads/2024/02/ECOMAN2.0_manual.pdf). 

# Software website

https://newtonproject.geoscienze.unipd.it/ecoman/


# Developers

[Manuele Faccenda](mailto:manuele.faccenda@unipd.it)

[Brandon VanderBeek](mailto:brandon.p.vanderbeek@gmail.com)

Albert de Montserrat

Jianfeng Yang
