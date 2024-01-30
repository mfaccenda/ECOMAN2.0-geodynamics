 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !!
 !!    Copyright (c) 2018-2023, Universita' di Padova, Manuele Faccenda
 !!    All rights reserved.
 !!
 !!    This software package was developed at:
 !!
 !!         Dipartimento di Geoscienze
 !!         Universita' di Padova, Padova         
 !!         via Gradenigo 6,            
 !!         35131 Padova, Italy 
 !!
 !!    project:    ECOMAN
 !!    funded by:  ERC StG 758199 - NEWTON
 !!
 !!    ECOMAN is free software package: you can redistribute it and/or modify
 !!    it under the terms of the GNU General Public License as published
 !!    by the Free Software Foundation, version 3 of the License.
 !!
 !!    ECOMAN is distributed WITHOUT ANY WARRANTY; without even the implied
 !!    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 !!    See the GNU General Public License for more details.
 !!
 !!    You should have received a copy of the GNU General Public License
 !!    along with ECOMAN. If not, see <http://www.gnu.org/licenses/>.
 !!
 !!
 !!    Contact:
 !!        Manuele Faccenda    [manuele.faccenda@unipd.it]
 !!        Brandon VanderBeek  [brandon.vanderbeek@unipd.it]
 !!
 !!
 !!    Main development team:
 !!        Manuele Faccenda    [manuele.faccenda@unipd.it]
 !!        Brandon VanderBeek  [brandon.vanderbeek@unipd.it]
 !!        Albert de Montserrat Navarro
 !!        Jianfeng Yang   
 !!
 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   Load elastic tensors for upper mantle and transition zone pahses !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE elastic_database(S0,dS0dp,dS0dp2,dS0dt,dS0dpt,dS0dp5,dS0dt5)

   IMPLICIT NONE
 
   INTEGER :: t
   DOUBLE PRECISION, DIMENSION(20,6,6) :: S0,dS0dp,dS0dp2,dS0dt,dS0dpt
   DOUBLE PRECISION, DIMENSION(6,6) :: dS0dp5,dS0dt5

!!! Stiffness matrices
!!! S0(1,:,:) : olivine, existing in D-Rex
!!! S0(2,:,:) : enstatite, existing in D-Rex
!!! S0(3,:,:) : wadsleyite, from Zha et al., 1997, Mg2SiO4 at ambient P-T
!!! S0(4,:,:) : wadsleyite, from Sinogeikin et al., 1998, (Mg0.92Fe0.08)2SiO4 at ambient P-T
!!! S0(5,:,:) : ringwoodite, from Sinogeikin et al., 1998, at ambient P-T
!!! S0(6,:,:) : pyrope garnet (Chai et al., 1998, at ambient P-T)
!!! S0(7,:,:) : Mg perovskite, Wentzcovitch et al., 2004, at ambient P-T (P shifted of 5GPa)
!!! S0(8,:,:) : Mg perovskite, Zhang et al, 2013, EPSL, LDA ab-initio, at 34 GPa and extraplated at 2000 K
!!! S0(9,:,:) : Mg perovskite, Zhang et al, 2013, EPSL, GGA ab-initio, at 48 GPa and extraplated at 2000 K
!!! S0(10,:,:) : MgO periclase, Karki et al., 2000, at ambient P-T
!!! S0(11,:,:) : post-Mg perovskite, Zhang et al., 2016, EPSL, at 99.95 GPa and 1500 K

   S0 = 0d0 ; dS0dp = 0d0 ; dS0dp2 = 0d0 ; dS0dt = 0d0 ; dS0dpt = 0d0

!!! Olivine ( Isaak et al., 1992, RAM sample, J. Geophys. Res.)
   
   t = 1

   !(GPa)
   S0(t,1,1) = 320.71d0 ; S0(t,1,2) = 69.84d0 ; S0(t,1,3) = 71.22d0
   S0(t,2,1) = S0(t,1,2)  ; S0(t,2,2) = 197.25d0 ; S0(t,2,3) = 74.8d0
   S0(t,3,1) = S0(t,1,3)  ; S0(t,3,2) = S0(t,2,3) ; S0(t,3,3) = 234.32d0
   S0(t,4,4) = 63.77d0  ; S0(t,5,5) = 77.67d0 ; S0(t,6,6) = 78.36d0

   ! Temperature derivative (-GPa/K) (Isaak et al., 1992, RAM sample, J. Geophys. Res.)
   dS0dt(t,1,1) = 4.02d-2 ; dS0dt(t,1,2) = 1.14d-2 ; dS0dt(t,1,3) = 0.96d-2 
   dS0dt(t,2,1) = dS0dt(t,1,2) ; dS0dt(t,2,2) = 3.10d-2 ; dS0dt(t,2,3) = 0.72d-2 
   dS0dt(t,3,1) = dS0dt(t,1,3) ; dS0dt(t,3,2) = dS0dt(t,2,3) ; dS0dt(t,3,3) = 3.53d-2 
   dS0dt(t,4,4) = 1.26d-2 ; dS0dt(t,5,5) = 1.30d-2 ; dS0dt(t,6,6) = 1.56d-2 

   ! Pressure derivative (Calculated for P <= 14 GPa from Zha et al.,1998, EPSL))
   dS0dp(t,1,1) = 5.41d0 ; dS0dp(t,1,2) = 1.88d0 ; dS0dp(t,1,3) = 1.53d0 
   dS0dp(t,2,1) = dS0dp(t,1,2) ; dS0dp(t,2,2) = 5.26d0 ; dS0dp(t,2,3) = 1.60d0 
   dS0dp(t,3,1) = dS0dp(t,1,3) ; dS0dp(t,3,2) = dS0dp(t,2,3) ; dS0dp(t,3,3) = 3.78d0 
   dS0dp(t,4,4) = 1.08d0 ; dS0dp(t,5,5) = 1.42 ; dS0dp(t,6,6) = 1.80d0 

!!! Enstatite (existing in D-Rex)
   
   t = 2
   
   !(GPa)
   S0(t,1,1) = 236.9d0 ; S0(t,1,2) = 79.6d0  ; S0(t,1,3) = 63.2d0
   S0(t,2,1) = S0(t,1,2)  ; S0(t,2,2) = 180.5d0 ; S0(t,2,3) = 56.8d0
   S0(t,3,1) = S0(t,1,3)  ; S0(t,3,2) = S0(t,2,3) ; S0(t,3,3) = 230.4d0
   S0(t,4,4) = 84.3d0  ; S0(t,5,5) = 79.4d0  ; S0(t,6,6) = 80.1d0

   ! Temperature derivative (-GPa/K) (Jackson et al., 2007, PEPI)
   dS0dt(t,1,1) = 3.64d-2 ; dS0dt(t,1,2) = 1.52d-2 ; dS0dt(t,1,3) = 2.29d-2 
   dS0dt(t,2,1) = dS0dt(t,1,2) ; dS0dt(t,2,2) = 3.43d-2 ; dS0dt(t,2,3) = 1.63d-2 
   dS0dt(t,3,1) = dS0dt(t,1,3) ; dS0dt(t,3,2) = dS0dt(t,2,3) ; dS0dt(t,3,3) = 5.7d-2 
   dS0dt(t,4,4) = 1.23d-2 ; dS0dt(t,5,5) = 1.58d-2 ; dS0dt(t,6,6) = 1.51d-2 

   ! Pressure derivative (Chai et al.,1997, JGR)
   dS0dp(t,1,1) = 10.27d0 ; dS0dp(t,1,2) = 6.22d0 ; dS0dp(t,1,3) = 6.63d0 
   dS0dp(t,2,1) = dS0dp(t,1,2) ; dS0dp(t,2,2) = 8.87d0 ; dS0dp(t,2,3) = 7.26d0 
   dS0dp(t,3,1) = dS0dp(t,1,3) ; dS0dp(t,3,2) = dS0dp(t,2,3) ; dS0dp(t,3,3) = 11.07d0 
   dS0dp(t,4,4) = 1.23d0 ; dS0dp(t,5,5) = 0.75d0 ; dS0dp(t,6,6) = 2.78d0 

   ! Second pressure derivatives (GPa^-1, Chai et al.,1997, JGR)
   dS0dp2(t,1,1) = -0.47d0 ; dS0dp2(t,1,2) = -0.33d0 ; dS0dp2(t,1,3) = -0.26d0 
   dS0dp2(t,2,1) = dS0dp2(t,1,2) ; dS0dp2(t,2,2) = -0.38d0 ; dS0dp2(t,2,3) = -0.31d0 
   dS0dp2(t,3,1) = dS0dp2(t,1,3) ; dS0dp2(t,3,2) = dS0dp2(t,2,3) ; dS0dp2(t,3,3) = -0.53d0 

!!! Wadsleyite (Zha et al., 1997) Mg2SiO4 at ambient P-T
   
   t = 3
   
   !(GPa)
   S0(t,1,1) = 370.5d0 ; S0(t,1,2) = 65.6d0  ; S0(t,1,3) = 95.2d0
   S0(t,2,1) = S0(t,1,2)  ; S0(t,2,2) = 367.7d0 ; S0(t,2,3) = 105.1d0
   S0(t,3,1) = S0(t,1,3)  ; S0(t,3,2) = S0(t,2,3) ; S0(t,3,3) = 272.4d0
   S0(t,4,4) = 111.2d0  ; S0(t,5,5) = 122.5d0  ; S0(t,6,6) = 103.1d0

   ! Temperature derivatives (for the moment, same as olivine)
   dS0dt(t,:,:) = dS0dt(1,:,:)
   
   ! Pressure derivative (Zha et al., 1997) Mg2SiO4
   dS0dp(t,1,1) = 5.21d0 ; dS0dp(t,1,2) = 4.06d0 ; dS0dp(t,1,3) = 3.3d0 
   dS0dp(t,2,1) = dS0dp(t,1,2) ; dS0dp(t,2,2) = 6.83d0 ; dS0dp(t,2,3) = 3.3d0 
   dS0dp(t,3,1) = dS0dp(t,1,3) ; dS0dp(t,3,2) = dS0dp(t,2,3) ; dS0dp(t,3,3) = 8.06d0 
   dS0dp(t,4,4) = 1.39d0 ; dS0dp(t,5,5) = 0d0 ; dS0dp(t,6,6) = 1.88d0 

!!! Wadsleyite (Zhou et al., 2022, Nat. Comm) Mg(0.906Fe0.094)2SiO4 H2O (0.15 wt%) at ambient P-T
   
   t = 4
   
   !(GPa)
   S0(t,1,1) = 341d0 ; S0(t,1,2) = 75d0  ; S0(t,1,3) = 99d0
   S0(t,2,1) = S0(t,1,2)  ; S0(t,2,2) = 358d0 ; S0(t,2,3) = 102d0
   S0(t,3,1) = S0(t,1,3)  ; S0(t,3,2) = S0(t,2,3) ; S0(t,3,3) = 244d0
   S0(t,4,4) = 106d0  ; S0(t,5,5) = 109d0  ; S0(t,6,6) = 90.8d0

   ! Temperature derivatives (-GPa/K)
   dS0dt(t,1,1) = 2.90d-2 ; dS0dt(t,1,2) = 0.60d-2 ; dS0dt(t,1,3) = 1.90d-2 
   dS0dt(t,2,1) = dS0dt(t,1,2) ; dS0dt(t,2,2) = 3.60d-2 ; dS0dt(t,2,3) = 1.20d-2 
   dS0dt(t,3,1) = dS0dt(t,1,3) ; dS0dt(t,3,2) = dS0dt(t,2,3) ; dS0dt(t,3,3) = 4.00d-2 
   dS0dt(t,4,4) = 1.80d-2 ; dS0dt(t,5,5) = 1.10d-2 ; dS0dt(t,6,6) = 1.30d-2 
   
   ! Pressure derivative
   dS0dp(t,1,1) = 7.80d0 ; dS0dp(t,1,2) = 4.60d0 ; dS0dp(t,1,3) = 2.33d0 
   dS0dp(t,2,1) = dS0dp(t,1,2) ; dS0dp(t,2,2) = 5.37d0 ; dS0dp(t,2,3) = 4.2d0 
   dS0dp(t,3,1) = dS0dp(t,1,3) ; dS0dp(t,3,2) = dS0dp(t,2,3) ; dS0dp(t,3,3) =10.5d0 
   dS0dp(t,4,4) = 0.87d0 ; dS0dp(t,5,5) = 0.84d0 ; dS0dp(t,6,6) = 2.43d0 

   ! Second pressure derivatives (GPa^-1)
   dS0dp2(t,1,1) = -0.26d0 ; dS0dp2(t,1,2) = -0.11d0 ; dS0dp2(t,1,3) = 0.00d0 
   dS0dp2(t,2,1) = dS0dp2(t,1,2) ; dS0dp2(t,2,2) = 0.00d0 ; dS0dp2(t,2,3) = -0.14d0 
   dS0dp2(t,3,1) = dS0dp2(t,1,3) ; dS0dp2(t,3,2) = dS0dp2(t,2,3) ; dS0dp2(t,3,3) = -0.46d0 
   dS0dp2(t,4,4) = 0.00d0 ; dS0dp2(t,5,5) = 0.00d0 ; dS0dp2(t,6,6) = -0.053d0 

!!! Ringwoodite (Sinogeikin et al., 2003) (Mg0.91Fe0.09)2SiO4 at ambient P-T
   
   t = 5
   
   !(GPa)
   !S0(t,1,1) = 332d0 ; S0(t,1,2) =117.6d0  ; S0(t,1,3) = S0(t,1,2)
   S0(t,1,1) = 329d0 ; S0(t,1,2) =118d0  ; S0(t,1,3) = S0(t,1,2)
   S0(t,2,1) = S0(t,1,2)  ; S0(t,2,2) = S0(t,1,1) ; S0(t,2,3) = S0(t,1,2)
   S0(t,3,1) = S0(t,1,2)  ; S0(t,3,2) = S0(t,1,2) ; S0(t,3,3) = S0(t,1,1)
   !S0(t,4,4) = 129.59d0 ; S0(t,5,5) = S0(t,4,4)  ; S0(t,6,6) = S0(t,4,4)
   S0(t,4,4) = 130d0 ; S0(t,5,5) = S0(t,4,4)  ; S0(t,6,6) = S0(t,4,4)

   ! Temperature derivative (-GPa/K) (Sinogeikin et al., 2003)          
   dS0dt(t,1,1) = 4.9d-2 ; dS0dt(t,1,2) = 0.7d-2 ; dS0dt(t,1,3) = dS0dt(t,1,2) 
   dS0dt(t,2,1) = dS0dt(t,1,2) ; dS0dt(t,2,2) = dS0dt(t,1,1) ; dS0dt(t,2,3) = dS0dt(t,1,2)  
   dS0dt(t,3,1) = dS0dt(t,1,2) ; dS0dt(t,3,2) = dS0dt(t,1,2) ; dS0dt(t,3,3) = dS0dt(t,1,1) 
   dS0dt(t,4,4) = 1.30d-2 ; dS0dt(t,5,5) = dS0dt(t,4,4) ; dS0dt(t,6,6) = dS0dt(t,4,4) 

   ! Pressure derivative (Sinogeikin et al., 2003,  ambient P-T)
  !dS0dp(t,1,1) = 6.06d0 ; dS0dp(t,1,2) = 0.82d0 ; dS0dp(t,1,3) = dS0dp(t,1,2)
   dS0dp(t,1,1) = 6.20d0 ; dS0dp(t,1,2) = 0.80d0 ; dS0dp(t,1,3) = dS0dp(t,1,2)
   dS0dp(t,2,1) = dS0dp(t,1,2) ; dS0dp(t,2,2) = dS0dp(t,1,1) ; dS0dp(t,2,3) = dS0dp(t,1,2) 
   dS0dp(t,3,1) = dS0dp(t,1,2) ; dS0dp(t,3,2) = dS0dp(t,1,2) ; dS0dp(t,3,3) = dS0dp(t,1,1) 
  !dS0dp(t,4,4) = 2.77d0 ; dS0dp(t,5,5) = dS0dp(t,4,4) ; dS0dp(t,6,6) = dS0dp(t,4,4) 
   dS0dp(t,4,4) = 2.8d0 ; dS0dp(t,5,5) = dS0dp(t,4,4) ; dS0dp(t,6,6) = dS0dp(t,4,4)

!!! Garnet (Chai et al., 1997, GRL at ambient P-T, Py-Gr-Alm)
   
   t = 6
   
   !(GPa)
   S0(t,1,1) = 299.1d0 ; S0(t,1,2) =106.7d0  ; S0(t,1,3) = S0(t,1,2)
   S0(t,2,1) = S0(t,1,2)  ; S0(t,2,2) = S0(t,1,1) ; S0(t,2,3) = S0(t,1,2)
   S0(t,3,1) = S0(t,1,2)  ; S0(t,3,2) = S0(t,1,2) ; S0(t,3,3) = S0(t,1,1)
   S0(t,4,4) =  93.7d0 ; S0(t,5,5) = S0(t,4,4)  ; S0(t,6,6) = S0(t,4,4)

   ! Temperature derivative (-GPa/K) (Sinogeikin & Bass, 2002, Pyrope 100 %)
   dS0dt(t,1,1) = 3.05d-2 ; dS0dt(t,1,2) = 0.58d-2 ; dS0dt(t,1,3) = dS0dt(t,1,2) 
   dS0dt(t,2,1) = dS0dt(t,1,2) ; dS0dt(t,2,2) = dS0dt(t,1,1) ; dS0dt(t,2,3) = dS0dt(t,1,2)  
   dS0dt(t,3,1) = dS0dt(t,1,2) ; dS0dt(t,3,2) = dS0dt(t,1,2) ; dS0dt(t,3,3) = dS0dt(t,1,1) 
   dS0dt(t,4,4) = 0.71d-2 ; dS0dt(t,5,5) = dS0dt(t,4,4) ; dS0dt(t,6,6) = dS0dt(t,4,4) 

   ! Pressure derivative (from Chai et al., 1997, GRL at ambient P-T)
   dS0dp(t,1,1) = 6.54d0 ; dS0dp(t,1,2) = 2.87d0 ; dS0dp(t,1,3) = dS0dp(t,1,2)
   dS0dp(t,2,1) = dS0dp(t,1,2) ; dS0dp(t,2,2) = dS0dp(t,1,1) ; dS0dp(t,2,3) = dS0dp(t,1,2) 
   dS0dp(t,3,1) = dS0dp(t,1,2) ; dS0dp(t,3,2) = dS0dp(t,1,2) ; dS0dp(t,3,3) = dS0dp(t,1,1) 
   dS0dp(t,4,4) = 1.72d0 ; dS0dp(t,5,5) = dS0dp(t,4,4) ; dS0dp(t,6,6) = dS0dp(t,4,4) 

!!! MgSiO3 perovskite (Wentzcovitch et al., 2004, Phys. Rev. Lett, 2004 at ambient T, P = 5 GPa)

   t = 7

   !(GPa)
   S0(t,1,1) = 484d0 ; S0(t,1,2) =146d0  ; S0(t,1,3) = S0(t,1,2)
   S0(t,2,1) = S0(t,1,2)  ; S0(t,2,2) = 542d0 ; S0(t,2,3) = 162d0
   S0(t,3,1) = S0(t,1,3)  ; S0(t,3,2) = S0(t,2,3) ; S0(t,3,3) = 477d0
   S0(t,4,4) = 195d0   ; S0(t,5,5) = 172d0 ; S0(t,6,6) = 151d0

   ! Temperature derivative (-GPa/K) at 100 GPa, 300 °K
   dS0dt(t,1,1) = 1.8d-2 ; dS0dt(t,1,2) = 0.4d-2 ; dS0dt(t,1,3) = 0.2d-2
   dS0dt(t,2,1) = dS0dt(t,1,2) ; dS0dt(t,2,2) = 3.3d-2 ; dS0dt(t,2,3) = 0.5d-2
   dS0dt(t,3,1) = dS0dt(t,1,3) ; dS0dt(t,3,2) = dS0dt(t,2,3) ; dS0dt(t,3,3) = 2.2d-2
   dS0dt(t,4,4) = 2.2d-2 ; dS0dt(t,5,5) = 0.6d-2 ; dS0dt(t,6,6) = 1.6d-2

   ! Temperature derivative at ambient conditions shifted by 5GPa
   dS0dt5 = 0d0
   dS0dt5(1,1) = 4.6d-2 ; dS0dt5(1,2) = 1.3d-2 ; dS0dt5(1,3) = 0.1d-2
   dS0dt5(2,1) = dS0dt5(1,2) ; dS0dt5(2,2) = 6.4d-2 ; dS0dt5(2,3) = 0.0d-2
   dS0dt5(3,1) = dS0dt5(1,3) ; dS0dt5(3,2) = dS0dt5(2,3) ; dS0dt5(3,3) = 5.4d-2
   dS0dt5(4,4) = 3.2d-2 ; dS0dt5(5,5) = 1.6d-2 ; dS0dt5(6,6) = 2.6d-2

   ! Pressure derivative at 100 GPa, 300 °K, shifted by 5GPa
   dS0dp(t,1,1) = 3.4d0 ; dS0dp(t,1,2) = 3.3d0 ; dS0dp(t,1,3) = 2.4d0
   dS0dp(t,2,1) = dS0dp(t,1,2) ; dS0dp(t,2,2) = 5.3d0 ; dS0dp(t,2,3) = 2.5d0
   dS0dp(t,3,1) = dS0dp(t,1,3) ; dS0dp(t,3,2) = dS0dp(t,2,3) ; dS0dp(t,3,3) = 5.2d0
   dS0dp(t,4,4) = 1.4d0 ; dS0dp(t,5,5) = 0.8d0 ; dS0dp(t,6,6) = 1.4

   ! Pressure derivative at ambient conditions shifted by 5GPa
   dS0dp5 = 0d0
   dS0dp5(1,1) = 6.3d0 ; dS0dp5(1,2) = 4.2d0 ; dS0dp5(1,3) = 3.0d0
   dS0dp5(2,1) = dS0dp5(1,2) ; dS0dp5(2,2) = 7.5d0 ; dS0dp5(2,3) = 3.3d0
   dS0dp5(3,1) = dS0dp5(1,3) ; dS0dp5(3,2) = dS0dp5(2,3) ; dS0dp5(3,3) = 8.2d0
   dS0dp5(4,4) = 2.2d0 ; dS0dp5(5,5) = 1.7d0 ; dS0dp5(6,6) = 2.4

   ! Pressure/Temperature second derivative (-GPa^2/GPa/K)
   dS0dpt(t,1,1) = -1.7d-5 ; dS0dpt(t,1,2) = -0.6d-5 ; dS0dpt(t,1,3) = -0.6d-5
   dS0dpt(t,2,1) = dS0dpt(t,1,2) ; dS0dpt(t,2,2) = -1.6d-5 ; dS0dpt(t,2,3) = -0.7d-5
   dS0dpt(t,3,1) = dS0dpt(t,1,3) ; dS0dpt(t,3,2) = dS0dpt(t,2,3) ; dS0dpt(t,3,3) = -1.6d-5
   dS0dpt(t,4,4) = -5.0d-7 ; dS0dpt(t,5,5) = -5.0d-7 ; dS0dpt(t,6,6) = -0.4d-5

!!! MgSiO3 perovskite (Zhang et al., 2013, EPSL, LDA ab-initio, interplation of data)
!!! T0 = 1500 K, P0 = 34 GPa

   t = 8 

   !(GPa)
   S0(t,1,1) = 581.42d0 ; S0(t,1,2) =249.02d0  ; S0(t,1,3) = 212.93d0
   S0(t,2,1) = S0(t,1,2)  ; S0(t,2,2) = 669.06d0 ; S0(t,2,3) = 239.01d0
   S0(t,3,1) = S0(t,1,3)  ; S0(t,3,2) = S0(t,2,3) ; S0(t,3,3) = 626.31d0
   S0(t,4,4) = 228.63d0   ; S0(t,5,5) = 200.26d0 ; S0(t,6,6) = 183.89d0

   ! Temperature derivative (-GPa/K) at 100 GPa, 300 °K
   dS0dt(t,1,1) = 3.99d-2 ; dS0dt(t,1,2) = 1.76d-2 ; dS0dt(t,1,3) = -0.64d-2
   dS0dt(t,2,1) = dS0dt(t,1,2) ; dS0dt(t,2,2) = 5.53d-2 ; dS0dt(t,2,3) = 0.26d-2
   dS0dt(t,3,1) = dS0dt(t,1,3) ; dS0dt(t,3,2) = dS0dt(t,2,3) ; dS0dt(t,3,3) = 4.32d-2
   dS0dt(t,4,4) = 2.12d-2 ; dS0dt(t,5,5) = 1.25d-2 ; dS0dt(t,6,6) = 1.96d-2

   ! Pressure derivative at 100 GPa, 300 °K, shifted by 5GPa
   dS0dp(t,1,1) = 3.73d0 ; dS0dp(t,1,2) = 3.3d0 ; dS0dp(t,1,3) = 2.48d0
   dS0dp(t,2,1) = dS0dp(t,1,2) ; dS0dp(t,2,2) = 5.54d0 ; dS0dp(t,2,3) = 2.53d0
   dS0dp(t,3,1) = dS0dp(t,1,3) ; dS0dp(t,3,2) = dS0dp(t,2,3) ; dS0dp(t,3,3) = 5.47d0
   dS0dp(t,4,4) = 1.38d0 ; dS0dp(t,5,5) = 0.85d0 ; dS0dp(t,6,6) = 1.56

   ! Pressure/Temperature second derivative (-GPa^2/GPa/K)
   dS0dpt(t,1,1) = 1.6d-4 ; dS0dpt(t,1,2) = 1.74d-4 ; dS0dpt(t,1,3) = -6.44d-5
   dS0dpt(t,2,1) = dS0dpt(t,1,2) ; dS0dpt(t,2,2) = 4.36d-5 ; dS0dpt(t,2,3) = 3.53d-5
   dS0dpt(t,3,1) = dS0dpt(t,1,3) ; dS0dpt(t,3,2) = dS0dpt(t,2,3) ; dS0dpt(t,3,3) = 5.11d-5
   dS0dpt(t,4,4) = -1.73d-5 ; dS0dpt(t,5,5) = 4.37d-5 ; dS0dpt(t,6,6) = 4.16d-5

!!! MgSiO3 perovskite (Zhang et al., 2013, EPSL, GGA ab-initio, at 48 GPa and extraplated at 2000 K)
!!! T0 = 1500 K, P0 = 48 GPa

   t = 9 
   
   !(GPa)
   S0(t,1,1) = 601.83d0 ; S0(t,1,2) =271.81d0  ; S0(t,1,3) = 236.00d0
   S0(t,2,1) = S0(t,1,2)  ; S0(t,2,2) = 698.17d0 ; S0(t,2,3) = 265.58d0
   S0(t,3,1) = S0(t,1,3)  ; S0(t,3,2) = S0(t,2,3) ; S0(t,3,3) = 652.27d0
   S0(t,4,4) = 233.53d0   ; S0(t,5,5) = 205.41d0 ; S0(t,6,6) = 190.43d0

   ! Temperature derivative (-GPa/K) at 100 GPa, 300 °K
   dS0dt(t,1,1) = 3.14d-2 ; dS0dt(t,1,2) = 1.94d-2 ; dS0dt(t,1,3) = -0.55d-2
   dS0dt(t,2,1) = dS0dt(t,1,2) ; dS0dt(t,2,2) = 5.94d-2 ; dS0dt(t,2,3) = 0.59d-2
   dS0dt(t,3,1) = dS0dt(t,1,3) ; dS0dt(t,3,2) = dS0dt(t,2,3) ; dS0dt(t,3,3) = 4.53d-2
   dS0dt(t,4,4) = 2.15d-2 ; dS0dt(t,5,5) = 1.26d-2 ; dS0dt(t,6,6) = 2.50d-2

   ! Pressure derivative at 100 GPa, 300 °K, shifted by 5GPa
   dS0dp(t,1,1) = 3.55d0 ; dS0dp(t,1,2) = 3.24d0 ; dS0dp(t,1,3) = 2.34d0
   dS0dp(t,2,1) = dS0dp(t,1,2) ; dS0dp(t,2,2) = 5.25d0 ; dS0dp(t,2,3) = 2.34d0
   dS0dp(t,3,1) = dS0dp(t,1,3) ; dS0dp(t,3,2) = dS0dp(t,2,3) ; dS0dp(t,3,3) = 5.31d0
   dS0dp(t,4,4) = 1.36d0 ; dS0dp(t,5,5) = 0.78d0 ; dS0dp(t,6,6) = 1.62

   ! Pressure/Temperature second derivative (-GPa^2/GPa/K)
   dS0dpt(t,1,1) = 0.66d-4 ; dS0dpt(t,1,2) = 1.58d-4 ; dS0dpt(t,1,3) = -7.53d-6
   dS0dpt(t,2,1) = dS0dpt(t,1,2) ; dS0dpt(t,2,2) = 2.55d-5 ; dS0dpt(t,2,3) = 1.60d-4
   dS0dpt(t,3,1) = dS0dpt(t,1,3) ; dS0dpt(t,3,2) = dS0dpt(t,2,3) ; dS0dpt(t,3,3) = -3.36d-5
   dS0dpt(t,4,4) =  9.61d-6 ; dS0dpt(t,5,5) = 3.59d-5 ; dS0dpt(t,6,6) = 1.85d-6

!!! MgO ferripericlase - magnesiowuestite ent P-T)

   t = 10

   !(GPa)
   S0(t,1,1) = 300d0 ; S0(t,1,2) =93.6d0  ; S0(t,1,3) = S0(t,1,2)
   S0(t,2,1) = S0(t,1,2)  ; S0(t,2,2) = S0(t,1,1) ; S0(t,2,3) = S0(t,1,2)
   S0(t,3,1) = S0(t,1,2)  ; S0(t,3,2) = S0(t,1,2) ; S0(t,3,3) = S0(t,1,1)
   S0(t,4,4) = 147d0 ; S0(t,5,5) = S0(t,4,4)  ; S0(t,6,6) = S0(t,4,4)

   ! Temperature derivative (from Karki et al., 2000, Phys. Rev. B;
   ! reversed sign of (1,2))
   dS0dt(t,1,1) = 5.98d-2 ; dS0dt(t,1,2) =  0.89d-2 ; dS0dt(t,1,3) = dS0dt(t,1,2)
   dS0dt(t,2,1) = dS0dt(t,1,2) ; dS0dt(t,2,2) = dS0dt(t,1,1) ; dS0dt(t,2,3) = dS0dt(t,1,2)
   dS0dt(t,3,1) = dS0dt(t,1,2) ; dS0dt(t,3,2) = dS0dt(t,1,2) ; dS0dt(t,3,3) = dS0dt(t,1,1)
   dS0dt(t,4,4) = 0.88d-2 ; dS0dt(t,5,5) = dS0dt(t,4,4) ; dS0dt(t,6,6) = dS0dt(t,4,4)

   ! Pressure derivative (from Karki et al., 2000, Phys. Rev. B)
   dS0dp(t,1,1) = 9.56d0 ; dS0dp(t,1,2) = 1.45d0 ; dS0dp(t,1,3) = dS0dp(t,1,2)
   dS0dp(t,2,1) = dS0dp(t,1,2) ; dS0dp(t,2,2) = dS0dp(t,1,1) ; dS0dp(t,2,3) = dS0dp(t,1,2)
   dS0dp(t,3,1) = dS0dp(t,1,2) ; dS0dp(t,3,2) = dS0dp(t,1,2) ; dS0dp(t,3,3) = dS0dp(t,1,1)
   dS0dp(t,4,4) = 1.03d0 ; dS0dp(t,5,5) = dS0dp(t,4,4) ; dS0dp(t,6,6) = dS0dp(t,4,4)

   ! P-T cross derivative (from Karki et al., 2000, Phys. Rev. B)
   ! reversed all signs
   dS0dpt(t,1,1) = -5.6d-4 ; dS0dpt(t,1,2) = 6d-5 ; dS0dpt(t,1,3) = dS0dpt(t,1,2)
   dS0dpt(t,2,1) = dS0dpt(t,1,2) ; dS0dpt(t,2,2) = dS0dpt(t,1,1) ; dS0dpt(t,2,3) = dS0dpt(t,1,2)
   dS0dpt(t,3,1) = dS0dpt(t,1,2) ; dS0dpt(t,3,2) = dS0dpt(t,1,2) ; dS0dpt(t,3,3) = dS0dpt(t,1,1)
   dS0dpt(t,4,4) = -2d-4 ; dS0dpt(t,5,5) = dS0dpt(t,4,4) ; dS0dpt(t,6,6) = dS0dpt(t,4,4)


!!! MgSiO3 Post-perovskite (Zhang et al., 2016, EPSL, interpolation of data)
!!! T0 = 1500 K, P0 = 99.95 GPa

   t = 11 

   !(GPa)
   S0(t,1,1) =1124.75d0 ; S0(t,1,2) =350.69d0  ; S0(t,1,3) = 287.53d0
   S0(t,2,1) = S0(t,1,2)  ; S0(t,2,2) = 884.22d0 ; S0(t,2,3) = 461.05d0
   S0(t,3,1) = S0(t,1,3)  ; S0(t,3,2) = S0(t,2,3) ; S0(t,3,3) =1124.32d0
   S0(t,4,4) = 254.57d0   ; S0(t,5,5) = 231.19d0 ; S0(t,6,6) = 344.05d0

   ! Temperature derivative (-GPa/K) at 99.95 GPa, 1500 K
   dS0dt(t,1,1) = 1.01d-2 ; dS0dt(t,1,2) =-0.915d-2 ; dS0dt(t,1,3) = -0.51d-2
   dS0dt(t,2,1) = dS0dt(t,1,2) ; dS0dt(t,2,2) = 5.24d-2 ; dS0dt(t,2,3) = 1.78d-2
   dS0dt(t,3,1) = dS0dt(t,1,3) ; dS0dt(t,3,2) = dS0dt(t,2,3) ; dS0dt(t,3,3) = 4.90d-2
   dS0dt(t,4,4) = 0.81d-2 ; dS0dt(t,5,5) = 2.33d-2 ; dS0dt(t,6,6) = 2.19d-2

   ! Pressure derivative at 99.95 GPa, 1500 K
   dS0dp(t,1,1) = 5.85d0 ; dS0dp(t,1,2) = 3.06d0; dS0dp(t,1,3) = 3.00d0
   dS0dp(t,2,1) = dS0dp(t,1,2) ; dS0dp(t,2,2) = 2.70d0 ; dS0dp(t,2,3) = 2.33d0
   dS0dp(t,3,1) = dS0dp(t,1,3) ; dS0dp(t,3,2) = dS0dp(t,2,3) ; dS0dp(t,3,3) = 5.44d0
   dS0dp(t,4,4) = 1.37d0 ; dS0dp(t,5,5) = 1.27d0 ; dS0dp(t,6,6) = 2.06

   ! Pressure/Temperature second derivative (GPa^2/GPa/K) 99.95 GPa, 1500 K
   dS0dpt(t,1,1) = -9.37d-4; dS0dpt(t,1,2) = 0.60d-4 ; dS0dpt(t,1,3) = -8.01d-4
   dS0dpt(t,2,1) = dS0dpt(t,1,2) ; dS0dpt(t,2,2) = 7.09d-4 ; dS0dpt(t,2,3) = 6.92d-4
   dS0dpt(t,3,1) = dS0dpt(t,1,3) ; dS0dpt(t,3,2) = dS0dpt(t,2,3) ; dS0dpt(t,3,3) = 7.97d-4
   dS0dpt(t,4,4) = -1.07d-4 ; dS0dpt(t,5,5) = 1.86d-4 ; dS0dpt(t,6,6) = -0.16d-4

   RETURN

   END SUBROUTINE elastic_database
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
