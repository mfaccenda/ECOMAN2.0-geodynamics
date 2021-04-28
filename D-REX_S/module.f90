 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !!
 !!    Copyright (c) 2018-2020, Universita' di Padova, Manuele Faccenda
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

   MODULE comvar

   INTEGER :: rocktype(1)
   DOUBLE PRECISION :: pressure,tkelv,rho(1),fractvoigt,pi,deg2rad
   DOUBLE PRECISION :: dt,timesum,max_strain,strain_max
   DOUBLE PRECISION, DIMENSION(3,3,1) :: Fij
   double precision, dimension(21) :: XEC
   ! finite strain tensor 

   CHARACTER(100) ::output_name 
   ! name of the imput file containing X1, X3, Ui and Dij

!!! LPO calculation

   DOUBLE PRECISION, DIMENSION(1,3,3) :: l,e
   ! velocity gradient tensor and strain rate tensor

   DOUBLE PRECISION :: epsnot(1),mtk(1),mpgpa(1),mx1(1),mx2(1),mx3(1)
   ! reference strain rate

!!! Rock properties

   DOUBLE PRECISION, DIMENSION(4) :: Xol,stressexp,lambda,Mob,chi,fractdislrock,top,bot
   DOUBLE PRECISION, DIMENSION(4,12) :: tau
   INTEGER, DIMENSION(4,2) :: single_crystal_elastic_db
   ! Xol = fraction of anisotropic phase in the aggregate
   ! stressexp = stress exponent for non-Newtonian behavior
   ! lambda = nucleation parameter
   ! Mob = grain mobility
   ! chi = threshold volume fraction for activation of grain boundary sliding
   ! tau = RSS for the slip systems of mineral phases

!!! Dynamic recrystallization

   INTEGER :: size3, size ! size = size3^3
   ! number of points in the (metric) Eulerian space

   DOUBLE PRECISION, DIMENSION(3,3,3) :: alt ! \epsilon_{ijk} tensor
   DOUBLE PRECISION, DIMENSION(3,3) :: del ! \delta_{ij} tensor
   INTEGER, DIMENSION(3,3) :: ijkl ! tensor of indices to form Cijkl from Sij
   INTEGER, DIMENSION(6) :: l1,l2 ! tensot of indices to form Sij from Cijkl

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: odf
   ! volume fraction of the olivine grains
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: odf_ens
   ! volume fraction of the enstatite grains

   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: acs,acs_ens
   DOUBLE PRECISION, DIMENSION(:,:,:),   ALLOCATABLE :: acs0
   !! matrix of direction cosine

   INTEGER degstp,nxy,nz
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phix,phiy,phiz
   ! min/max Vp, dVs

!!! Output

   DOUBLE PRECISION, DIMENSION(3) :: phi_fse
   ! orientation of the long axis of the FSE

   DOUBLE PRECISION :: ln_fse
   ! finite strain = ln(longaxis/shortaxis)

   DOUBLE PRECISION, DIMENSION(3) :: phi_a  
   ! average orientation of a-axis

   DOUBLE PRECISION :: perc_a
   ! percentage of S wave anisotropy

!!! Elastic tensor 

   DOUBLE PRECISION, DIMENSION(6,6) :: Voigt,Reuss,Mixed
   DOUBLE PRECISION, DIMENSION(20,6,6) :: S0,dS0dp,dS0dp2,dS0dt,dS0dpt
   DOUBLE PRECISION, DIMENSION(6,6) :: dS0dp5,dS0dt5,mandel_scale
   DOUBLE PRECISION, DIMENSION(6,6,1) :: Sav
   ! stiffness matrix

!!! SPO parameters

   INTEGER :: spomod,ptmod,eosmod,meltspomod,spograinmod,sporockmod,numminermax,numminer(5),rocknum=5,marknum,fsemod=0
   CHARACTER :: meltfilename              
   DOUBLE PRECISION :: volfractrock(5),ro_back,ro_incl,phi_spo
   DOUBLE PRECISION :: x10size,x20size,X1n(1),X2n(1)
   INTEGER :: thetanum,phinum
   DOUBLE PRECISION :: Vstp,Vmax,stp,esa(3)
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: thetaspo,phispo,meltmark 
   DOUBLE PRECISION, DIMENSION(3,3) :: kron ! tensor of indices to form Cijkl from Sij
   DOUBLE PRECISION, DIMENSION(6,6) :: Cback,Cinc
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: Is,Id 
   !Parameters set in inputspo.dat
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AA,BB,CC,FF,LLL,NN,RRho
   !TI moduli f(P,T) and melt content

!!! Density,Vp,Vs thermodynamic databases
   
   INTEGER :: tnum,pnum
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: rhorock,vsrock,vprock,srock
   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: rhominer,vsminer,vpminer,volfractminer
   !HeFESTo
      
   INTEGER :: tknum,pbnum,ptnum
   DOUBLE PRECISION :: tkmin,pbmin,tkstp,pbstp
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: td_rho,td_vp,td_vs

   !PERPLE_X

   END MODULE comvar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

