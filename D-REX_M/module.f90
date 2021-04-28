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

   !! interval of calculation of the LPO in time      
   INTEGER :: Tinit, Tstep, Tend, OutputStep, tnum, numdt, nt
   ! time step and file number
   DOUBLE PRECISION :: dt,timesum,dt0,timesum0,timemax,pi,deg2rad
   CHARACTER (3) :: dt_to_str
   ! Operating modes
   INTEGER :: fsemod,ptmod,fractdislmod,fabrictransformmod,uppermantlemod,eosmod,fossilfabric

!!! Density thermodynamic database
   
   INTEGER :: tknum,pbnum
   DOUBLE PRECISION :: tkmin,pbmin,tkstp,pbstp
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: td_rho,td_vp,td_vs

!!! Eulerian grid points - Velocity field - Temperature - Pressure - Fraction of dislocation creep

   INTEGER :: invdepthaxis,axisorder(2),dimensions,cartspher,basicstag,x1periodic,x2periodic,x3periodic,yinyang
   INTEGER :: nx1, nx2, nx3 !! grid size in X1 and X2 direction
   INTEGER :: nodenum !! number of Eulerian grid nodes       
   DOUBLE PRECISION :: x1min,x2min,x3min     
   DOUBLE PRECISION :: x1max,x2max,x3max     
   DOUBLE PRECISION :: x1stp,x2stp,x3stp     

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X1,X2,X3
   ! X,Y grid points coordinate

   DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE :: Ui,Ux
   ! velocity vector - U1(x1,x2,x3)=Ui(1,x1,x2,x3)

   DOUBLE PRECISION, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: Dij
   ! velocity gradient tensor - D11(x1,x2,x3)=Dij(1,1,x1,x2,x3)

   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Tk,Pa,Fd
   ! temperature, pressure, fraction of dislocation creep

!!! Lagrangian grid - crystal aggregates   

   INTEGER :: mnx1, mnx2, mnx3 !! number of aggregates in X1 and X2 direction
   INTEGER :: marknum,marknum1 !! number of aggregates
   DOUBLE PRECISION :: mx1min,mx2min,mx3min     
   DOUBLE PRECISION :: mx1max,mx2max,mx3max     
   DOUBLE PRECISION :: mx1stp,mx2stp,mx3stp 
   DOUBLE PRECISION :: mx1minfab,mx1maxfab,mx2minfab,mx2maxfab,mx3minfab,mx3maxfab !! X1-X2-X3 range for setting fossil fabric

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mx1,mx2,mx3
   ! X,Y grid points coordinates 

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rho,max_strain,time
   ! density, cumulative maximum strain

   INTEGER, DIMENSION(:), ALLOCATABLE :: rocktype,mYY
   ! rocktype        

   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Fij
   ! finite strain tensor - F11(marknum)=Fij(1,1,marknum)

!!! Rock properties

   DOUBLE PRECISION, DIMENSION(4) :: Xol,stressexp,lambda,Mob,chi,fractdislrock,minx2,maxx2
   DOUBLE PRECISION, DIMENSION(4,12) :: tau
   INTEGER, DIMENSION(4,2) :: single_crystal_elastic_db
   ! Xol = fraction of anisotropic phase in the aggregate
   ! stressexp = stress exponent for non-Newtonian behavior
   ! lambda = nucleation parameter
   ! Mob = grain mobility
   ! threshold volume fraction for activation of grain boundary sliding
   ! tau = RSS for the slip systems of mineral phases

!!! LPO calculation

   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: l,e
   ! local velocity gradient tensor and strain rate tensor

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: epsnot
   ! reference strain rate

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

!!! Output

!!! Elastic tensor 

   DOUBLE PRECISION, DIMENSION(20,6,6) :: S0,dS0dp,dS0dp2,dS0dt,dS0dpt
   DOUBLE PRECISION, DIMENSION(6,6) :: dS0dp5,dS0dt5,mandel_scale
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Sav
   DOUBLE PRECISION :: fractvoigt
   ! stiffness matrix

   END MODULE comvar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

