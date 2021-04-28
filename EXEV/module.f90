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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Module of common variables                                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   MODULE comvar

   INTEGER :: Tinit,Tend,Tstep,dimensions,cartspher,basicstag,x1periodic,x2periodic,x3periodic,yinyang
   INTEGER :: Lagrangian,uppermantlemod,rocktypemod,fsemod,fse3Dmod,fseminmod,fsemaxmod,TIaxismod,vpmaxmod,dvsmaxmod
   INTEGER :: Eulerian,vpvsmod,dvpvsmod,zoeppritzmod,radialmod,azimod,reflectmod,replicateZmod,specfem3Dmod,syntomomod 
   DOUBLE PRECISION :: ln_fse_min,koef,cosx1,cosx2,cosx3,qhat(3),Incangle,pi,deg2rad
   DOUBLE PRECISION :: aziscalex1,aziscalex2,aziscalex3
   DOUBLE PRECISION :: dt,timesum
   ! File number and domain geometry

!!! Grid points of the domain where to interpolate the elastic tensors

   INTEGER :: nodenum,cellnum !! number of nodes,cells  
   INTEGER :: nx1ref,nx3ref !! horizontal plane indexes for reference Vp, Vs, rho profile
   INTEGER :: nx11, nx21, nx31 !! Number of nodes of the velocity model
   DOUBLE PRECISION :: nstepx1, nstepx2, nstepx3
   DOUBLE PRECISION :: n1first,n2first,n3first,n1last,n2last,n3last
   !New Tomographic grid
   INTEGER :: nx1, nx2, nx3 !! Number of nodes of the Drex model 
   DOUBLE PRECISION :: x10size, x20size, x30size
   DOUBLE PRECISION :: i1first,i2first,i3first,i1last,i2last,i3last
   DOUBLE PRECISION :: x1max,x2max,x3max
   !Old D-REX_S grid

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X1,X1n
   ! X grid points coordinates 
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X2,X2n
   ! Y grid points coordinates 
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X3,X3n
   ! Z grid points coordinates 

!!! D-REX_S aggregates distribution and properties   

   INTEGER :: marknum,marknum0 !! number of markers
   DOUBLE PRECISION :: mx1min,mx2min,mx3min     
   DOUBLE PRECISION :: mx1max,mx2max,mx3max     
   DOUBLE PRECISION :: mx1stp,mx2stp,mx3stp
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mx1
   ! X grid points coordinates 
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mx2
   ! Y grid points coordinates 
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mx3
   ! Z grid points coordinates 
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ln_fse,amin,amed,amax
   ! finite strain = ln(longaxis/shortaxis)

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rho,mtk,mpgpa
   !Density, temperature, pressure,dummy rocktype
   INTEGER, DIMENSION(:), ALLOCATABLE :: rocktype,rt1,mYY
   !Rocktype        

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: epsnot
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: l,e
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Fij,Rotm
   DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE :: Ui,Ux
   DOUBLE PRECISION, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: Dij

!!! Drex Parameters

   DOUBLE PRECISION, DIMENSION(5) :: top,bot
   DOUBLE PRECISION, DIMENSION(5) :: Xol,stressexp,lambda,Mob,chi,fractdislrock
   DOUBLE PRECISION, DIMENSION(5,12) :: tau
   ! Xol = fraction of anisotropic phase in the aggregate
   ! stressexp = stress exponent for non-Newtonian behavior
   ! lambda = nucleation parameter
   ! Mob = grain mobility
   ! threshold volume fraction for activation of grain boundary sliding
   ! tau = RSS for the slip systems of mineral phases
 
!!! Min/max Vp,dVs   

   INTEGER degstp,nxy,nz
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phix,phiy,phiz

!!! Elastic tensor 

   DOUBLE PRECISION, DIMENSION(20,6,6) :: S0,dS0dp,dS0dp2,dS0dt,dS0dpt
   DOUBLE PRECISION, DIMENSION(6,6) :: dS0dp5,dS0dt5,mandel_scale
   INTEGER, DIMENSION(3,3) :: ijkl ! tensor of indices to form Cijkl from Sij
   INTEGER, DIMENSION(6) :: l1,l2 ! tensot of indices to form Sij from Cijkl
   DOUBLE PRECISION, DIMENSION (:,:,:), ALLOCATABLE :: Sav,Savn,Sav0
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: Rhon,fsen       
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: Tk, Pa
   DOUBLE PRECISION ,DIMENSION(:,:), ALLOCATABLE :: X1save
   ! stiffness matrix
   INTEGER :: size3, size ! size = size3^3

!!! SPO parameters

   INTEGER :: spomod,ptmod,meltspomod,spograinmod,sporockmod,numminermax,numminer(5),rocknum=5
   DOUBLE PRECISION :: volfractrock(5),ro_back,ro_incl,phi_spo,etac1,volinc1
   CHARACTER(200) :: meltfilename
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

   INTEGER :: tnum,pnum,ptnum
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: rhorock,vsrock,vprock,srock
   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: rhominer,vsminer,vpminer,volfractminer
   !HeFESTo

   INTEGER :: tknum,pbnum
   DOUBLE PRECISION :: tkmin,pbmin,tkstp,pbstp
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: td_rho
   !PERPLE_X

   END MODULE comvar

