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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine eulerian_output, saving Vp,Vs,dVp,dVs,radial & azimuthal anisi, etc. !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE eulerian_output(dt_str4,cijkl_dir,output_dir)

USE comvar
USE hdf5   

IMPLICIT NONE

INTEGER :: i,j,k,i1,i2,i3,m,m1,gi,giref,gijk,ijk,x,y,z,t,zn,yyy
INTEGER :: yc,zc,dum_int(3),azinum,flag
DOUBLE PRECISION, DIMENSION(21) :: XE
DOUBLE PRECISION, DIMENSION (3) :: n1,n2,n3   
DOUBLE PRECISION, DIMENSION (3,3) :: acs     
DOUBLE PRECISION, DIMENSION (3,3,3,3) :: cijkl
DOUBLE PRECISION :: PwaveMod,SwaveMod,Nparam,Lparam,Gc,Gs,G3D,phi,lr,cr
DOUBLE PRECISION :: xx,yy,zz,phi1,theta,phi2,azi(3)
DOUBLE PRECISION :: Vp01,Vp02,Vs01,Vs02,ro1,ro2,reflPP,reflSS,reflPS
DOUBLE PRECISION :: Pcritangle
DOUBLE PRECISION :: a,b,c,d,E1,F,G1,H,D1,p,k1,k2,j1,j2
!Buffers
INTEGER, DIMENSION(:), ALLOCATABLE :: idx
INTEGER, DIMENSION(:,:), ALLOCATABLE :: connectivity 
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dVp,dVs,Vp,Vs,E_chi,G,Vpref,Vsref,Rhoref,Rhon0,fsen0,depth
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Rpp,Rps,Tpp,Tps,ERpp,ERps,ETpp,ETps
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: vertices,azimuthal,xyz 
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Savref

CHARACTER (500) :: filename,filenamexmf
CHARACTER (len=*) :: cijkl_dir,output_dir
CHARACTER (len(trim(cijkl_dir))) :: str1
CHARACTER (len(trim(output_dir))) :: str2
CHARACTER (4) :: dt_str4
INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl,memtype ! Handles
INTEGER     ::   error  ! Error flag
   
str1=trim(cijkl_dir)
str2=trim(output_dir)

yyy=yinyang

IF(nx11 < 2) nx11 = 2 
IF(nx21 < 2) nx21 = 2 
IF(dimensions == 3 .AND. nx31 < 2) nx31 = 2 

! Define velocity model grid where to interpolate Sav and rho
! The grid spacing should be >= to initial aggregate spacing set in DREX
nstepx1=(n1last-n1first)/(nx11 - 1)
nstepx2=(n2last-n2first)/(nx21 - 1)
nstepx3=(n3last-n3first)/(nx31 - 1)
    
nodenum=nx11*nx21*nx31
cellnum=(nx11-1)*(nx21-1)*(nx31-1)
!Check for the 3rd dimension in 2D models
IF(dimensions == 2 .AND. replicateZmod == 0) THEN
   nx31 = 1 ; nx3ref = 1
   nstepx3 = 0d0
   n3first = pi/2.0 !Set 2D model on the equatorial plane
   nodenum=nx11*nx21*nx31
   cellnum=(nx11-1)*(nx21-1)
END IF
IF(dimensions == 2 .AND. replicateZmod > 0 .AND. nx31 == 1) THEN
   print *,'ReplicateZmod active, but nx31 is 1. Increase nx31'
   stop
END IF

IF(nodenum > 16777216) write(*,"(a,i10,a)"),' WARNING : nodenum = ',nodenum,' > 2^24 = 16777216. There is an issue with Paraview reading grids with such a large number of nodes from HDF5 files. Decrease the grid resolution' 

!Allocate arrays
ALLOCATE(Savn(6,6,nodenum*yyy),Rhon(nodenum*yyy),fsen(nodenum*yyy))
ALLOCATE(X1(nx11),X2(nx21),X3(nx31))

!Grid coordinates
DO i=1,nx11
   X1(i)=n1first + (i-1)*nstepx1
END DO
n1last = X1(nx11)

DO i=1,nx21
   X2(i)=n2first + (i-1)*nstepx2
END DO
n2last = X2(nx21)

DO i=1,nx31
   X3(i)=n3first + (i-1)*nstepx3
END DO
n3last = X3(nx31)

!Check Eulerian grid if compatible with DREX_S model
IF(cartspher == 1 .AND. nstepx1 < mx1stp)  WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : nstepx1 = ',nstepx1,' < mx1stp = ',mx1stp,'. Decrease nx11'
IF(cartspher == 2 .AND. nstepx1*X2(1)*sin(MINVAL(X3)) < mx1stp)  WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : nstepx1 = ',nstepx1*X2(1)*sin(MINVAL(X3)),' < mx1stp = ',mx1stp,'. Decrease nx11'
IF(nstepx2 < mx2stp) WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : nstepx2 = ',nstepx2,' < mx2stp = ',mx2stp,'. Decrease nx21'
IF(reflectmod /= -1 .AND. n1first < i1first*0.9999) WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : n1first = ',n1first,' < i1first = ',i1first,'. Set n1first >= i1first'
IF(reflectmod /=  1 .AND. n1last  > i1last*1.0001 ) WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : n1last  = ',n1last ,' > i1last  = ',i1last ,'. Set n1last  >= i1last '
IF(n2first < i2first*0.9999) WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : n2first = ',n2first,' < i2first = ',i2first,'. Set n2first >= i2first'
IF(n2last  > i2last*1.0001 ) WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : n2last  = ',n2last ,' > i2last  = ',i2last ,'. Set n2last  >= i2last '

IF(dimensions == 3) THEN
   IF(cartspher == 1 .AND. nstepx3 < mx3stp*0.9999)  WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : nstepx3 = ',nstepx3,' < mx3stp = ',mx3stp,'. Decrease nx31'
   IF(cartspher == 2 .AND. nstepx3*X2(1) < mx3stp*0.9999)  WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : nstepx3 = ',nstepx3*X2(1),' < mx3stp = ',mx3stp,'. Decrease nx31'
   IF(reflectmod /= -3 .AND. n3first < i3first*0.9999) WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : n3first = ',n3first,' < i3first = ',i3first,'. Set n3first >= i3first'
   IF(reflectmod /=  3 .AND. n3last  > i3last*1.0001 ) WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : n3last  = ',n3last ,' > i3last  = ',i3last ,'. Set n3last  >= i3last '
   IF(reflectmod /=  0 .AND. yinyang == 2) THEN 
      WRITE(*,"(a)"),' WARNING : YIN-YANG grid ACTIVE!. Set reflectmod = 0'
      reflectmod = 0
   END IF
END IF

!Print grid infos
if(cartspher==1) then
   write(*,'(a,i0,a,i0,a,i0)',advance='no') ' Tomographic grid with number of nodes: X = ',nx11,', Y = ',nx21,', Z = ',nx31
   write(*,'(a,i0)') ', Number of cells = ',cellnum 
   print *
   write(*,'(a)') ' Domain size is:'
   write(*,'(a,f12.4,a,f12.4,a)') ' X min = ',n1first/1d3,'; X max = ',n1last/1d3,' (km)'
   write(*,'(a,f12.4,a,f12.4,a)') ' Y min = ',n2first/1d3,'; Y max = ',n2last/1d3,' (km)'
   write(*,'(a,f12.4,a,f12.4,a)') ' Z min = ',n3first/1d3,'; Z max = ',n3last/1d3,' (km)'
   print *
end if
if(cartspher==2) then
   write(*,'(a,i0,a,i0,a,i0)',advance='no') ' Tomographic grid with number of nodes: Long = ',nx11,', Depth = ',nx21,', Colat = ',nx31
   write(*,'(a,i0)') ', Number of cells = ',cellnum 
   print *
   write(*,'(a)') ' Domain size is:'
   write(*,'(a,f10.2,a,f10.2,a)') ' Long   min = ',n1first/deg2rad,'; Long   max = ',n1last/deg2rad,' (deg)'
   write(*,'(a,f10.2,a,f10.2,a)') ' Radius min = ',n2first/1d3,    '; Radius max = ',n2last/1d3,' (km)'
   write(*,'(a,f10.2,a,f10.2,a)') ' Colat  min = ',n3first/deg2rad,'; Colat  max = ',n3last/deg2rad,' (deg)'
   print *

   if(yinyang==2) then
      print *,'YIN-YANG grids ACTIVE!' 
      print *
   end if

end if
      
!!! Interpolate aggregates elastic tensor to the tomographic model grid by (bi-/tri-)linear interpolation
CALL mark2node

x=1
y=nx11
z=nx11*nx21
yc=nx11-1
zc=(nx11-1)*(nx21-1)

!Replicate along Z (Colat) direction
IF(dimensions == 2 .AND. replicateZmod > 0 .AND. nx31 > 1) THEN

   ALLOCATE(Sav0(6,6,z),Rhon0(z),fsen0(z))

   !Save 2D reference plane
   DO i2=1,nx21
   DO i1=1,nx11
      gi   = i1 + (i2-1)*y 
      Sav0(:,:,gi) = Savn(:,:,gi)
      Rhon0(gi) = Rhon(gi)
      fsen0(gi) = fsen(gi)
   END DO;END DO

   !$omp parallel & 
   !$omp shared(Savn,Rhon,fsen,Sav0,Rhon0,fsen0,X3) &
   !$omp private(i1,i2,i3,gi,gijk,phi1,theta,phi2,acs) &    
   !$omp firstprivate(cartspher,pi,x,y,z,nx11,nx21,nx31,i1first,i1last,i3first,i3last)
   !$omp do schedule(guided,8)
   DO i2=1,nx21
   DO i1=1,nx11
      gi   = i1 + (i2-1)*y 
      DO i3=1,nx31

         gijk = i1 + (i2-1)*y + (i3-1)*z

         Savn(:,:,gijk) = Sav0(:,:,gi)
         Rhon(gijk) = Rhon0(gi)
         fsen(gijk) = fsen0(gi)

         IF(cartspher == 2) THEN

            !Rotate the elastic tensor back with respect to colat (X3(i3)) around X axis
            phi1 = 0d0 ; theta = +(X3(i3) - pi/2.0); phi2 = 0d0

            !Transform Euler angles into direction cosine matrix
            CALL rotmatrixZXZ(phi1,theta,phi2,acs)

            CALL rotvoigt(Savn(:,:,gijk),acs,Savn(:,:,gijk))

         END IF

      END DO
   END DO;END DO
   !$omp end do
   !$omp end parallel
  
   DEALLOCATE(Sav0,Rhon0,fsen0)

END IF
!END Replicate along Z (Colat) direction

!!! Reflect aggregates
IF(reflectmod /= 0) THEN

   acs = 0 
   acs(1,1) = 1; acs(2,2) = 1; acs(3,3) = 1 ; 
   kron = acs
   IF(cartspher==1) THEN

      IF(ABS(reflectmod) == 1) acs(1,1) =-1 !Reflect with respect to plane normal to x1
      IF(ABS(reflectmod) == 3) acs(3,3) =-1 !Reflect with respect to plane normal to x3

   ELSE

      !Reflect with respect to plane normal to x1 
      IF(ABS(reflectmod) == 1) THEN

         !Find normal to reflection plane
         IF(reflectmod==-1) THEN
            n3(1) = cos(i1first+pi/2.0); n3(2) = sin(i1first+pi/2.0); n3(3) = 0
         END IF
         IF(reflectmod== 1) THEN
            n3(1) = cos(i1last +pi/2.0); n3(2) = sin(i1last +pi/2.0); n3(3) = 0 
         END IF
      
      END IF

      !Reflect with respect to plane normal to x3 
      IF(dimensions == 3 .AND. ABS(reflectmod) == 3) THEN

         !Find normal to reflection plane
         IF(reflectmod==-3) THEN
            n3(1) = sin(i3first+pi/2.0)*cos((i1first+i1last)/2.0); n3(2) = sin(i3first+pi/2.0)*sin((i1first+i1last)/2.0); n3(3) = cos(i3first+pi/2.0)
         END IF
         IF(reflectmod== 3) THEN
            n3(1) = sin(i3last +pi/2.0)*cos((i1first+i1last)/2.0); n3(2) = sin(i3last +pi/2.0)*sin((i1first+i1last)/2.0); n3(3) = cos(i3last +pi/2.0)
         END IF
      
      END IF

      ! R = I - 2*n*n, where n is the normal to the reflection surface
      DO i1=1,3 ; DO i2=1,3
         acs(i1,i2) = kron(i1,i2) - 2d0*n3(i2)*n3(i1)
      END DO; END DO

   END IF

   !$omp parallel & 
   !$omp shared(Savn,Rhon,fsen,X1,X3) &
   !$omp private(i1,i2,i3,gi,gijk,flag) &    
   !$omp firstprivate(reflectmod,acs,x,y,z,nx11,nx21,nx31,i1first,i1last,i3first,i3last)
   !$omp do schedule(guided,8)
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      flag = 0
      IF      (reflectmod ==-1 .AND. X1(i1) < i1first) THEN
         gi   = i1 + (i2-1)*y + (i3-1)*z 
         gijk = nx11 - (i1-1) + (i2-1)*y + (i3-1)*z
         flag = 1
      ELSE IF (reflectmod == 1 .AND. X1(i1) < i1last ) THEN
         gijk = i1 + (i2-1)*y + (i3-1)*z 
         gi   = nx11 - (i1-1) + (i2-1)*y + (i3-1)*z
         flag = 1
      ELSE IF (reflectmod ==-3 .AND. X3(i3) < i3first) THEN
         gi   = i1 + (i2-1)*y + (i3-1)*z 
         gijk = i1 + (i2-1)*y + (nx31-i3)*z
         flag = 1
      ELSE IF (reflectmod == 3 .AND. X3(i3) < i3last ) THEN
         gijk = i1 + (i2-1)*y + (i3-1)*z 
         gi   = i1 + (i2-1)*y + (nx31-i3)*z
         flag = 1
      END IF          
 
      IF(flag /=0 ) THEN
         fsen(gi) = fsen(gijk)
         Rhon(gi) = Rhon(gijk)
         !Mirror Cijkl tensor
         Savn(:,:,gi) = Savn(:,:,gijk)
         CALL rotvoigt(Savn(:,:,gi),acs,Savn(:,:,gi))
      END IF

   END DO;END DO;END DO
   !$omp end do
   !$omp end parallel

END IF
!END Reflect aggregates

!Find the average density and elastic tensor at each depth
ALLOCATE(Rhoref(nx21),Savref(6,6,nx21))

DO i2=1,nx21
   Rhoref(i2)=0d0
   Savref(:,:,i2)=0d0
   zn=0
   DO i1=1,nx11
   DO i3=1,nx31
      gi = i1 + (i2-1)*y + (i3-1)*z
      IF(Rhon(gi) > 0) THEN
         Rhoref(i2)= Rhoref(i2) + Rhon(gi)
         Savref(:,:,i2)= Savref(:,:,i2) + Savn(:,:,gi)
         zn=zn+1
      END IF
   END DO
   END DO
   Rhoref(i2)=Rhoref(i2)/REAL(zn)
   Savref(:,:,i2)=Savref(:,:,i2)/REAL(zn)
END DO

!Fill empty nodes
!$omp parallel & 
!$omp shared(Savn,Rhon,Savref,Rhoref) &
!$omp private(yyy,i1,i2,i3,gi,PwaveMod,SwaveMod,XE) &    
!$omp firstprivate(x,y,z,nx11,nx21,nx31)
!$omp do schedule(guided,8)
DO yyy=1,yinyang
DO i3=1,nx31
DO i2=1,nx21
DO i1=1,nx11

   gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z

   IF(Rhon(gi) == 0) THEN     

      Rhon(gi) = Rhoref(i2)
      Savn(:,:,gi) = Savref(:,:,i2)

      PwaveMod = 0d0; SwaveMod = 0d0
      !Set isotropic elastic tensor
      XE = 0d0
      call V21D(Savn(:,:,gi),XE)
      PwaveMod = 3D0/15D0*(XE(1)+XE(2)+XE(3))+sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+2D0/15D0*(XE(7)+XE(8)+XE(9))
      SwaveMod = (2D0/15D0*(XE(1)+XE(2)+XE(3))-sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+1D0/5D0*(XE(7)+XE(8)+XE(9)))/2
      Savn(:,:,gi) = 0d0
      Savn(1,1,gi) = PwaveMod; Savn(2,2,gi) = Savn(1,1,gi) ; Savn(3,3,gi) = Savn(1,1,gi)  
      Savn(1,2,gi) = PwaveMod - 2*SwaveMod ; Savn(2,1,gi) = Savn(1,2,gi)   
      Savn(1,3,gi) = Savn(1,2,gi) ; Savn(3,1,gi) = Savn(1,2,gi)   
      Savn(2,3,gi) = Savn(1,2,gi) ; Savn(3,2,gi) = Savn(1,2,gi)   
      Savn(4,4,gi) = SwaveMod; Savn(5,5,gi) = Savn(4,4,gi) ; Savn(6,6,gi) = Savn(4,4,gi) 

   END IF
   
END DO;END DO;END DO;END DO
!$omp end do
!$omp end parallel

!END Fill empty nodes

IF(specfem3Dmod > 0) THEN

   print *,' Save input file for specfem3D globe'
   print *,''
   
   !Depth must be 0 at the surface and increasing upward
   ALLOCATE(depth(nx21))
   DO i2=1,nx21
      depth(i2) = n2last - X2(i2)
   END DO

   filenamexmf = str2//'drex_model.xyz'

   OPEN(20,file=filenamexmf,status='replace')
   WRITE(20,'(a)') '#    NX      NY      NZ'
   WRITE(20,'(3i)'),nx11,nx31,nx21
   WRITE(20,'(a)') '#    Lon'
   WRITE(20,'(3000f12.3)'),X1
   WRITE(20,'(a)') '#    Lat'
   WRITE(20,'(3000f12.3)'),X3
   WRITE(20,'(a)') '#    Depth'
   WRITE(20,'(3000f12.3)'),depth
   WRITE(20,'(a)') '#    RHO          C11          C12          C13          C14          C15          C16          C22          C23          C24          C25          C26          C33          C34          C35          C36          C44          C45          C46          C55          C56          C66'

   !filenamexmf = str2//'drex_model.bin'

   !OPEN(20,file=filenamexmf,status='replace',form='unformatted')
   !WRITE(20),nx11,nx31,nx21
   !WRITE(20),X1
   !WRITE(20),X3
   !WRITE(20),depth

   DO yyy=1,yinyang
   DO i2=1,nx21
   DO i3=1,nx31
   DO i1=1,nx11

      gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z

      !In specfem3D_globe, X is east (long), Y is North (colat) and Z is radial distance positive upward
      !As a consequence, need to swap Y and Z axes
      WRITE(20,'(1f8.1,21f12.3)'),Rhon(gi),&
      Savn(1,1,gi),Savn(1,3,gi),Savn(1,2,gi),Savn(1,4,gi),Savn(1,6,gi),Savn(1,5,gi),&
                   Savn(3,3,gi),Savn(3,2,gi),Savn(3,4,gi),Savn(3,6,gi),Savn(3,5,gi),&
                                Savn(2,2,gi),Savn(2,4,gi),Savn(2,6,gi),Savn(2,5,gi),&
                                             Savn(4,4,gi),Savn(4,6,gi),Savn(4,5,gi),&
                                                          Savn(6,6,gi),Savn(6,5,gi),&
                                                                       Savn(5,5,gi)  

   END DO;END DO;END DO; END DO

   CLOSE(20)

   DEALLOCATE(depth)

END IF

IF(syntomomod > 0) THEN

   print *,' Save input file for SYNTOMO'            
   print *,''
   
   ! Create a new file using default properties.
   filename = str2//'syntomo'//dt_str4//'.h5'
  
   !Initialize FORTRAN interface.
   CALL H5open_f (error)

   CALL H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

   CALL loadsave_double(0,1,file_id,nx11,H5T_NATIVE_DOUBLE, X1,'X1',1) ! x axis
   CALL loadsave_double(0,1,file_id,nx21,H5T_NATIVE_DOUBLE, X2,'X2',1) ! y axis
   CALL loadsave_double(0,1,file_id,nx31,H5T_NATIVE_DOUBLE, X3,'X3',1) ! z axis
   CALL loadsave_double(0,1,file_id,nodenum*yyy,H5T_NATIVE_DOUBLE,Rhon,'Rho',1)

   ALLOCATE(X1save(21,nodenum*yyy))
   DO m=1,nodenum*yyy
      !Save elastic tensor: only 21 elastic moduli
       X1save(1,m)  = Savn(1,1,m)
       X1save(2,m)  = Savn(2,2,m)
       X1save(3,m)  = Savn(3,3,m)
       X1save(4,m)  = Savn(2,3,m)
       X1save(5,m)  = Savn(1,3,m)
       X1save(6,m)  = Savn(1,2,m)
       X1save(7,m)  = Savn(4,4,m)
       X1save(8,m)  = Savn(5,5,m)
       X1save(9,m)  = Savn(6,6,m)
       X1save(10,m) = Savn(1,4,m)
       X1save(11,m) = Savn(2,5,m)
       X1save(12,m) = Savn(3,6,m)
       X1save(13,m) = Savn(3,4,m)
       X1save(14,m) = Savn(1,5,m)
       X1save(15,m) = Savn(2,6,m)
       X1save(16,m) = Savn(2,4,m)
       X1save(17,m) = Savn(3,5,m)
       X1save(18,m) = Savn(1,6,m)
       X1save(19,m) = Savn(5,6,m)
       X1save(20,m) = Savn(4,6,m)
       X1save(21,m) = Savn(4,5,m)
   END DO

   dum_int(1)=21
   dum_int(2)=nodenum*yyy
   CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,X1save,'Sav',1)
   DEALLOCATE(X1save)

   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   CALL H5close_f(error)
  
END IF

!!! Rotate tensor toward Z axis direction if polar/spherical coordinates
IF(cartspher==2 .AND. (vpvsmod == 2 .OR. radialmod > 0 .OR. azimod > 0 .OR. zoeppritzmod > 0)) THEN

   !$omp parallel & 
   !$omp shared(Savn,Rhon,fsen,X1,X3) &
   !$omp private(yyy,i1,i2,i3,gi,phi1,theta,phi2,acs,lr,cr) &    
   !$omp firstprivate(pi,ln_fse_min,x,y,z,nx11,nx21,nx31)
   !$omp do schedule(guided,8)
   DO yyy=1,yinyang
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z
      
      IF(Rhon(gi) > 0 .AND. fsen(gi) >= ln_fse_min) THEN

         !Rotate the elastic tensor to Z axis 
         phi1 = 0d0; theta = -X3(i3); phi2 = -X1(i1) + pi*1.5d0  
         IF(yyy == 2) THEN
            CALL yin2yang(X1(i1),X3(i3),lr,cr)
            theta = -cr ; phi2 = -lr + pi*1.5d0
         END IF

         !Transform Euler angles into direction cosine matrix
         CALL rotmatrixZXZ(phi1,theta,phi2,acs)

         CALL rotvoigt(Savn(:,:,gi),acs,Savn(:,:,gi))

      END IF

   END DO;END DO;END DO;END DO
   !$omp end do
   !$omp end parallel
     
END IF
     
!!! Calculate Vp and Vs
IF(vpvsmod > 0 ) THEN
      
   ALLOCATE(Vp(nodenum*yyy),Vs(nodenum*yyy))

   Vp = 0d0 ; Vs = 0d0 
   
   !$omp parallel & 
   !$omp shared(Savn,Rhon,fsen,Vp,Vs,X1,X3) &
   !$omp private(yyy,i1,i2,i3,gi,PwaveMod,SwaveMod,XE) &    
   !$omp firstprivate(vpvsmod,pi,ln_fse_min,x,y,z,nx11,nx21,nx31)
   !$omp do schedule(guided,8)
   DO yyy=1,yinyang
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z
      
      IF(Rhon(gi) > 0) THEN     
     
         IF(vpvsmod == 2 .AND. fsen(gi) >= ln_fse_min) THEN
     
            !Calculate Vp,Vs for a given direction
            CALL SLOWNESS(gi,PwaveMod,SwaveMod)

         END IF

         IF(vpvsmod == 1 .OR. fsen(gi) < ln_fse_min) THEN
     
            PwaveMod = 0d0; SwaveMod = 0d0 
            !Isotropic Vp,Vs, useful to see only thermal anomalies
            !Calculate isotropic Vp and Vs by finding isotropic part of the elastic tensor
            XE = 0d0 
            call V21D(Savn(:,:,gi),XE)
            PwaveMod = 3D0/15D0*(XE(1)+XE(2)+XE(3))+sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+2D0/15D0*(XE(7)+XE(8)+XE(9))
            SwaveMod = (2D0/15D0*(XE(1)+XE(2)+XE(3))-sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+1D0/5D0*(XE(7)+XE(8)+XE(9)))/2d0
     
         END IF
     
         Vp(gi) = (PwaveMod*1e+9/Rhon(gi))**0.5d0
         Vs(gi) = (SwaveMod*1e+9/Rhon(gi))**0.5d0

      END IF
   
   END DO;END DO;END DO;END DO
   !$omp end do
   !$omp end parallel

   !Calculate dVp and dVs
   IF(dvpvsmod > 0) THEN

   !Find the average Vp and Vs velocities at each depth
   ALLOCATE(Vpref(nx21),Vsref(nx21))

   DO i2=1,nx21
      Vpref(i2)=0d0
      Vsref(i2)=0d0
      zn=0
      DO i1=1,nx11
      DO i3=1,nx31
         gi = i1 + (i2-1)*y + (i3-1)*z
         IF(Rhon(gi) > 0) THEN
            Vpref(i2)= Vpref(i2) + Vp(gi)
            Vsref(i2)= Vsref(i2) + Vs(gi)
            zn=zn+1
         END IF
      END DO
      END DO
      Vpref(i2)=Vpref(i2)/REAL(zn)
      Vsref(i2)=Vsref(i2)/REAL(zn)
   END DO

   ALLOCATE(dVp(nodenum*yyy),dVs(nodenum*yyy))

   dVp = 0d0 ; dVs = 0d0
   
   !$omp parallel & 
   !$omp shared(Savn,Rhon,fsen,Vp,Vs,dVp,dVs,Vpref,Vsref,Rhoref,X1,X3) &
   !$omp private(yyy,i1,i2,i3,gi,giref) &    
   !$omp firstprivate(dvpvsmod,nx1ref,nx3ref,x,y,z,nx11,nx21,nx31)
   !$omp do schedule(guided,8)
   DO yyy=1,yinyang
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z
      giref = nx1ref + (i2-1)*y + (nx3ref-1)*z

      IF(Rhon(gi) > 0) THEN     
     
         !Compare with velocities from reference vertical profile
         IF(dvpvsmod > 1) THEN
            IF(Vp(gi) > 0 .AND. Vp(giref) > 0) dVp(gi) = (Vp(gi) - Vp(giref))/Vp(giref)*100d0
            IF(Vs(gi) > 0 .AND. Vs(giref) > 0) dVs(gi) = (Vs(gi) - Vs(giref))/Vs(giref)*100d0
         !Compare with mean velocities at given depth
         ELSE
            IF(Vp(gi) > 0 .AND. Vpref(i2) > 0) dVp(gi) = (Vp(gi) - Vpref(i2))/Vpref(i2)*100d0
            IF(Vs(gi) > 0 .AND. Vsref(i2) > 0) dVs(gi) = (Vs(gi) - Vsref(i2))/Vsref(i2)*100d0
         END IF

      !If empty, assign reference values 
      ELSE

         dVp(gi) = 0d0
         dVs(gi) = 0d0

      END IF
   
   END DO;END DO;END DO;END DO
   !$omp end do
   !$omp end parallel

   DEALLOCATE(Vpref,Vsref)

   END IF
   !End dvpvsmod

END IF
!End vpvsmod

DEALLOCATE(Rhoref,Savref)

IF(zoeppritzmod .NE. 0) THEN
   
   ALLOCATE(Rpp(nodenum*yyy),Rps(nodenum*yyy),Tpp(nodenum*yyy),Tps(nodenum*yyy))
   ALLOCATE(ERpp(nodenum*yyy),ERps(nodenum*yyy),ETpp(nodenum*yyy),ETps(nodenum*yyy))

   Rpp = 0d0 ; Rps = 0d0 ; Tpp = 0d0 ; Tps = 0d0  
   ERpp = 0d0 ; ERps = 0d0 ; ETpp = 0d0 ; ETps = 0d0
   
   !$omp parallel & 
   !$omp shared(Savn,Rhon,fsen,Rpp,Rps,Tpp,Tps,ERpp,ERps,ETpp,ETps) &
   !$omp private(yyy,i1,i2,i3,gi,PwaveMod,SwaveMod,XE,ro1,ro2,Vp01,Vp02,Vs01,Vs02,Pcritangle,p,k1,k2,j1,j2,a,b,c,d,E1,F,G1,H,D1) &    
   !$omp firstprivate(zoeppritzmod,ln_fse_min,x,y,z,nx11,nx21,nx31,Incangle)
   !$omp do schedule(guided,8)
   DO yyy=1,yinyang
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z
    
      !Set number of nodes to move along vertical direction 
      IF(Rhon(gi) > 0) THEN     
     
         PwaveMod = 0d0; SwaveMod = 0d0 ; ro1 = 0d0 ; ro2 = 0d0
         Vp01 = 0d0 ; Vs01 = 0d0 ; Vp02 = 0d0 ; Vs02 = 0d0
         !Set axis of vertical direction 
         if(i2 > 1 .and. Rhon(gi) > 0d0 .and. Rhon(gi-y) > 0d0) then

            !Calculate Vp, Vs
            if(zoeppritzmod == 2 .AND. fsen(gi) >= ln_fse_min) then
               CALL SLOWNESS(gi-y,PwaveMod,SwaveMod)
            end if
            if(zoeppritzmod == 1 .OR. fsen(gi) < ln_fse_min) then
               XE = 0d0
               call V21D(Savn(:,:,gi-y),XE)
               PwaveMod = 3D0/15D0*(XE(1)+XE(2)+XE(3))+sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+2D0/15D0*(XE(7)+XE(8)+XE(9))
               SwaveMod = (2D0/15D0*(XE(1)+XE(2)+XE(3))-sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+1D0/5D0*(XE(7)+XE(8)+XE(9)))/2
            end if

            Vp01 = (PwaveMod*1e+9/Rhon(gi-y))**0.5d0
            Vs01 = (SwaveMod*1e+9/Rhon(gi-y))**0.5d0
            ro1 = Rhon(gi-y)

            !Calculate Vp, Vs
            if(zoeppritzmod == 2 .AND. fsen(gi) >= ln_fse_min) then
               CALL SLOWNESS(gi,PwaveMod,SwaveMod)
            end if
            if(zoeppritzmod == 1 .OR. fsen(gi) < ln_fse_min) then
               XE = 0d0
               call V21D(Savn(:,:,gi),XE)
               PwaveMod = 3D0/15D0*(XE(1)+XE(2)+XE(3))+sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+2D0/15D0*(XE(7)+XE(8)+XE(9))
               SwaveMod = (2D0/15D0*(XE(1)+XE(2)+XE(3))-sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+1D0/5D0*(XE(7)+XE(8)+XE(9)))/2d0
            end if

            Vp02 = (PwaveMod*1e+9/Rhon(gi))**0.5d0
            Vs02 = (SwaveMod*1e+9/Rhon(gi))**0.5d0
            ro2 = Rhon(gi)

            !Zoeppritz's equation
            !Exact solution from Aki and Richards, Quantitative seismology, 2002,pp. 144-145
            !Valid till critical angle of P wave
            Pcritangle = asin(Vp01/Vp02)
            if(Pcritangle .le. Incangle) then
               print *,'Imposed incidence angle higher than P critical angle'
               print *,Vp01,Vp02,Pcritangle/deg2rad,Incangle/deg2rad
               stop
            end if

            !Ray parameter from Snell's law
            p=sin(Incangle)/Vp01;
            !k1 = reflected P wave incidence angle
            k1=asin(p*Vp01);
            !k2 = transmitted P wave incidence angle
            k2=asin(p*Vp02);
            !j1 = reflected S wave incidence angle
            j1=asin(p*Vs01);
            !j2 = transmitted S wave incidence angle
            j2=asin(p*Vs02);
  
            a=ro2*(1d0-2d0*(Vs02*p)**2d0)-ro1*(1d0-2d0*(Vs01*p)**2d0);! = (ro2-ro1) -p^2*d = ro2 - c = b - ro1 
            b=ro2*(1d0-2d0*(Vs02*p)**2d0)+2d0*ro1*(Vs01*p)**2d0; ! = ro2 -p^2*d
            c=ro1*(1d0-2d0*(Vs01*p)**2d0)+2d0*ro2*(Vs02*p)**2d0; ! = ro1 +p^2*d
            d=2d0*(ro2*Vs02**2d0-ro1*Vs01**2);! = 2*(mu2-mu1) = 2*dmu
            E1=b*cos(k1)/Vp01+c*cos(k2)/Vp02;
            F=b*cos(j1)/Vs01+c*cos(j2)/Vs02;
            G1=a-d*cos(k1)/Vp01*cos(j2)/Vs02;
            H=a-d*cos(k2)/Vp02*cos(j1)/Vs01;
            D1=E1*F+G1*H*p**2d0;
            !Reflected P
            Rpp(gi)=((b*cos(k1)/Vp01-c*cos(k2)/Vp02)*F-(a+d*cos(k1)/Vp01*cos(j2)/Vs02)*H*p**2d0)/D1;
            !Reflected S
            Rps(gi)=-2d0*cos(k1)/Vp01*(a*b+c*d*cos(k2)/Vp02*cos(j2)/Vs02)*p*Vp01/Vs01/D1;
            !Transmitted P
            Tpp(gi)=2d0*ro1*cos(k1)/Vp01*F*Vp01/(Vp02*D1);
            !Transmitted S
            Tps(gi)=2d0*ro1*cos(k1)/Vp01*H*p*Vp01/(Vs02*D1)
   
            !!! Compute energy fluxex normalized with respect to
            !!! energy flux of the incident P-wave (eg. 5.42) 
            ERpp(gi)=ro1*Vp01*cos(k1)*Rpp(gi)**2d0/(ro1*Vp01*cos(k1))
            ERps(gi)=ro1*Vs01*cos(j1)*Rps(gi)**2d0/(ro1*Vp01*cos(k1))
            ETpp(gi)=ro2*Vp02*cos(k2)*Tpp(gi)**2d0/(ro1*Vp01*cos(k1))
            ETps(gi)=ro2*Vs02*cos(j2)*Tps(gi)**2d0/(ro1*Vp01*cos(k1))
            !Energy(gi)=(ERpp(gi)+ERps(gi)+ETpp(gi)+ETps(gi)) ! Sum of normalized energies must be 1.0

         END IF

      END IF

   END DO;END DO;END DO;END DO
   !$omp end do
   !$omp end parallel

   !Make the first row = second row
   DO yyy=1,yinyang
   DO i3=1,nx31
   DO i1=1,nx11

      gi = nodenum*(yyy-1) + i1 + (i3-1)*z

      Rpp(gi) = Rpp(gi+y)
      Rps(gi) = Rps(gi+y)
      Tpp(gi) = Tpp(gi+y)
      Tps(gi) = Tps(gi+y)
      ERpp(gi) = ERpp(gi+y)
      ERps(gi) = ERps(gi+y)
      ETpp(gi) = ETpp(gi+y)
      ETps(gi) = ETps(gi+y)

   END DO;END DO;END DO

END IF
!End zoeppritzmod

IF(radialmod .NE. 0) THEN
      
   ALLOCATE(E_chi(nodenum*yyy))

   E_chi = 1d0
   
   !$omp parallel & 
   !$omp shared(Savn,Rhon,fsen,E_chi) &
   !$omp private(yyy,i1,i2,i3,gi,Nparam,Lparam) &    
   !$omp firstprivate(radialmod,ln_fse_min,x,y,z,nx11,nx21,nx31)
   !$omp do schedule(guided,8)
   DO yyy=1,yinyang
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z
      
      IF(Rhon(gi) > 0 .AND. fsen(gi) >= ln_fse_min) THEN     
     
         IF(cartspher == 1) THEN
                Nparam = 1D0/8D0*(Savn(1,1,gi)+Savn(3,3,gi))-0.25D0*Savn(1,3,gi)+0.5D0*Savn(5,5,gi)
                Lparam = 0.5D0*(Savn(4,4,gi)+Savn(6,6,gi))
         ELSE
                Nparam = 1D0/8D0*(Savn(1,1,gi)+Savn(2,2,gi))-0.25D0*Savn(1,2,gi)+0.5D0*Savn(6,6,gi)
                Lparam = 0.5D0*(Savn(4,4,gi)+Savn(5,5,gi))
         END IF
         
         IF(radialmod > 1) THEN

            E_chi(gi) = 100d0*((Nparam/Lparam)**0.5d0 - 1d0)

         ELSE

            E_chi(gi) = Nparam/Lparam

            IF(E_chi(gi) < 0) THEN
               print *,'Negative radial anisotropy : ',i1,i2,i3,E_chi(gi),Nparam,Lparam,Savn(1,1,gi),Savn(3,3,gi),Savn(1,3,gi),Savn(4,4,gi),Savn(5,5,gi),Savn(6,6,gi)
               stop
            END IF

         END IF

      END IF
   
   END DO;END DO;END DO;END DO
   !$omp end do
   !$omp end parallel

END IF
!End radialmod

IF(azimod .NE. 0) THEN
      
   !!! 2D
   IF(dimensions == 2) THEN

   ALLOCATE(G(nodenum))

   G = 0d0

   !$omp parallel & 
   !$omp shared(Savn,Rhon,fsen,G) &
   !$omp private(i1,i2,i3,gi,Gc,Gs) &    
   !$omp firstprivate(azimod,ln_fse_min,x,y,z,nx11,nx21,nx31)
   !$omp do schedule(guided,8)
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      gi = i1 + (i2-1)*y + (i3-1)*z
      
      IF(Rhon(gi) > 0 .AND. fsen(gi) >= ln_fse_min) THEN     
     
         IF(azimod == 1) THEN

            IF(cartspher == 1) THEN
               Gc = 0.5D0*(Savn(6,6,gi)-Savn(4,4,gi))
               Gs = Savn(4,6,gi)
            ELSE 
               Gc = 0.5D0*(Savn(5,5,gi)-Savn(4,4,gi))
               Gs = Savn(4,5,gi)
            END IF

            G(gi) = (Gc**2 + Gs**2)**0.5d0

         ELSE

            IF(cartspher == 1) G(gi) = 100d0*((Savn(6,6,gi)/Savn(4,4,gi))**0.5 - 1d0)
            IF(cartspher == 2) G(gi) = 100d0*((Savn(5,5,gi)/Savn(4,4,gi))**0.5 - 1d0)

         END IF

      END IF
   
   END DO;END DO;END DO
   !$omp end do
   !$omp end parallel

   !!! 3D
   ELSE
   
   azinum = 0
   DO yyy=1,yinyang
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z

      IF(yyy == 2) CALL yin2yang(X1(i1),X3(i3),lr,cr)

      IF(yyy == 2 .AND. lr>X1(1) .AND. lr<X1(nx11) .AND. cr>X3(1) .AND. cr<X3(nx31)) THEN
      ELSE
      IF(Rhon(gi) > 0 .AND. fsen(gi) >= ln_fse_min) THEN
         !Remove nodes of Yang grid falling in Yin grid domain

         azinum = azinum + 1

      END IF; END IF

   END DO;END DO;END DO;END DO

   !Add scale
   azinum = azinum + 1

   ALLOCATE(azimuthal(3,azinum),idx(azinum))

   azimuthal = 0d0 ; idx = 0

   azinum = 0
   DO yyy=1,yinyang
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z

      IF(yyy == 2) CALL yin2yang(X1(i1),X3(i3),lr,cr)

      IF(yyy == 2 .AND. lr>X1(1) .AND. lr<X1(nx11) .AND. cr>X3(1) .AND. cr<X3(nx31)) THEN
      ELSE
      IF(Rhon(gi) > 0 .AND. fsen(gi) >= ln_fse_min) THEN

         IF(cartspher == 1) THEN
            Gc = 0.5D0*(Savn(6,6,gi)-Savn(4,4,gi))
            Gs = Savn(4,6,gi)
         ELSE
            Gc = 0.5D0*(Savn(5,5,gi)-Savn(4,4,gi))
            Gs = Savn(4,5,gi)
         END IF

         G3D  = (Gc**2 + Gs**2)**0.5d0

         azinum = azinum + 1
         idx(azinum) = gi
         phi =  atan2(Gs,Gc)*0.5D0
         !Set angle from 0 to 180 degrees
         if(phi > pi) phi = phi - pi
         if(phi < 0) phi = phi + pi
         !Find horizontal length of the vector
         IF(cartspher == 1) THEN
            azimuthal(1,azinum) = cos(phi)*G3D
            azimuthal(2,azinum) = 0d0
            azimuthal(3,azinum) = sin(phi)*G3D
         ELSE
            azimuthal(1,azinum) = cos(phi)*G3D
            azimuthal(2,azinum) = sin(phi)*G3D
            azimuthal(3,azinum) = 0d0

            !Rotate toward the Z axis in spherical coordinates
            !Rotate azimuth according to long and colat of the station
            !This rotation is inverse to the previous one
            phi1 = X1(i1) - pi*1.5d0; theta = X3(i3) ; phi2 = 0d0
            IF(yyy == 2) THEN
               theta = cr ; phi1 = lr - pi*1.5d0
            END IF

            CALL rotmatrixZXZ(phi1,theta,phi2,acs)

            azi(:) = azimuthal(:,azinum)
            azimuthal(:,azinum) = 0d0
            DO j=1,3; DO k=1,3
               azimuthal(j,azinum) = azimuthal(j,azinum) + acs(j,k)*azi(k)
            END DO; END DO

         END IF

      END IF; END IF

   END DO;END DO;END DO;END DO

   !Add scale
   azinum = azinum + 1

   END IF
   
END IF
!End azimod

!!! Write infos in hdf5 format
IF(vpvsmod > 0 .OR. zoeppritzmod > 0 .OR. radialmod > 0 .OR. azimod > 0) THEN

   ! Create a new file using default properties.
   filename = str2//'eulerian'//dt_str4//'.h5'
  
   !Initialize FORTRAN interface.
   CALL H5open_f (error)

   CALL H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

   !Make grid
   !3D grid 
   IF(dimensions == 3 .OR. (replicateZmod > 0 .AND. nx31 > 1)) THEN

   ALLOCATE(vertices(3,nodenum),connectivity(8,cellnum))

   !$omp parallel & 
   !$omp shared(vertices,X1,X2,X3) &
   !$omp private(i1,i2,i3,gi) &    
   !$omp firstprivate(cartspher,x,y,z,nx11,nx21,nx31)
   !$omp do schedule(guided,8)
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      gi = i1 + (i2-1)*y + (i3-1)*z

      IF(cartspher == 1) THEN
         vertices(1,gi)=X1(i1);
         vertices(2,gi)=X2(i2);
         vertices(3,gi)=X3(i3);
      ELSE IF(cartspher == 2) THEN
         vertices(1,gi)=X2(i2)*sin(X3(i3))*cos(X1(i1))
         vertices(2,gi)=X2(i2)*sin(X3(i3))*sin(X1(i1))
         vertices(3,gi)=X2(i2)*cos(X3(i3))
      END IF

   END DO;END DO;END DO
   !$omp end do
   !$omp end parallel

   !Cells' connectivity     
   !$omp parallel & 
   !$omp shared(connectivity) &
   !$omp private(i1,i2,i3,gi) &    
   !$omp firstprivate(x,y,z,yc,zc,nx11,nx21,nx31)
   !$omp do schedule(guided,8)
   DO i3=1,nx31-1
   DO i2=1,nx21-1
   DO i1=1,nx11-1

      gi = i1 + (i2-1)*yc + (i3-1)*zc
      connectivity(1,gi)=(i1-1)+(i2-1)*y+(i3-1)*z;
      connectivity(2,gi)=connectivity(1,gi)+1;
      connectivity(3,gi)=connectivity(1,gi)+y+1;
      connectivity(4,gi)=connectivity(1,gi)+y;
      connectivity(5,gi)=connectivity(1,gi)+z;
      connectivity(6,gi)=connectivity(1,gi)+z+1;
      connectivity(7,gi)=connectivity(1,gi)+z+y+1;
      connectivity(8,gi)=connectivity(1,gi)+z+y;

   END DO;END DO;END DO
   !$omp end do
   !$omp end parallel

   dum_int(1)=3
   dum_int(2)=nodenum
   CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,vertices,'vertices',1)
   dum_int(1)=8
   dum_int(2)=cellnum
   CALL loadsave_integer(0,2,file_id,dum_int(1:2),H5T_NATIVE_INTEGER,connectivity,'connectivity',1)

   !2D grid
   ELSE

   ALLOCATE(vertices(2,nodenum),connectivity(4,cellnum))
   DO i2=1,nx21
   DO i1=1,nx11

      gi = i1 + (i2-1)*y
      
      IF(cartspher == 1) THEN
         vertices(1,gi)=X1(i1);
         vertices(2,gi)=X2(i2);
      ELSE IF(cartspher == 2) THEN
         xx = cos(X1(i1))
         yy = sin(X1(i1))
         vertices(1,gi)=X2(i2)*xx;
         vertices(2,gi)=X2(i2)*yy;
      END IF

   END DO;END DO
   
   !Cells' connectivity     
   DO i2=1,nx21-1
   DO i1=1,nx11-1

      gi = i1 + (i2-1)*yc
      connectivity(1,gi)=(i1-1)+(i2-1)*nx11
      connectivity(2,gi)=connectivity(1,gi)+1;
      connectivity(3,gi)=connectivity(1,gi)+nx11+1;
      connectivity(4,gi)=connectivity(1,gi)+nx11;
   
   END DO;END DO

   dum_int(1)=2
   dum_int(2)=nodenum
   CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,vertices,'vertices',1)
   dum_int(1)=4
   dum_int(2)=cellnum
   CALL loadsave_integer(0,2,file_id,dum_int(1:2),H5T_NATIVE_INTEGER,connectivity,'connectivity',1)

   END IF
   !End make grid
        
   IF(vpvsmod > 0) THEN

      print *
      IF(vpvsmod == 1) print *,'Save isotropic seismic velocities and anomalies'
      IF(vpvsmod == 2) THEN
         write(*,'(a)') ' Save seismic velocities and anomalies along direction with angular distances from:'
         write(*,'(a,1f6.2,a,1f6.2,a,1f6.2,a)') ' axis 1 =',cosx1,' (deg) , axis 2 = ',cosx2,' (deg), axis 3 = ',cosx3,' (deg)'
      END IF
      print *

      CALL loadsave_double(0,1,file_id,nx11,H5T_NATIVE_DOUBLE, X1,'X1',1) ! x axis
      CALL loadsave_double(0,1,file_id,nx21,H5T_NATIVE_DOUBLE, X2,'X2',1) ! y axis
      CALL loadsave_double(0,1,file_id,nx31,H5T_NATIVE_DOUBLE, X3,'X3',1) ! z axis
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Rhon(1:nodenum),'Rho',1)
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Vp(1:nodenum),'Vp',1)
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Vs(1:nodenum),'Vs',1)
      IF(yinyang /= 2) DEALLOCATE(Vp,Vs)
      IF(dvpvsmod > 0) THEN
         CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,dVp(1:nodenum),'dVp',1)
         CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,dVs(1:nodenum),'dVs',1)
         IF(yinyang /= 2) DEALLOCATE(dVp,dVs)
      END IF

      !dum_int(1)=6
      !dum_int(2)=6
      !dum_int(3)=nodenum
      !CALL loadsave_double(0,3,file_id,dum_int(1:3),H5T_NATIVE_DOUBLE,Savn,'Sav',1)

   END IF

   IF(zoeppritzmod > 0) THEN
      write(*,"(a,f4.1,a)"),' Save reflected/transmitted P-S energy with incidence angle = ',Incangle/deg2rad, '(deg)'
      print *,''
      !CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Rpp,'Rpp',1)
      !CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Rpp,'Rps',1)
      !CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Rpp,'Tpp',1)
      !CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Rpp,'Tps',1)
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,ERpp(1:nodenum),'ERpp',1)
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,ERps(1:nodenum),'ERps',1)
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,ETpp(1:nodenum),'ETpp',1)
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,ETps(1:nodenum),'ETps',1)
      IF(yinyang /= 2) DEALLOCATE(Rpp,Rps,Tpp,Tps,ERpp,ERps,ETpp,ETps)
   END IF

   IF(radialmod > 0) THEN
      print *,'Save radial anisotropy'
      print *,''
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,E_chi(1:nodenum),'RadialAnisotropy',1)
      IF(yinyang /= 2) DEALLOCATE(E_chi)
   END IF

   IF(azimod > 0 .AND. dimensions == 2) THEN

      print *,'Save azimuthal anisotropy'
      print *,''
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,G,'AzimuthalAnisotropy',1)
      DEALLOCATE(G)

   END IF

   IF(yinyang /= 2) DEALLOCATE(vertices,connectivity)

   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   CALL H5close_f(error)
  
END IF
!END  Write infos in hdf5 format

!Make xmf file for visualization in Paraview
IF(vpvsmod > 0 .OR. zoeppritzmod > 0 .OR. radialmod > 0 .OR. (dimensions == 2 .AND. azimod > 0)) THEN
   
   filenamexmf = str2//'XDMF.eulerian'//dt_str4//'.xmf'

   OPEN(19,file=filenamexmf,status='replace')
   WRITE(19,'(a)') '<?xml version="1.0" ?>'
   WRITE(19,'(a)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') ' <Domain>'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '   <Grid Name="Eulerian Grid" >'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a,1PE12.4,a)') '      <Time Value="',timesum,'" />'
   WRITE(19,'(a)') ' '

   !3D grid 
   IF(dimensions == 3 .OR. (dimensions == 2 .AND. replicateZmod > 0 .AND. nx31 > 1)) THEN

   WRITE(19,'(a,i0,a)') '         <Topology Type="Hexahedron" NumberOfElements="',cellnum,'">'
   WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" DataType="Int"  Dimensions="',cellnum,' 8">eulerian',dt_str4,'.h5:/connectivity</DataItem>'
   WRITE(19,'(a)') '         </Topology>'

   WRITE(19,'(a)') '         <Geometry Type="XYZ">'
   WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 3">eulerian',dt_str4,'.h5:/vertices </DataItem>'
   WRITE(19,'(a)') '         </Geometry>'
   WRITE(19,'(a)') ' '

   !2D grid
   ELSE
   
   WRITE(19,'(a,i0,a)') '         <Topology Type="Quadrilateral" NumberOfElements="',cellnum,'">' 
   WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" DataType="Int"  Dimensions="',cellnum,' 4">eulerian',dt_str4,'.h5:/connectivity</DataItem>'
   WRITE(19,'(a)') '         </Topology>'

   WRITE(19,'(a)') '         <Geometry Type="XY">'
   WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 2">eulerian',dt_str4,'.h5:/vertices </DataItem>'
   WRITE(19,'(a)') '         </Geometry>'
   WRITE(19,'(a)') ' '
   
   END IF

   IF(vpvsmod > 0) THEN
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="Rho">'
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerian',dt_str4,'.h5:/Rho</DataItem>'
      WRITE(19,'(a)') '            </DataItem>'
      WRITE(19,'(a)') '         </Attribute>'
      WRITE(19,'(a)') ' '
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="Vp">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerian',dt_str4,'.h5:/Vp</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="Vs">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerian',dt_str4,'.h5:/Vs</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
     
      IF(dvpvsmod > 0 ) THEN

      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="dVp">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerian',dt_str4,'.h5:/dVp</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="dVs">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerian',dt_str4,'.h5:/dVs</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '

      END IF

   END IF
   
   IF(radialmod > 0) THEN
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="Radial Anisotropy">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerian',dt_str4,'.h5:/RadialAnisotropy</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
   END IF
   
   IF(dimensions == 2 .AND. azimod > 0) THEN
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="Azimuthal Anisotropy">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerian',dt_str4,'.h5:/AzimuthalAnisotropy</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
   END IF
   
   IF(zoeppritzmod > 0) THEN
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="ERpp">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerian',dt_str4,'.h5:/ERpp</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="ERps">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerian',dt_str4,'.h5:/ERps</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="ETpp">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerian',dt_str4,'.h5:/ETpp</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="ETps">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerian',dt_str4,'.h5:/ETps</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
   END IF

   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '   </Grid>'
   WRITE(19,'(a)') ' '

   IF(yinyang == 2) THEN

   WRITE(19,'(a)') '   <Grid Name="Eulerian Grid" >'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a,1PE12.4,a)') '      <Time Value="',timesum,'" />'
   WRITE(19,'(a)') ' '

   WRITE(19,'(a,i0,a)') '         <Topology Type="Hexahedron" NumberOfElements="',cellnum,'">'
   WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" DataType="Int"  Dimensions="',cellnum,' 8">eulerianY',dt_str4,'.h5:/connectivity</DataItem>'
   WRITE(19,'(a)') '         </Topology>'

   WRITE(19,'(a)') '         <Geometry Type="XYZ">'
   WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 3">eulerianY',dt_str4,'.h5:/vertices </DataItem>'
   WRITE(19,'(a)') '         </Geometry>'
   WRITE(19,'(a)') ' '

   IF(vpvsmod > 0) THEN
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="Rho">'
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerianY',dt_str4,'.h5:/Rho</DataItem>'
      WRITE(19,'(a)') '            </DataItem>'
      WRITE(19,'(a)') '         </Attribute>'
      WRITE(19,'(a)') ' '
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="Vp">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerianY',dt_str4,'.h5:/Vp</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="Vs">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerianY',dt_str4,'.h5:/Vs</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
     
      IF(dvpvsmod > 0 ) THEN

      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="dVp">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerianY',dt_str4,'.h5:/dVp</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="dVs">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerianY',dt_str4,'.h5:/dVs</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '

      END IF

   END IF
   
   IF(radialmod > 0) THEN
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="Radial Anisotropy">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerianY',dt_str4,'.h5:/RadialAnisotropy</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
   END IF
   
   IF(zoeppritzmod > 0) THEN
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="ERpp">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerianY',dt_str4,'.h5:/ERpp</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="ERps">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerianY',dt_str4,'.h5:/ERps</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="ETpp">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerianY',dt_str4,'.h5:/ETpp</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="ETps">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerianY',dt_str4,'.h5:/ETps</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
   END IF

   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '   </Grid>'
   WRITE(19,'(a)') ' '

   END IF

   WRITE(19,'(a)') ' </Domain>'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '</Xdmf>'

   CLOSE(19)

END IF
!END Make xmf file for visualization in Paraview


!Save data of YANG grid    
IF(yinyang == 2) THEN

IF(vpvsmod > 0 .OR. zoeppritzmod > 0 .OR. radialmod > 0 .OR. azimod > 0) THEN

   ! Create a new file using default properties.
   filename = str2//'eulerianY'//dt_str4//'.h5'
  
   !Initialize FORTRAN interface.
   CALL H5open_f (error)

   CALL H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

   !Make grid
   !3D grid 

   !$omp parallel & 
   !$omp shared(vertices) &
   !$omp private(i1,yy) &    
   !$omp firstprivate(nodenum)
   !$omp do schedule(guided,8)
   DO i1=1,nodenum

      vertices(1,i1) = -vertices(1,i1)
      yy = vertices(2,i1)
      vertices(2,i1) = vertices(3,i1)
      vertices(3,i1) = yy

   END DO
   !$omp end do
   !$omp end parallel

   !Cells' connectivity     
   !$omp parallel & 
   !$omp shared(connectivity) &
   !$omp private(i1,i2,i3,gi) &    
   !$omp firstprivate(x,y,z,yc,zc,nx11,nx21,nx31)
   !$omp do schedule(guided,8)
   DO i3=1,nx31-1
   DO i2=1,nx21-1
   DO i1=1,nx11-1

      gi = i1 + (i2-1)*yc + (i3-1)*zc
      connectivity(1,gi)=(i1-1)+(i2-1)*y+(i3-1)*z;
      connectivity(2,gi)=connectivity(1,gi)+1;
      connectivity(3,gi)=connectivity(1,gi)+y+1;
      connectivity(4,gi)=connectivity(1,gi)+y;
      connectivity(5,gi)=connectivity(1,gi)+z;
      connectivity(6,gi)=connectivity(1,gi)+z+1;
      connectivity(7,gi)=connectivity(1,gi)+z+y+1;
      connectivity(8,gi)=connectivity(1,gi)+z+y;

   END DO;END DO;END DO
   !$omp end do
   !$omp end parallel

   dum_int(1)=3
   dum_int(2)=nodenum
   CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,vertices,'vertices',1)
   dum_int(1)=8
   dum_int(2)=cellnum
   CALL loadsave_integer(0,2,file_id,dum_int(1:2),H5T_NATIVE_INTEGER,connectivity,'connectivity',1)

   !End make grid
        
   i1 = nodenum + 1
   i2 = nodenum*yinyang

   IF(vpvsmod > 0) THEN

      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Rhon(i1:i2),'Rho',1)
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Vp(i1:i2),'Vp',1)
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Vs(i1:i2),'Vs',1)
      DEALLOCATE(Vp,Vs)
      IF(dvpvsmod > 0) THEN
         CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,dVp(i1:i2),'dVp',1)
         CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,dVs(i1:i2),'dVs',1)
         DEALLOCATE(dVp,dVs)
      END IF

      !dum_int(1)=6
      !dum_int(2)=6
      !dum_int(3)=nodenum
      !CALL loadsave_double(0,3,file_id,dum_int(1:3),H5T_NATIVE_DOUBLE,Savn,'Sav',1)

   END IF

   IF(zoeppritzmod > 0) THEN
      !CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Rpp(i1:i2),'Rpp',1)
      !CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Rpp(i1:i2),'Rps',1)
      !CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Rpp(i1:i2),'Tpp',1)
      !CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,Rpp(i1:i2),'Tps',1)
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,ERpp(i1:i2),'ERpp',1)
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,ERps(i1:i2),'ERps',1)
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,ETpp(i1:i2),'ETpp',1)
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,ETps(i1:i2),'ETps',1)
      DEALLOCATE(Rpp,Rps,Tpp,Tps,ERpp,ERps,ETpp,ETps)
   END IF

   IF(radialmod > 0) THEN
      CALL loadsave_double(0,1,file_id,nodenum,H5T_NATIVE_DOUBLE,E_chi(i1:i2),'RadialAnisotropy',1)
      DEALLOCATE(E_chi)
   END IF

   DEALLOCATE(vertices,connectivity)

   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   CALL H5close_f(error)
  
END IF

END IF

!End save data of YANG grid

IF(dimensions == 3 .AND. azimod > 0) THEN

   IF(azinum == 0) THEN

      print *,'No azimuthal anisotropy for 3D model'
      GOTO 10

   END IF


   print *,'Save azimuthal anisotropy'
   print *,''
   filename = str2//'azianis'//dt_str4//'.h5'
   !Initialize FORTRAN interface.
   CALL H5open_f (error)

   ! Create a new file using default properties.
   CALL H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

   !Save grid
   ALLOCATE(xyz(3,azinum))
   k = 1
   DO yyy=1,yinyang
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z

      IF(gi == idx(k)) THEN

         IF(cartspher == 1) THEN
            xyz(1,k) = X1(i1)
            xyz(2,k) = X2(i2)
            xyz(3,k) = X3(i3)
         ELSE IF(cartspher == 2) THEN
            IF(yyy == 1) THEN
               xyz(1,k)=X2(i2)*sin(X3(i3))*cos(X1(i1))
               xyz(2,k)=X2(i2)*sin(X3(i3))*sin(X1(i1))
               xyz(3,k)=X2(i2)*cos(X3(i3))
            ELSE
               CALL yin2yang(X1(i1),X3(i3),lr,cr)
               xyz(1,k)=X2(i2)*sin(cr)*cos(lr)
               xyz(2,k)=X2(i2)*sin(cr)*sin(lr)
               xyz(3,k)=X2(i2)*cos(cr)
            END IF
         END IF

         k = k + 1

      END IF
   END DO; END DO; END DO; END DO

   !Set scale
   G3D = 1.0
   phi = 0.0 !parallel to x1 axis
   azimuthal(1,azinum) = cos(phi)*G3D
   azimuthal(2,azinum) = sin(phi)*G3D
   azimuthal(3,azinum) = 0d0

   !Rotate toward the Z axis in spherical coordinates
   IF(cartspher == 2) THEN

      aziscalex1 = aziscalex1*deg2rad
      aziscalex3 = aziscalex3*deg2rad

      !Rotate azimuth according to long and colat of the station
      !This rotation is inverse to the previous one
      phi1 = aziscalex1 - pi*1.5d0; theta = aziscalex3; phi2 = 0d0
      CALL rotmatrixZXZ(phi1,theta,phi2,acs)

      azi(:) = azimuthal(:,azinum)
      azimuthal(:,azinum) = 0d0
      DO j=1,3; DO k=1,3
         azimuthal(j,azinum) = azimuthal(j,azinum) + acs(j,k)*azi(k)
      END DO; END DO

   END IF

   k = azinum
   IF(cartspher == 1) THEN
      xyz(1,k) = aziscalex1
      xyz(2,k) = aziscalex2
      xyz(3,k) = aziscalex3
   ELSE IF(cartspher == 2) THEN
      xyz(1,k)=aziscalex2*sin(aziscalex3)*cos(aziscalex1)
      xyz(2,k)=aziscalex2*sin(aziscalex3)*sin(aziscalex1)
      xyz(3,k)=aziscalex2*cos(aziscalex3)
   END IF
   dum_int(1)=3
   dum_int(2)=azinum
   CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,xyz,'xyz',1)
   CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,azimuthal,'Azimuthal Anisotropy',1)
   DEALLOCATE(azimuthal,xyz,idx)

   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   CALL H5close_f(error)
  
   filenamexmf = str2//'XDMF.azianis'//dt_str4//'.xmf'

   OPEN(19,file=filenamexmf,status='replace')
   WRITE(19,'(a)') '<?xml version="1.0" ?>'
   WRITE(19,'(a)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') ' <Domain>'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '   <Grid Name="materialSwarm" >'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a,1f10.8,a)') '      <Time Value="',timesum,'" />'
   WRITE(19,*) ' '
   WRITE(19,'(a,i0,a)') '          <Topology Type="POLYVERTEX" NodesPerElement="',azinum,'"> </Topology>'
   WRITE(19,'(a)') '          <Geometry Type="XYZ">'
   WRITE(19,'(a,i0,a,a,a)') '             <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',azinum,' 3">azianis',dt_str4,'.h5:/xyz </DataItem>'
   WRITE(19,'(a)') '          </Geometry>'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '          <Attribute Type="Vector" Center="Node" Name="Azimuthal Anisotropy">'
   WRITE(19,'(a,i0,a,a,a)') '             <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',azinum,' 3">azianis',dt_str4,'.h5:/Azimuthal Anisotropy</DataItem>'
   WRITE(19,'(a)') '          </Attribute>'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '   </Grid>'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') ' </Domain>'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '</Xdmf>'

   CLOSE(19)
 
10 END IF

END SUBROUTINE eulerian_output 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine UPPERLEFT2D, calculation of upper left node                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE upperleft2D(x10,x20,i1,i2)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: x10,x20

   INTEGER :: i1,i2

   i1 = 1; i2 = 1

   DO WHILE (X1(i1) .GT. x10 .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
   DO WHILE (X1(i1+1) .LT. x10 .AND. i1 .LT. nx11-1) ; i1 = i1+1 ; ENDDO

   DO WHILE (X2(i2) .GT. x20 .AND. i2 .GT. 1) ; i2 = i2-1 ; ENDDO
   DO WHILE (X2(i2+1) .LT. x20 .AND. i2 .LT. nx21-1) ; i2 = i2+1 ; ENDDO

   RETURN

   END SUBROUTINE upperleft2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine UPPERLEFT3D, calculation of upper left node                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE upperleft3D(x10,x20,x30,i1,i2,i3)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: x10,x20,x30

   INTEGER :: i1,i2,i3

   i1 = 1; i2 = 1; i3 = 1

   DO WHILE (X1(i1) .GT. x10 .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
   DO WHILE (X1(i1+1) .LT. x10 .AND. i1 .LT. nx11-1) ; i1 = i1+1 ; ENDDO

   DO WHILE (X2(i2) .GT. x20 .AND. i2 .GT. 1) ; i2 = i2-1 ; ENDDO
   DO WHILE (X2(i2+1) .LT. x20 .AND. i2 .LT. nx21-1) ; i2 = i2+1 ; ENDDO

   DO WHILE (X3(i3) .GT. x30 .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
   DO WHILE (X3(i3+1) .LT. x30 .AND. i3 .LT. nx31-1) ; i3 = i3+1 ; ENDDO

   RETURN

   END SUBROUTINE upperleft3D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine interp, bilinear interpolation of a function                !!!
!!! following the method given in Numerical Recipies, Press et al., p96    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE mark2node

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: i1,i2,i3,i,j,k,m,gi,gijk,x,y,z,tid,nrot,yyy,yyn,yydum
   INTEGER, DIMENSION(26) :: surrnodes
   ! indices of the UP-LEFT grid point closest to the considered point

   DOUBLE PRECISION :: xx,yy,zz,w,w1,w2,w3,w4,lr,rr,cr,lr1,cr1
   DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: wt
   ! finite strain = ln(longaxis/shortaxis)

!	                   
!	                   
!	                   
!	          3----4   
!       +Y	  |    |   
!	|         |    |  
!	|         |    | 
!	|         1----2
!	|  
!	----- +X

!Interpolate Sav to nodes
yyy = yinyang
ALLOCATE(wt(nodenum*yyy))
Savn = 0 ; Rhon = 0; fsen = 0; wt = 0

!!! 2D
IF(dimensions == 2) THEN

x=1
y=nx11
DO m=1,marknum0

      IF(Sav(1,1,m) > 0d0) THEN
      IF(mx1(m) >= n1first .AND. mx1(m) <= n1last) THEN
      IF(mx2(m) >= n2first .AND. mx2(m) <= n2last) THEN

      !Find nearest upper left node
      CALL upperleft2D(mx1(m),mx2(m),i1,i2)

      xx = (mx1(m)-X1(i1))/(X1(i1+1)-X1(i1))
      IF(mx1(m) >= X1(nx11)) xx = 1d0
      IF(xx > 1d0 .or. xx < 0d0) THEN
      print *,'X',xx,mx1(m),X1(i1),X1(i1+1)-X1(i1)
      stop
      END IF

      yy = (mx2(m)-X2(i2))/(X2(i2+1)-X2(i2))
      IF(mx2(m) >= X2(nx21)) yy = 1d0
      IF(yy > 1d0 .or. yy < 0d0) THEN
      print *,'Y',yy,mx2(m),X2(i2),X2(i2+1)-X2(i2)
      stop
      END IF

      gi = i1 + (i2-1)*y

      !Node 1
      w = (1d0-xx)*(1d0-yy)
      Savn(:,:,gi      ) = Savn(:,:,gi      ) + Sav(:,:,m)*w;
      Rhon(gi      ) = Rhon(gi      ) + rho(m)*w;
      fsen(gi      ) = fsen(gi      ) + ln_fse(m)*w;
      wt(gi) = wt(gi) + w
      
      !Node 2
      w = xx*(1d0-yy)
      Savn(:,:,gi+x    ) = Savn(:,:,gi+x    ) + Sav(:,:,m)*w;
      Rhon(gi+x    ) = Rhon(gi+x    ) + rho(m)*w;
      fsen(gi+x    ) = fsen(gi+x    ) + ln_fse(m)*w;
      wt(gi+x) = wt(gi+x) + w

      !Node 3
      w = (1d0-xx)*yy
      Savn(:,:,gi+y    ) = Savn(:,:,gi+y    ) + Sav(:,:,m)*w;
      Rhon(gi+y    ) = Rhon(gi+y    ) + rho(m)*w;
      fsen(gi+y    ) = fsen(gi+y    ) + ln_fse(m)*w;
      wt(gi+y) = wt(gi+y) + w

      !Node 4
      w = xx*yy
      Savn(:,:,gi+x+y  ) = Savn(:,:,gi+x+y  ) + Sav(:,:,m)*w;
      Rhon(gi+x+y  ) = Rhon(gi+x+y  ) + rho(m)*w;
      fsen(gi+x+y  ) = fsen(gi+x+y  ) + ln_fse(m)*w;
      wt(gi+x+y) = wt(gi+x+y) + w

      IF(wt(gi) < 0d0) THEN  
      print *,gi,i1,i2,wt(gi)!,Sav(:,:,m)
      END IF

      END IF;END IF;END IF

END DO

!Normalization
gi=0
DO i=1,nodenum
  IF(wt(i) > 0d0) THEN
     Savn(:,:,i)=Savn(:,:,i)/wt(i)
     Rhon(i)=Rhon(i)/wt(i)
     fsen(i)=fsen(i)/wt(i)
     gi=gi+1
  ELSE
     print *,'Empty node:',i
  END IF
END DO

IF(0 == 1) THEN
!For empty nodes, interpolate Savn, Rhon from surrounding nodes
print *
print *,'Interpolation for empty nodes is active !!!'

DO i2=1,nx21
DO i1=1,nx11
     
  gi = i1 + (i2-1)*y
  IF(wt(gi) .eq. 0d0) THEN
     !Find surrounding nodes
     surrnodes = 0d0
     !Node 1
     IF(i1 < nx11 .AND. wt(gi+x) > 0) surrnodes(1) = gi + x
     !Node 2
     IF(i1 > 1 .AND. wt(gi-x) > 0) surrnodes(2) = gi - x
     !Node 3
     IF(i2 < nx21 .AND. wt(gi+y) > 0) surrnodes(3) = gi + y
     !Node 4
     IF(i2 > 1 .AND. wt(gi-y) > 0) surrnodes(4) = gi - y
     !Node 5
     IF(i1 < nx11 .AND. i2 > 1 .AND. wt(gi+x-y) > 0) surrnodes(7) = gi + x - y
     !Node 6
     IF(i1 < nx11 .AND. i2 < nx21 .AND. wt(gi+x+y) > 0) surrnodes(8) = gi + x + y
     !Node 7 
     IF(i1 > 1 .AND. i2 > 1 .AND. wt(gi-x-y) > 0) surrnodes(9) = gi - x - y
     !Node 8 
     IF(i1 > 1 .AND. i2 < nx21 .AND. wt(gi-x+y) > 0) surrnodes(10) = gi - x + y
        
     i = 0
     DO j = 1,8  
        IF(surrnodes(j) > 0) THEN
        IF(wt(surrnodes(j)) > 0) THEN
           i = i + 1
           Savn(:,:,gi)=Savn(:,:,gi) + Savn(:,:,surrnodes(j))
           Rhon(gi)=Rhon(gi) + Rhon(surrnodes(j))
           fsen(gi)=fsen(gi) + fsen(surrnodes(j))
        END IF 
        END IF 
    END DO

    IF(i > 0) THEN
       Savn(:,:,gi) = Savn(:,:,gi)/i
       Rhon(gi) = Rhon(gi)/i
       fsen(gi) = fsen(gi)/i
    END IF
    
  END IF

END DO
END DO
    
END IF

END IF

!!! 3D
IF(dimensions == 3) THEN

x=1
y=nx11
z=nx11*nx21

DO m=1,marknum

      lr = mx1(m) ; rr = mx2(m) ; cr = mx3(m)

      IF(Sav(1,1,m) > 0d0 .AND. Sav(1,1,m) < 1d4) THEN
      IF(lr >= n1first .AND. lr <= n1last) THEN
      IF(rr >= n2first .AND. rr <= n2last) THEN
      IF(cr >= n3first .AND. cr <= n3last) THEN

      yydum = 1;
      yyy = mYY(m)
      
      DO yyn = 1,yinyang

      !Check if possible to interpolate to YANG grid
      IF(yyn == 2) THEN
         yydum = 0
         CALL yin2yang(lr,cr,lr1,cr1)
         IF(lr1 >= n1first .AND. lr1 <= n1last) THEN
         IF(cr1 >= n3first .AND. cr1 <= n3last) THEN
            yydum = 1
            lr = lr1
            cr = cr1
            IF(yyy == 1) THEN
               yyy = 2
            ELSE
               yyy = 1
            END IF
         END IF;END IF
      END IF
      
      IF(yydum > 0) THEN

      !Find nearest upper left node
      CALL upperleft3D(lr,rr,cr,i1,i2,i3)

      xx = (lr-X1(i1))/(X1(i1+1)-X1(i1))
      if(lr >= X1(nx11)) xx = 1d0
      if(xx > 1d0 .or. xx < 0d0) then
      print *,'X',xx,lr,X1(i1),X1(i1+1)-X1(i1)
      stop
      end if

      yy = (rr-X2(i2))/(X2(i2+1)-X2(i2))
      if(rr >= X2(nx21)) yy = 1d0
      if(yy > 1d0 .or. yy < 0d0) then
      print *,'Y',yy,rr,X2(i2),X2(i2+1)-X2(i2)
      stop
      end if

      zz = (cr-X3(i3))/(X3(i3+1)-X3(i3))
      if(cr >= X3(nx31)) zz = 1d0
      if(zz > 1d0 .or. zz < 0d0) then
      print *,'Z',zz,cr,X3(i3),X3(i3+1)-X3(i3)
      stop
      end if
      
      gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z

      !Node 1
      w = (1d0-xx)*(1d0-yy)*(1d0-zz)
      Savn(:,:,gi      ) = Savn(:,:,gi      ) + Sav(:,:,m)*w;
      Rhon(gi      ) = Rhon(gi      ) + rho(m)*w;
      fsen(gi      ) = fsen(gi      ) + ln_fse(m)*w;
      wt(gi) = wt(gi) + w
      
      !Node 2
      w = xx*(1d0-yy)*(1d0-zz)
      Savn(:,:,gi+x    ) = Savn(:,:,gi+x    ) + Sav(:,:,m)*w;
      Rhon(gi+x    ) = Rhon(gi+x    ) + rho(m)*w;
      fsen(gi+x    ) = fsen(gi+x    ) + ln_fse(m)*w;
      wt(gi+x) = wt(gi+x) + w

      !Node 3
      w = xx*(1d0-yy)*zz
      Savn(:,:,gi+x+z  ) = Savn(:,:,gi+x+z  ) + Sav(:,:,m)*w;
      Rhon(gi+x+z  ) = Rhon(gi+x+z  ) + rho(m)*w;
      fsen(gi+x+z  ) = fsen(gi+x+z  ) + ln_fse(m)*w;
      wt(gi+x+z) = wt(gi+x+z) + w

      !Node 4
      w = (1d0-xx)*(1d0-yy)*zz
      Savn(:,:,gi+z    ) = Savn(:,:,gi+z    ) + Sav(:,:,m)*w;
      Rhon(gi+z    ) = Rhon(gi+z    ) + rho(m)*w;
      fsen(gi+z    ) = fsen(gi+z    ) + ln_fse(m)*w;
      wt(gi+z) = wt(gi+z) + w

      !Node 5
      w = (1d0-xx)*yy*(1d0-zz)
      Savn(:,:,gi+y    ) = Savn(:,:,gi+y    ) + Sav(:,:,m)*w;
      Rhon(gi+y    ) = Rhon(gi+y    ) + rho(m)*w;
      fsen(gi+y    ) = fsen(gi+y    ) + ln_fse(m)*w;
      wt(gi+y) = wt(gi+y) + w

      !Node 6
      w = xx*yy*(1d0-zz)
      Savn(:,:,gi+x+y  ) = Savn(:,:,gi+x+y  ) + Sav(:,:,m)*w;
      Rhon(gi+x+y  ) = Rhon(gi+x+y  ) + rho(m)*w;
      fsen(gi+x+y  ) = fsen(gi+x+y  ) + ln_fse(m)*w;
      wt(gi+x+y) = wt(gi+x+y) + w

      !Node 7
      w = xx*yy*zz
      Savn(:,:,gi+x+y+z) = Savn(:,:,gi+x+y+z) + Sav(:,:,m)*w;
      Rhon(gi+x+y+z) = Rhon(gi+x+y+z) + rho(m)*w;
      fsen(gi+x+y+z) = fsen(gi+x+y+z) + ln_fse(m)*w;
      wt(gi+x+y+z) = wt(gi+x+y+z) + w

      !Node 8
      w = (1d0-xx)*yy*zz
      Savn(:,:,gi+y+z  ) = Savn(:,:,gi+y+z  ) + Sav(:,:,m)*w;
      Rhon(gi+y+z  ) = Rhon(gi+y+z  ) + rho(m)*w;
      fsen(gi+y+z  ) = fsen(gi+y+z  ) + ln_fse(m)*w;
      wt(gi+y+z) = wt(gi+y+z) + w

      IF(wt(gi) < 0d0 .OR. Savn(1,1,gi)/wt(gi) > 1500 .OR. Rhon(gi)/wt(gi) < 100) THEN  
      print *,gi,i1,i2,i3,wt(gi),Sav(1,1,m),Savn(1,1,gi)/wt(gi),Rhon(gi)/wt(gi),X1(i1),X2(i2),X3(i3),lr,rr,cr,w
      END IF

      END IF

      END DO

      END IF;END IF;END IF;END IF

END DO

!Normalization
gi=0
do i=1,nodenum*yinyang
  if(wt(i) .gt. 0d0) THEN
     Savn(:,:,i)=Savn(:,:,i)/wt(i)
     Rhon(i)=Rhon(i)/wt(i)
     fsen(i)=fsen(i)/wt(i)
     gi=gi+1
  ELSE
     print *,'Empty node:',i
  end if
end do

if(0 == 1) then
!For empty nodes, interpolate Savn, Rhon from surrounding nodes
print *
print *,'Interpolation for empty nodes is active !!!'
print *

do yyy=1,yinyang
do i3=1,nx31
do i2=1,nx21
do i1=1,nx11
     
  gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z

  if(wt(gi) .eq. 0d0) THEN
     !Find surrounding nodes
     surrnodes = 0d0
     !Node 1
     if(i1 .lt. nx11 .and. wt(gi+x)>0) surrnodes(1) = gi + x
     !Node 2
     if(i1 .gt. 1 .and. wt(gi-x)>0) surrnodes(2) = gi - x
     !Node 3
     if(i2 .lt. nx21 .and. wt(gi+y)>0) surrnodes(3) = gi + y
     !Node 4
     if(i2 .gt. 1 .and. wt(gi-y)>0) surrnodes(4) = gi - y
     !Node 5
     if(i3 .lt. nx31 .and. wt(gi+z)>0) surrnodes(5) = gi + z
     !Node 6
     if(i3 .gt. 1 .and. wt(gi-z)>0) surrnodes(6) = gi - z
     !Node 7
     if(i1 .lt. nx11 .and. i2 .gt. 1 .and. wt(gi+x-y)>0) surrnodes(7) = gi + x - y
     !Node 8
     if(i1 .lt. nx11 .and. i2 .lt. nx21 .and. wt(gi+x+y)>0) surrnodes(8) = gi + x + y
     !Node 9 
     if(i1 .gt. 1 .and. i2 .gt. 1 .and. wt(gi-x-y)>0) surrnodes(9) = gi - x - y
     !Node 10
     if(i1 .gt. 1 .and. i2 .lt. nx21 .and. wt(gi-x+y)>0) surrnodes(10) = gi - x + y
     !Node 11
     if(i1 .lt. nx11 .and. i2 .gt. 1 .and. i3 .lt. nx31 .and. wt(gi+x-y+z)>0) surrnodes(11) = gi + x - y + z
     !Node 12
     if(i1 .lt. nx11 .and. i2 .lt. nx21 .and. i3 .lt. nx31 .and. wt(gi+x+y+z)>0) surrnodes(12) = gi + x + y + z
     !Node 13
     if(i1 .gt. 1 .and. i2 .gt. 1 .and. i3 .lt. nx31 .and. wt(gi-x-y+z)>0) surrnodes(13) = gi - x - y + z
     !Node 14
     if(i1 .gt. 1 .and. i2 .lt. nx21 .and. i3 .lt. nx31 .and. wt(gi-x+y+z)>0) surrnodes(14) = gi - x + y + z
     !Node 15
     if(i1 .lt. nx11 .and. i3 .gt. 1 .and. wt(gi+x-z)>0) surrnodes(15) = gi + x - z
     !Node 16
     if(i1 .lt. nx11 .and. i3 .lt. nx31 .and. wt(gi+x+z)>0) surrnodes(16) = gi + x + z
     !Node 17
     if(i1 .gt. 1 .and. i3 .gt. 1 .and. wt(gi-x-z)>0) surrnodes(17) = gi - x - z
     !Node 18
     if(i1 .gt. 1 .and. i3 .lt. nx31 .and. wt(gi-x+z)>0) surrnodes(18) = gi - x + z
     !Node 19
     if(i2 .lt. nx21 .and. i3 .gt. 1 .and. wt(gi+y-z)>0) surrnodes(19) = gi + y - z
     !Node 20
     if(i2 .lt. nx21 .and. i3 .lt. nx31 .and. wt(gi+y+z)>0) surrnodes(20) = gi + y + z
     !Node 21
     if(i2 .gt. 1 .and. i3 .gt. 1 .and. wt(gi-y-z)>0) surrnodes(21) = gi - y - z
     !Node 22
     if(i2 .gt. 1 .and. i3 .lt. nx31 .and. wt(gi-y+z)>0) surrnodes(22) = gi - y + z
     !Node 23
     if(i1 .lt. nx11 .and. i2 .gt. 1 .and. i3 .gt. 1 .and. wt(gi+x-y-z)>0) surrnodes(23) = gi + x - y - z
     !Node 24
     if(i1 .lt. nx11 .and. i2 .lt. nx21 .and. i3 .gt. 1 .and. wt(gi+x+y-z)>0) surrnodes(24) = gi + x + y - z
     !Node 25
     if(i1 .gt. 1 .and. i2 .gt. 1 .and. i3 .gt. 1 .and. wt(gi-x-y-z)>0) surrnodes(25) = gi - x - y - z
     !Node 26
     if(i1 .gt. 1 .and. i2 .lt. nx21 .and. i3 .gt. 1 .and. wt(gi-x+y-z)>0) surrnodes(26) = gi - x + y - z
        
     i = 0
     do j = 1,26 
        if(surrnodes(j) .gt. 0) then
        if(wt(surrnodes(j)) .gt. 0) then
           i = i + 1
           Savn(:,:,gi)=Savn(:,:,gi) + Savn(:,:,surrnodes(j))
           Rhon(gi)=Rhon(gi) + Rhon(surrnodes(j))
           fsen(gi)=fsen(gi) + fsen(surrnodes(j))
        end if 
        end if 
    end do

    if(i .gt. 0) then
       Savn(:,:,gi) = Savn(:,:,gi)/i
       Rhon(gi) = Rhon(gi)/i
       fsen(gi) = fsen(gi)/i
    end if
    
  end if

end do
end do
end do
end do
    
end if

!Boundary nodes of Yin-Yang grid
if(yinyang > 1) THEN

do yyy=1,yinyang
do i3=1,nx31
do i2=1,nx21
do i1=1,nx11
    
   if(i1 == 1 .OR. i1 == nx11 .OR. i3 == 1 .OR. i3 == nx31) THEN

      gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z

      CALL yin2yang(X1(i1),X3(i3),lr,cr)
      IF(yyy == 1) THEN
         yy = 2
      ELSE
         yy = 1
      END IF
     
      !Find nearest upper left node
      CALL upperleft3D(lr,X2(i2),cr,i,j,k)

      xx = (lr-X1(i))/(X1(i+1)-X1(i))
      if(lr >= X1(nx11)) xx = 1d0
      if(xx > 1d0 .or. xx < 0d0) then
      print *,'X',xx,lr,X1(i),X1(i+1)-X1(i)
      stop
      end if

      zz = (cr-X3(k))/(X3(k+1)-X3(k))
      if(cr >= X3(nx31)) zz = 1d0
      if(zz > 1d0 .or. zz < 0d0) then
      print *,'Z',zz,cr,X3(k),X3(k+1)-X3(k)
      stop
      end if
      
      j = i2
      gijk = nodenum*(yy-1) + i + (j-1)*y + (k-1)*z

      !Node 1
      w1 = (1d0-xx)*(1d0-zz)
      w2 =      xx *(1d0-zz)
      w3 = (1d0-xx)* zz 
      w4 =      xx * zz 

      Savn(:,:,gi) = Savn(:,:,gijk  )*w1 + Savn(:,:,gijk+x  )*w2 + &
                     Savn(:,:,gijk+z)*w3 + Savn(:,:,gijk+x+z)*w4     
      Rhon(    gi) = Rhon(    gijk  )*w1 + Rhon(    gijk+x  )*w2 + &
                     Rhon(    gijk+z)*w3 + Rhon(    gijk+x+z)*w4     
      fsen(    gi) = fsen(    gijk  )*w1 + fsen(    gijk+x  )*w2 + &
                     fsen(    gijk+z)*w3 + fsen(    gijk+x+z)*w4     

  end if

end do
end do
end do
end do

end if
!End Boundary nodes of Yin-Yang grid

DEALLOCATE(wt)

END IF
!END 3D

RETURN

END SUBROUTINE mark2node

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine rotmatrix, build rotation matrix                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rotmatrixZXZ(phi1,theta,phi2,acs)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: phi1,theta,phi2
   DOUBLE PRECISION, DIMENSION(3,3) :: acs

   !Transform Euler angles into direction cosine matrix
   !Z-X-Z: Euler angles in radians
   acs(1,1)=COS(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*SIN(phi2)
   acs(2,1)=COS(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*SIN(phi2)
   acs(3,1)=SIN(phi2)*SIN(theta)

   acs(1,2)=-SIN(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*COS(phi2)
   acs(2,2)=-SIN(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*COS(phi2)
   acs(3,2)=COS(phi2)*SIN(theta)

   acs(1,3)=SIN(theta)*SIN(phi1)
   acs(2,3)=-SIN(theta)*COS(phi1)
   acs(3,3)=COS(theta)
            
   RETURN

   END SUBROUTINE rotmatrixZXZ      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine rotmatrix, build rotation matrix                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rotmatrixZYZ(phi1,theta,phi2,acs)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: phi1,theta,phi2
   DOUBLE PRECISION, DIMENSION(3,3) :: acs

   !Transform Euler angles into direction cosine matrix
   !Z-X-Z: Euler angles in radians
   acs(1,1)=COS(phi2)*COS(phi1)*COS(theta)-SIN(phi1)*SIN(phi2)
   acs(2,1)=COS(phi1)*SIN(phi2)+COS(theta)*COS(phi2)*SIN(phi1)
   acs(3,1)=-COS(phi2)*SIN(theta)

   acs(1,2)=-SIN(phi1)*COS(phi2)-COS(theta)*SIN(phi2)*COS(phi1)
   acs(2,2)=COS(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*SIN(phi2)
   acs(3,2)=SIN(phi2)*SIN(theta)

   acs(1,3)=SIN(theta)*COS(phi1)
   acs(2,3)=-SIN(theta)*SIN(phi1)
   acs(3,3)=COS(theta)
            
   RETURN

   END SUBROUTINE rotmatrixZYZ      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SLOWNESS(gi,PwaveMod,SwaveMod)
   
USE comvar

IMPLICIT NONE
  
INTEGER :: ti(1),gi,i,j,k,ll
DOUBLE PRECISION :: PwaveMod,SwaveMod
DOUBLE PRECISION, DIMENSION(3) :: n,evals,evals0 
DOUBLE PRECISION, DIMENSION(3,3) :: mtens,evects,evects0
DOUBLE PRECISION, DIMENSION(3,3,3,3) :: c

CALL TENS4(Savn(:,:,gi),c)

!Christoffel tensor
mtens = 0d0
DO i = 1 , 3
   DO j = 1 , 3
      DO k = 1 , 3
         DO ll = 1 , 3
            mtens(i,k) = mtens(i,k) + c(i,j,k,ll)*qhat(j)*qhat(ll)
         END DO
      END DO
   END DO
END DO
!   print *,mtens

!Find 3 wave moduli
CALL DSYEVQ3(mtens,evects0,evals0)

!Sort from smallest to biggest
DO j = 1,3
   ti = MINLOC(evals0)
   evects(:,j) = evects0(:,ti(1))
   evals(j) = evals0(ti(1))
   evals0(ti(1))= 1d60
END DO

PwaveMod = evals(3)
SwaveMod = evals(2)
!SwaveMod = 0.5d0*(evals(2) + evals(1))

RETURN

END SUBROUTINE SLOWNESS
