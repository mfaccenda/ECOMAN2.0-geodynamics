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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine eulerian_output, saving Vp,Vs,dVp,dVs,radial & azimuthal anisi, etc. !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE eulerian_output(dt_str4,cijkl_dir,output_dir)

USE comvar
USE hdf5   

IMPLICIT NONE

INTEGER :: i,j,k,i1,i2,i3,m,m1,gi,giref,gijk,ijk,x,y,z,t,zn,yyy
INTEGER :: yc,zc,dum_int(3),azinum,flag
DOUBLE PRECISION, DIMENSION (3) :: n1,n2,n3   
DOUBLE PRECISION, DIMENSION (3,3) :: acs     
DOUBLE PRECISION, DIMENSION (3,3,3,3) :: cijkl
DOUBLE PRECISION :: Nparam,Lparam,Gc,Gs,G3D,phi,lr,cr
DOUBLE PRECISION :: xx,yy,zz,phi1,theta,phi2,azi(3)
!Buffers
INTEGER, DIMENSION(:), ALLOCATABLE :: idx
INTEGER, DIMENSION(:,:), ALLOCATABLE :: connectivity 
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: E_chi,G
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: vertices,azimuthal,xyz 

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

IF(nodenum > 16777216) write(*,"(a,i10,a)"),' WARNING : nodenum = ',nodenum,' > 2^24 = 16777216. There is an issue with Paraview reading grids with such a large number of nodes from HDF5 files. Decrease the grid resolution' 

!Allocate arrays
ALLOCATE(Savn(6,6,nodenum*yyy),fsen(nodenum*yyy))
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
IF(n2first < i2first*0.9999) WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : n2first = ',n2first,' < i2first = ',i2first,'. Set n2first >= i2first'
IF(n2last  > i2last*1.0001 ) WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : n2last  = ',n2last ,' > i2last  = ',i2last ,'. Set n2last  >= i2last '

IF(dimensions == 3) THEN
   IF(cartspher == 1 .AND. nstepx3 < mx3stp*0.9999)  WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : nstepx3 = ',nstepx3,' < mx3stp = ',mx3stp,'. Decrease nx31'
   IF(cartspher == 2 .AND. nstepx3*X2(1) < mx3stp*0.9999)  WRITE(*,"(a,f12.4,a,f12.4,a)"),' WARNING : nstepx3 = ',nstepx3*X2(1),' < mx3stp = ',mx3stp,'. Decrease nx31'
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

!!! Rotate tensor toward Z axis direction if polar/spherical coordinates
!!! WARNINNG: this is not needed when you want to measure anisotropy in the FSE
!!! reference frame as DEM tensors are already oriented like this
IF(1==0 .AND. cartspher==2 .AND. (radialmod > 0 .OR. azimod > 0)) THEN

   !$omp parallel & 
   !$omp shared(Savn,fsen,X1,X3) &
   !$omp private(yyy,i1,i2,i3,gi,phi1,theta,phi2,acs,lr,cr) &    
   !$omp firstprivate(pi,ln_fse_min,x,y,z,nx11,nx21,nx31)
   !$omp do schedule(guided,8)
   DO yyy=1,yinyang
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z
      
      IF(fsen(gi) >= ln_fse_min) THEN

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
     
IF(radialmod .NE. 0) THEN
      
   ALLOCATE(E_chi(nodenum*yyy))

   E_chi = 1d0
   
   !$omp parallel & 
   !$omp shared(Savn,fsen,E_chi) &
   !$omp private(yyy,i1,i2,i3,gi,Nparam,Lparam) &    
   !$omp firstprivate(radialmod,ln_fse_min,x,y,z,nx11,nx21,nx31)
   !$omp do schedule(guided,8)
   DO yyy=1,yinyang
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      gi = nodenum*(yyy-1) + i1 + (i2-1)*y + (i3-1)*z
      
      IF(fsen(gi) >= ln_fse_min) THEN     
     
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
   !$omp shared(Savn,fsen,G) &
   !$omp private(i1,i2,i3,gi,Gc,Gs) &    
   !$omp firstprivate(azimod,ln_fse_min,x,y,z,nx11,nx21,nx31)
   !$omp do schedule(guided,8)
   DO i3=1,nx31
   DO i2=1,nx21
   DO i1=1,nx11

      gi = i1 + (i2-1)*y + (i3-1)*z
      
      IF(fsen(gi) >= ln_fse_min) THEN     
     
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
      IF(fsen(gi) >= ln_fse_min) THEN
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
      IF(fsen(gi) >= ln_fse_min) THEN

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
IF(radialmod > 0 .OR. azimod > 0) THEN

   ! Create a new file using default properties.
   filename = str2//'eulerian'//dt_str4//'.h5'
  
   !Initialize FORTRAN interface.
   CALL H5open_f (error)

   CALL H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

   !Make grid
   !3D grid 
   IF(dimensions == 3) THEN

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
        
   IF(radialmod > 0) THEN
      print *,''
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
IF(radialmod > 0 .OR. (dimensions == 2 .AND. azimod > 0)) THEN
   
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

   IF(radialmod > 0) THEN
      WRITE(19,'(a)') '         <Attribute Type="Scalar" Center="Node" Name="Radial Anisotropy">' 
      WRITE(19,'(a,i0,a)') '            <DataItem ItemType="HyperSlab" Dimensions="',nodenum,' 1">'
      WRITE(19,'(a,i0,a)') '            <DataItem Dimensions="3 2" Format="XML"> 0 0 1 1 ',nodenum,' 1 </DataItem>'
      WRITE(19,'(a,i0,a,a,a)') '            <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',nodenum,' 1">eulerianY',dt_str4,'.h5:/RadialAnisotropy</DataItem>'
      WRITE(19,'(a)') '            </DataItem>' 
      WRITE(19,'(a)') '         </Attribute>' 
      WRITE(19,'(a)') ' '
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

IF(radialmod > 0 .OR. azimod > 0) THEN

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
Savn = 0 ; fsen = 0; wt = 0

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
      fsen(gi      ) = fsen(gi      ) + ln_fse(m)*w;
      wt(gi) = wt(gi) + w
      
      !Node 2
      w = xx*(1d0-yy)
      Savn(:,:,gi+x    ) = Savn(:,:,gi+x    ) + Sav(:,:,m)*w;
      fsen(gi+x    ) = fsen(gi+x    ) + ln_fse(m)*w;
      wt(gi+x) = wt(gi+x) + w

      !Node 3
      w = (1d0-xx)*yy
      Savn(:,:,gi+y    ) = Savn(:,:,gi+y    ) + Sav(:,:,m)*w;
      fsen(gi+y    ) = fsen(gi+y    ) + ln_fse(m)*w;
      wt(gi+y) = wt(gi+y) + w

      !Node 4
      w = xx*yy
      Savn(:,:,gi+x+y  ) = Savn(:,:,gi+x+y  ) + Sav(:,:,m)*w;
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
     fsen(i)=fsen(i)/wt(i)
     gi=gi+1
  ELSE
     print *,'Empty node:',i
  END IF
END DO

IF(0 == 1) THEN
!For empty nodes, interpolate Savn from surrounding nodes
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
           fsen(gi)=fsen(gi) + fsen(surrnodes(j))
        END IF 
        END IF 
    END DO

    IF(i > 0) THEN
       Savn(:,:,gi) = Savn(:,:,gi)/i
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
      fsen(gi      ) = fsen(gi      ) + ln_fse(m)*w;
      wt(gi) = wt(gi) + w
      
      !Node 2
      w = xx*(1d0-yy)*(1d0-zz)
      Savn(:,:,gi+x    ) = Savn(:,:,gi+x    ) + Sav(:,:,m)*w;
      fsen(gi+x    ) = fsen(gi+x    ) + ln_fse(m)*w;
      wt(gi+x) = wt(gi+x) + w

      !Node 3
      w = xx*(1d0-yy)*zz
      Savn(:,:,gi+x+z  ) = Savn(:,:,gi+x+z  ) + Sav(:,:,m)*w;
      fsen(gi+x+z  ) = fsen(gi+x+z  ) + ln_fse(m)*w;
      wt(gi+x+z) = wt(gi+x+z) + w

      !Node 4
      w = (1d0-xx)*(1d0-yy)*zz
      Savn(:,:,gi+z    ) = Savn(:,:,gi+z    ) + Sav(:,:,m)*w;
      fsen(gi+z    ) = fsen(gi+z    ) + ln_fse(m)*w;
      wt(gi+z) = wt(gi+z) + w

      !Node 5
      w = (1d0-xx)*yy*(1d0-zz)
      Savn(:,:,gi+y    ) = Savn(:,:,gi+y    ) + Sav(:,:,m)*w;
      fsen(gi+y    ) = fsen(gi+y    ) + ln_fse(m)*w;
      wt(gi+y) = wt(gi+y) + w

      !Node 6
      w = xx*yy*(1d0-zz)
      Savn(:,:,gi+x+y  ) = Savn(:,:,gi+x+y  ) + Sav(:,:,m)*w;
      fsen(gi+x+y  ) = fsen(gi+x+y  ) + ln_fse(m)*w;
      wt(gi+x+y) = wt(gi+x+y) + w

      !Node 7
      w = xx*yy*zz
      Savn(:,:,gi+x+y+z) = Savn(:,:,gi+x+y+z) + Sav(:,:,m)*w;
      fsen(gi+x+y+z) = fsen(gi+x+y+z) + ln_fse(m)*w;
      wt(gi+x+y+z) = wt(gi+x+y+z) + w

      !Node 8
      w = (1d0-xx)*yy*zz
      Savn(:,:,gi+y+z  ) = Savn(:,:,gi+y+z  ) + Sav(:,:,m)*w;
      fsen(gi+y+z  ) = fsen(gi+y+z  ) + ln_fse(m)*w;
      wt(gi+y+z) = wt(gi+y+z) + w

      IF(wt(gi) < 0d0 .OR. Savn(1,1,gi)/wt(gi) > 1500000) THEN  
      print *,gi,i1,i2,i3,wt(gi),Sav(1,1,m),Savn(1,1,gi)/wt(gi),X1(i1),X2(i2),X3(i3),lr,rr,cr,w
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
     fsen(i)=fsen(i)/wt(i)
     gi=gi+1
  else
     print *,'Empty node:',i
  end if
end do

if(0 == 1) then
!For empty nodes, interpolate Savn from surrounding nodes
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
           fsen(gi)=fsen(gi) + fsen(surrnodes(j))
        end if 
        end if 
    end do

    if(i .gt. 0) then
       Savn(:,:,gi) = Savn(:,:,gi)/i
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

