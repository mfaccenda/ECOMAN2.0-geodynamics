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

PROGRAM VIZVISC 

USE comvarvis
USE comvar   

IMPLICIT NONE

CHARACTER (4) :: dt_str4
CHARACTER (500) :: cijkl_dir,output_dir
INTEGER :: m,t
DOUBLE PRECISION :: etac,volinc,r1b,r2b,r1,r2
DOUBLE PRECISION, DIMENSION(3,3) :: acsV

!!! Read database of viscous tensors generated with DEM modelling
CALL loadtensordem

!!! Read initial aggregates distribution and tomographic model grid from input file        
CALL read_infos(cijkl_dir,output_dir)

etac = etac1
volinc = volinc1

DO t = Tinit , Tend, Tstep

   if(t > Tend) stop

   write(dt_str4,'(i4.4)') t

   !!! Read crystal aggregates Fij, Cij, etc.
   write(*,*)
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)
   write(*,'(a)'),' READ AGGREGATE INFOS FROM FILE:'
   write(*,*)

   CALL readcijkl(dt_str4,cijkl_dir)
   
   write(*,'(a)'),' READ AGGREGATE INFOS FROM FILE OK!'
   write(*,*)
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)

   !!! Read crystal aggregates Fij, Cij, etc.
   write(*,*)
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)
   write(*,'(a)'),' COMPUTE FSE SHAPE AND ORIENTATION'
   write(*,*)

   CALL fsecalc
   
   write(*,'(a)'),' COMPUTE FSE SHAPE AND ORIENTATION OK!'
   write(*,*)
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)

   !!! Read crystal aggregates Fij, Cij, etc.
   write(*,*)
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)
   write(*,'(a)'),' INTERPOLATE VISCOUS TENSOR FROM DEM DATABASE'
   write(*,*)

   DO m = 1, marknum    

      !!! Bulk rock aspect ratio
      r1b = log10(amax(m)/amed(m))
      r2b = log10(amed(m)/amin(m))

      !!! Find mean inclusions aspect ratio
      CALL fseinclusions(etac,volinc,r1b,r2b,r1,r2)
      
      !!! Interpolate viscous tesnor from DEM database
      CALL interpolate_tensor(etac,volinc,r1,r2,Sav(:,:,m))

      !!! Limit minimum shear viscosity for weak inclusions
      IF(etac < 1.0d0) CALL viscositycutoff(etac,volinc,Sav(:,:,m))

      !!! Rotate tensor toward fse semiaxes
      !!! WARNINNG: this is not needed when you want to measure anisotropy in the FSE
      !!! reference frame as DEM tensors are already oriented like this
      !acsV(:,1) = Rotm(:,3,m)
      !acsV(:,2) = Rotm(:,2,m)
      !acsV(:,3) = Rotm(:,1,m)
      !CALL rotvoigt(Sav(:,:,m),acsV,Sav(:,:,m)) 

   END DO

   write(*,'(a)'),' INTERPOLATE VISCOUS TENSOR FROM DEM DATABASE OK!'
   write(*,*)
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)

   !!! Calculate FSE semiaxes length and orientation
   
   IF(Lagrangian .NE. 0) THEN

      write(*,*)
      write(*,"(a)"),'********************************************************'
      write(*,"(a)"),'********************************************************'
      write(*,*)
      write(*,"(a)"),' SAVE LAGRANGIAN AGGREGATES PROPERTIES  '
      write(*,*)
      write(*,'(a,i0)'),' Number of aggregates = ',marknum
      write(*,*)
      
      CALL lagrangian_output(dt_str4,cijkl_dir,output_dir)
   
      write(*,"(a)"),' SAVE LAGRANGIAN AGGREGATES PROPERTIES OK!'
      write(*,*)
      write(*,"(a)"),'********************************************************'
      write(*,"(a)"),'********************************************************'
      write(*,*)

   END IF

   IF(Eulerian .NE. 0) THEN

      write(*,*)
      write(*,"(a)"),'********************************************************'
      write(*,"(a)"),'********************************************************'
      write(*,*)
      write(*,"(a)"),' SAVE EULERIAN GRID PROPERTIES  '
      write(*,*)

      CALL eulerian_output(dt_str4,cijkl_dir,output_dir)

      write(*,"(a)"),' SAVE EULERIAN GRID PROPERTIES OK!'
      write(*,*)
      write(*,"(a)"),'********************************************************'
      write(*,"(a)"),'********************************************************'
      write(*,*)

   END IF

DEALLOCATE(mx1,mx2,mx3,mYY,ln_fse,Sav,Fij,rocktype,Rotm,amin,amed,amax)
DEALLOCATE(X1n,X2n,X3n)

END DO

STOP

END PROGRAM VIZVISC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!   Read grid infos    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE read_infos(cijkl_dir,output_dir)

USE comvar

IMPLICIT NONE

CHARACTER (500) :: cijkl_dir,output_dir,command
CHARACTER (50) :: arg,inputfilename
INTEGER :: i,j1,j2,dum(3)

!Read viztomoinput.dat
call read_input_file(cijkl_dir,output_dir)

pi = 3.141592653589793238462643383279
deg2rad = pi/180.0

IF(vpvsmod == 2 .OR. dvpvsmod == 2 .OR. zoeppritzmod > 0) THEN

   qhat(1)= cos(cosx1*deg2rad)
   qhat(2)= cos(cosx2*deg2rad)
   qhat(3)= cos(cosx3*deg2rad)
   write(*,"(a,f6.2,a,f6.2,a,f6.2,a)"),'Unit vector for incoming wave: cosx1 = ',cosx1,', cosx2 = ',cosx2,' cosx3 = ',cosx3,' (deg)'
   IF((SUM(qhat**2.0))**0.5 < 0.9999 .OR. (SUM(qhat**2.0))**0.5 >1.00001) THEN
      write(*,"(a)"),'ERROR: unit vector of incoming wave direction is not 1'
      write(*,"(a)"),'Change cosx1, cosx2, cosx3 such that (SUM(cosx1^2 + cosx2^2 + cosx3^2))^0.5 = 1.0'
      stop
   END IF
      
   Incangle = atan((qhat(1)**2.0 + qhat(3)**2.0)**0.5,qhat(2))
   
END IF

!!! Tensors of indices
ijkl(1,1) = 1 ; ijkl(1,2) = 6 ; ijkl(1,3) = 5
ijkl(2,1) = 6 ; ijkl(2,2) = 2 ; ijkl(2,3) = 4
ijkl(3,1) = 5 ; ijkl(3,2) = 4 ; ijkl(3,3) = 3

l1(1) = 1 ; l1(2) = 2 ; l1(3) = 3
l1(4) = 2 ; l1(5) = 3 ; l1(6) = 1
l2(1) = 1 ; l2(2) = 2 ; l2(3) = 3
l2(4) = 3 ; l2(5) = 1 ; l2(6) = 2
   
mandel_scale = 1d0

DO j1=1,6 ; DO j2=1,6
   IF(j1 .GT. 3 .AND. j2 .GT. 3) THEN
      mandel_scale(j1,j2) = 2
   ELSE IF(j1 .GT. 3 .OR. j2 .GT. 3) THEN
      mandel_scale(j1,j2) = 2**0.5d0
   END IF
END DO ; END DO
   
END SUBROUTINE read_infos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!  Read crystal aggregates elastic tensor !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE readcijkl(dt_str4,cijkl_dir)

   USE comvar
   USE hdf5

   IMPLICIT NONE

   CHARACTER(500) :: filename
   CHARACTER(4) :: dt_str4
   CHARACTER (len=*) :: cijkl_dir
   CHARACTER (len(trim(cijkl_dir))) :: i_dir
   INTEGER :: i,m,m1,dum_int(3)
   DOUBLE PRECISION :: datadb(3)
   DOUBLE PRECISION ,DIMENSION(:,:), ALLOCATABLE :: Xsave
   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl,memtype !Handles
   INTEGER     ::   error  ! Error flag
   LOGICAL     ::   link_exists

   i_dir=trim(cijkl_dir)

   filename=i_dir//'Cijkl'//dt_str4//'.h5'

   write(*,'(a,a)'),' ',trim(filename)
   write(*,*)
         
   !Initialize FORTRAN interface.

   CALL H5open_f (error)

   ! Create a new file using default properties.

   CALL H5Fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

   !Load domain infos      
   CALL loadsave_integer(1,1,file_id,2,H5T_NATIVE_INTEGER,dum_int(1:2),'CoordinateSystem',0)
   dimensions = dum_int(1)
   cartspher  = dum_int(2)

   !Swap axes 2 and 3 in spherical coordinates as tensors are rotated such that
   !the vertical axis is parallel to axis 3.
   IF(cartspher == 2) THEN
      qhat(2)= cos(cosx3*deg2rad)
      qhat(3)= cos(cosx2*deg2rad)
   END IF   

   CALL loadsave_double(1,1,file_id,2,H5T_NATIVE_DOUBLE,datadb,'Time',0)
   dt=datadb(1)
   timesum=datadb(2)

   CALL loadsave_integer(1,1,file_id,3,H5T_NATIVE_INTEGER,dum_int,'Periodic',0)
   x1periodic = dum_int(1)
   x2periodic = dum_int(2)
   x3periodic = dum_int(3)
   yinyang = 1
   IF(x1periodic > 0 .AND. x3periodic > 0 .AND. cartspher == 2 .AND. dimensions == 3) yinyang = 2

   !Convert time from sec to Myr
   !timesum = timesum /365.25/86400/1e+6

   !Define model grid
   write(*,"(a)"),' MODEL GRID'
   write(*,*)
   if(dimensions==2) write(*,"(a,i5,a)"),' dimensions:   ',dimensions,' --> 2D model'
   if(dimensions==3) write(*,"(a,i5,a)"),' dimensions:   ',dimensions,' --> 3D model'
   if(cartspher==1)  write(*,"(a,i5,a)"),' cartspher :   ',cartspher, ' --> cartesian coordinates system'
   if(cartspher==2 .AND. dimensions == 2) write(*,"(a,i5,a)"),' cartspher :   ',cartspher,' --> polar coordinates system'
   if(cartspher==2 .AND. dimensions == 3) write(*,"(a,i5,a)"),' cartspher :   ',cartspher,' --> spherical coordinates system'

   !Open group with Eulerian nodes infos
   CALL H5Gopen_f(file_id, "/Nodes", group_id, error)
   !Read Attribute
   CALL loadsave_integer(1,1,group_id,dimensions,H5T_NATIVE_INTEGER,dum_int,'nxyz',0)
   nx1=dum_int(1)
   nx2=dum_int(2)
   nx3=dum_int(3)

   ALLOCATE(X1n(nx1),X2n(nx2),X3n(nx3))


   !Read Datasets  
   CALL loadsave_double(0,1,group_id,nx1,H5T_NATIVE_DOUBLE,X1n,'gx',0)
   CALL loadsave_double(0,1,group_id,nx2,H5T_NATIVE_DOUBLE,X2n,'gy',0)
   CALL loadsave_double(0,1,group_id,nx3,H5T_NATIVE_DOUBLE,X3n,'gz',0)

   !IF polar coordinates, i1first, i1last, and x10size are in radians
   i1first = X1n(1)
   x1max = i1last
   i2first = X2n(1)
   x2max = i2last
   i1last  = X1n(nx1)
   i2last  = X2n(nx2)
   x10size = i1last-i1first
   x20size = i2last-i2first

   IF(dimensions == 3) THEN 
      i3first = X3n(1)
      i3last  = X3n(nx3)
      x3max = i3last
      x30size = i3last-i3first
   END IF

   CALL H5Gclose_f(group_id, error)

   IF(cartspher==2) THEN
      n1first = n1first*deg2rad ; n1last = n1last*deg2rad
      n3first = n3first*deg2rad ; n3last = n3last*deg2rad
   END IF

   CALL H5Gopen_f(file_id, "/Particles", group_id, error)

   CALL loadsave_integer(1,1,group_id,1,H5T_NATIVE_INTEGER,dum_int(1:1),'marknum',0)
   marknum = dum_int(1)
   marknum0 = marknum
   IF (reflectmod /= 0) marknum0 = marknum * 2

   !Lagrangian aggregates domain in DREX_S
   CALL loadsave_double(0,1,group_id,3,H5T_NATIVE_DOUBLE,datadb(1:3),'mxstp',0)
   mx1stp = datadb(1)
   mx2stp = datadb(2)
   mx3stp = datadb(3)
   CALL loadsave_double(0,1,group_id,3,H5T_NATIVE_DOUBLE,datadb(1:3),'mxmin',0)
   mx1min = datadb(1)
   mx2min = datadb(2)
   mx3min = datadb(3)
   CALL loadsave_double(0,1,group_id,3,H5T_NATIVE_DOUBLE,datadb(1:3),'mxmax',0)
   mx1max = datadb(1)
   mx2max = datadb(2)
   mx3max = datadb(3)

   ALLOCATE(rocktype(marknum0),mx1(marknum0),mx2(marknum0),mx3(marknum0),Fij(3,3,marknum0))
   ALLOCATE(ln_fse(marknum0),Rotm(3,3,marknum0),amin(marknum0),amed(marknum0),amax(marknum0))
   ALLOCATE(mYY(marknum0),Sav(6,6,marknum0))
   
   CALL loadsave_integer(0,1,group_id,marknum,H5T_NATIVE_INTEGER,rocktype(1:marknum),'rocktype',0)
   CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mx1,'mx1',0)
   CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mx2,'mx2',0)
   IF(dimensions == 3) THEN
      CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mx3,'mx3',0)
   ELSE 
      IF(cartspher == 1) mx3 = 0d0
      IF(cartspher == 2) mx3 = pi/2.0
   END IF

   mYY = 1
   IF(yinyang == 2) CALL loadsave_integer(0,1,group_id,marknum,H5T_NATIVE_INTEGER,mYY(1:marknum),'mYY',0)

   !FSE
   dum_int(1)=3
   dum_int(2)=3
   dum_int(3)=marknum
   CALL loadsave_double(0,3,group_id,dum_int,H5T_NATIVE_DOUBLE,Fij,'fse',0)

   CALL H5Gclose_f(group_id, error)

   !Terminate access to the file.

   CALL H5Fclose_f(file_id, error)

   CALL H5close_f(error)
   !Close FORTRAN interface.
  
   END SUBROUTINE readcijkl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine fsecalc, calculate fse semiaxes and their orientation       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE fsecalc   

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: j,m,ti(1)
   DOUBLE PRECISION, DIMENSION(3) :: evals,c2
   DOUBLE PRECISION, DIMENSION (3,3) :: LSij,evects,fseacs

   ln_fse = -1d0

!$omp parallel & 
!$omp shared(rocktype,Fij,ln_fse,Rotm,amin,amed,amax) &
!$omp private(j,m,LSij,evals,evects,fseacs,c2,ti) &    
!$omp firstprivate(marknum)
!$omp do schedule(guided,8)
   DO m = 1,marknum

      IF(rocktype(m) > 0 .AND. rocktype(m) < 10) THEN

!!! Left-stretch tensor for FSE calculation
      LSij = MATMUL(Fij(:,:,m),TRANSPOSE(Fij(:,:,m)))
      CALL DSYEVQ3(LSij,evects,evals)

      !Order from smallest to largest semiaxis
      DO j = 1,3
         ti = MINLOC(evals) ;fseacs(:,j) = evects(:,ti(1)) ; c2(j) = evals(ti(1))**0.5 ; evals(ti(1))= 1d60
      END DO

      !Save semiaxes orientation
      Rotm(:,:,m) = fseacs(:,:)

      !Save semiaxes lengths
      amin(m) = c2(1) ; amed(m) =c2(2) ; amax(m) = c2(3)
      IF(amax(m) > 1e+3) THEN
          ln_fse(m) = 0.0
      ELSE
!!! natural strain = ln(a/c) where a is the long axis = maxval(evals)**0.5
          ln_fse(m) = LOG(amax(m)/amin(m))
      END IF

      END IF

   END DO

!$omp end do
!$omp end parallel

   END SUBROUTINE fsecalc   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
