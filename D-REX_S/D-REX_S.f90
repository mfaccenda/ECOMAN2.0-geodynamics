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

   PROGRAM DREX_S  

   USE comvar
   USE hdf5    

   IMPLICIT NONE

   INTEGER :: step1 !! number of points on a streamline

   INTEGER :: i,j,m,n,nx(3),ti(1)

   DOUBLE PRECISION, DIMENSION(3) :: evals,c2
   DOUBLE PRECISION, DIMENSION (3,3) ::LSij,evects,fseacs,acs1,acs2,Rotm,ee,Rotm1,lrot
   DOUBLE PRECISION :: fractdisl,phi1,theta,phi2,a0,w3
   DOUBLE PRECISION :: pphi1,ttheta,pphi2
   DOUBLE PRECISION, DIMENSION(6,6) :: Cstilwe
   DOUBLE PRECISION, DIMENSION(1000,6,6) :: Cdem

   CHARACTER (3) :: dt_str3
   CHARACTER (4) :: dt_str4
   CHARACTER (100) :: fname,str

   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl,memtype ! Handles
   INTEGER     ::   error  ! Error flag

!!! initialization

   CALL init0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  LPO  and FSE calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Constant timestep   
   dt = 1d-2/epsnot(1)

   write(*,"(a,f8.2)"),' Timestep = ',dt
   write(*,*)

   fractdisl = fractdislrock(rocktype(1))

   n = 0

   DO WHILE (max_strain <= strain_max)

      n = n + 1

!!! FSE and LPO calculation

      IF(strain_max > 0) THEN

         CALL strain(1,1,fractdisl)

         max_strain = max_strain + dt*epsnot(1)!*fractdisl

         timesum = timesum + dt

      END IF

      IF(n==10 .OR. max_strain == strain_max) THEN

      write(*,"(a,f8.2)"),' Finite strain:',max_strain

!!! Cijkl tensor (using Voigt average)
      CALL stifftenz(1)

!!! Percentage of anisotropy and orientation of axis of hexagonal symmetry
      !CALL DECSYM(Voigt,perc_a,phi_a)

      !OUTPUT      

      !CALL output
      
      n = 0

      !!! Write infos in hdf5 format

      !Initialize FORTRAN interface.

      CALL H5open_f (error)
   
      ! Create a new file using default properties.
      IF(max_strain < 9.999) THEN 
         write(dt_str3,'(1f3.1)') max_strain
         fname = trim(output_name)//dt_str3//'.h5'
      END IF
      IF(max_strain>= 9.999 .AND. max_strain < 100) THEN 
         write(dt_str4,'(1f4.1)') max_strain
         fname = trim(output_name)//dt_str4//'.h5'
      END IF

      print *
      print *,' Save to ',trim(fname)
      print *

      CALL H5Fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)

      nx(1)=size
      nx(2)=3
      nx(3)=3
      CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs(:,:,:,1),'acs',1)
      CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,odf(1,:),'odf',1)
      IF(rocktype(1) == 1) THEN
         CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs_ens(:,:,:,1),'acs_ens',1)
         CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,odf_ens(1,:),'odf_ens',1)
      END IF
      nx(1)=3
      nx(2)=3
      CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Fij(:,:,1),'Fij',1)
      nx(1)=6
      nx(2)=6
      CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Voigt,'Voigt',1)
      CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Reuss,'Reuss',1)
      CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Mixed,'Mixed',1)
      !Density
      CALL loadsave_double(0,1,file_id,1,H5T_NATIVE_DOUBLE,rho,'Density',1)
      !Rocktype
      CALL loadsave_integer(0,1,file_id,1,H5T_NATIVE_INTEGER,rocktype,'Rocktype',1)

      !Terminate access to the file.
      CALL H5Fclose_f(file_id, error)

      !Close FORTRAN interface.
      CALL H5close_f(error)

      IF(strain_max == 0) max_strain = 1.0

      END IF
      !END OUTPUT

   END DO

   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)
   print *,'Voigt average'
   write(*,'(6f10.3)') Voigt
   print *
   print *,'Reuss average'
   write(*,'(6f10.3)') Reuss
   print *
   print *,'Mixed Voigt/Ruess average'
   write(*,'(6f10.3)') Mixed
   print *
   write(*,"(a)"),'--------------------------------------------------------'
   print *

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  SPO calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   IF(spomod > 0) THEN
  
      pi = 3.141592653589793238462643383279
      deg2rad = pi/180.0

      !!! Read input file infos for DEM model and allocate TD database matrices
      CALL initspo

      i = 1       

      !Find principal axes of FSE
      IF(spomod == 1) THEN

         LSij = MATMUL(Fij(:,:,1),TRANSPOSE(Fij(:,:,1)))
         CALL DSYEVQ3(LSij,evects,evals)

         !Order from smallest to largest semiaxis
         DO j = 1,3
            ti = MINLOC(evals) ;fseacs(:,j) = evects(:,ti(1)) ; c2(j) = evals(ti(1))**0.5 ; evals(ti(1))= 1d60
         END DO

         !Save semiaxes orientation
         !1st column: orientaton mininum axis
         !2nd column: orientaton medium  axis
         !3rd column: orientaton maximum axis
         Rotm1 = fseacs
         IF(phi_spo /=0) CALL rot3Dspo(fseacs(:,:),Rotm1(:,:),fseacs(:,2),phi_spo)

         !Rotate tensor parallel to max semiaxis of FSE
         phi1  = atan2(Rotm1(1,1),-Rotm1(2,1)) !Angle of FSE max semiaxis with x1 axis 
         theta =  acos(Rotm1(3,1)) !Angle of FSE min semiaxis with x3 axis
         phi2  = atan2(Rotm1(3,3),Rotm1(3,2))
         CALL rotmatrixZXZ(phi1,theta,phi2,acs2)

         !STILWE (Backus, 1962, JGR)
         CALL stilwe(1,Sav(:,:,i)) 
         
         IF(ptmod == 0) rho(1) = ro_back*(1.0-Vmax) + ro_incl*Vmax

      END IF
    
      IF(spomod > 1) THEN

         ee = e(1,:,:)
         CALL DSYEVQ3(ee,evects,evals)
         
         !Order from smallest to largest semiaxis
         DO j = 1,3
            ti = MINLOC(evals) ;fseacs(:,j) = evects(:,ti(1)) ; c2(j) = evals(ti(1)) ; evals(ti(1))= 1d60
         END DO

         !Save semiaxes orientation
         !1st column: orientaton compressive  axis
         !2nd column: orientaton intermediate axis
         !3rd column: orientaton extensionsal axis
         Rotm1 = fseacs

         !Rotate tensor parallel to compressive (minimum) strain rate axis
         phi1  = atan2(Rotm1(1,3),-Rotm1(2,3)) !Angle of compressive strain rate axis with x1 axis 
         theta =  acos(Rotm1(3,3)) !Angle of extensional strain rate axis with x3 axis
         phi2  = atan2(Rotm1(3,1),Rotm1(3,2))

         a0 = phi_spo
         IF(phi_spo /=0) THEN
            !Rotate velocity gradient tensor parallel to intermediate strain rate axis
            pphi1  = atan2(Rotm1(1,2),-Rotm1(2,2)) 
            ttheta =  acos(Rotm1(3,2))
            pphi2  = 0d0
            CALL rotmatrixZXZ(pphi1,ttheta,pphi2,acs1)
            !L' = R'*L*R
            lrot = MATMUL(TRANSPOSE(acs1),l(1,:,:))
            lrot = MATMUL(lrot,acs1)
            !Vorticity
            w3 = 0.5*(lrot(2,1) - lrot(1,2))
            IF(w3 < 0) a0 = -phi_spo
            !CALL rot3Dspo(fseacs(:,:),Rotm1(:,:),fseacs(:,2),a0)

         END IF

         IF(spomod == 3) Sav(:,:,i) = Mixed

         CALL rotmatrixZXZ(phi1,theta,phi2,acs2)
         CALL tensorrot_aggr(i,TRANSPOSE(acs2)) 
   
         !Now rotate around Y axis if phi_spo different from 0
         IF(phi_spo /= 0) THEN
            CALL rotmatrixZYZ(0d0,a0,0d0,acs1)
            CALL tensorrot_aggr(i,TRANSPOSE(acs1)) !Reverse rotation with transpose rot. matrix
         END IF

         !DEM (Mainprice, 1997, EPSL)
         CALL DEM(1,Vmax,Cdem) 

         !Find DEM elastic tensor according to porosity of the geodynamic model
         j = INT(Vmax/Vstp) + 1
         Sav(:,:,i) = Cdem(j,:,:)

         IF(ptmod == 0 .OR. spomod == 2) rho(1) = ro_back*(1.0-Vmax) + ro_incl*Vmax
         IF(ptmod >  0) rho(1) = rho(1)*(1.0-Vmax)  + ro_incl*Vmax

         !Now rotate around Y axis if phi_spo different from 0
         IF(phi_spo /= 0) THEN
            CALL tensorrot_aggr(i,acs1) !Reverse rotation with transpose rot. matrix
         END IF

      END IF

      !Rotate tensor parallel to:
      !spomod = 1: max semiaxis of FSE 
      !spomod > 1: compressive axis
      CALL tensorrot_aggr(i,acs2) 
    
      Mixed = Sav(:,:,i)

      fname = trim(output_name)//'SPO.h5'
      print *
      print *,'Save to ',trim(fname)
      print *

      !!! Write infos in hdf5 format

      !Initialize FORTRAN interface.

      CALL H5open_f (error)
   
      CALL H5Fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)

      nx(1)=3
      nx(2)=3
      IF(spomod == 1) CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Fij(:,:,1),'Fij',1)
      IF(spomod >  1) CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,ee,'Eij',1)
      nx(1)=6
      nx(2)=6
      CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Mixed,'Mixed',1)
      !Density
      CALL loadsave_double(0,1,file_id,1,H5T_NATIVE_DOUBLE,rho,'Density',1)
      !Rocktype
      CALL loadsave_integer(0,1,file_id,1,H5T_NATIVE_INTEGER,spomod,'spomod',1)

      !Terminate access to the file.
      CALL H5Fclose_f(file_id, error)

      !Close FORTRAN interface.
      CALL H5close_f(error)

      print *,'Mixed Voigt/Ruess average + SPO'
      IF(spomod == 1) print *,'Tensor rotated parallel to FSE max semiaxis'
      IF(spomod >  1) print *,'Tensor rotated parallel to compressive axis'
      write(*,'(6f10.3)') Mixed 
      print *
      write(*,"(a)"),'--------------------------------------------------------'
      print *
 
   END IF

   str='mkdir '//trim(output_name)
   CALL system(str)
   str='mv '//trim(output_name)//'*.h5 '//trim(output_name)
   CALL system(str)

   !!! Write infos in hdf5 format

   !Initialize FORTRAN interface.

   CALL H5open_f (error)

   CALL H5Fcreate_f('fossilfabric.h5', H5F_ACC_TRUNC_F, file_id, error)

   nx(1)=size
   nx(2)=3
   nx(3)=3
   CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs(:,:,:,1),'acs',1)
   CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,odf(1,:),'odf',1)
   CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs_ens(:,:,:,1),'acs_ens',1)
   CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,odf_ens(1,:),'odf_ens',1)
   nx(1)=3
   nx(2)=3
   CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,Fij(:,:,1),'Fij',1)

   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   !Close FORTRAN interface.
   CALL H5close_f(error)

   STOP

   END PROGRAM DREX_S

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine INIT0, initialization of variables, arrays and parameters   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE init0

   USE comvar
   USE hdf5  

   IMPLICIT NONE

   INTEGER :: gi,j,i1,i3,i,j1,j2,j3,nrot ! loop counters

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ran0
   ! matrix of random numbers used to generate initial random LPO

   DOUBLE PRECISION :: phi1,theta,phi2
   ! eulerian angles

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xe1,xe2,xe3
   ! matrixes of initial random eulerian angles

   DOUBLE PRECISION, DIMENSION(3) :: evals
   DOUBLE PRECISION, DIMENSION(3,3) :: evects,ee
   ! eigen values and vectors in jacobi

   DOUBLE PRECISION phi
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phix0,phiy0,phiz0,phiz1
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Rhodum,Vpdum,Vsdum

   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl,memtype ! Handles
   INTEGER     ::   error  ! Error flag

   l = 0d0 ; e = 0d0

!!! Name of input files from file input_fabric0.dat

   CALL read_input_file

!!! initial size
   size = size3**3

   Xol = Xol /100.0

!!! strain rate tensor
   e(1,1,1) = l(1,1,1) ; e(1,3,3) = l(1,3,3) ; e(1,2,2) = l(1,2,2)
   e(1,3,1) = (l(1,1,3)+l(1,3,1))/2d0 ; e(1,1,3) = e(1,3,1)
   e(1,1,2) = (l(1,2,1)+l(1,1,2))/2d0 ; e(1,2,1) = e(1,1,2)
   e(1,2,3) = (l(1,3,2)+l(1,2,3))/2d0 ; e(1,3,2) = e(1,2,3)

!! reference strain rate
   ee = e(1,:,:)
   CALL DSYEVQ3(ee,evects,evals)
   epsnot(1) = MAXVAL(ABS(evals))

   write(*,"(a)"),' VELOCITY GRADIENT TENSOR Dij'
   write(*,*)
   write(*,"(3f8.2)"),(l(1,1,:))
   write(*,"(3f8.2)"),(l(1,2,:))
   write(*,"(3f8.2)"),(l(1,3,:))
   write(*,*)
   write(*,"(a)"),' STRAIN RATE tensor'
   write(*,*)
   write(*,"(3f8.2)"),(e(1,1,:))
   write(*,"(3f8.2)"),(e(1,2,:))
   write(*,"(3f8.2)"),(e(1,3,:))
   write(*,*)
   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)

   max_strain = 0d0
        
!!! Initial deformation gradient tensor
   Fij = 0d0 ; Fij(1,1,:) = 1d0 ; Fij(2,2,:) = 1d0 ; Fij(3,3,:) = 1d0

!!! tensor \epsilon_{ijk}

   alt=0d0
   alt(1,2,3) = 1d0 ; alt(2,3,1) = 1d0 ; alt(3,1,2) = 1d0
   alt(1,3,2) = -1d0 ; alt(2,1,3) = -1d0 ; alt(3,2,1) = -1d0

!!! tensor \delta_{ij}

   del=0d0
   del(1,1) = 1d0 ; del(2,2) = 1d0 ; del(3,3) = 1d0

!!! tensors of indices

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

!!! Loading stiffness tensors (GPa)

   CALL elastic_database(S0,dS0dp,dS0dp2,dS0dt,dS0dpt,dS0dp5,dS0dt5)

!!! allocation of the dimensions of the arrays

   ALLOCATE(xe1(size),xe2(size),xe3(size))

   ALLOCATE(odf(1,size))
   ALLOCATE(odf_ens(1,size))

   ALLOCATE(ran0(3*size))

   ALLOCATE(acs(size,3,3,1),acs0(size,3,3))
   ALLOCATE(acs_ens(size,3,3,1))

!!! initialization of orientations - uniformally random distribution
!!! Rmq cos(theta) used to sample the metric Eulerian space

   CALL RANDOM_NUMBER(ran0)

   i = 1

   DO j1 =1, size3 ; DO j2 =1, size3 ; DO j3 =1, size3
      xe1(i) = (REAL(j1)-ran0(i))/REAL(size3)*ACOS(-1d0)
      xe2(i) = ACOS(-1d0 + (REAL(j2)-ran0(size+i))/REAL(size3)*2d0)
      xe3(i) = (REAL(j3)-ran0(i+2*size))/REAL(size3)*ACOS(-1d0)
      i = i + 1
   END DO ; END DO ; END DO

   DO i = 1 , size

      phi1 = xe1(i) ; theta = xe2(i) ; phi2 = xe3(i)

!!! Direction cosine matrix
!!! acs(k,j) = cosine of the angle between the kth crystallographic axes and the jth external axes 

      acs0(i,1,1)=COS(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*SIN(phi2)
      acs0(i,1,2)=COS(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*SIN(phi2)
      acs0(i,1,3)=SIN(phi2)*SIN(theta)

      acs0(i,2,1)=-SIN(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*COS(phi2)
      acs0(i,2,2)=-SIN(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*COS(phi2)
      acs0(i,2,3)=COS(phi2)*SIN(theta)

      acs0(i,3,1)=SIN(theta)*SIN(phi1)
      acs0(i,3,2)=-SIN(theta)*COS(phi1)
      acs0(i,3,3)=COS(theta)

   END DO
   
!!! Set random initial LPO and same grain size

   DO i = 1 , 1        

      odf(i,:) = 1d0/REAL(size3**3)
      odf_ens(i,:) = odf(i,:)
      acs(:,:,:,i) = acs0
      acs_ens(:,:,:,i) = acs0

   END DO

   DEALLOCATE(xe1,xe2,xe3,ran0)

!!! Cartesian coordinates of a spherical object    

   degstp=15 !!! step in degrees
   nxy=360/degstp
   nz=90/degstp + 1
    
   ALLOCATE(phix0(nxy),phiy0(nxy),phiz0(nz),phiz1(nz))
   ALLOCATE(phix(nz,nxy),phiy(nz,nxy),phiz(nz,nxy))

   DO i = 1 , nz 
      phi = degstp*(i-1)*acos(-1.0)/180
      phiz0(i) = cos(phi); 
      phiz1(i) = sin(phi); 
   END DO 

   DO i = 1 , nxy
      phi = degstp*(i-1)*acos(-1.0)/180
      phix0(i) = cos(phi); 
      phiy0(i) = sin(phi); 
   END DO 
  
   DO i = 1 , nz
      DO j = 1 , nxy
      phix(i,j) = phix0(j) * phiz1(i)
      phiy(i,j) = phiy0(j) * phiz1(i)
      phiz(i,j) = phiz0(i)           
      END DO
   END DO

   DEALLOCATE(phix0,phiy0,phiz0,phiz1)

   !!! Load mantle density database made with PERPLE_X
   rho(1) = 3353d0
   IF(ptmod > 0) THEN

      !P-T domain grid
      tknum = 85 ; tkmin = 300; tkstp = 50
      pbnum = 1401 ; pbmin = 0d0; pbstp = 1d3
      ptnum = tknum*pbnum

      ALLOCATE(Rhodum(ptnum),Vpdum(ptnum),Vsdum(ptnum))

      CALL H5open_f (error)

      IF(eosmod == 1) CALL H5Fopen_f("../DATABASES/MMA-EoS/dunite.h5", H5F_ACC_RDONLY_F, file_id, error)
      IF(eosmod == 2) CALL H5Fopen_f("../DATABASES/MMA-EoS/hartzburgite.h5", H5F_ACC_RDONLY_F, file_id, error)
      IF(eosmod == 3) CALL H5Fopen_f("../DATABASES/MMA-EoS/pyrolite.h5", H5F_ACC_RDONLY_F, file_id, error)
      IF(eosmod == 4) CALL H5Fopen_f("../DATABASES/MMA-EoS/morb.h5", H5F_ACC_RDONLY_F, file_id, error)
      IF(eosmod == 5) CALL H5Fopen_f("../DATABASES/MMA-EoS/pyroxenite.h5", H5F_ACC_RDONLY_F, file_id, error)

      CALL H5Gopen_f(file_id, "/rock", group_id, error)

      CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Rhodum,'Rho',0)
      CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vpdum,'Vp',0)
      CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vsdum,'Vs',0)

      CALL H5Gclose_f(group_id, error)
      CALL H5Fclose_f(file_id, error)
      CALL H5close_f(error)

      ALLOCATE(td_rho(tknum,pbnum),td_vp(tknum,pbnum),td_vs(tknum,pbnum))

      DO j = 1 , tknum ; DO i = 1 , pbnum
 
         gi = j + (i-1)*tknum
         td_rho(j,i) = Rhodum(gi) !kg/m^3
         td_vp(j,i)  = Vpdum(gi)  !km/s
         td_vs(j,i)  = Vsdum(gi)  !km/s

      END DO ; END DO

      DEALLOCATE(Rhodum,Vpdum,Vsdum)

   END IF

   RETURN

   END SUBROUTINE init0
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Rotate tensor parallel to shortes axis of FSE  !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   SUBROUTINE tensorrot_aggr(m,acs1)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i,j,k,ll,p,q,r,ss,m
   DOUBLE PRECISION, DIMENSION (3,3,3,3) :: C0,Cav
   DOUBLE PRECISION, DIMENSION(3,3) :: acs1

   !Angles are commonly defined according to the right-hand rule. Namely, they
   !have positive values when they represent a rotation that appears clockwise
   !when looking in the positive direction of the rotating axis, and negative
   !values when
   !the rotation appears counter-clockwise. 

   Cav = 0d0 ; C0 = 0d0

   !Convert to Cijkl 4th order tensor
   DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
      C0(i,j,k,ll) = Sav(ijkl(i,j),ijkl(k,ll),m)
   END DO ; END DO ; END DO ; END DO

   !Rotate 4th order tensor
   DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
      DO p = 1 , 3 ; DO q = 1 , 3 ; DO r = 1 , 3 ; DO ss = 1 , 3
          Cav(i,j,k,ll) = Cav(i,j,k,ll) + acs1(i,p)*acs1(j,q)*acs1(k,r)*acs1(ll,ss)*C0(p,q,r,ss)
      END DO ; END DO ; END DO ; END DO
   END DO ; END DO ; END DO ; END DO

   !Convert to Voigt notation
   DO i = 1 , 6 ; DO j = 1 , 6
      Sav(i,j,m) = Cav(l1(i),l2(i),l1(j),l2(j))
   END DO ; END DO

   RETURN

   END SUBROUTINE tensorrot_aggr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine rotmatrix, build rotation matrix                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rotmatrixZXZ(phi1,theta,phi2,acs1)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: phi1,theta,phi2
   DOUBLE PRECISION, DIMENSION(3,3) :: acs1

   !Transform Euler angles into direction cosine matrix
   !Z-X-Z: Euler angles in radians
   acs1(1,1)=COS(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*SIN(phi2)
   acs1(2,1)=COS(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*SIN(phi2)
   acs1(3,1)=SIN(phi2)*SIN(theta)

   acs1(1,2)=-SIN(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*COS(phi2)
   acs1(2,2)=-SIN(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*COS(phi2)
   acs1(3,2)=COS(phi2)*SIN(theta)

   acs1(1,3)=SIN(theta)*SIN(phi1)
   acs1(2,3)=-SIN(theta)*COS(phi1)
   acs1(3,3)=COS(theta)

   RETURN

   END SUBROUTINE rotmatrixZXZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine rotmatrix, build rotation matrix                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rotmatrixZYZ(phi1,theta,phi2,acs2)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: phi1,theta,phi2
   DOUBLE PRECISION, DIMENSION(3,3) :: acs2

   !Transform Euler angles into direction cosine matrix
   !Z-Y-Z: Euler angles in radians
   acs2(1,1)=COS(phi2)*COS(phi1)*COS(theta)-SIN(phi1)*SIN(phi2)
   acs2(2,1)=COS(phi1)*SIN(phi2)+COS(theta)*COS(phi2)*SIN(phi1)
   acs2(3,1)=-COS(phi2)*SIN(theta)

   acs2(1,2)=-SIN(phi1)*COS(phi2)-COS(theta)*SIN(phi2)*COS(phi1)
   acs2(2,2)=COS(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*SIN(phi2)
   acs2(3,2)=SIN(phi2)*SIN(theta)

   acs2(1,3)=SIN(theta)*COS(phi1)
   acs2(2,3)=-SIN(theta)*SIN(phi1)
   acs2(3,3)=COS(theta)

   RETURN

   END SUBROUTINE rotmatrixZYZ

 
