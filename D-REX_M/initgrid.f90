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
!!! subroutine INIT0, initialization of variables, arrays and parameters   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE init0(input_dir,output_dir)

   USE comvar
   USE hdf5
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: i1,i2,i3,i,j,gi,j1,j2,j3,t,m,mrnum,mcolatnum,mlongnum,nx(3),ptmod0,ptnum,yy,marknum00 ! loop counters

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ran0
   ! matrix of random numbers used to generate initial random LPO

   DOUBLE PRECISION :: phi1,theta,phi2,phi,memnodes,memmark,dummy,lr,cr
   ! eulerian angles

   DOUBLE PRECISION :: mradiusmin,mcolatmin,mlongmin
   DOUBLE PRECISION :: mradiusmax,mcolatmax,mlongmax
   DOUBLE PRECISION :: colatlen,longlen,mR,mcolat,mlong
   DOUBLE PRECISION :: mcolatstp,mlongstp,mradiusstp
   DOUBLE PRECISION :: mcolatradiansstp,mlongradiansstp                      

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xe1,xe2,xe3
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Rhodum,Vpdum,Vsdum
   ! matrixes of initial random eulerian angles

   CHARACTER (len=500) :: input_dir,output_dir

   ! input and output directories
   CHARACTER(50) :: arg,inputfilename
   CHARACTER(500) :: str

   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl,memtype !Handles
   INTEGER     ::   error  ! Error flag

   !$omp parallel    
   nt = OMP_GET_NUM_THREADS()
   !$omp end parallel

   print *,' '
   write(*,'(a,i0)'),' Number of threads = ',nt
   print *,' '

   ALLOCATE(l(nt,3,3),e(nt,3,3),epsnot(nt))

   pi = 3.141592653589793238462643383279
   deg2rad = pi/180.0

   tau = 1d60

   !Read input file infos
   CALL read_input_file(input_dir,output_dir)

   !Number of files to be processed
   tnum=(Tend-Tinit)/Tstep+1

   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)
   write(*,'(a)'),' INITIALIZE EULERIAN GRID, AGGREGATES, ETC.'     
   write(*,*)
   
   !Scale coordinates 
   IF(cartspher>1) THEN
      !Convert longitude to radians
      x1min = x1min*deg2rad ; x1max = x1max*deg2rad 
      mx1min = mx1min*deg2rad  ; mx1max = mx1max*deg2rad 
      mx1minfab = mx1minfab*deg2rad  ; mx1maxfab = mx1maxfab*deg2rad 
      x3min = x3min*deg2rad  ; x3max = x3max*deg2rad 
      mx3min = mx3min*deg2rad  ; mx3max = mx3max*deg2rad 
      mx3minfab = mx3minfab*deg2rad  ; mx3maxfab = mx3maxfab*deg2rad 
   END IF

   !Check if domain is defined correctly
   IF(x1max < x1min) THEN ; print *,' Wrong Eulerian domain because x1max < x1min '; stop; END IF
   IF(x2max < x2min) THEN ; print *,' Wrong Eulerian domain because x2max < x2min '; stop; END IF
   IF(x3max < x3min) THEN ; print *,' Wrong Eulerian domain because x3max < x3min '; stop; END IF
   IF(mx1max < mx1min) THEN ; print *,' Wrong Lagrangian domain because mx1max < mx1min '; stop; END IF
   IF(mx2max < mx2min) THEN ; print *,' Wrong Lagrangian domain because mx2max < mx2min '; stop; END IF
   IF(mx3max < mx3min) THEN ; print *,' Wrong Lagrangian domain because mx3max < mx3min '; stop; END IF
   IF(mx1min < x1min) THEN ; print *,' Markers defined outside the domain because mx1min < x1min '; stop; END IF
   IF(mx2min < x2min) THEN ; print *,' Markers defined outside the domain because mx2min < x2min '; stop; END IF
   IF(mx3min < x3min) THEN ; print *,' Markers defined outside the domain because mx3min < x3min '; stop; END IF
   IF(mx1max > x1max) THEN ; print *,' Markers defined outside the domain because mx1max > x1max '; stop; END IF
   IF(mx2max > x2max) THEN ; print *,' Markers defined outside the domain because mx2max > x2max '; stop; END IF
   IF(mx3max > x3max) THEN ; print *,' Markers defined outside the domain because mx3max > x3max '; stop; END IF

   IF(cartspher==2) THEN
      IF(x1max <  0d0) THEN; print *,'Wrong Eulerian domain because Longitude  <    0'; stop; END IF
      IF(x1max > 2*pi) THEN; print *,'Wrong Eulerian domain because Longitude  > 2*pi'; stop; END IF
      IF(mx1max <  0d0) THEN; print *,'Wrong Lagrangian domain because Longitude  <    0'; stop; END IF
      IF(mx1max > 2*pi) THEN; print *,'Wrong Lagrangian domain because Longitude  > 2*pi'; stop; END IF
      IF(x3max <  0d0) THEN; print *,'Wrong Eulerian domain because Colatitude  <    0'; stop; END IF
      IF(x3max > pi) THEN; print *,'Wrong Eulerian domain because Colatitude  > pi'; stop; END IF
      IF(mx3max <  0d0) THEN; print *,'Wrong Lagrangian domain because Colatitude  <    0'; stop; END IF
      IF(mx3max > pi) THEN; print *,'Wrong Lagrangian domain because Colatitude  > pi'; stop; END IF
   END IF

   IF(fossilfabric > 0) THEN
      IF(mx1maxfab < mx1minfab) THEN ; print *,' Wrong fossil fabric domain because mx1maxfab < mx1minfab '; stop; END IF
      IF(mx2maxfab < mx2minfab) THEN ; print *,' Wrong fossil fabric domain because mx2maxfab < mx2minfab '; stop; END IF
      IF(mx3maxfab < mx3minfab) THEN ; print *,' Wrong fossil fabric domain because mx3maxfab < mx3minfab '; stop; END IF
      IF(mx1minfab < x1min) THEN ; print *,' Fossil fabric defined outside the domain because mx1minfab < x1min '; stop; END IF
      IF(mx2minfab < x2min) THEN ; print *,' Fossil fabric defined outside the domain because mx2minfab < x2min '; stop; END IF
      IF(mx3minfab < x3min) THEN ; print *,' Fossil fabric defined outside the domain because mx3minfab < x3min '; stop; END IF
      IF(mx1maxfab > x1max) THEN ; print *,' Fossil fabric defined outside the domain because mx1maxfab > x1max '; stop; END IF
      IF(mx2maxfab > x2max) THEN ; print *,' Fossil fabric defined outside the domain because mx2maxfab > x2max '; stop; END IF
      IF(mx3maxfab > x3max) THEN ; print *,' Fossil fabric defined outside the domain because mx3maxfab > x3max '; stop; END IF
   END IF

   !Check if LPO parameters are defined correctly
   DO i=1,4

      IF(Xol(i) < 0 .OR. Xol(i) > 100d0) THEN ; print *,'Wrong volume fraction of main phase for rocktype ',i,', it should be >= 0 and <= 1 '; stop; END IF 
      IF(chi(i) < 0 .OR. chi(i) > 1d0) THEN ; print *,'Wrong GBS volume threshold for rocktype ',i,', it should be >= 0 and <= 1 '; stop; END IF   
      IF(fractdislrock(i) < 0 .OR. fractdislrock(i) > 1d0) THEN ; print *,'Wrong fraction of deformation accommodated by dislocation creep for rocktype ',i,', it should be >= 0 and <= 1 '; stop; END IF   
      IF(single_crystal_elastic_db(i,1) < 1 .OR.  single_crystal_elastic_db(i,1) > 20) THEN ; print *,'Wrong single crystal elastic tensor for main phase for rocktype ',i,', it should be > 0 and < 20'; stop; END IF   
      IF(single_crystal_elastic_db(i,2) < 1 .OR.  single_crystal_elastic_db(i,2) > 20) THEN ; print *,'Wrong single crystal elastic tensor for minor phase for rocktype ',i,', it should be > 0 and < 20'; stop; END IF   
      IF(stressexp(i) < 0) THEN ; print *,'Negative power-law exponent for rocktype ',i,', it should be >= 1'; stop; END IF   
      IF(Mob(i) < 0)       THEN ; print *,'Negative Mob parameter for rocktype ',i,', it should be > 0'; stop; END IF   
      IF(lambda(i) < 0)    THEN ; print *,'Negative lambda parameter for rocktype ',i,', it should be > 0'; stop; END IF    
      DO j=1,12
         IF(tau(i,j) < 1)  THEN ; print *,'nCRSS < 1 for rocktype ',i, ', it should be >= 1'; stop; END IF 
      END DO

   END DO
  
   !Create ouput directory
   str='mkdir '//trim(output_dir)
   CALL SYSTEM(str)

!!! initial size
   size = size3**3

!!! Xol given as a percentage
   Xol = Xol/1d2

!!! Initialization of the Eulerian fields: velocity, velocity gradient, FSE, P, T, etc.
   yy = yinyang
   nodenum = nx1*nx2*nx3

   ALLOCATE(X1(nx1),X2(nx2),X3(nx3),Ui(yy,3,nx1,nx2,nx3),Ux(yy,3,nx1,nx2,nx3),Dij(yy,3,3,nx1,nx2,nx3))
   IF(ptmod > 0) ALLOCATE(Tk(yy,nx1,nx2,nx3),Pa(yy,nx1,nx2,nx3))
   IF(fractdislmod > 0) ALLOCATE(Fd(yy,nx1,nx2,nx3))

   x1stp = (x1max-x1min)/(nx1-1)
   x2stp = (x2max-x2min)/(nx2-1)
   x3stp = (x3max-x3min)/(nx3-1)

   !Define Eulerian grid coordinates
   !X
   X1(1) = x1min
   DO i1 = 2, nx1
      X1(i1) = X1(i1-1) + x1stp     
   END DO   
   
   !Y
   X2(1) = x2min
   DO i2 = 2, nx2
      X2(i2) = X2(i2-1) + x2stp   
   END DO  
   
   !Z
   X3(1) = x3min
   DO i3 = 2, nx3
      X3(i3) = X3(i3-1) + x3stp   
   END DO  
   
!!! Initialization of the Lagrangian grid with crystal aggregates

   IF(cartspher==1) THEN
  
   mnx1=NINT((mx1max-mx1min)/mx1stp)
   mnx2=NINT((mx2max-mx2min)/mx2stp)
   mnx3=NINT((mx3max-mx3min)/mx3stp)
   IF(dimensions == 2) mnx3 = 1

   marknum=mnx1*mnx2*mnx3

   ELSE

   !Define quantities appropriate for polar coordinates 
   
   !X is Colatitude
   mlongmin = mx1min !Radians 
   mlongmax = mx1max !Radians
   !Y is Longitude
   mradiusmin = mx2min !Distance
   mradiusmax = mx2max !Distance
   !Z is Colatitude
   mcolatmin = mx3min !Radians 
   mcolatmax = mx3max !Radians

   mlongstp   = mx1stp !Distance
   mradiusstp = mx2stp !Distance 
   mcolatstp  = mx3stp !Distance
   
   !Calculate number of aggregates for array/matrix memory allocation
   marknum=0
   
   !Number of vertical layers
   mrnum=floor((mradiusmax-mradiusmin)/mradiusstp)

   IF(dimensions == 2) THEN

   DO i2 = 1, mrnum
      !LAYER i2
      !Find radius
      mR = mradiusmin + mradiusstp*(i2-0.5)
      !Find lenght of latitude arc (great circle arc)
      longlen = (mlongmax-mlongmin)*mR
      !Number of aggregates along latitude arc
      mlongnum = floor(longlen/mlongstp)
      !Increase number of markers
      marknum = marknum + mlongnum
   END DO

   END IF

   IF(dimensions == 3) THEN
   
   DO i2 = 1, mrnum
      !LAYER i2
      !Find radius
      mR = mradiusmin + mradiusstp*(i2-0.5)
      !Find lenght of latitude arc (great circle arc)
      colatlen = (mcolatmax-mcolatmin)*mR
      !Number of aggregates along latitude arc
      mcolatnum = floor(colatlen/mcolatstp)
      !Latitude step in radians
      mcolatradiansstp = (mcolatmax-mcolatmin)/REAL(mcolatnum)
      DO i3 = 1, mcolatnum
         !Colatitude in radians
         mcolat = mcolatmin + mcolatradiansstp*(i3-0.5)
         !Find lenght of longitude arc 
         longlen = (mlongmax-mlongmin)*mR*sin(mcolat)
         !Number of markers along longitude    
         mlongnum = floor(longlen/mlongstp)
         !Increase number of markers
         marknum = marknum + mlongnum
      END DO
   END DO

   !Double number of markers when yinyang grids are active
   IF(yinyang == 2) marknum = marknum*2

   END IF

   END IF

 
   write(*,*)
   write(*,'(a,i0,a,i0)'),' Number of aggregates in the domain = ',marknum
   write(*,*)

   memnodes=nx1+nx2+nx3+nx1*nx2*nx3*27*yy
   write(*,'(a,1f6.2,a)'),' Memory needed for nodes is ',memnodes*8/1d9,' GB'
   print*
  
   memmark = marknum*(5 + 9)

   !Account for Sav, P, T
   memmark = memmark + marknum*(36+2)

   dummy=marknum*size
   dummy=dummy*(9 + 9 + 1 + 1) 
   memmark = memmark + dummy           

   write(*,'(a,1f6.2,a)'),' Memory needed for aggregates is ',((memmark-1)*8+1)/1d9,' GB'
   print*
   write(*,'(a,1f6.2,a)'),' Total memory needed is about ',((memnodes+memmark)*8/1d9)*1.05d0, ' GB' !Add 5% to account for other variables declared
   print*

   ALLOCATE(mx1(marknum),mx2(marknum),mx3(marknum),mYY(marknum),rocktype(marknum),rho(marknum))
   ALLOCATE(Fij(3,3,marknum))

   mx1 = 0d0 ; mx2 = 0d0 ; mx3 = 0d0 ; mYY = 1 ; rocktype = 1 ; rho = 0d0

!!! Initial deformation gradient tensor
   Fij = 0d0 ; Fij(1,1,:) = 1d0 ; Fij(2,2,:) = 1d0 ; Fij(3,3,:) = 1d0

   m = 0
  
   IF(cartspher == 1) THEN

   DO i1 = 1, mnx1
      DO i2 = 1, mnx2
         DO i3 = 1, mnx3
            !Marker Index
            m = m + 1
            mx1(m)=mx1min+mx1stp*(i1 - 0.5d0)
            mx2(m)=mx2min+mx2stp*(i2 - 0.5d0)
            mx3(m)=mx3min+mx3stp*(i3 - 0.5d0)
         END DO
      END DO
   END DO

   !Set Z coordinate equal to 0
   IF(dimensions == 2) mx3 = 0d0

   ELSE

   IF(dimensions == 2) THEN

   DO i2 = 1, mrnum
      !LAYER i2
      !Find radius
      mR = mradiusmin + mradiusstp*(i2-0.5)
      !Find lenght of latitude arc (great circle arc)
      longlen = (mlongmax-mlongmin)*mR
      !Number of aggregates along latitude arc
      mlongnum = floor(longlen/mlongstp)
      !Latitude step in radians
      mlongradiansstp = (mlongmax-mlongmin)/REAL(mlongnum)
      DO i1 = 1, mlongnum
         !Colatitude in radians
         mlong = mlongmin + mlongradiansstp*(i1-0.5)
         m = m + 1
         !Polar coordinates
         !Longitude
         mx1(m) = mlong  ! mR * cos(mlong)
         !Radial
         mx2(m) = mR     ! mR * sin(mlong)
      END DO
   END DO

   !Set colatitude on the equatorial plane
   mx3 = pi/2.0

   END IF
   
   IF(dimensions == 3) THEN

   DO i2 = 1, mrnum
      !LAYER i3
      !Find radius
      mR = mradiusmin + mradiusstp*(i2-0.5)
      !Find lenght of latitude arc (great circle arc)
      colatlen = (mcolatmax-mcolatmin)*mR
      !Number of aggregates along latitude arc
      mcolatnum = floor(colatlen/mcolatstp)
      !Latitude step in radians
      mcolatradiansstp = (mcolatmax-mcolatmin)/REAL(mcolatnum)
      DO i3 = 1, mcolatnum
         !Colatitude in radians
         mcolat = mcolatmin + mcolatradiansstp*(i3-0.5)
         !Find lenght of longitude arc 
         longlen = (mlongmax-mlongmin)*mR*sin(mcolat)
         !Number of aggregates along longitude   
         mlongnum = floor(longlen/mlongstp)
         !Longitude step in radians            
         mlongradiansstp = (mlongmax-mlongmin)/REAL(mlongnum)
         DO i1 = 1, mlongnum
            !Longitude in radians
            mlong = mlongmin + mlongradiansstp*(i1-0.5)
            m = m + 1
            !Spherical coordinates
            !Longitude
            mx1(m) = mlong  !mR * sin(mcolat) * cos(mlong)
            !Radial
            mx2(m) = mR     !mR * sin(mcolat) * sin(mlong)
            !Colatitude
            mx3(m) = mcolat !mR * cos(mcolat)
         END DO
      END DO
   END DO

   !Add markers for Yang grid
   IF(yy == 2) THEN

     marknum00 = marknum/2
     DO i1=marknum00+1,marknum
        mx1(i1) = mx1(i1-marknum00)
        mx2(i1) = mx2(i1-marknum00)
        mx3(i1) = mx3(i1-marknum00)
        mYY(i1) = 2                
     END DO   

   END IF
   
   END IF
   
   END IF
   
   IF(fsemod == 0) THEN

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

!!! scaling factor for elastic tensor inversion

   mandel_scale = 1d0

   DO j1=1,6 ; DO j2=1,6
      IF(j1 > 3 .AND. j2 > 3) THEN
         mandel_scale(j1,j2) = 2d0
      ELSE IF(j1 > 3 .OR. j2 > 3) THEN
         mandel_scale(j1,j2) = 2**0.5d0
      END IF
   END DO ; END DO

!!! Loading stiffness tensors (GPa)

   CALL elastic_database(S0,dS0dp,dS0dp2,dS0dt,dS0dpt,dS0dp5,dS0dt5)

!!! allocation of the dimensions of the arrays

   ALLOCATE(xe1(size),xe2(size),xe3(size))

   ALLOCATE(odf(marknum,size))
   ALLOCATE(odf_ens(marknum,size))

   ALLOCATE(ran0(3*size))

   ALLOCATE(acs(size,3,3,marknum),acs0(size,3,3))
   ALLOCATE(acs_ens(size,3,3,marknum))

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

   DO i = 1 , marknum
   
      odf(i,:) = 1d0/REAL(size3**3)
      odf_ens(i,:) = odf(i,:)
      acs(:,:,:,i) = acs0
      acs_ens(:,:,:,i) = acs0
   
   END DO

   DEALLOCATE(xe1,xe2,xe3,ran0)

!!! Pre-existing mantle fabric

!!! Load mantle density database made with PERPLE_X
   rho = 3353d0
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

!!! Set initial rocktype
   ptmod0 = ptmod
   ptmod = 1
   DO m = 1, marknum

      CALL rocktypecheck(m)

   END DO
   ptmod = ptmod0

   END IF

   !Delete markers of YANG gridi overlapping with YIN grid
   IF(yy == 2) THEN

     DO m=1,marknum
        !IF(mYY(m) == 2 .AND. (mx1(m)<-pi/2.0 .OR. mx1(m)>pi/2.0 .OR. mx3(m)<0.25*pi .OR. mx3(m)>0.75d0*pi)) THEN
        IF(mYY(m) == 2) THEN
           CALL yin2yang(mx1(m),mx3(m),lr,cr)
           IF(lr>-pi/2.0 .AND. lr<pi/2.0 .AND. cr>0.25*pi .AND. cr<0.75*pi) rocktype(m)= rocktype(m) + 100
        END IF
     END DO   

   END IF

   write(*,"(a)"),' INITIALIZE EULERIAN GRID, AGGREGATES, ETC. OK!'
   write(*,*)
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)

   RETURN

   1000 FORMAT(101(1pe14.6))

   END SUBROUTINE init0
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE set_fossil_fabric 

   USE comvar
   USE omp_lib
   USE hdf5    

   IMPLICIT NONE

   INTEGER :: m,nx(3)
   DOUBLE PRECISION :: lr,rr,cr
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: odf0,odf_ens0
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: acs00,acs_ens00
   DOUBLE PRECISION, DIMENSION(3,3) :: fse0
   ! matrixes of pre-existing crystal volume fraction, orientation and FSE

   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl, memtype ! Handles
   INTEGER     ::   error  ! Error flag

   print *,'Impose preexisting fabric'
   print *

   ALLOCATE(acs00(size,3,3),acs_ens00(size,3,3),odf0(size),odf_ens0(size))

   !!! Read infos from hdf5 file format

   !Initialize FORTRAN interface.

   CALL H5open_f (error)

   !  Open file using default properties.

   CALL H5Fopen_f('../D-REX_S/fossilfabric.h5', H5F_ACC_RDONLY_F, file_id, error)

   nx(1)=size
   nx(2)=3
   nx(3)=3
   CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs00,'acs',0)
   CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,odf0,'odf',0)
   CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,acs_ens00,'acs_ens',0)
   CALL loadsave_double(0,1,file_id,nx(1:1),H5T_NATIVE_DOUBLE,odf_ens0,'odf_ens',0)
   nx(1)=3
   nx(2)=3
   CALL loadsave_double(0,2,file_id,nx(1:2),H5T_NATIVE_DOUBLE,fse0,'Fij',0)

   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   !Close FORTRAN interface.
   CALL H5close_f(error)

   DO m = 1 , marknum

      lr = mx1(m) ; rr = mx2(m) ; cr = mx3(m)
      IF(mYY(m) == 2) CALL yin2yang(mx1(m),mx3(m),lr,cr)

      IF(lr >= mx1minfab .AND. lr <= mx1maxfab) THEN
      IF(rr >= mx2minfab .AND. rr <= mx2maxfab) THEN
      IF(cr >= mx3minfab .AND. cr <= mx3maxfab) THEN

         IF(fsemod == 0) THEN

         odf(m,:) = odf0
         acs(:,:,:,m) = acs00
         IF(rocktype(m) == 1) THEN
            acs_ens(:,:,:,m) = acs_ens00
            odf_ens(m,:) = odf_ens0
         ELSE
            acs_ens(:,:,:,m) = acs0
            odf_ens(m,:) = 1d0/REAL(size3**3)
         END IF

         END IF

         Fij(:,:,m) = fse0

         IF(cartspher == 2) CALL rottensor(m)

      END IF; END IF; END IF

   END DO

   DEALLOCATE(acs00,acs_ens00,odf0,odf_ens0)

   END SUBROUTINE set_fossil_fabric

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rottensor(m)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: m,i,j,k,p,q 

   DOUBLE PRECISION :: phi1,theta,phi2,angle0,angle01,lr,cr
   DOUBLE PRECISION, DIMENSION(3) :: a,b
   DOUBLE PRECISION, DIMENSION(3,3) :: acs1,acs2,acs_ens2,fse1
   ! eulerian angles

   !In D-REX_S, the vertical direction is axis 2, which coincide with the
   !Y axis in cartesian coordinates, so remove 90° from colatitude 
   phi1 = mx1(m) - pi*3.0/2.0 ; theta = mx3(m) - pi/2.0 ; phi2 = 0d0
   IF(mYY(m) == 2) THEN
      CALL yin2yang(mx1(m),mx3(m),lr,cr)
      phi1 = lr - pi*3.0/2.0 ; theta = cr - pi/2.0
   END IF

   !Rotate fse with a single rotation 
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

   fse1(:,:) = Fij(:,:,m)
   Fij(:,:,m) = 0d0

   DO j=1,3; DO k=1,3; DO p=1,3; DO q=1,3
      Fij(j,k,m) = Fij(j,k,m) + acs1(j,p)*acs1(k,q)*fse1(p,q)
   END DO; END DO; END DO; END DO

   !Rotate cosine direction matrix with 2 separate rotations
   IF(fsemod == 0) THEN

   !the minus sign is because we need counter-clockwise rotation
   angle0=-phi1/deg2rad
   angle01=-theta/deg2rad
   !First axis of rotation: z-axis
   a=0d0
   a(3)=1;
   !Second axis of rotation: the new x-axis
   acs2=0d0 
   acs2(1,1)=1
   acs2(2,2)=1
   acs2(3,3)=1
   CALL rot3D(acs2,acs1,a,-angle0)
   b(:)=acs1(:,1)
  
   DO i=1,size

      ! First rotate around the z-axis by long
      acs2(:,:) = acs(i,:,:,m)
      acs_ens2(:,:) = acs_ens(i,:,:,m)
      CALL rot3D(acs2,acs(i,:,:,m),a,angle0)
      CALL rot3D(acs_ens2,acs_ens(i,:,:,m),a,angle0)

      ! Second rotate around the new x-axis by colat
      acs2(:,:) = acs(i,:,:,m)
      acs_ens2(:,:) = acs_ens(i,:,:,m)
      CALL rot3D(acs2,acs(i,:,:,m),b,angle01)
      CALL rot3D(acs_ens2,acs_ens(i,:,:,m),b,angle01)
   
   END DO

   END IF

   END SUBROUTINE rottensor

