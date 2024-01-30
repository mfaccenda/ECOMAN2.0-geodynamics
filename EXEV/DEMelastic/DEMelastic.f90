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
!!! Module of common variables                                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   MODULE comvardem

   INTEGER :: phinum,thetanum,volnum,esanum(3)
   DOUBLE PRECISION :: pi,stp,Vstp
   DOUBLE PRECISION, DIMENSION(6,6) :: Cback,Cinc,mandel_scale
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: Is,Id
   INTEGER, DIMENSION(3,3) :: ijkl ! tensor of indices to form Cijkl from Sij
   DOUBLE PRECISION, DIMENSION(3,3) :: kron ! tensor of indices to form Cijkl from Sij
   INTEGER, DIMENSION(6) :: l1,l2 ! tensot of indices to form Sij from Cijkl
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: thetaspo,phispo,esa1,esa2,esa3,Vol,ro
   ! input file name before file number
   CHARACTER (300) :: output_filename

   END MODULE comvardem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Read input file spo_input.dat and allocate matrices    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE initspo

   USE comvardem
   USE omp_lib
   
   INTEGER :: ii,jj,kk,ll,j1,j2,rt,nt,axislogscalemod
   DOUBLE PRECISION :: Vmax,ro_back,ro_incl
   DOUBLE PRECISION :: esamin(3),esamax(3),esastp(3)

!!! Read input values for SPO model       
   OPEN(15,file='dem_input.dat')

   !Output filename                
   !call read_par_char(15,output_filename)
   output_filename = 'elastictensordem.h5'

   !Constant Cij
   call read_par_double6(15,Cback(1,1:6))
   call read_par_double6(15,Cback(2,1:6))
   call read_par_double6(15,Cback(3,1:6))
   call read_par_double6(15,Cback(4,1:6))
   call read_par_double6(15,Cback(5,1:6))
   call read_par_double6(15,Cback(6,1:6))
   call read_par_double(15,ro_back)
   call read_par_double6(15,Cinc(1,1:6))
   call read_par_double6(15,Cinc(2,1:6))
   call read_par_double6(15,Cinc(3,1:6))
   call read_par_double6(15,Cinc(4,1:6))
   call read_par_double6(15,Cinc(5,1:6))
   call read_par_double6(15,Cinc(6,1:6))
   call read_par_double(15,ro_incl)

   !Inclusions shape
   call read_par_int(15,axislogscalemod)
   call read_par_double3(15,esamin)  !!Ellipsoid min semiaxes
   call read_par_double3(15,esamax)  !!Ellipsoid max semiaxes
   call read_par_double3(15,esastp)  !!Ellipsoid semiaxes increment

   !DEM volume increment
   call read_par_double(15,Vstp)

   !Maximum inclusions volme fraction
   call read_par_double(15,Vmax)

   !!! Make semiaxes arrays
   esanum = NINT((esamax - esamin)/esastp) + 1
   ALLOCATE(esa1(esanum(1)),esa2(esanum(2)),esa3(esanum(3)))
   DO ii=1,esanum(1)
      esa1(ii) = esamin(1) + esastp(1)*(ii-1)
      IF(axislogscalemod > 0) esa1(ii) = 10d0**esa1(ii)
   END DO
   DO ii=1,esanum(2)
      esa2(ii) = esamin(2) + esastp(2)*(ii-1)
      IF(axislogscalemod > 0) esa2(ii) = 10d0**esa2(ii)
   END DO
   DO ii=1,esanum(3)
      esa3(ii) = esamin(3) + esastp(3)*(ii-1)
      IF(axislogscalemod > 0) esa3(ii) = 10d0**esa3(ii)
   END DO

   print *,'Major  semiaxis: ',esa1
   print *,'Medium semiaxis: ',esa2
   print *,'Minor  semiaxis: ',esa3
   print *

   ! Volume increment array

   volnum = NINT(Vmax/Vstp)+1
   ALLOCATE(Vol(volnum))
   ALLOCATE(ro(volnum))
   DO ii = 1,volnum
      Vol(ii) = Vstp*(ii-1)
      ro(ii) = ro_back*(1d0-Vol(ii)) + ro_incl*Vol(ii)
   END DO

   print *,'Computing volume steps',Vol     
   print *
   print *,'Two-phase medium density',ro_back,ro_incl,ro    
   print *

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
         mandel_scale(j1,j2) = 2**0.5
      END IF
   END DO ; END DO

   Is = 0d0 ; kron = 0d0

   !Kronecker delta
   kron(1,1)=1;
   kron(2,2)=1;
   kron(3,3)=1;
   
   !Symmetricic symmetric fourth-rank unit tensor
   DO ii = 1,3; DO jj = 1,3; DO kk = 1,3; DO ll=1,3
      Is(ii,jj,kk,ll) = Is(ii,jj,kk,ll) + 0.5*(kron(ii,kk)*kron(jj,ll)+kron(ii,ll)*kron(jj,kk));
   END DO; END DO; END DO; END DO

   !Make spherical grid
   pi = 3.141592653589793238462643383279
   stp=pi/100
   thetanum = INT(pi/stp)+1
   ALLOCATE(thetaspo(thetanum))
   DO ii = 1, thetanum
      thetaspo(ii) = stp*(ii)
   END DO

   phinum = INT(2*pi/stp)+1
   ALLOCATE(phispo(phinum))
   DO ii = 1, phinum
      phispo(ii) = stp*(ii)
   END DO

   END SUBROUTINE initspo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  SPO calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   PROGRAM DEMelastic

   USE comvardem
   USE omp_lib
   USE hdf5

   IMPLICIT NONE

   INTEGER :: ii,jj,kk,ll,m,pp,qq,rr,ss,tt,uu,cyc,vcyc,nu,tid,nt
   INTEGER :: j1,j2,cyc_esa1,cyc_esa2,cyc_esa,gi,nx(3)
   DOUBLE PRECISION :: N,L,E_chi,Gc,Gs,G,phiazi,azimuthal(3),esa(3),dum_db(3),wtime
   DOUBLE PRECISION :: vol_fract
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: esa12,esa23
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: C11,C22,C33,C44,C55,C66,C12,C13,C23
   DOUBLE PRECISION, DIMENSION(6,6) :: Cdem
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: c0dem,Cback3,Cinc3,Ctemp,dC1,dC2,dC3,dC4
   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl,memtype ! Handles
   INTEGER     ::   error  ! Error flag

   !$omp parallel   
   nt = OMP_GET_NUM_THREADS()
   !$omp end parallel   
   print *,' '
   write(*,'(a,i0)'),' Number of threads = ',nt
   print *,' '

   wtime = OMP_GET_WTIME()

   !!! Differential effective medium (Mainprice, 1997, EPSL)              
   CALL initspo

   !Convert from Voigt notation to 4th order tensor
   DO ii = 1,3; DO jj = 1,3; DO kk = 1,3; DO ll=1,3
      Cback3(ii,jj,kk,ll) = Cback(ijkl(ii,jj),ijkl(kk,ll))
      Cinc3(ii,jj,kk,ll)  = Cinc(ijkl(ii,jj),ijkl(kk,ll))
   END DO; END DO; END DO; END DO

   !Multiple ellipsoids
   if(esanum(1)*esanum(2) > 1) THEN

   !Allocate memory for output 3D arrays
   ALLOCATE(C11(volnum,esanum(1),esanum(2)))
   ALLOCATE(C22(volnum,esanum(1),esanum(2)))
   ALLOCATE(C33(volnum,esanum(1),esanum(2)))
   ALLOCATE(C44(volnum,esanum(1),esanum(2)))
   ALLOCATE(C55(volnum,esanum(1),esanum(2)))
   ALLOCATE(C66(volnum,esanum(1),esanum(2)))
   ALLOCATE(C12(volnum,esanum(1),esanum(2)))
   ALLOCATE(C13(volnum,esanum(1),esanum(2)))
   ALLOCATE(C23(volnum,esanum(1),esanum(2)))
   ALLOCATE(esa12(esanum(1),esanum(2)))
   ALLOCATE(esa23(esanum(1),esanum(2)))

   !Inclusions aspect ratios
   DO cyc_esa2 = 1, esanum(2)
   DO cyc_esa1 = 1, esanum(1)

     cyc_esa = cyc_esa1 - cyc_esa2 + 1

     esa(1)=esa1(cyc_esa1)
     esa(2)=esa2(cyc_esa2)
     esa(3)=1
     IF(esa(2) <= esa(1)) THEN
        esa12(cyc_esa1,cyc_esa2)= log10(esa(1)/esa(2))
        esa23(cyc_esa1,cyc_esa2)= log10(esa(2)/esa(3))
     END IF

   END DO; END DO


   !$omp parallel & 
   !$omp shared(esa1,esa2,esa3,Vol,C11,C22,C33,C44,C55,C66,C12,C13,C23,esa12,esa23) &
   !$omp private(tid,esa,cyc_esa1,cyc_esa2,cyc_esa,cyc,vol_fract,c0dem,Ctemp,Cdem,dC1,dC2,dC3,dC4,ii,jj) &
   !$omp firstprivate(esanum,volnum,nu,Cback3,Cinc3,Vstp,l1,l2)   
   !$omp do schedule(static) collapse(2) 
   DO cyc_esa2 = 1, esanum(2)
   DO cyc_esa1 = 1, esanum(1)

   tid = OMP_GET_THREAD_NUM()

   cyc_esa = cyc_esa1 - cyc_esa2 + 1

   esa(1)=esa1(cyc_esa1)
   esa(2)=esa2(cyc_esa2)
   esa(3)=1

   !Skip ellipsoids when a2 > a1
   IF(esa(2) > esa(1)) GOTO 10

   DO cyc=1,volnum

      IF(Vol(cyc) .EQ. 0d0) THEN

         c0dem = Cback3

      ELSE IF(Vol(cyc) .EQ. 1d0) THEN

         c0dem = Cinc3

      ELSE

         vol_fract = log(1d0+Vstp/(1d0-Vol(cyc))) 

         !Calculate symmetrical tensor Green function
         !4th order Runge-Kutta integration scheme    
         Ctemp = c0dem 
         CALL greenfunction(Ctemp,Cinc3,vol_fract,esa,dC1)
         Ctemp = c0dem + 0.5d0*dC1
         CALL greenfunction(Ctemp,Cinc3,vol_fract,esa,dC2)
         Ctemp = c0dem + 0.5d0*dC2
         CALL greenfunction(Ctemp,Cinc3,vol_fract,esa,dC3)
         Ctemp = c0dem + dC3
         CALL greenfunction(Ctemp,Cinc3,vol_fract,esa,dC4)
         
         c0dem = c0dem + (dC1 + 2d0*dC2 + 2d0*dC3 + dC4)/6d0
 
      END IF
        
      !Convert to Voigt notation
      DO jj = 1, 6 ; DO ii = 1 , 6
         Cdem(ii,jj) = c0dem(l1(ii),l2(ii),l1(jj),l2(jj))
      END DO ; END DO

      !Compute global index in the output matrix
      C11(cyc,cyc_esa1,cyc_esa2)= Cdem(1,1)
      C22(cyc,cyc_esa1,cyc_esa2)= Cdem(2,2)
      C33(cyc,cyc_esa1,cyc_esa2)= Cdem(3,3)
      C44(cyc,cyc_esa1,cyc_esa2)= Cdem(4,4)
      C55(cyc,cyc_esa1,cyc_esa2)= Cdem(5,5)
      C66(cyc,cyc_esa1,cyc_esa2)= Cdem(6,6)
      C12(cyc,cyc_esa1,cyc_esa2)= Cdem(1,2)
      C13(cyc,cyc_esa1,cyc_esa2)= Cdem(1,3)
      C23(cyc,cyc_esa1,cyc_esa2)= Cdem(2,3)

   END DO 

10 END DO; END DO
   !$omp end do
   !$omp end parallel

   !!! Write infos in hdf5 format

   !Initialize FORTRAN interface.

   CALL H5open_f (error)

   ! Create a new file using default properties.

   CALL H5Fcreate_f(output_filename, H5F_ACC_TRUNC_F, file_id, error)

   !Save output matrix
   nx(1)=volnum
   nx(2)=esanum(1)
   nx(3)=esanum(2)
   CALL loadsave_integer(0,1,file_id,3,H5T_NATIVE_INTEGER,nx,'gridnum',1)
   CALL loadsave_double(0,1,file_id,nx(1),H5T_NATIVE_DOUBLE,Vol,'Vol',1)
   CALL loadsave_double(0,1,file_id,nx(1),H5T_NATIVE_DOUBLE,ro,'ro',1)
   CALL loadsave_double(0,2,file_id,nx(2:3),H5T_NATIVE_DOUBLE,esa12,'esa12',1)
   CALL loadsave_double(0,2,file_id,nx(2:3),H5T_NATIVE_DOUBLE,esa23,'esa23',1)
   CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,C11,'C11',1)
   CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,C22,'C22',1)
   CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,C33,'C33',1)
   CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,C44,'C44',1)
   CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,C55,'C55',1)
   CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,C66,'C66',1)
   CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,C12,'C12',1)
   CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,C13,'C13',1)
   CALL loadsave_double(0,3,file_id,nx(1:3),H5T_NATIVE_DOUBLE,C23,'C23',1)

   CALL H5Fclose_f(file_id, error)

   CALL H5close_f(error)
   !Close FORTRAN interface.

   DEALLOCATE(Vol,C11,C22,C33,C44,C55,C66,C12,C13,C23,esa12,esa23)

   !Single ellipsoid shape
   ELSE

   esa(1)=esa1(1)
   esa(2)=esa2(1)
   esa(3)=1

   print *,''
   print *,'Matrix elastic tensor:'
   write(*,'(6f12.5)') Cback
   print *,''
   print *,'Inclusion elastic tensor:'
   write(*,'(6f12.5)') Cinc
   print *,''
   write(*,'(a,3f14.5)') 'Inclusion semiaxes:',esa
   print *,''
   write(*,'(a)'),'      Cxxxx       Cyyyy       Czzzz       Cyzyz       Cxzxz        Cxyxy       Cxxyy       Cxxzz       Cyyzz      Vol fract'

   DO cyc=1,volnum
      IF(Vol(cyc) .EQ. 0d0) THEN
         c0dem = Cback3
      ELSE IF(Vol(cyc) .EQ. 1d0) THEN
         c0dem = Cinc3
      ELSE
         !vol_fract = Vstp/(1d0-Vol(cyc)) 
         vol_fract = log(1d0+Vstp/(1d0-Vol(cyc)))
         !Calculate symmetrical tensor Green function
         !4th order Runge-Kutta integration scheme    
         Ctemp = c0dem
         !write(*,'(9f)') Ctemp
         CALL greenfunction(Ctemp,Cinc3,vol_fract,esa,dC1)
         Ctemp = c0dem + 0.5d0*dC1
         CALL greenfunction(Ctemp,Cinc3,vol_fract,esa,dC2)
         Ctemp = c0dem + 0.5d0*dC2
         CALL greenfunction(Ctemp,Cinc3,vol_fract,esa,dC3)
         Ctemp = c0dem + dC3
         CALL greenfunction(Ctemp,Cinc3,vol_fract,esa,dC4)

         c0dem = c0dem + (dC1 + 2d0*dC2 + 2d0*dC3 + dC4)/6d0

      END IF

      !Convert to Voigt notation
      DO ii = 1, 6 ; DO jj = 1 , 6
         Cdem(ii,jj) = c0dem(l1(ii),l2(ii),l1(jj),l2(jj))
      END DO ; END DO

      write(* ,'(10f12.5)') Cdem(1,1),Cdem(2,2),Cdem(3,3),Cdem(4,4),Cdem(5,5),Cdem(6,6),Cdem(1,2),Cdem(1,3),Cdem(2,3),Vol(cyc)

   END DO;

   END IF

   print *
   write(*,'(a,1f6.2,a)'),' Time to calculate dem = ',OMP_GET_WTIME() - wtime,' sec'
   print *

   DEALLOCATE(thetaspo,phispo)
   DEALLOCATE(esa1,esa2,esa3)

   END PROGRAM DEMelastic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Green function 4th order tensor                                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE greenfunction(c0dem,Cinc3,vol_fract,esa,dc0dem)

   USE comvardem

   IMPLICIT NONE

   INTEGER :: ii,jj,kk,ll
   DOUBLE PRECISION :: vol_fract,esa(3) 
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: c0dem,Cinc3,dc0dem,Ai,Ainv,dC,Gs,G,SdC
   DOUBLE PRECISION, DIMENSION(6,6) :: Ai_mandel,Ainv_mandel 
  
   !Calculate interaction tensor
   CALL interactiontensor(c0dem,esa,G)

   ! Symmetric interaction tensor
   Gs = 0d0
   DO ii = 1,3; DO jj = 1,3; DO kk = 1,3; DO ll=1,3
      !Mainprice
      !Gs(ii,kk,jj,ll) = 0.5d0*(G(ii,kk,jj,ll)+G(jj,kk,ii,ll))
      !Hornby
      Gs(ii,jj,kk,ll) = 0.5d0*(G(ii,kk,jj,ll)+G(jj,kk,ii,ll))
   END DO; END DO; END DO; END DO
   
   !Difference in viscous tensor between inclusion and matrix 
   dC = Cinc3 - c0dem

   !Eshelby tensor = S = Jd*T*Cm = Ts*Cm 
   CALL tensorcontraction(Gs,dC,SdC)
   
   !Add symmetric identity 4th-order tensor
   Ai =  Is + SdC

   !Convert Ai to Mandel notation
   DO ii = 1, 6 ; DO jj = 1 , 6
      Ai_mandel(ii,jj) = Ai(l1(ii),l2(ii),l1(jj),l2(jj))*mandel_scale(ii,jj)
   END DO ; END DO
   
   !Invert Ai
   CALL inverse(Ai_mandel,Ainv_mandel,6)

   !Convert from Mandel notation to 4th order tensor
   DO ii = 1 , 3 ; DO jj = 1 , 3 ; DO kk = 1 , 3 ; DO ll = 1 , 3
      Ainv(ii,jj,kk,ll) = Ainv_mandel(ijkl(ii,jj),ijkl(kk,ll))/mandel_scale(ijkl(ii,jj),ijkl(kk,ll))
   END DO ; END DO ; END DO ; END DO
         
   !New composite stiffness tensor
   CALL tensorcontraction(dC,Ainv,dc0dem)
   
   !Multiply by volume fraction increment
   dc0dem = vol_fract*dc0dem

   END SUBROUTINE greenfunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE interactiontensor(c0dem,esa,G)
   
   USE comvardem

   IMPLICIT NONE

   INTEGER :: ii,jj,kk,ll,pp,qq,rr,ss,tt,uu
   DOUBLE PRECISION :: SurfArea,x(3),esa(3)  
   DOUBLE PRECISION, DIMENSION(3,3) :: Kinv,K    
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: c0dem,G
  
   G = 0d0
   SurfArea = 0d0
   DO pp=1,thetanum; DO qq=1,phinum
      
      !Directions
      x(1)=sin(thetaspo(pp))*cos(phispo(qq))/esa(1)
      x(2)=sin(thetaspo(pp))*sin(phispo(qq))/esa(2)
      x(3)=cos(thetaspo(pp))/esa(3)
      
      !Christoffel stiffness tensor K
      K = 0d0
      DO rr = 1,3; DO ss = 1,3; DO tt = 1,3; DO uu=1,3
         K(rr,tt)=K(rr,tt)+c0dem(rr,ss,tt,uu)*x(ss)*x(uu)
      END DO; END DO; END DO; END DO
       
      !Find inverse of K
      CALL inverse(K,Kinv,3)

      DO ii = 1,3; DO kk = 1,3; DO jj = 1,3; DO ll=1,3
         !Hornby, 1994, Geophysics
         G(ii,jj,kk,ll) = G(ii,jj,kk,ll) + Kinv(ii,jj)*x(kk)*x(ll)*sin(thetaspo(pp))!*stp*stp
         !Mainprice, 2015, Treat. Geophys.
         !G(ii,kk,jj,ll) = G(ii,kk,jj,ll) + Kinv(ii,jj)*x(kk)*x(ll)*sin(thetaspo(pp))!*stp*stp
      END DO; END DO; END DO; END DO
      SurfArea = SurfArea + sin(thetaspo(pp))!*stp*stp
   END DO; END DO

   G = G/SurfArea

   END SUBROUTINE interactiontensor
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Double index 4th order tensor contraction                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE tensorcontraction(A,B,C)

   IMPLICIT NONE

   INTEGER :: ii,jj,kk,ll,mm,nn
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: A,B,C

   C = 0d0
   DO ii = 1,3; DO jj = 1,3; DO kk = 1,3; DO ll=1,3; DO mm = 1,3; DO nn = 1,3
      C(ii,jj,kk,ll) = C(ii,jj,kk,ll) + A(ii,jj,mm,nn)*B(mm,nn,kk,ll)
   END DO; END DO; END DO; END DO; END DO; END DO

   END SUBROUTINE tensorcontraction
