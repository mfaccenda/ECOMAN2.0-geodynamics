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

   MODULE comvardem

   INTEGER :: phinum,thetanum,etacontrnum,volnum,Volsavenum,esanum(3),devmod,isotropicmod
   DOUBLE PRECISION :: pi,stp,Vstpdem,Volsavemin,Volsavemax,Volsavestp
   DOUBLE PRECISION, DIMENSION(6,6) :: Cmat,Cinc,mandel_scale
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: Is,Id
   INTEGER, DIMENSION(3,3) :: ijkl ! tensor of indices to form Cijkl from Sij
   DOUBLE PRECISION, DIMENSION(3,3) :: kron ! tensor of indices to form Cijkl from Sij
   INTEGER, DIMENSION(6) :: l1,l2 ! tensot of indices to form Sij from Cijkl
   INTEGER, DIMENSION(:), ALLOCATABLE :: Volsave
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: theta,phi,esa1,esa2,esa3,etacontr,Vol,Volsavedb
   ! input file name before file number
   CHARACTER (300) :: output_filename

   END MODULE comvardem


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Read input file inputspo.dat and allocate matrices    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE initspo

   USE comvardem
   USE omp_lib
   
   IMPLICIT NONE

   INTEGER :: ii,jj,kk,ll,j1,j2,etacontrlogscalemod,axislogscalemod
   DOUBLE PRECISION :: Vmin,Vmax,nu_inc,nu_back
   DOUBLE PRECISION :: esamin(3),esamax(3),esastp(3)    
   DOUBLE PRECISION :: etacontrmin,etacontrmax,etacontrstp    

!!! Read input values for SPO model       
   OPEN(15,file='dem_input.dat')

   !Output filename                
   !call read_par_char(15,output_filename)
   output_filename='viscoustensordem.h5'

   call read_par_int(15,devmod)

   !Set viscosity 
   call read_par_int(15,isotropicmod)
   call read_par_int(15,etacontrlogscalemod)
   call read_par_double(15,etacontrmin)
   call read_par_double(15,etacontrmax)
   call read_par_double(15,etacontrstp)

   call read_par_double6(15,Cmat(1,1:6))
   call read_par_double6(15,Cmat(2,1:6))
   call read_par_double6(15,Cmat(3,1:6))
   call read_par_double6(15,Cmat(4,1:6))
   call read_par_double6(15,Cmat(5,1:6))
   call read_par_double6(15,Cmat(6,1:6))

   call read_par_double6(15,Cinc(1,1:6))
   call read_par_double6(15,Cinc(2,1:6))
   call read_par_double6(15,Cinc(3,1:6))
   call read_par_double6(15,Cinc(4,1:6))
   call read_par_double6(15,Cinc(5,1:6))
   call read_par_double6(15,Cinc(6,1:6))

   !Inclusions shape
   call read_par_int(15,axislogscalemod)
   call read_par_double3(15,esamin)  !!Ellipsoid min semiaxes
   call read_par_double3(15,esamax)  !!Ellipsoid max semiaxes
   call read_par_double3(15,esastp)  !!Ellipsoid semiaxes increment

   !Inclusions volme fractions to be saved in the database
   call read_par_double(15,Volsavemin)   
   call read_par_double(15,Volsavemax)   
   call read_par_double(15,Volsavestp)   

   !DEM volume increment
   call read_par_double(15,Vstpdem)   

   if(devmod > 0) then

      write(*,*)
      write(*,"(a)"),' COMPUTE DEVIATORIC VISCOUS TENSOR '
      write(*,*)

   else

      write(*,*)
      write(*,"(a)"),' COMPUTE TOTAL VISCOUS TENSOR '
      write(*,*)

   end if

   write(*,"(a,a)"),' SAVE VISCOUS TENSOR TO ',trim(output_filename)

   if(isotropicmod > 0) then

      write(*,*)
      write(*,"(a)"),' INPUT ISOTROPIC VISCOUS TENSOR '
      write(*,*)

   else

      write(*,*)
      write(*,"(a)"),' INPUT VISCOUS TENSOR FROM CIJ COMPONENTS '
      write(*,*)

   end if

   if(etacontrlogscalemod > 0) then

      write(*,*)
      write(*,"(a)"),' DEFINE VISCOSITY CONTRAST ARRAY IN LOG10 SCALE'
      write(*,*)

   end if

   etacontrnum = NINT((etacontrmax - etacontrmin)/etacontrstp) + 1
   ALLOCATE(etacontr(etacontrnum))
   DO ii=1,etacontrnum
      etacontr(ii) = etacontrmin + etacontrstp*(ii-1)
      IF(etacontrlogscalemod > 0) etacontr(ii) = 10d0**etacontr(ii)
   END DO

   print *,'Viscosity contrast: ',etacontr

   if(axislogscalemod > 0) then

      write(*,*)
      write(*,"(a)"),' DEFINE SEMIAXES ARRAYS IN LOG10 SCALE'
      write(*,*)

   end if

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
   Vmin = 0d0
   Vmax = Volsavemax

   volsavenum = NINT((Volsavemax-Volsavemin)/Volsavestp)+1
   ALLOCATE(Volsavedb(volsavenum))
   DO ii = 1,volsavenum
      Volsavedb(ii) = Volsavemin + Volsavestp*(ii-1)
   END DO

   volnum = NINT((Vmax-Vmin)/Vstpdem)+1
   ALLOCATE(Vol(volnum),Volsave(volnum))
   kk = 0
   DO ii = 1,volnum
      Vol(ii) = Vstpdem*(ii-1)
      Volsave(ii) = 0
      DO jj = 1,volsavenum
         IF(Vol(ii) == Volsavedb(jj)) THEN
            kk = kk + 1
            Volsave(ii) = kk
         END IF
      END DO
   END DO
       
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
  
   Is = 0d0 ; Id= 0d0 ; kron = 0d0

   !Kronecker delta
   kron(1,1)=1;
   kron(2,2)=1;
   kron(3,3)=1;
   
   !Deviatoric symmetric fourth-rank unit tensor
   DO ii = 1,3; DO jj = 1,3; DO kk = 1,3; DO ll=1,3
      Is(ii,jj,kk,ll) = Is(ii,jj,kk,ll) + 0.5*(kron(ii,kk)*kron(jj,ll)+kron(ii,ll)*kron(jj,kk));
      Id(ii,jj,kk,ll) = Id(ii,jj,kk,ll) + 0.5*(kron(ii,kk)*kron(jj,ll)+kron(ii,ll)*kron(jj,kk)) - kron(ii,jj)*kron(kk,ll)/3;
   END DO; END DO; END DO; END DO

   pi = 3.141592653589793238462643383279
   stp=pi/100d0   
   thetanum = INT(pi/stp)+1
   ALLOCATE(theta(thetanum))
   DO ii = 1, thetanum
      theta(ii) = stp*(ii)
   END DO  
   phinum = INT(2*pi/stp)+1
   ALLOCATE(phi(phinum))
   DO ii = 1, phinum
      phi(ii) = stp*(ii)
   END DO  

   END SUBROUTINE initspo  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  SPO calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   PROGRAM DEM              

   USE comvardem
   USE omp_lib    
   USE hdf5       

   IMPLICIT NONE

   INTEGER :: ii,jj,kk,ll,m,pp,qq,rr,ss,tt,uu,cyc,vcyc,nu,tid,nt
   INTEGER :: j1,j2,cyc_esa1,cyc_esa2,cyc_esa,gi,nx(4) 
   DOUBLE PRECISION :: N,L,E_chi,Gc,Gs,G,phiazi,azimuthal(3),esa(3),dum_db(3),wtime
   DOUBLE PRECISION :: vol_fract
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: esa12,esa23
   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: C11,C22,C33,C44,C55,C66,C12,C13,C23
   DOUBLE PRECISION, DIMENSION(6,6) :: Cdem
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: c0dem,Cmat3,Cinc3,Ctemp,dC1,dC2,dC3,dC4
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
   
   !Multiple viscosity contrasts/ellipsoids
   if((isotropicmod > 0 .AND. etacontrnum > 1) .OR. esanum(1)*esanum(2) > 1 .OR. volsavenum > 1) THEN

   !Allocate memory for output 3D arrays
   ALLOCATE(C11(etacontrnum,volsavenum,esanum(1),esanum(2)))
   ALLOCATE(C22(etacontrnum,volsavenum,esanum(1),esanum(2)))
   ALLOCATE(C33(etacontrnum,volsavenum,esanum(1),esanum(2)))
   ALLOCATE(C44(etacontrnum,volsavenum,esanum(1),esanum(2)))
   ALLOCATE(C55(etacontrnum,volsavenum,esanum(1),esanum(2)))
   ALLOCATE(C66(etacontrnum,volsavenum,esanum(1),esanum(2)))
   ALLOCATE(C12(etacontrnum,volsavenum,esanum(1),esanum(2)))
   ALLOCATE(C13(etacontrnum,volsavenum,esanum(1),esanum(2)))
   ALLOCATE(C23(etacontrnum,volsavenum,esanum(1),esanum(2)))
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
        esa12(cyc_esa,cyc_esa2)= log10(esa(1)/esa(2))   
        esa23(cyc_esa,cyc_esa2)= log10(esa(2)/esa(3))    
!print *,cyc_esa2,cyc_esa1,cyc_esa
!print *,esa
!print *,esa23(cyc_esa2,cyc_esa),esa12(cyc_esa2,cyc_esa)
!print *
!read(*,*)
     END IF

   END DO; END DO
!stop
   !Viscosity contrast loop
   DO nu = 1, etacontrnum

   Cmat = 0 ; Cinc = 0
   DO jj=1,6; DO ii=1,6
      IF(ii == jj) THEN
         Cmat(ii,jj)=1d0     
         Cinc(ii,jj)=etacontr(nu)
         IF(ii < 4) Cmat(ii,jj) = Cmat(ii,jj)*2
         IF(ii < 4) Cinc(ii,jj) = Cinc(ii,jj)*2
      END IF
   END DO; END DO

   !Convert from Voigt notation to 4th order tensor
   DO ll = 1,3; DO kk = 1,3; DO jj = 1,3; DO ii=1,3
      Cmat3(ii,jj,kk,ll) = Cmat(ijkl(ii,jj),ijkl(kk,ll))
      Cinc3(ii,jj,kk,ll) = Cinc(ijkl(ii,jj),ijkl(kk,ll))
   END DO; END DO; END DO; END DO

   !Compute deviatoric viscous tensor if devmod is active
   IF(devmod) THEN

      c0dem=Cmat3
      Ctemp=Cinc3
      CALL tensorcontraction(c0dem,Id,Cmat3)
      CALL tensorcontraction(Ctemp,Id,Cinc3)

      DO ll = 1,3; DO kk = 1,3; DO jj = 1,3; DO ii=1,3
         Cmat(ijkl(ii,jj),ijkl(kk,ll)) = Cmat3(ii,jj,kk,ll)
         Cinc(ijkl(ii,jj),ijkl(kk,ll)) = Cinc3(ii,jj,kk,ll)
      END DO; END DO; END DO; END DO

   END IF

   write(*,"(a)"),'Matrix'
   write(*,"(6f12.5)") Cmat 
   write(*,*)
   write(*,"(a)"),'Inclusions'
   write(*,"(6f12.5)") Cinc
   write(*,*)

   !Skip null viscosity contrast
   IF(etacontr(nu) == 1d0) THEN
      C11(nu,:,:,:)= Cmat(1,1)
      C22(nu,:,:,:)= Cmat(2,2)
      C33(nu,:,:,:)= Cmat(3,3)
      C44(nu,:,:,:)= Cmat(4,4)
      C55(nu,:,:,:)= Cmat(5,5)
      C66(nu,:,:,:)= Cmat(6,6)
      C12(nu,:,:,:)= Cmat(1,2)
      C13(nu,:,:,:)= Cmat(1,3)
      C23(nu,:,:,:)= Cmat(2,3)
      GOTO 20
   END IF

   !$omp parallel & 
   !$omp shared(esa1,esa2,esa3,Vol,etacontr,C11,C22,C33,C44,C55,C66,C12,C13,C23,esa12,esa23)  &
   !$omp private(tid,esa,cyc_esa1,cyc_esa2,cyc_esa,cyc,vol_fract,c0dem,Ctemp,Cdem,dC1,dC2,dC3,dC4,ii,jj) &
   !$omp firstprivate(esanum,volnum,nu,Cmat3,Cinc3,Vstpdem,l1,l2)   
   !$omp do schedule(static) collapse(2) 
   DO cyc_esa2 = 1, esanum(2)
   DO cyc_esa1 = 1, esanum(1)
  
   tid = OMP_GET_THREAD_NUM()

   cyc_esa = cyc_esa1 - cyc_esa2 + 1

   esa(1)=esa1(cyc_esa1)
   esa(2)=esa2(cyc_esa2)
   esa(3)=1

!print *,tid,cyc_esa1,cyc_esa2,esa(1),esa(2)

   !Skip ellipsoids when a2 > a1
   IF(esa(2) > esa(1)) GOTO 10

   DO cyc=1,volnum
      IF(Vol(cyc) .EQ. 0d0) THEN
         c0dem = Cmat3
      ELSE IF(Vol(cyc) .EQ. 1d0) THEN
         c0dem = Cinc3
      ELSE
         !vol_fract = Vstpdem/(1d0-Vol(cyc)) 
         vol_fract = log(1d0+Vstpdem/(1d0-Vol(cyc))) 
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
        
      IF(Volsave(cyc) > 0) THEN

         !Convert to Voigt notation
         DO jj = 1, 6 ; DO ii = 1 , 6
            Cdem(ii,jj) = c0dem(l1(ii),l2(ii),l1(jj),l2(jj))
         END DO ; END DO
        
         vcyc = Volsave(cyc)
         !Compute global index in the output matrix
         C11(nu,vcyc,cyc_esa,cyc_esa2)= Cdem(1,1)
         C22(nu,vcyc,cyc_esa,cyc_esa2)= Cdem(2,2)
         C33(nu,vcyc,cyc_esa,cyc_esa2)= Cdem(3,3)
         C44(nu,vcyc,cyc_esa,cyc_esa2)= Cdem(4,4)
         C55(nu,vcyc,cyc_esa,cyc_esa2)= Cdem(5,5)
         C66(nu,vcyc,cyc_esa,cyc_esa2)= Cdem(6,6)
         C12(nu,vcyc,cyc_esa,cyc_esa2)= Cdem(1,2)
         C13(nu,vcyc,cyc_esa,cyc_esa2)= Cdem(1,3)
         C23(nu,vcyc,cyc_esa,cyc_esa2)= Cdem(2,3)
      END IF

   END DO; 

10 END DO; END DO
   !$omp end do
   !$omp end parallel

20 END DO 
   !End Viscosity contrast loop

   !!! Write infos in hdf5 format

   !Initialize FORTRAN interface.

   CALL H5open_f (error)

   ! Create a new file using default properties.

   CALL H5Fcreate_f(output_filename, H5F_ACC_TRUNC_F, file_id, error)

   !Save output matrix
   nx(1)=etacontrnum   
   nx(2)=volsavenum
   nx(3)=esanum(1)
   nx(4)=esanum(2)
   CALL loadsave_integer(0,1,file_id,4,H5T_NATIVE_INTEGER,nx,'gridnum',1)
   CALL loadsave_double(0,1,file_id,nx(1),H5T_NATIVE_DOUBLE,etacontr,'etacontr',1)
   CALL loadsave_double(0,1,file_id,nx(2),H5T_NATIVE_DOUBLE,Volsavedb,'Vol',1)
   CALL loadsave_double(0,2,file_id,nx(3:4),H5T_NATIVE_DOUBLE,esa12,'esa12',1)
   CALL loadsave_double(0,2,file_id,nx(3:4),H5T_NATIVE_DOUBLE,esa23,'esa23',1)
   CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C11,'C11',1)
   CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C22,'C22',1)
   CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C33,'C33',1)
   CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C44,'C44',1)
   CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C55,'C55',1)
   CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C66,'C66',1)
   CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C12,'C12',1)
   CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C13,'C13',1)
   CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C23,'C23',1)

   CALL H5Fclose_f(file_id, error)

   CALL H5close_f(error)
   !Close FORTRAN interface.

   DEALLOCATE(Vol,C11,C22,C33,C44,C55,C66,C12,C13,C23,esa12,esa23,etacontr)
   
   !Single ellipsoid
   ELSE

   IF(isotropicmod > 0) THEN

   DO jj=1,6; DO ii=1,6
      IF(ii == jj) THEN
         Cmat(ii,jj)=1d0     
         Cinc(ii,jj)=etacontr(1)
         IF(ii < 4) Cmat(ii,jj) = Cmat(ii,jj)*2
         IF(ii < 4) Cinc(ii,jj) = Cinc(ii,jj)*2
      END IF
   END DO; END DO

   END IF

   !Convert from Voigt notation to 4th order tensor
   DO ii = 1,3; DO jj = 1,3; DO kk = 1,3; DO ll=1,3
      Cmat3(ii,jj,kk,ll) = Cmat(ijkl(ii,jj),ijkl(kk,ll))
      Cinc3(ii,jj,kk,ll) = Cinc(ijkl(ii,jj),ijkl(kk,ll))
   END DO; END DO; END DO; END DO

   !Compute deviatoric viscous tensor if devmod is active
   IF(devmod) THEN

      c0dem=Cmat3
      Ctemp=Cinc3
      CALL tensorcontraction(c0dem,Id,Cmat3)
      CALL tensorcontraction(Ctemp,Id,Cinc3)

      DO ii = 1,3; DO jj = 1,3; DO kk = 1,3; DO ll=1,3
         Cmat(ijkl(ii,jj),ijkl(kk,ll)) = Cmat3(ii,jj,kk,ll)
         Cinc(ijkl(ii,jj),ijkl(kk,ll)) = Cinc3(ii,jj,kk,ll)
      END DO; END DO; END DO; END DO

   END IF

   esa(1)=esa1(1)
   esa(2)=esa2(1)
   esa(3)=1

   print *,''
   print *,'Matrix viscous tensor:'
   write(*,'(6f12.5)') Cmat
   print *,''
   print *,'Inclusion viscous tensor:'
   write(*,'(6f12.5)') Cinc
   print *,''
   write(*,'(a,3f14.5)') 'Inclusion semiaxes:',esa
   print *,''
   write(*,'(a)'),'      Nxxxx       Nyyyy       Nzzzz       Nyzyz       Nxzxz       Nxyxy       Nxxyy       Nxxzz       Nyyzz      Vol fract'

   DO cyc=1,volnum
      IF(Vol(cyc) .EQ. 0d0) THEN
         c0dem = Cmat3
      ELSE IF(Vol(cyc) .EQ. 1d0) THEN
         c0dem = Cinc3
      ELSE
         !vol_fract = Vstpdem/(1d0-Vol(cyc)) 
         vol_fract = log(1d0+Vstpdem/(1d0-Vol(cyc))) 
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
 
   DEALLOCATE(theta,phi)
   DEALLOCATE(esa1,esa2,esa3)
  
   END PROGRAM DEM    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Double index 4th order tensor contraction                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE tensorcontraction(A,B,C)

   USE comvardem

   IMPLICIT NONE

   INTEGER :: ii,jj,kk,ll,mm,nn
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: A,B,C

   C = 0d0
   DO ii = 1,3; DO jj = 1,3; DO kk = 1,3; DO ll=1,3; DO mm = 1,3; DO nn = 1,3
      C(ii,jj,kk,ll) = C(ii,jj,kk,ll) + A(ii,jj,mm,nn)*B(mm,nn,kk,ll)
   END DO; END DO; END DO; END DO; END DO; END DO

   END SUBROUTINE tensorcontraction
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Green function 4th order tensor                                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE greenfunction(c0dem,Cinc3,vol_fract,esa,dc0dem)

   USE comvardem

   IMPLICIT NONE

   INTEGER :: ii,jj,kk,ll
   DOUBLE PRECISION :: vol_fract,esa(3)
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: c0dem,Cinc3,dc0dem,Ai,Ainv,dC,Ts,T,SdC
   DOUBLE PRECISION, DIMENSION(6,6) :: Ai_mandel,Ainv_mandel 
  
   !Calculate interaction tensor
   CALL interactiontensor(c0dem,esa,T)

   ! Symmetric interaction tensor
   Ts = 0d0
   DO ii = 1,3; DO jj = 1,3; DO kk = 1,3; DO ll=1,3
      Ts(ii,jj,kk,ll) = Ts(ii,jj,kk,ll) + 0.25*(T(ii,jj,kk,ll) + T(jj,ii,kk,ll) + T(ii,jj,ll,kk) + T(jj,ii,ll,kk))
   END DO; END DO; END DO; END DO
   
   !Difference in viscous tensor between inclusion and matrix 
   dC = Cinc3 - c0dem

   !Eshelby tensor = S = Jd*T*Cm = Ts*Cm 
   CALL tensorcontraction(Ts,dC,SdC)
   
   !Add symmetric identity 4th-order tensor
   !In reality it should be Ai = Id + SdC, which is however not invertible in the
   !cartesian space (Lebenson et al., 1998)
   Ai =  Is + SdC

   !Convert Ai to Voigt notation
   DO ii = 1, 6 ; DO jj = 1 , 6
      Ai_mandel(ii,jj) = Ai(l1(ii),l2(ii),l1(jj),l2(jj))*mandel_scale(ii,jj)
   END DO ; END DO
   
   !Invert Ai
   CALL inverse(Ai_mandel,Ainv_mandel,6)

   !Convert from Voigt notation to 4th order tensor
   Ai = 0d0
   DO ii = 1 , 3 ; DO jj = 1 , 3 ; DO kk = 1 , 3 ; DO ll = 1 , 3
      Ai(ii,jj,kk,ll) = Ainv_mandel(ijkl(ii,jj),ijkl(kk,ll))/mandel_scale(ijkl(ii,jj),ijkl(kk,ll))
   END DO ; END DO ; END DO ; END DO
         
   !Take deviatoric part of Ai     
   CALL tensorcontraction(Id,Ai,Ainv)

   !New composite stiffness tensor
   CALL tensorcontraction(dC,Ainv,dc0dem)
   
   !Multiply by volume fraction increment
   dc0dem = vol_fract*dc0dem

   END SUBROUTINE greenfunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE interactiontensor(c0dem,esa,T)
   
   USE comvardem

   IMPLICIT NONE

   INTEGER :: ii,jj,kk,ll,pp,qq,rr,ss,tt,uu
   DOUBLE PRECISION :: SurfArea,x(3),esa(3) 
   DOUBLE PRECISION, DIMENSION(4,4) :: K4inv,K4    
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: c0dem,T
  
   T = 0d0
   SurfArea = 0d0
  
   DO pp=1,thetanum; DO qq=1,phinum/2

      !Directions
      x(1)=sin(theta(pp))*cos(phi(qq))/esa(1)
      x(2)=sin(theta(pp))*sin(phi(qq))/esa(2)
      x(3)=cos(theta(pp))/esa(3)

      !Christoffel viscosity tensor
      K4 = 0d0
      DO rr = 1,3; DO ss = 1,3; DO tt = 1,3; DO uu=1,3
         K4(rr,tt)=K4(rr,tt)+c0dem(rr,ss,tt,uu)*x(ss)*x(uu)
      END DO; END DO; END DO; END DO

      !Viscosity stiffness tensor for incompressible materials (Jiang, 2014,
      !JSG)
      K4(4,1:3) = x(:)
      K4(1:3,4) = x(:)
      K4(4,4) = 0d0

      !Find inverse of K matrix
      CALL inv4x4(K4,K4inv)

      DO ii = 1,3; DO kk = 1,3; DO jj = 1,3; DO ll=1,3
         !Jiang, 2014, JSG
         T(ii,jj,kk,ll) = T(ii,jj,kk,ll) + K4inv(ii,kk)*x(jj)*x(ll)*sin(theta(pp))
      END DO; END DO; END DO; END DO
      SurfArea = SurfArea + sin(theta(pp))!*stp*stp
   END DO; END DO

   T = T/SurfArea

   END SUBROUTINE interactiontensor

   !============================================================
   subroutine det4x4(A,detA)
   !============================================================
   ! determinant of 4 by 4 matrix
   !===========================================================
   ! implicit none
   
   double precision, dimension(4,4) :: A
   double precision :: detA
   
   detA = A(1,1)*A(2,2)*A(3,3)*A(4,4) + & 
          A(1,1)*A(2,3)*A(3,4)*A(4,2) + & 
          A(1,1)*A(2,4)*A(3,2)*A(4,3) - & 
          A(1,1)*A(2,4)*A(3,3)*A(4,2) - & 
          A(1,1)*A(2,3)*A(3,2)*A(4,4) - & 
          A(1,1)*A(2,2)*A(3,4)*A(4,3) - & 
          A(1,2)*A(2,1)*A(3,3)*A(4,4) - & 
          A(1,3)*A(2,1)*A(3,4)*A(4,2) - & 
          A(1,3)*A(2,1)*A(3,2)*A(4,3) + & 
          A(1,4)*A(2,1)*A(3,3)*A(4,2) + & 
          A(1,3)*A(2,1)*A(3,2)*A(4,4) + & 
          A(1,2)*A(2,1)*A(3,4)*A(4,3) + &
          A(1,2)*A(2,3)*A(3,1)*A(4,4) + & 
          A(1,3)*A(2,4)*A(3,1)*A(4,2) + & 
          A(1,4)*A(2,2)*A(3,1)*A(4,3) - &
          A(1,4)*A(2,3)*A(3,1)*A(4,2) - & 
          A(1,3)*A(2,2)*A(3,1)*A(4,4) - & 
          A(1,2)*A(2,4)*A(3,1)*A(4,3) - &
          A(1,2)*A(2,3)*A(3,4)*A(4,1) - & 
          A(1,3)*A(2,4)*A(3,2)*A(4,1) - & 
          A(1,4)*A(2,2)*A(3,3)*A(4,1) + &
          A(1,4)*A(2,3)*A(3,2)*A(4,1) + &
          A(1,3)*A(2,2)*A(3,4)*A(4,1) + &
          A(1,2)*A(2,4)*A(3,3)*A(4,1)
   
   end subroutine det4x4

   !============================================================
   subroutine adjointA(A,Aprime)
   !============================================================
   ! adjoint of 4 by 4 matrix
   !===========================================================
   ! implicit none
   double precision, dimension(4,4) :: A,Aprime
   
   Aprime(1,1) = (A(2,2)*A(3,3)*A(4,4) + A(2,3)*A(3,4)*A(4,2) + A(2,4)*A(3,2)*A(4,3) - &
                     A(2,4)*A(3,3)*A(4,2) - A(2,3)*A(3,2)*A(4,4) - A(2,2)*A(3,4)*A(4,3) )

   Aprime(1,2) = (-A(1,2)*A(3,3)*A(4,4) - A(1,3)*A(3,4)*A(4,2) - A(1,4)*A(3,2)*A(4,3) + &
                      A(1,4)*A(3,3)*A(4,2) + A(1,3)*A(3,2)*A(4,4) + A(1,2)*A(3,4)*A(4,3) )

   Aprime(1,3) = (A(1,2)*A(2,3)*A(4,4) + A(1,3)*A(2,4)*A(4,2) + A(1,4)*A(2,2)*A(4,3) - &
                     A(1,4)*A(2,3)*A(4,2) - A(1,3)*A(2,2)*A(4,4) - A(1,2)*A(2,4)*A(4,3) )

   Aprime(1,4) = (-A(1,2)*A(2,3)*A(3,4) - A(1,3)*A(2,4)*A(3,2) - A(1,4)*A(2,2)*A(3,3) + &
                      A(1,4)*A(2,3)*A(3,2) + A(1,3)*A(2,2)*A(3,4) + A(1,2)*A(2,4)*A(3,3) )

   Aprime(2,1) = (-A(2,1)*A(3,3)*A(4,4) - A(2,3)*A(3,4)*A(4,1) - A(2,4)*A(3,1)*A(4,3) + &
                      A(2,4)*A(3,3)*A(4,1) + A(2,3)*A(3,1)*A(4,4) + A(2,1)*A(3,4)*A(4,3) )

   Aprime(2,2) = (A(1,1)*A(3,3)*A(4,4) + A(1,3)*A(3,4)*A(4,2) + A(1,4)*A(3,1)*A(4,3) - &
                     A(1,4)*A(3,3)*A(4,1) - A(1,3)*A(3,1)*A(4,4) - A(1,1)*A(3,4)*A(4,3) )

   Aprime(2,3) = (-A(1,1)*A(2,3)*A(4,4) - A(1,3)*A(2,4)*A(4,1) - A(1,4)*A(2,1)*A(4,3) + &
                      A(1,4)*A(2,3)*A(4,1) + A(1,3)*A(2,1)*A(4,4) + A(1,1)*A(2,4)*A(4,3) )

   Aprime(2,4) = (A(1,1)*A(2,3)*A(3,4) + A(1,3)*A(2,4)*A(3,1) + A(1,4)*A(2,1)*A(3,3) - &
                     A(1,4)*A(2,3)*A(3,1) - A(1,3)*A(2,1)*A(3,4) - A(1,1)*A(2,4)*A(3,3) )

   Aprime(3,1) = (A(2,1)*A(3,2)*A(4,4) + A(2,2)*A(3,4)*A(4,1) + A(2,4)*A(3,1)*A(4,2) - &
                     A(2,4)*A(3,2)*A(4,1) - A(2,2)*A(3,1)*A(4,4) - A(2,1)*A(3,4)*A(4,2) )

   Aprime(3,2) = (-A(1,1)*A(3,2)*A(4,4) - A(1,2)*A(3,4)*A(4,1) - A(1,4)*A(3,1)*A(4,2) + &
                      A(1,4)*A(3,2)*A(4,1) + A(1,2)*A(3,1)*A(4,4) + A(1,1)*A(3,4)*A(4,2) )

   Aprime(3,3) = (A(1,1)*A(2,2)*A(4,4) + A(1,2)*A(2,4)*A(4,1) + A(1,4)*A(2,1)*A(4,2) - &
                     A(1,4)*A(2,2)*A(4,1) - A(1,2)*A(2,1)*A(4,4) - A(1,1)*A(2,4)*A(4,2) )

   Aprime(3,4) = (-A(1,1)*A(2,2)*A(3,4) - A(1,2)*A(2,4)*A(3,1) - A(1,4)*A(2,1)*A(3,2) + &
                      A(1,4)*A(2,2)*A(3,1) + A(1,2)*A(2,1)*A(3,4) + A(1,1)*A(2,4)*A(3,2) )

   Aprime(4,1) = (-A(2,1)*A(3,2)*A(4,3) - A(2,2)*A(3,3)*A(4,1) - A(2,3)*A(3,1)*A(4,2) + &
                      A(2,3)*A(3,2)*A(4,1) + A(2,2)*A(3,1)*A(4,3) + A(2,1)*A(3,3)*A(4,2) )

   Aprime(4,2) = (A(1,1)*A(3,2)*A(4,3) + A(1,2)*A(3,3)*A(4,1) + A(1,3)*A(3,1)*A(4,2) - &
                     A(1,3)*A(3,2)*A(4,1) - A(1,2)*A(3,1)*A(4,3) - A(1,1)*A(3,3)*A(4,2) )

   Aprime(4,3) = (-A(1,1)*A(2,2)*A(4,3) - A(1,2)*A(2,3)*A(4,1) - A(1,3)*A(2,1)*A(4,2) + &
                      A(1,3)*A(2,2)*A(4,1) + A(1,2)*A(2,1)*A(4,3) + A(1,1)*A(2,3)*A(4,2) )

   Aprime(4,4) = (A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) -  &
                     A(1,3)*A(2,2)*A(3,1) - A(1,2)*A(2,1)*A(3,3) - A(1,1)*A(2,3)*A(3,2) ) 
   
   end subroutine adjointA

   !============================================================
   subroutine inv4x4(A,Aprime)
   !============================================================
   ! inverse of 4 by 4 matrix
   !===========================================================
   implicit none
   integer :: i,j
   double precision, dimension(4,4) :: A,Aprime
   double precision :: d
   
   CALL adjointA(A,Aprime)
   CALL det4x4(A,d)

   do i=1,4
      do  j =1,4
         Aprime(i,j) = Aprime(i,j) / d
      end do 
   end do

   end subroutine inv4x4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE inverse(a,c,n)
   !============================================================
   ! Inverse matrix
   ! Method: Based on Doolittle LU factorization for Ax=b
   ! Alex G. December 2009
   !-----------------------------------------------------------
   ! input ...
   ! a(n,n) - array of coefficients for matrix A
   ! n      - dimension
   ! output ...
   ! c(n,n) - inverse matrix of A
   ! comments ...
   ! the original matrix a(n,n) will be destroyed 
   ! during the calculation
   !===========================================================
   implicit none
   integer n
   double precision a(n,n),c(n,n)
   double precision L(n,n),U(n,n), b(n), d(n), x(n)
   double precision coeff
   integer i, j, k

   ! step 0: initialization for matrices L and U and b
   ! Fortran 90/95 aloows such operations on matrices
   L=0.0d0
   U=0.0d0
   b=0.0d0

   ! step 1: forward elimination
   do k=1, n-1
      do i=k+1,n
         coeff=a(i,k)/a(k,k)
         L(i,k) = coeff
         do j=k+1,n
            a(i,j) = a(i,j)-coeff*a(k,j)
         end do
       end do
   end do

   ! Step 2: prepare L and U matrices 
   ! L matrix is a matrix of the elimination coefficient
   ! + the diagonal elements are 1.0
   do i=1,n
      L(i,i) = 1.0d0
   end do
   ! U matrix is the upper triangular part of A
   do j=1,n
      do i=1,j
         U(i,j) = a(i,j)
      end do
   end do

   ! Step 3: compute columns of the inverse matrix C
   do k=1,n
      b(k)=1.0d0
      d(1) = b(1)
   ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
         d(i)=b(i)
         do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
         end do
      end do
   ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
         x(i) = d(i)
         do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
         end do
         x(i) = x(i)/u(i,i)
      end do
   ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
         c(i,k) = x(i)
      end do
      b(k)=0.0d0
   end do

   END SUBROUTINE inverse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE read_par_char(unitt,str)

   IMPLICIT NONE

   INTEGER :: i,ier,unitt
   CHARACTER(500) :: str

10 READ(unitt,"(a)",advance='yes',IOSTAT=ier) str
   i=index(str,'#')
   if(i==1 .OR. LEN_TRIM(str)==0) then
      goto 10
   end if

   END SUBROUTINE read_par_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE read_par_int(unitt,intg)

   IMPLICIT NONE

   INTEGER ::i,ier,intg,unitt
   CHARACTER(500) :: str

10 READ(unitt,"(a)",advance='yes',IOSTAT=ier) str  
   i=index(str,'#')
   if(i==1 .OR. LEN_TRIM(str)==0) then
      goto 10
   else
      read(str,*) intg
   end if

   END SUBROUTINE read_par_int 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE read_par_double(unitt,dbl)

   IMPLICIT NONE

   INTEGER ::i,ier,unitt
   DOUBLE PRECISION :: dbl
   CHARACTER(500) :: str

10 READ(unitt,"(a)",advance='yes',IOSTAT=ier) str  
   i=index(str,'#')
   if(i==1 .OR. LEN_TRIM(str)==0) then
      goto 10
   else
      read(str,*) dbl
   end if

   END SUBROUTINE read_par_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE read_par_double3(unitt,dbl)

   IMPLICIT NONE

   INTEGER ::i,ier,unitt
   DOUBLE PRECISION :: dbl(3)
   CHARACTER(500) :: str

10 READ(unitt,"(a)",advance='yes',IOSTAT=ier) str  
   i=index(str,'#')
   if(i==1 .OR. LEN_TRIM(str)==0) then
      goto 10
   else
      read(str,*) dbl(:)
   end if

   END SUBROUTINE read_par_double3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE read_par_double6(unitt,dbl)

   IMPLICIT NONE

   INTEGER ::i,ier,unitt
   DOUBLE PRECISION :: dbl(6)
   CHARACTER(500) :: str

10 READ(unitt,"(a)",advance='yes',IOSTAT=ier) str  
   i=index(str,'#')
   if(i==1 .OR. LEN_TRIM(str)==0) then
      goto 10
   else
      read(str,*) dbl(:)
   end if

   END SUBROUTINE read_par_double6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Load/Save dataset, format double precision                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   SUBROUTINE loadsave_double(data_attr,rank,loc_id,dims,memtype,buf1,path,mode)

   USE comvardem
   USE hdf5

   IMPLICIT NONE

   CHARACTER(LEN=*), INTENT(IN) :: path
   INTEGER     ::   rank,dims(rank),data_attr,mode
   INTEGER(HID_T)  :: loc_id, dataspace_id, dataset_id, attr_id, memtype !Handles
   INTEGER(HSIZE_T), DIMENSION(1:rank) :: dims1
   INTEGER     ::   error  ! Error flag
   DOUBLE PRECISION, DIMENSION(*) :: buf1

   dims1=dims
   !Dataset
   IF(data_attr == 0) THEN
      IF(mode == 0) THEN
         CALL H5Dopen_f(loc_id, path, dataset_id, error)
         CALL H5Dread_f(dataset_id, memtype, buf1, dims1, error)
         CALL H5Dclose_f(dataset_id,error)
      ELSE
         CALL H5Screate_simple_f(rank,dims1,dataspace_id,error)
         CALL H5Dcreate_f(loc_id, path, memtype, dataspace_id, dataset_id,error)
         CALL H5Dwrite_f (dataset_id, memtype, buf1, dims1, error)
         CALL H5Dclose_f (dataset_id, error)
         CALL H5Sclose_f (dataspace_id, error)
      END IF
   END IF
   !Attribute
   IF(data_attr == 1) THEN
      IF(mode == 0) THEN
         CALL H5Aopen_f(loc_id, path, attr_id, error)
         CALL H5Aread_f(attr_id, memtype, buf1, dims1, error)
         CALL H5Aclose_f(attr_id, error)
      ELSE
         CALL H5Screate_simple_f(rank, dims1, dataspace_id, error)
         CALL H5Acreate_f(loc_id, path, memtype, dataspace_id, attr_id, error)
         CALL H5Awrite_f(attr_id, memtype, buf1, dims1, error)
         CALL H5Aclose_f(attr_id, error)
         CALL H5Sclose_f(dataspace_id, error)
      END IF
   END IF

   RETURN
   
   END SUBROUTINE loadsave_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Load/Save dataset, format integer                                      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   SUBROUTINE loadsave_integer(data_attr,rank,loc_id,dims,memtype,buf1,path,mode)

   USE comvardem
   USE hdf5

   IMPLICIT NONE

   CHARACTER(LEN=*), INTENT(IN) :: path
   INTEGER     ::   rank,dims(rank),data_attr,mode
   INTEGER(HID_T)  :: loc_id, dataspace_id, dataset_id, attr_id, memtype !Handles
   INTEGER(HSIZE_T), DIMENSION(1:rank) :: dims1
   INTEGER     ::   error  ! Error flag
   INTEGER, DIMENSION(*) :: buf1

   dims1=dims
   !Dataset
   IF(data_attr .EQ. 0) THEN
      IF(mode .eq. 0) THEN
         CALL H5Dopen_f(loc_id, path, dataset_id, error)
         CALL H5Dread_f(dataset_id, memtype, buf1, dims1, error)
         CALL H5Dclose_f(dataset_id,error)
      ELSE
         CALL H5Screate_simple_f(rank,dims1,dataspace_id,error)
         CALL H5Dcreate_f(loc_id, path, memtype, dataspace_id, dataset_id,error)
         CALL H5Dwrite_f (dataset_id, memtype, buf1, dims1, error)
         CALL H5Dclose_f (dataset_id, error)
         CALL H5Sclose_f (dataspace_id, error)
      END IF
   END IF
   !Attribute
   IF(data_attr .EQ. 1) THEN
      IF(mode .eq. 0) THEN
         CALL H5Aopen_f(loc_id, path, attr_id, error)
         CALL H5Aread_f(attr_id, memtype, buf1, dims1, error)
         CALL H5Aclose_f(attr_id, error)
      ELSE
         CALL H5Screate_simple_f(rank, dims1, dataspace_id, error)
         CALL H5Acreate_f(loc_id, path, memtype, dataspace_id, attr_id, error)
         CALL H5Awrite_f(attr_id, memtype, buf1, dims1, error)
         CALL H5Aclose_f(attr_id, error)
         CALL H5Sclose_f(dataspace_id, error)
      END IF
   END IF

   RETURN

   END SUBROUTINE loadsave_integer

