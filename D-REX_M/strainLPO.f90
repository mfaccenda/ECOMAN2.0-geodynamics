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
!!!                                                                        !!!
!!! Most of the subroutines in this file are present in the original D-REX !!!
!!!                                                                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine STRAIN - Calculation of strain along pathlines              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE strain(tid,m,fractdisl)

   USE comvar 
   USE omp_lib

   IMPLICIT NONE
   INTEGER :: m,i,j,k,ll,j1,j2,nnum,N_strain,tid,n,n1,n2,nrot
   INTEGER , DIMENSION(1) :: ti ! reordering array
   ! loop counters

   DOUBLE PRECISION :: dt_strain,fractdisl,dt_straintot,dt_straindiff,theta
   ! time spent on a given point of the streamline and DRex time step

   DOUBLE PRECISION, DIMENSION(3,3) :: fse,fsei,fse0,fsediff,Q,fseacs,fseacs0,lx,ex
   DOUBLE PRECISION, DIMENSION(3) :: evals,evals0
   ! local finite deformation tensor

   DOUBLE PRECISION, DIMENSION(3,3) :: kfse1,kfse2,kfse3,kfse4
   ! RGK intermediates for Finite Strain Ellipsoid

   DOUBLE PRECISION, DIMENSION(size) :: dotodf,dotodf_ens,odfi,odfi_ens
   DOUBLE PRECISION, DIMENSION(size) :: kodf1,kodf2,kodf3,kodf4
   DOUBLE PRECISION, DIMENSION(size) :: kodf1_ens,kodf2_ens,kodf3_ens,kodf4_ens
   ! intermediates for odf

   DOUBLE PRECISION, DIMENSION(size,3,3) :: dotacs,dotacs_ens,acsi,acsi_ens
   DOUBLE PRECISION, DIMENSION(size,3,3) :: kac1,kac2,kac3,kac4
   DOUBLE PRECISION, DIMENSION(size,3,3) :: kac1_ens,kac2_ens,kac3_ens,kac4_ens
   ! intermediates for matrix of direction cosine

!! time stepping for LPO calculation
   dt_straintot = MIN(dt,1d-2/epsnot(tid))

!!! WARNING
!   dt_straintot = dt

!! number of iterations in the LPO loop
   N_strain = NINT(dt/dt_straintot)

!!! Account for time left over
   IF(N_strain > 1 .AND. dt - REAL(N_strain)*dt_straintot > 0) N_strain = N_strain + 1

!!! Dummy fse
   fse0= Fij(:,:,m) 

!!! Main strain cycle
   DO nnum = 1 , N_strain

   IF(N_strain > 1 .AND. nnum == N_strain) dt_straintot = dt - REAL(N_strain - 1)*dt_straintot

   IF(dt_straintot < 0 .OR. dt_straintot > dt ) THEN
      print *,dt_straintot, timesum,N_strain,dt,REAL(N_strain)*dt_straintot
      stop
   END IF

!!! Calculate new fse
   fsei = fse0
   kfse1 = MATMUL(l(tid,:,:),fsei)*dt_straintot
   fsei = fse0 + 0.5d0*kfse1
   kfse2 = MATMUL(l(tid,:,:),fsei)*dt_straintot
   fsei = fse0 + 0.5d0*kfse2
   kfse3 = MATMUL(l(tid,:,:),fsei)*dt_straintot
   fsei = fse0 + kfse3
   kfse4 = MATMUL(l(tid,:,:),fsei)*dt_straintot
   fse = fse0 + (kfse1/2d0+kfse2+kfse3+kfse4/2d0)/3d0
 
!!! Compute fabric for two-phase mantle aggregates, except lower transition zone as ringwoodite and majoritic garnet are almost isotropic
!!! Compute fabric for both phases for upper mantle aggregates   

   IF(fsemod == 0 .AND. rocktype(m) /= 3) THEN

!!! Strain-induced LPO when dislocation creep is active           
   IF(fractdisl > 0d0) THEN

      !!! Scale time stepping for amount of deformation accommodated by disl. creep/anisotropic phases 
      dt_strain = dt_straintot*fractdisl

      !!! 1st increment
   
      !!! Main phase
      odfi = odf(m,:) ; acsi = acs(:,:,:,m)
      !!! Minor phase
      odfi_ens = odf_ens(m,:) ; acsi_ens = acs_ens(:,:,:,m)

      CALL deriv(tid,m,dotodf,dotodf_ens,odfi,odfi_ens,dotacs,dotacs_ens,acsi,acsi_ens)

      !!! Main phase
      kodf1 = dotodf*dt_strain*epsnot(tid)
      kac1 = dotacs*dt_strain*epsnot(tid)
      odfi = odf(m,:) + 0.5d0*kodf1
      acsi = acs(:,:,:,m) + 0.5d0*kac1
     
      DO j = 1 , size ; DO j1 = 1 , 3 ; DO j2 = 1, 3
         IF (acsi(j,j1,j2) .GT. 1d0) acsi(j,j1,j2) = 1d0
         IF (acsi(j,j1,j2) .LT. -1d0) acsi(j,j1,j2) = -1d0
      END DO ; END DO ; END DO
      DO j = 1 , size
         IF (odfi(j) .LE. 0 ) odfi(j) = 0d0
      END DO
      odfi = odfi/SUM(odfi)
      
      !!! Minor phase
      IF(rocktype(m) == 1) THEN 

         kodf1_ens = dotodf_ens*dt_strain*epsnot(tid)
         kac1_ens = dotacs_ens*dt_strain*epsnot(tid)

         odfi_ens = odf_ens(m,:) + 0.5d0*kodf1_ens
         acsi_ens = acs_ens(:,:,:,m) + 0.5d0*kac1_ens

         DO j = 1 , size ; DO j1 = 1 , 3 ; DO j2 = 1, 3
            IF (acsi_ens(j,j1,j2) .GT. 1d0) acsi_ens(j,j1,j2) = 1d0
            IF (acsi_ens(j,j1,j2) .LT. -1d0) acsi_ens(j,j1,j2) = -1d0
         END DO ; END DO ; END DO
         DO j = 1 , size
            IF (odfi_ens(j) .LE. 0 ) odfi_ens(j) = 0d0
         END DO
         odfi_ens = odfi_ens/SUM(odfi_ens) 

      END IF

      !!! 2nd increment
   
      CALL deriv(tid,m,dotodf,dotodf_ens,odfi,odfi_ens,dotacs,dotacs_ens,acsi,acsi_ens)

      !!! Main phase
      kodf2 = dotodf*dt_strain*epsnot(tid)
      kac2 = dotacs*dt_strain*epsnot(tid)
      odfi = odf(m,:) + 0.5d0*kodf2
      acsi = acs(:,:,:,m) + 0.5d0*kac2

      DO j = 1 , size ; DO j1 = 1 , 3 ; DO j2 = 1, 3
         IF (acsi(j,j1,j2) .GT. 1d0) acsi(j,j1,j2) = 1d0
         IF (acsi(j,j1,j2) .LT. -1d0) acsi(j,j1,j2) = -1d0
      END DO ; END DO ; END DO
      DO j = 1 , size
         IF (odfi(j) .LE. 0 ) odfi(j) = 0d0
      END DO
      odfi = odfi/SUM(odfi)

      !!! Minor phase
      IF(rocktype(m) == 1) THEN 

         kodf2_ens = dotodf_ens*dt_strain*epsnot(tid)
         kac2_ens = dotacs_ens*dt_strain*epsnot(tid)

         odfi_ens = odf_ens(m,:) + 0.5d0*kodf2_ens
         acsi_ens = acs_ens(:,:,:,m) + 0.5d0*kac2_ens

         DO j = 1 , size ; DO j1 = 1 , 3 ; DO j2 = 1, 3
            IF (acsi_ens(j,j1,j2) .GT. 1d0) acsi_ens(j,j1,j2) = 1d0
            IF (acsi_ens(j,j1,j2) .LT. -1d0) acsi_ens(j,j1,j2) = -1d0
         END DO ; END DO ; END DO
         DO j = 1 , size
            IF (odfi_ens(j) .LE. 0 ) odfi_ens(j) = 0d0
         END DO
         odfi_ens = odfi_ens/SUM(odfi_ens) 

      END IF

      !!! 3rd increment
   
      CALL deriv(tid,m,dotodf,dotodf_ens,odfi,odfi_ens,dotacs,dotacs_ens,acsi,acsi_ens)

      !!! Main phase
      kodf3 = dotodf*dt_strain*epsnot(tid)
      kac3 = dotacs*dt_strain*epsnot(tid)

      odfi = odf(m,:) + kodf3
      acsi = acs(:,:,:,m) + kac3

      DO j = 1 , size ; DO j1 = 1 , 3 ; DO j2 = 1, 3
         IF (acsi(j,j1,j2) .GT. 1d0) acsi(j,j1,j2) = 1d0
         IF (acsi(j,j1,j2) .LT. -1d0) acsi(j,j1,j2) = -1d0
      END DO ; END DO ; END DO
      DO j = 1 , size
         IF (odfi(j) .LE. 0 ) odfi(j) = 0d0
      END DO
      odfi = odfi/SUM(odfi)

      !!! Minor phase
      IF(rocktype(m) == 1) THEN 

         kodf3_ens = dotodf_ens*dt_strain*epsnot(tid)
         kac3_ens = dotacs_ens*dt_strain*epsnot(tid)

         odfi_ens = odf_ens(m,:) + kodf3_ens
         acsi_ens = acs_ens(:,:,:,m) + kac3_ens

         DO j = 1 , size ; DO j1 = 1 , 3 ; DO j2 = 1, 3
            IF (acsi_ens(j,j1,j2) .GT. 1d0) acsi_ens(j,j1,j2) = 1d0
            IF (acsi_ens(j,j1,j2) .LT. -1d0) acsi_ens(j,j1,j2) = -1d0
         END DO ; END DO ; END DO
         DO j = 1 , size
            IF (odfi_ens(j) .LE. 0 ) odfi_ens(j) = 0d0
         END DO
         odfi_ens = odfi_ens/SUM(odfi_ens) 

      END IF

      !!! 4th increment
   
      CALL deriv(tid,m,dotodf,dotodf_ens,odfi,odfi_ens,dotacs,dotacs_ens,acsi,acsi_ens)

      !!! Main phase
      kodf4 = dotodf*dt_strain*epsnot(tid)
      kac4 = dotacs*dt_strain*epsnot(tid)

      acs(:,:,:,m) = acs(:,:,:,m) + (kac1/2d0+kac2+kac3+kac4/2d0)/3d0
      odf(m,:) = odf(m,:) + (kodf1/2d0+kodf2+kodf3+kodf4/2d0)/3d0
      
      DO j = 1 , size ; DO j1 = 1 , 3 ; DO j2 = 1, 3
         IF (acs(j,j1,j2,m) .GT. 1d0) acs(j,j1,j2,m) = 1d0
         IF (acs(j,j1,j2,m) .LT. -1d0) acs(j,j1,j2,m) = -1d0
      END DO ; END DO ; END DO

      odf(m,:) = odf(m,:)/SUM(odf(m,:))

      !!! Minor phase
      IF(rocktype(m) == 1) THEN 

         kodf4_ens = dotodf_ens*dt_strain*epsnot(tid)
         kac4_ens = dotacs_ens*dt_strain*epsnot(tid)

         acs_ens(:,:,:,m) = acs_ens(:,:,:,m) + (kac1_ens/2d0+kac2_ens+kac3_ens+kac4_ens/2d0)/3d0
         odf_ens(m,:) = odf_ens(m,:) + (kodf1_ens/2d0+kodf2_ens+kodf3_ens+kodf4_ens/2d0)/3d0

         DO j = 1 , size ; DO j1 = 1 , 3 ; DO j2 = 1, 3
            IF (acs_ens(j,j1,j2,m) .GT. 1d0) acs_ens(j,j1,j2,m) = 1d0
            IF (acs_ens(j,j1,j2,m) .LT. -1d0) acs_ens(j,j1,j2,m) = -1d0
         END DO ; END DO ; END DO

         odf_ens(m,:) = odf_ens(m,:)/SUM(odf_ens(m,:))

      END IF

   END IF
!!! End Strain-induced LPO when dislocation creep is active           

!!! Fluid Deformation Rotation when creep mechanisms other than dislocation creep are active
   IF(fractdisl < 1d0) THEN

      !Skip undeformed FSE
      IF(fse0(1,1) == 1d0 .AND. fse0(2,2) == 1d0 .AND. fse0(3,3) == 1d0) GOTO 50

      !Find eigenvalues and eigenvectors of new fse for fluid rotation
      CALL eigen(m,fse0,fseacs0,evals0,1)
      IF(rocktype(m) > 100) GOTO 60
 
      !Timestep for diffusion deformation
      dt_straindiff = dt_straintot*(1d0 -fractdisl)

      !!! Calculate new fse for fluid rotation
      fsei = fse0
      kfse1 = MATMUL(l(tid,:,:),fsei)*dt_straindiff
      fsei = fse0 + 0.5d0*kfse1
      kfse2 = MATMUL(l(tid,:,:),fsei)*dt_straindiff
      fsei = fse0 + 0.5d0*kfse2
      kfse3 = MATMUL(l(tid,:,:),fsei)*dt_straindiff
      fsei = fse0 + kfse3
      kfse4 = MATMUL(l(tid,:,:),fsei)*dt_straindiff
      fsediff = fse0 + (kfse1/2d0+kfse2+kfse3+kfse4/2d0)/3d0

      !Find eigenvalues and eigenvectors of new fse for fluid rotation
      CALL eigen(m,fsediff,fseacs,evals,1)
      IF(rocktype(m) > 100) GOTO 60
 
      !Cosine direction matrix between pair of eigenvectors of fseacs0 and fseacs 
      Q = MATMUL(TRANSPOSE(fseacs0),fseacs)

      !Check if semiaxis is too far, then choose its other half
      DO j = 1, 3 
         IF(Q(j,j) < 0d0) THEN
            fseacs(:,j) = -fseacs(:,j)
         END IF
      END DO
         
      !Check if too much rotation, which can be caused by swapping of semiaxes
      !after deformation (an ellipsoidal fse which becomes spherical and then
      !ellipsoidal again)
      Q = MATMUL(fseacs0,TRANSPOSE(fseacs))
      IF(Q(1,1)<0.5d0 .OR. Q(2,2)<0.5d0 .OR. Q(3,3)<0.5d0) GOTO 50

     !Rotate crystals
     DO j = 1 , size  
        !Q = f1*f2' --> First rotate forward with f1, then rotate backward with f2
        acs(j,:,:,m) = MATMUL(acs(j,:,:,m),Q)
        DO j1 = 1 , 3 ; DO j2 = 1, 3
         IF (acs(j,j1,j2,m) .GT. 1d0) acs(j,j1,j2,m) = 1d0
         IF (acs(j,j1,j2,m) .LT. -1d0) acs(j,j1,j2,m) = -1d0
        END DO ; END DO
     END DO

     !!! Minor phase
     IF(rocktype(m) == 1) THEN

        DO j = 1 , size  
           !Q = f1*f2' --> First rotate forward with f1, then rotate backward with f2
           acs_ens(j,:,:,m) = MATMUL(acs_ens(j,:,:,m),Q)
           DO j1 = 1 , 3 ; DO j2 = 1, 3
            IF (acs_ens(j,j1,j2,m) .GT. 1d0) acs_ens(j,j1,j2,m) = 1d0
            IF (acs_ens(j,j1,j2,m) .LT. -1d0) acs_ens(j,j1,j2,m) = -1d0
           END DO ; END DO
        END DO

     END IF

50 END IF
!!! End Fluid Deformation Rotation

   END IF
!!! End Compute fabric for two-phase mantle aggregates, except lower transition zone as ringwoodite and majoritic garnet are almost isotropic
!!! End Compute fabric for both phases for upper mantle aggregates   

   fse0 = fse

   END DO
!!! End Main strain cycle

   IF(fsemod == 0 .AND. (isnan(acs(1,1,1,m)) .OR. isnan(odf(m,1)) .OR. isnan(acs_ens(1,1,1,m)) .OR. isnan(odf_ens(m,1)))) THEN
      print *,"acs or odf Nan at",m,rocktype(m),mx1(m),mx2(m),acs(1,1,1,m),odf(m,1)
      acs(:,:,:,m) = acs0
      acs_ens(:,:,:,m) = acs0
      odf(m,:) = 1d0/REAL(size3**3)
      odf_ens(m,:) = odf(m,:)
      fse = 0d0; fse(1,1) = 1d0; fse(2,2) = 1d0 ; fse(3,3) = 1d0
   END IF

!!! Update aggregate fse   
60 Fij(:,:,m) = fse
   
   RETURN

   END SUBROUTINE strain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE eigen(m,fse,evects,evals,stretchmod)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: j,m,ti(1),stretchmod ! counters

   DOUBLE PRECISION, DIMENSION(3,3) :: fse,evects,dum,V2
   DOUBLE PRECISION, DIMENSION(3) :: evals2,evals
   
   !Left stretch tensor
   IF(stretchmod == 1) V2 = MATMUL(fse,TRANSPOSE(fse))
   !Right stretch tensor
   IF(stretchmod == 2) V2 = MATMUL(TRANSPOSE(fse),fse)
   CALL DSYEVQ3(V2,dum,evals2)

   IF(MINVAL(evals2)<0d0 .OR. MAXVAL(evals2)>1e+8) THEN

     !print *,'Too big FSE for marker ',m,' X = ',mx1(m),' Y = ',mx2(m)
     !print *,'V2'
     !write(*,'(3e17.8)'),V2
     !print *,'fse'
     !write(*,'(3e17.8)'),fse
     !print *,'evals2'
     !write(*,'(3e17.8)'),evals2
     !print *,'evects'
     !write(*,'(3e17.8)'),evects

     !!! Reset LPO and FSE of aggregate
     acs(:,:,:,m) = acs0
     acs_ens(:,:,:,m) = acs0
     odf(m,:) = 1d0/REAL(size3**3)
     odf_ens(m,:) = odf(m,:)
     fse = 0d0; fse(1,1) = 1d0; fse(2,2) = 1d0 ; fse(3,3) = 1d0

     !Remove aggregate
     !rocktype(m) = rocktype(m) + 100

     RETURN

   END IF

   !Sort semiaxes from smallest to biggest
   DO j = 1,3
      ti = MINLOC(evals2) 
      evects(:,j) = dum(:,ti(1))
      evals(j) = evals2(ti(1))**0.5d0
      evals2(ti(1))= 1d60
   END DO

   RETURN

   END SUBROUTINE eigen 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine DERIV, calculation of the rotation vector and slip rate     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE deriv(tid,m,dotodf,dotodf_ens,odfi,odfi_ens,dotacs,dotacs_ens,acsi,acsi_ens)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: i,i1,i2,i3,i4,j,k,tid,m,maxsp ! counters
   INTEGER , DIMENSION(1) :: ti ! reordering array
   INTEGER , DIMENSION(12) :: idx ! list of indices

   DOUBLE PRECISION, DIMENSION(size) :: dotodf,dotodf_ens,odfi,odfi_ens,rt,rt_ens
   DOUBLE PRECISION, DIMENSION(size,3,3) :: dotacs,dotacs_ens,acsi,acsi_ens
   DOUBLE PRECISION :: Emean,rt1,rt2,rt3,rt4,rt5,Emean_ens,rt0_ens
   !! surface averaged aggregate NRJ
   !! dislocation density for each slip system

   DOUBLE PRECISION :: gam0,alpha_x,alpha_z
   ! slip rate on the softest slip system

   DOUBLE PRECISION :: R1,R2
   DOUBLE PRECISION :: sn1,rat,rti
   !!! dummies

   DOUBLE PRECISION, DIMENSION(12) :: bigi,q,qab,qi ! intermediates for G calc

   DOUBLE PRECISION, DIMENSION(12) :: gam
   ! ratios of strain between softest slip system and slip system s for Olivine

   DOUBLE PRECISION, DIMENSION(3) :: rot
   !! rotation rate vector

   DOUBLE PRECISION, DIMENSION(3,3) :: g
   ! slip tensor

   DOUBLE PRECISION, DIMENSION(3,3) :: lx,ex
   ! dimensionless velocity gradient and strain rate tensors

   DOUBLE PRECISION, DIMENSION(3,3) :: acsnsp,acssd,acssd0
   DOUBLE PRECISION, DIMENSION(3,3) :: acsnsp011,acsnsp021,acsnsp101,acsnsp110,acsnsp10_1,acsnsp_110
   DOUBLE PRECISION, DIMENSION(3,3) :: acssd110,acssd111,acssd11_1,acssd1_11,acssd1_1_1,acssd1_10
   DOUBLE PRECISION, DIMENSION(3,3) :: acssd_2110,acssd1_210,acssd11_20
   ! rotated direction cosine matrix for normal to slip plane and slip direction

!!! Dimensionless strain rate and velocity gradient tensors
   lx = l(tid,:,:)/epsnot(tid) ; ex = e(tid,:,:)/epsnot(tid)

   rt = 0d0 ; rt_ens = rt
   
!!! Plastic deformation + dynamic recrystallization

   DO i=1,size

!!! Calculate invariants e_{pr} T_{pr} for the four slip systems of olivine

   bigi=0d0 ; gam = 0d0 ; g = 0d0 ; idx = 0d0 ; qi = 0d0

!!! Olivine
!!! orthorombic
IF(rocktype(m) == 1) THEN
   DO i1 = 1,3 ; DO i2 = 1,3
      ![100](010)
      bigi(1) = bigi(1)+ex(i1,i2)*acsi(i,1,i1)*acsi(i,2,i2)
      ![100](001)
      bigi(2) = bigi(2)+ex(i1,i2)*acsi(i,1,i1)*acsi(i,3,i2)
      ![001](010)
      bigi(3) = bigi(3)+ex(i1,i2)*acsi(i,3,i1)*acsi(i,2,i2)
      ![001](100)
      bigi(4) = bigi(4)+ex(i1,i2)*acsi(i,3,i1)*acsi(i,1,i2)
   ENDDO ; ENDDO

   maxsp = 4

!!! Wadsleyite
!!! a = 5.67 | b = 11.58 | c = 8.26 |
!!! orthorombic
ELSE IF(rocktype(m) == 2) THEN

   alpha_x = 63.9119d0
   alpha_z = 32.6447d0

   !Find direction of normal to slip plane (011): rotation about the x-axis of 35.5
   !degrees, nsp orientation is [001]
   acsnsp011 = 0d0
   CALL rot3D(acsi(i,:,:),acsnsp011,acsi(i,1,:),35.5d0)


   !Find direction of normal to slip plane (021): rotation about the x-axis of 19.6287
   !degrees, nsp orientation is [001]
   acsnsp021 = 0d0
   CALL rot3D(acsi(i,:,:),acsnsp021,acsi(i,1,:),19.6287d0)


   !Find direction of normal to slip plane (101): rotation about the y-axis of -55.5
   !degrees, nsp orientation is [001]
   acsnsp101 = 0d0
   CALL rot3D(acsi(i,:,:),acsnsp101,acsi(i,2,:), -55.5d0)

   !Find direction of normal to slip plane (10_1): rotation about the y-axis of 55.5
   !degrees, nsp orientation is [001]
   acsnsp10_1 = 0d0
   CALL rot3D(acsi(i,:,:),acsnsp10_1,acsi(i,2,:),55.5d0)

   !Rotate slip direction into [111], sd orientation is [100]
   acssd110 = 0d0
   acssd111 = 0d0
   CALL rot3D(acsi(i,:,:),acssd110,acsi(i,3,:),-alpha_x)
   CALL rot3D(acssd110,acssd111,acssd110(2,:),alpha_z)

   !Rotate slip direction into [11_1], sd orientation is [100]
   acssd11_1 = 0d0
   CALL rot3D(acssd110,acssd11_1,acssd110(2,:),-alpha_z)

   !Rotate slip direction into [1_11], sd orientation is [100]
   acssd1_10 = 0d0
   acssd1_11 = 0d0
   CALL rot3D(acsi(i,:,:),acssd1_10,acsi(i,3,:),alpha_x)
   CALL rot3D(acssd1_10,acssd1_11,acssd1_10(2,:),alpha_z)

   !Rotate slip direction into [1_1_1], sd orientation is [100]
   acssd1_1_1 = 0d0
   CALL rot3D(acssd1_10,acssd1_1_1,acssd1_10(2,:),-alpha_z)

   DO i1 = 1,3 ; DO i2 = 1,3
      !!! [100](001)
      bigi(1) = bigi(1)+ex(i1,i2)*acsi(i,1,i1)*acsi(i,3,i2)
      !!! [100](010)
      bigi(2) = bigi(2)+ex(i1,i2)*acsi(i,1,i1)*acsi(i,2,i2)
      !!! [100](011)
      bigi(3) = bigi(3)+ex(i1,i2)*acsi(i,1,i1)*acsnsp011(3,i2)
      !!! [100](021)
      bigi(4) = bigi(4)+ex(i1,i2)*acsi(i,1,i1)*acsnsp021(3,i2)
      !!! [111](10-1)
      bigi(5) = bigi(5)+ex(i1,i2)*acssd111(1,i1)*acsnsp10_1(3,i2)
      !!! [11_1](101)
      bigi(6) = bigi(6)+ex(i1,i2)*acssd11_1(1,i1)*acsnsp101(3,i2)
      !!! [1_11](10-1)
      bigi(7) = bigi(7)+ex(i1,i2)*acssd1_11(1,i1)*acsnsp10_1(3,i2)
      !!! [1_1_1](101)
      bigi(8) = bigi(8)+ex(i1,i2)*acssd1_1_1(1,i1)*acsnsp101(3,i2)
      !!! [001](010)
      bigi(9) = bigi(9)+ex(i1,i2)*acsi(i,3,i1)*acsi(i,2,i2)
   ENDDO ; ENDDO

   maxsp = 9

!!! Bridgmanite
!!! a = 4.675 | b = 4.814 | c = 6.703 (at 38 GPa, 2500 K)
!!! orthorombic
ELSE IF(rocktype(m) == 4) THEN

   alpha_x = 45.8392d0

   !Find direction of normal to slip plane (_110): rotation about the z-axis of 45
   !degrees, nsp orientation is [100]
   acsnsp_110 = 0d0
   CALL rot3D(acsi(i,:,:),acsnsp_110,acsi(i,3,:),alpha_x)

   !Find direction of normal to slip plane (110): rotation about the z-axis of 45
   !degrees, nsp orientation is [100]
   acsnsp110 = 0d0
   CALL rot3D(acsi(i,:,:),acsnsp110,acsi(i,3,:),-alpha_x)

   !!! From Mainprice et al., 2008. EPSL
   DO i1 = 1,3 ; DO i2 = 1,3
      !!! [100](010)
      bigi(1) = bigi(1)+ex(i1,i2)*acsi(i,1,i1)*acsi(i,2,i2)
      !!! [100](001)
      bigi(2) = bigi(2)+ex(i1,i2)*acsi(i,1,i1)*acsi(i,3,i2)
      !!! [010](100)
      bigi(3) = bigi(3)+ex(i1,i2)*acsi(i,2,i1)*acsi(i,1,i2)
      !!! [010](001)
      bigi(4) = bigi(4)+ex(i1,i2)*acsi(i,2,i1)*acsi(i,3,i2)
      !!! [001](100)
      bigi(5) = bigi(5)+ex(i1,i2)*acsi(i,3,i1)*acsi(i,1,i2)
      !!! [001](010)
      bigi(6) = bigi(6)+ex(i1,i2)*acsi(i,3,i1)*acsi(i,2,i2)
      !!! [001](110)
      bigi(7) = bigi(7)+ex(i1,i2)*acsi(i,3,i1)*acsnsp110(1,i2)
      !!! [001](-110)
      bigi(8) = bigi(8)+ex(i1,i2)*acsi(i,3,i1)*acsnsp_110(1,i2)
      !!! [110](001)
      bigi(9) = bigi(9)+ex(i1,i2)*acsnsp110(1,i1)*acsi(i,3,i2)
      !!! [-110](001)
      bigi(10) = bigi(10)+ex(i1,i2)*acsnsp_110(1,i1)*acsi(i,3,i2)
      !!! [110](-110)
      bigi(11) = bigi(11)+ex(i1,i2)*acsnsp110(1,i1)*acsnsp_110(1,i2)
      !!! [-110](110)
      bigi(12) = bigi(12)+ex(i1,i2)*acsnsp_110(1,i1)*acsnsp110(1,i2)
   ENDDO ; ENDDO

   maxsp = 12

END IF

!!! Quotients I/tau
   q = bigi/tau(rocktype(m),:)

!!! Reorder quotients I/tau according to absolute magnitude
   qab = ABS(q)
   DO j=1,maxsp     
      ti = MAXLOC(qab) ; idx(j) = ti(1) ; qab(idx(j))=-1d0
   END DO

!!! Calculate weighting factors gam_s relative to value gam_i for which
!!! I/tau is largest

   gam(idx(1))=1d0

   rat = tau(rocktype(m),idx(1))/bigi(idx(1))
   DO j=2,maxsp    
      qi(idx(j)) = rat*bigi(idx(j))/tau(rocktype(m),idx(j))
   END DO
   sn1 = stressexp(rocktype(m))-1d0

   DO j=2,maxsp    
      gam(idx(j))=qi(idx(j))*(abs(qi(idx(j))))**sn1
   END DO

!!! calculation of G tensor
IF(rocktype(m) == 1) THEN

   DO i1 = 1,3 ; DO i2 = 1,3
      g(i1,i2)=2d0*(gam(1)*acsi(i,1,i1)*acsi(i,2,i2) + &
                    gam(2)*acsi(i,1,i1)*acsi(i,3,i2) + &
                    gam(3)*acsi(i,3,i1)*acsi(i,2,i2) + &
                    gam(4)*acsi(i,3,i1)*acsi(i,1,i2))
   END DO ; END DO

ELSE IF(rocktype(m) == 2) THEN

   DO i1 = 1,3 ; DO i2 = 1,3
      g(i1,i2)=2d0*(gam(1)*acsi(i,1,i1)*acsi(i,3,i2) + & !!! [100](001)
                    gam(2)*acsi(i,1,i1)*acsi(i,2,i2) + & !!! [100](010)
                    gam(3)*acsi(i,1,i1)*acsnsp011(3,i2) + & !!! [100](011)
                    gam(4)*acsi(i,1,i1)*acsnsp021(3,i2) + & !!! [100](021)
                    gam(5)*acssd111(1,i1)*acsnsp10_1(3,i2) +&  !!! [111](10_1)
                    gam(6)*acssd11_1(1,i1)*acsnsp101(3,i2) +&  !!! [11_1](101)
                    gam(7)*acssd1_11(1,i1)*acsnsp10_1(3,i2) +& !!! [1_11](10_1)
                    gam(8)*acssd1_1_1(1,i1)*acsnsp101(3,i2) +&   !!! [1_1_1](101)
                    gam(9)*acsi(i,3,i1)*acsi(i,2,i2))   !!! [001](010)
   END DO ; END DO

ELSE IF(rocktype(m) == 4) THEN

   DO i1 = 1,3 ; DO i2 = 1,3
      g(i1,i2)=2d0*(gam(1)*acsi(i,1,i1)*acsi(i,2,i2) + &          !!! [100](010)
                    gam(2)*acsi(i,1,i1)*acsi(i,3,i2) + &          !!! [100](001)
                    gam(3)*acsi(i,2,i1)*acsi(i,1,i2) + &          !!! [010](100)
                    gam(4)*acsi(i,2,i1)*acsi(i,3,i2) + &          !!! [010](001)
                    gam(5)*acsi(i,3,i1)*acsi(i,1,i2) + &          !!! [001](100)
                    gam(6)*acsi(i,3,i1)*acsi(i,2,i2) + &          !!! [001](010)
                    gam(7)*acsi(i,3,i1)*acsnsp110(1,i2) + &       !!! [001](110)
                    gam(8)*acsi(i,3,i1)*acsnsp_110(1,i2) + &      !!! [001](-110)
                    gam(9)*acsnsp110(1,i1)*acsi(i,3,i2) + &       !!! [110](001)
                    gam(10)*acsnsp_110(1,i1)*acsi(i,3,i2) + &     !!! [-110](001)
                    gam(11)*acsnsp110(1,i1)*acsnsp_110(1,i2) + &  !!! [110](-110)
                    gam(12)*acsnsp_110(1,i1)*acsnsp110(1,i2))     !!! [-110](110)
   END DO ; END DO

END IF

!!! calculation of strain rate on the softest slip system

   R1 = 0d0 ; R2 = 0d0

   DO j= 1 , 3
      i2 = j + 2
      IF (i2 > 3) i2 = i2 - 3

      R1 = R1 - (g(j,i2)-g(i2,j))*(g(j,i2)-g(i2,j))
      R2 = R2 - (g(j,i2)-g(i2,j))*(lx(j,i2)-lx(i2,j))

      DO k = 1 , 3

         R1 = R1 + 2d0*g(j,k)*g(j,k)
         R2 = R2 + 2d0*lx(j,k)*g(j,k)

      END DO
   END DO

   gam0 = R2/R1
   
!!! dislocation density calculation

   DO j=1,maxsp

      rti=tau(rocktype(m),j)**(1.5d0-stressexp(rocktype(m)))*ABS(gam(j)*gam0)**(1.5d0/stressexp(rocktype(m)))

      rt(i) = rt(i) + rti*exp(-lambda(rocktype(m))*rti**2)

   END DO

!!! calculation of the rotation rate:
   
   rot(3) = (lx(2,1)-lx(1,2))/2d0-(g(2,1)-g(1,2))/2d0*gam0

   rot(2) = (lx(1,3)-lx(3,1))/2d0-(g(1,3)-g(3,1))/2d0*gam0
   
   rot(1) = (lx(3,2)-lx(2,3))/2d0-(g(3,2)-g(2,3))/2d0*gam0

!!! derivative of the matrix of direction cosine

   dotacs(i,:,:) = 0d0 

   DO i1 = 1 , 3 ; DO i2 = 1 , 3 ; DO i3 = 1 , 3 ; DO i4 = 1 , 3
      dotacs(i,i1,i2)=dotacs(i,i1,i2)+alt(i2,i3,i4)*acsi(i,i1,i4)*rot(i3)
   END DO ; END DO ; END DO ; END DO
   
!!! grain boundary sliding for small grains
   IF (odfi(i) < chi(rocktype(m))/REAL(size)) THEN
      dotacs(i,:,:) = 0d0
      rt(i) = 0d0
   END IF

   END DO

!!! Volume averaged energy
   Emean = SUM(odfi*rt)
   
!!! Change of volume fraction by grain boundary migration
   DO i = 1 , size
      dotodf(i) = Xol(rocktype(m)) * Mob(rocktype(m)) * odfi(i) * (Emean-rt(i))
   END DO

!!! ENSTATITE
!!! Minor phase only for upper mantle (ol + ens)
   IF(rocktype(m) == 1) THEN

   DO i=1,size

!!! Calculate slip tensor

   g = 0d0

   DO i1 = 1,3 ; DO i2 = 1,3

      g(i1,i2) = 2d0*acsi_ens(i,3,i1)*acsi_ens(i,1,i2)

   ENDDO ; ENDDO

!!! calculation of strain rate

   R1 = 0d0 ; R2 = 0d0

   DO j= 1 , 3
      i2 = j + 2
      IF (i2 > 3) i2 = i2 - 3

      R1 = R1 - (g(j,i2)-g(i2,j))*(g(j,i2)-g(i2,j))
      R2 = R2 - (g(j,i2)-g(i2,j))*(lx(j,i2)-lx(i2,j))

      DO k = 1 , 3

         R1 = R1 + 2d0*g(j,k)*g(j,k)
         R2 = R2 + 2d0*lx(j,k)*g(j,k)

      END DO
   END DO

   gam0 = R2/R1

! weight factor between olivine and enstatite

   gam0 = gam0*(1d0/tau(1,5))**stressexp(rocktype(m))

!!! dislocation density calculation

   rt0_ens=tau(1,5)**(1.5d0-stressexp(rocktype(m)))*ABS(gam0)**(1.5d0/stressexp(rocktype(m)))

   rt_ens(i) = rt0_ens*exp(-lambda(rocktype(m))*rt0_ens**2)
   
!!! calculation of the rotation rate: 
   
   rot(3) = (lx(2,1)-lx(1,2))/2d0-(g(2,1)-g(1,2))/2d0*gam0

   rot(2) = (lx(1,3)-lx(3,1))/2d0-(g(1,3)-g(3,1))/2d0*gam0

   rot(1) = (lx(3,2)-lx(2,3))/2d0-(g(3,2)-g(2,3))/2d0*gam0

   dotacs_ens(i,:,:) = 0d0

   DO i1 = 1 , 3 ; DO i2 = 1 , 3 ; DO i3 = 1 , 3 ; DO i4 = 1 , 3
      dotacs_ens(i,i1,i2)=dotacs_ens(i,i1,i2)+alt(i2,i3,i4)*acsi_ens(i,i1,i4)*rot(i3)
   END DO ; END DO ; END DO ; END DO

!!! grain boundary sliding for small grains
   IF (odfi_ens(i) < chi(rocktype(m))/REAL(size)) THEN
      dotacs_ens(i,:,:) = 0d0
      rt_ens(i) = 0d0
   END IF

   END DO

!!! Volume averaged energy
   Emean_ens = SUM(odfi_ens*rt_ens)

!!! Change of volume fraction by grain boundary migration
   DO i = 1 , size
      dotodf_ens(i) = (1d0-Xol(rocktype(m))) * Mob(rocktype(m)) * odfi_ens(i) * (Emean_ens-rt_ens(i))
   END DO

   END IF

   RETURN

   END SUBROUTINE deriv
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine rot 3D - direction cosine matrix 3D rotation around given axis!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rot3D(acsold,acsnew,a,angle0)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: angle,angle0
   DOUBLE PRECISION, DIMENSION(3) :: a
   DOUBLE PRECISION, DIMENSION(3,3) :: acsold,acsnew,acsrot

   !!! angle must be in radians
   angle = angle0/180d0*pi          

   acsnew = 0d0 ; acsrot = 0d0

   !!! Rotations are clockwise when the axis points toward the observer,
   !!! right handed coordinate system. 
   !!! Euler-Rodrigues's formula for rotation of a vector around an axis

   acsrot(1,1)= cos(angle) + (1-cos(angle))*a(1)*a(1)
   acsrot(1,2)= (1-cos(angle))*a(1)*a(2) - sin(angle)*a(3)
   acsrot(1,3)= (1-cos(angle))*a(1)*a(3) + sin(angle)*a(2)
   acsrot(2,1)= (1-cos(angle))*a(2)*a(1) + sin(angle)*a(3)
   acsrot(2,2)= cos(angle) + (1-cos(angle))*a(2)*a(2)
   acsrot(2,3)= (1-cos(angle))*a(2)*a(3) - sin(angle)*a(1)
   acsrot(3,1)= (1-cos(angle))*a(3)*a(1) - sin(angle)*a(2)
   acsrot(3,2)= (1-cos(angle))*a(3)*a(2) + sin(angle)*a(1)
   acsrot(3,3)= cos(angle) + (1-cos(angle))*a(3)*a(3)

   acsnew = MATMUL(acsold,acsrot)

   END SUBROUTINE rot3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
