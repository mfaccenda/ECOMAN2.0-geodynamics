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

 !!! Input
 ! x1m = in polar coordinates, longitude (double)
 ! x3m = in polar coordinates, colatitude (double). In 2D cylindrical/annulus, set to pi/2.0 
 ! cartspher = cartesian (1) or polar (2) domain (integer)
 ! yinyang = Yin-Yang grid for global simulations, active when = 2 (integer)
 ! fse0 = particle Fij (3x3 double) in cartesian coordinates
 ! l = particle velocity gradient (3x3 double) in cartesian coordinates
 ! dt = timestep (double)                
 ! etac = viscosity contrast (etaincl/etamatrix) (double). For hard inclusions
 !        (etac>1), allowed values are etac = 10, 100, 1000
 ! volinc = inclusions volume fraction (double). Allowed values are volinc = 0.1, 0.2, 0.3.

 !!! Output
 ! fse  = updated particle Fij (3x3 double) in cartesian coordinates
 ! etavoigt = mapped viscous tensor (6x6 double)

SUBROUTINE viscoustensor(x1m,x3m,cartspher,yinyang,fse0,l,dt,etac,volinc,fse,etavoigt)

USE comvarvis

IMPLICIT NONE

INTEGER :: j,ti(1),cartspher,yinyang ! counters
DOUBLE PRECISION :: x1m,x3m,dt,etac,volinc,r1b,r2b,r1,r2
DOUBLE PRECISION :: phi1,theta,phi2,lr,cr,pi
DOUBLE PRECISION, DIMENSION(3) :: evals2,evals
DOUBLE PRECISION, DIMENSION(3,3) :: l,dum,evects,acs,fse0,fse,fsei
DOUBLE PRECISION, DIMENSION(3,3) :: kfse1,kfse2,kfse3,kfse4
DOUBLE PRECISION, DIMENSION(6,6) :: etavoigt

pi = 3.141592653589793238462643383279

! A) Update fse
fsei  = fse0
kfse1 = MATMUL(l(:,:),fsei)*dt
fsei  = fse0 + 0.5d0*kfse1
kfse2 = MATMUL(l(:,:),fsei)*dt
fsei  = fse0 + 0.5d0*kfse2
kfse3 = MATMUL(l(:,:),fsei)*dt
fsei  = fse0 + kfse3
kfse4 = MATMUL(l(:,:),fsei)*dt
fse   = fse0 + (kfse1/2d0+kfse2+kfse3+kfse4/2d0)/3d0
 
! B) Compute eigenvalues and eigen vectors

!Left stretch tensor
fsei = MATMUL(fse,TRANSPOSE(fse))
CALL DSYEVQ3(fsei,dum,evals2)

!Sort semiaxes from biggest to smallest
DO j = 1,3
   ti             = MAXLOC(evals2) 
   evects(:,j)    = dum(:,ti(1))
   evals(j)       = evals2(ti(1))**0.5d0
   evals2(ti(1))  = -1d60
END DO

!Check for range of fse semiaxes
IF(MINVAL(evals)<0d0 .OR. MAXVAL(evals)>1e+8) THEN

  !!! Reset FSE of aggregate
   fse      = 0d0
   fse(1,1) = 1d0
   fse(2,2) = 1d0
   fse(3,3) = 1d0
  
  CALL interpolate_tensor(etac,volinc,1.0,1.0,etavoigt)

  RETURN

END IF

! B.1) Rotate eigenvector from global to local coordinates in spherical coordinates
! When not active, tensor rotated with respect to cartesian reference frame
IF(1==0 .AND. cartspher == 2) THEN

   !Rotate the viscous tensor to Z axis 
   phi1 = 0d0; theta = -x3m + pi/2.d0 ; phi2 = -x1m !+ pi*1.5d0  
   IF(yinyang == 2) THEN
      CALL yin2yang(x1m,x3m,lr,cr)
      theta = -cr + pi/2.d0 ; phi2 = -lr !+ pi*1.5d0
   END IF

   !Transform Euler angles into direction cosine matrix
   CALL rotmatrixZXZ(phi1,theta,phi2,acs)

   evects = MATMUL(acs,evects)
   !CALL rotvoigt(etavoigt,acs,etavoigt) 

END IF

! C) Find mean inclusions aspect ratio
! Bulk rock aspect ratio
r1b = evals(1)/evals(2)
r2b = evals(2)/evals(3)
CALL fseinclusions(etac,volinc,r1b,r2b,r1,r2)

! D.1) Interpolate DEM tensor
CALL interpolate_tensor(etac,volinc,r1,r2,etavoigt)

! D.2) Limit minimum shear viscosity
CALL viscositycutoff(etac,volinc,etavoigt)

! E) Rotate tensor toward fse semiaxes
CALL rotvoigt(etavoigt,evects,etavoigt) 

RETURN

END SUBROUTINE viscoustensor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE fseinclusions(etac,volinc,r1b,r2b,r1,r2)
 
   IMPLICIT NONE

   DOUBLE PRECISION :: etac, volinc, r1b, r2b, r1, r2, zeta, xi, chi, theta, psi

   if (etac < 1d0) then

      ! -------------------------------------------------------------------------
      ! Parametrisation of average inclusion shape for WEAK inclusions 
      !   (independent of viscosity contrast and volume fraction)
      ! -------------------------------------------------------------------------
      !
      ! r = zeta + xi * A + chi * B 
      !
      ! Where :
      !    (*) r1  = log10(a1 / a2) (INCLUSION FSE)
      !    (*) r2  = log10(a2 / a3) (INCLUSION FSE)
      !    (*) A = r1b = log10(a1b / a2b) (BULK FSE)
      !    (*) B = r2b = log10(a2b / a3b) (BULK FSE)  
      !   
      ! Input :
      !    (*) etac = viscosity contrast (etaincl/etamatrix)
      !    (*) volinc = inclusions volume fraction
      !    (*) A = r1b = log10(a1b / a2b) (BULK FSE)
      !    (*) B = r2b = log10(a2b / a3b) (BULK FSE)  
      !
      ! Output :
      !    (*) r1 = log10(a1 / a2) (INCLUSION FSE)
      !    (*) r1 = log10(a2 / a3) (INCLUSION FSE)
      ! -------------------------------------------------------------------------

      ! -- Calculate: r1 = log10(a1 / a2)
      zeta    =  0.015159 
      xi      =  1.1013
      chi     = -0.093104
      r1      = zeta + xi * r1b + chi * r2b
         
      ! -- Calculate: r2 = log10(a2 / a3)
      zeta    = 0.0028906
      xi      = 1.0533
      chi     = 0.23141
      r2      = zeta + xi * r1b + chi * r2b

   elseif (etac > 1d0) then

      ! -------------------------------------------------------------------------
      ! Parametrisation of average inclusion shape for STRONG inclusions 
      ! -------------------------------------------------------------------------
      !
      ! r = zeta + xi * A + chi * A^2 + theta * B + psi * B^2 
      !
      ! Where :
      !    (*) r1  = log10(a1 / a2) (INCLUSION FSE)
      !    (*) r2  = log10(a2 / a3) (INCLUSION FSE)
      !    (*) A = r1b = log10(a1b / a2b) (BULK FSE)
      !    (*) B = r2b = log10(a2b / a3b) (BULK FSE)  
      !    (*) C = etac = eta_inclusion / eta_matrix (BULK FSE)  
      !   
      ! Input :
      !    (*) etac = viscosity contrast (etaincl/etamatrix)
      !    (*) volinc = inclusions volume fraction
      !    (*) r1b = log10(a1b / a2b) (BULK FSE)
      !    (*) r2b = log10(a2b / a3b) (BULK FSE)  
      !
      ! Output :
      !    (*) r1 = log10(a1 / a2) (INCLUSION FSE)
      !    (*) r1 = log10(a2 / a3) (INCLUSION FSE)
      ! -------------------------------------------------------------------------

      if (etac == 10d0 .and. volinc == 0.1d0) then
         ! -- Calculate: r1 = log10(a1 / a2)
         zeta    = -0.0099044
         xi      = -0.81743
         chi     =  0.44954
         theta   =  1.3274
         psi     = -0.69978

         r1      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b
            
         ! -- Calculate: r2 = log10(a2 / a3)
         zeta    = -0.011368
         xi      = -3.6146
         chi     = 1.5061
         theta   = 4.1358
         psi     = -1.8329

         r2      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b

      elseif (etac == 10d0 .and. volinc == 0.2d0) then
         ! -- Calculate: r1 = log10(a1 / a2)
         zeta    =  0.00060967
         xi      = -0.96877
         chi     =  0.54536
         theta   =  1.3214
         psi     = -0.47325
      
         r1      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b
            
         ! -- Calculate: r2 = log10(a2 / a3)
         zeta    = -0.0024842
         xi      =  0.98501
         chi     = -0.1894
         theta   = -0.56291
         psi     =  0.098311

         r2      =  zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b

      elseif (etac == 10d0 .and. volinc == 0.3d0) then
         ! -- Calculate: r1 = log10(a1 / a2)
         zeta    = -0.0036649
         xi      = -8.0076
         chi     =  4.7426
         theta   =  8.4612
         psi     = -4.6643

         r1      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b 
            
         ! -- Calculate: r2 = log10(a2 / a3)
         zeta    =-0.00044919 
         xi      = 0.69164
         chi     = 0.297
         theta   =-0.32593
         psi     =-0.32645

         r2      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b 

      elseif (etac == 100d0 .and. volinc == 0.1d0) then
         ! -- Calculate: r1 = log10(a1 / a2)
         zeta    =  0.0040817
         xi      =  0.27039
         chi     = -0.12083
         theta   = -0.14853 
         psi     =  0.045963

         r1      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b
            
         ! -- Calculate: r2 = log10(a2 / a3)
         zeta    =  0.0027909
         xi      =  0.23834
         chi     = -0.10328
         theta   = -0.11141
         psi     =  0.037942

         r2      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b

      elseif (etac == 100d0 .and. volinc == 0.2d0) then
         ! -- Calculate: r1 = log10(a1 / a2)
         zeta    = 0.018488
         xi      = 0.71462
         chi     =-0.10733
         theta   =-0.51373
         psi     = 0.13718

         r1      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b
            
         ! -- Calculate: r2 = log10(a2 / a3)  
         zeta    = -0.022818
         xi      =  0.81919
         chi     = -0.34697
         theta   = -0.3002
         psi     =  0.076509

         r2      =  zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b

      elseif (etac == 100d0 .and. volinc == 0.3d0) then
         ! -- Calculate: r1 = log10(a1 / a2)
         zeta    =  0.010126
         xi      =  0.042156
         chi     =  0.10276
         theta   = -0.065527
         psi     =  0.017561

         r1      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b 

         ! -- Calculate: r2 = log10(a2 / a3)
         zeta    =  0.0019115
         xi      =  0.094268
         chi     =  0.042405
         theta   = -0.023625
         psi     =  0.006159

         r2      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b 

      elseif (etac == 1000d0 .and. volinc == 0.1d0) then
         ! -- Calculate: r1 = log10(a1 / a2)
         zeta    = -5.3533d-5
         xi      =  0.0069203
         chi     = -0.0022227
         theta   = -0.00058388
         psi     =  0.00015146

         r1      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b 

         ! -- Calculate: r2 = log10(a2 / a3)
         zeta    =-0.00030846
         xi      = 0.024959
         chi     =-7.2365d-5
         theta   =-0.017061
         psi     = 0.0036481

         r2      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b 

      elseif (etac == 1000d0 .and. volinc == 0.2d0) then
         ! -- Calculate: r1 = log10(a1 / a2)
         zeta    =  0.0026912
         xi      =  0.095508
         chi     = -0.016754
         theta   = -0.030502
         psi     =  0.0088581

         r1      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b
            
         ! -- Calculate: r2 = log10(a2 / a3)  
         zeta    = -0.0012058
         xi      =  0.11378
         chi     = -0.041931 
         theta   = -0.010006
         psi     =  0.00062852

         r2      =  zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b

      elseif (etac == 1000d0 .and. volinc == 0.3d0) then
         ! -- Calculate: r1 = log10(a1 / a2)
         zeta    = -0.00051641
         xi      = -0.082515
         chi     =  0.041083 
         theta   =  0.10888
         psi     = -0.052419

         r1      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b 
            
         ! -- Calculate: r2 = log10(a2 / a3)
         zeta    = -0.00090885
         xi      =  0.00080317
         chi     = -0.082507
         theta   =  0.035072
         psi     =  0.08577

         r2      = zeta + xi*r1b + chi*r1b*r1b + theta*r2b + psi*r2b*r2b 
      endif

   endif



RETURN

END SUBROUTINE fseinclusions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE viscositycutoff(etac,volinc,etavoigt)
 
   IMPLICIT NONE

   INTEGER :: i
   DOUBLE PRECISION :: etac, volinc, r1b, r2b, r1, r2, zeta, xi, chi, theta, psi, lambda, weakening
   DOUBLE PRECISION, DIMENSION (6,6) :: etavoigt  

   ! -------------------------------------------------------------------------
   ! Parametrisation of maximum weakening observed in the models
   ! -------------------------------------------------------------------------
   !
   ! omega = zeta*etac^psi + xi*etac^lambda + chi*etac + theta
   !
   ! Where :
   !    (*) omega := weakening defined as the ratio of
   !                 minimum effective composite and matrix viscosity    
   !   
   ! Input :
   !    (*) etac      = viscosity contrast (etaincl/etamatrix)
   !    (*) volinc   = inclusions volume fraction
   !    (*) etavoigt  = viscous tensor in Voigt notation
   !
   ! Output :   
   !    (*) etavoigt  = limited viscous tensor in Voigt notation 
   ! -------------------------------------------------------------------------

   if (volinc==0.1d0) then 
      zeta     =  0.903221
      xi       =  1.20365
      chi      = -1.148182
      theta    =  0.3735
      psi      =  0.699
      lambda   =  0.6619

   elseif (volinc == 0.2d0) then 
      zeta    =  43.646559
      xi      = -42.907958
      chi     =  0.0
      theta   =  0.261095
      psi     =  0.984129
      lambda  =  1.0

   elseif (volinc == 0.3d0) then 
      zeta    =  40.628095
      xi      = -39.823767
      chi     =  0.0
      theta   =  0.196188
      psi     =  0.985677
      lambda  =  1.0

   else
      ! THROW ERROR: VOLUME FRACTION NOT SUPPORTED
      print *,"ERROR: only phi = 0.1 , 0.2, 0.3 are supported"

   endif

   weakening   = zeta*etac**psi + xi*etac**lambda + chi*etac + theta
   
   do i = 4,6
      etavoigt(i,i) = max(etavoigt(i,i),weakening)
      ! NOTE: second argument is directly 'weakening' IF eta_matrix = 1
   end do

RETURN

END SUBROUTINE viscositycutoff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE interpolate_tensor(etac,volinc,a1a2,a2a3,etavoigt)

USE comvarvis 

IMPLICIT NONE

INTEGER :: nn,vv,aa1,aa2
DOUBLE PRECISION :: etac,volinc,a1a2,a2a3
DOUBLE PRECISION, DIMENSION(6,6) :: etavoigt
!DOUBLE PRECISION :: dnn,dvv,daa1,daa2

!Find nearest etacontr node
nn = etacontrnum
IF(etacontrnum > 1) THEN
   DO WHILE(etacontr(nn) > etac); nn = nn - 1 ; END DO
END IF
IF(nn < 1) nn = 1
!dnn = (etac - etacontr(nn))/(etacontr(nn+1)-etacontr(nn))

!Find nearest volume node
vv = volnum
IF(volnum > 1) THEN
   DO WHILE(Vol(vv) > volinc); vv = vv - 1 ; END DO
END IF
IF(vv < 1) vv = 1
!dvv = (volinc - vol(vv))/(vol(vv+1)-vol(vv))

!Find nearest a1a2 node
aa1 = esanum(1)
IF(esanum(1) > 1) THEN
   DO WHILE(esa12(aa1,1) > a1a2); aa1 = aa1 - 1 ; END DO
END IF
IF(aa1 < 1) aa1 = 1
!daa1 = (a1a2 - esa12(aa1))/(esa12(aa1+1)-esa12(aa1))

!Find nearest a2a3 node
aa2 = esanum(2)
IF(esanum(2) > 1) THEN
   DO WHILE(esa23(1,aa2) > a2a3); aa2 = aa2 - 1 ; END DO
END IF
IF(aa2 < 1) aa2 = 1
!daa2 = (a2a3 - esa23(aa2))/(esa23(aa2+1)-esa23(aa2))

etavoigt(1,1) = Cij(nn,vv,aa1,aa2,1)   
etavoigt(2,2) = Cij(nn,vv,aa1,aa2,2)   
etavoigt(3,3) = Cij(nn,vv,aa1,aa2,3)   
etavoigt(4,4) = Cij(nn,vv,aa1,aa2,4)   
etavoigt(5,5) = Cij(nn,vv,aa1,aa2,5)   
etavoigt(6,6) = Cij(nn,vv,aa1,aa2,6)   
etavoigt(1,2) = Cij(nn,vv,aa1,aa2,7)   
etavoigt(1,3) = Cij(nn,vv,aa1,aa2,8)   
etavoigt(2,3) = Cij(nn,vv,aa1,aa2,9)   
etavoigt(2,1) = etavoigt(1,2)
etavoigt(3,1) = etavoigt(1,3)
etavoigt(3,2) = etavoigt(2,3)

if(0==1) then
print *,etacontr
print *
print *,Vol
print *
print *,esa12(:,1)
print *
print *,esa23(1,:)
print *
print *,etac,volinc,a1a2,a2a3
print *
print *,nn,vv,aa1,aa2
!read(*,*)
end if

RETURN

END SUBROUTINE interpolate_tensor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine rotvoigt, rotate 4th order tensor in Voigt notation         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rotvoigt(C0,acs,C6)

IMPLICIT NONE

INTEGER :: i,j
DOUBLE PRECISION, DIMENSION (3,3) :: acs
DOUBLE PRECISION, DIMENSION (6,6) :: C0,C6,K6    

!C0  =  input stiffness tensor in Voigt notation (6x6)
!C6  = output stiffness tensor in Voigt notation (6x6)
!acs = rotation matrix (3x3) 

!Rotate stiffness tensor in Voigt notation
!Build 6x6 rotation matrix (Bond, 1943; Auld, 1973)              
K6 = 0d0

!Upper left quadrant
DO j = 1 , 3; DO i = 1, 3
   K6(i,j) = acs(i,j)**2d0
END DO; END DO
!Upper rigth quadrant
K6(1,4) = 2d0*acs(1,2)*acs(1,3); K6(1,5)= 2d0*acs(1,3)*acs(1,1); K6(1,6)= 2d0*acs(1,1)*acs(1,2);
K6(2,4) = 2d0*acs(2,2)*acs(2,3); K6(2,5)= 2d0*acs(2,3)*acs(2,1); K6(2,6)= 2d0*acs(2,1)*acs(2,2);
K6(3,4) = 2d0*acs(3,2)*acs(3,3); K6(3,5)= 2d0*acs(3,3)*acs(3,1); K6(3,6)= 2d0*acs(3,1)*acs(3,2);
!Lower left quadrant
K6(4,1) = acs(2,1)*acs(3,1); K6(4,2)= acs(2,2)*acs(3,2); K6(4,3)= acs(2,3)*acs(3,3);
K6(5,1) = acs(3,1)*acs(1,1); K6(5,2)= acs(3,2)*acs(1,2); K6(5,3)= acs(3,3)*acs(1,3);
K6(6,1) = acs(1,1)*acs(2,1); K6(6,2)= acs(1,2)*acs(2,2); K6(6,3)= acs(1,3)*acs(2,3);
!Lower right quadrant
K6(4,4) = acs(2,2)*acs(3,3)+acs(2,3)*acs(3,2); K6(4,5)= acs(2,3)*acs(3,1)+acs(2,1)*acs(3,3); K6(4,6)= acs(2,1)*acs(3,2)+acs(2,2)*acs(3,1);
K6(5,4) = acs(3,2)*acs(1,3)+acs(3,3)*acs(1,2); K6(5,5)= acs(3,3)*acs(1,1)+acs(3,1)*acs(1,3); K6(5,6)= acs(3,1)*acs(1,2)+acs(3,2)*acs(1,1);
K6(6,4) = acs(1,2)*acs(2,3)+acs(1,3)*acs(2,2); K6(6,5)= acs(1,3)*acs(2,1)+acs(1,1)*acs(2,3); K6(6,6)= acs(1,1)*acs(2,2)+acs(1,2)*acs(2,1);

C6=MATMUL(MATMUL(K6,C0),TRANSPOSE(K6))

END SUBROUTINE rotvoigt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine rotvoigt, rotate 4th order tensor in Voigt notation         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rotvoigt2(C0,R,C6)

IMPLICIT NONE

INTEGER :: i,j
DOUBLE PRECISION :: a,b,c,d,e,f,g,h,l,x1,x2,x3,y1,y2,y3,z1,z2,z3
DOUBLE PRECISION, DIMENSION (3,3) :: R
DOUBLE PRECISION, DIMENSION (6,6) :: C0,C6,K6    

!C0  =  input stiffness tensor in Voigt notation (6x6)
!C6  = output stiffness tensor in Voigt notation (6x6)
!R = rotation matrix (3x3) 

! Tensor coefficients
a = C0(1,1)
b = C0(2,2)
c = C0(3,3)
d = C0(4,4)
e = C0(5,5)
f = C0(6,6)
g = C0(1,2)
h = C0(1,3)
l = C0(2,3)
! Rotation matrix
x1 = R(1,1)
x2 = R(1,2)
x3 = R(1,3)
y1 = R(2,1)
y2 = R(2,2)
y3 = R(2,3)
z1 = R(3,1)
z2 = R(3,2)
z3 = R(3,3)

! Rotated components
C6(1,1) = 4 * ((e * x3 **2 + f * x2 **2) * x1 **2 + d * x2 **2 * x3 **2) + &
         (h * x1 **2 + l * x2 **2 + c * x3 **2) * x3 **2 + (g * x2 **2 + h * x3 **2 + a * x1 **2) * x1 **2 + &
         (g * x1 **2 + l * x3 **2 + b * x2 **2) * x2 **2

C6(1,2) = 4 * ((e * x3 * y3 + f * x2 * y2) * x1 * y1 + d * x2 * x3 * y2 * y3) + &
         (h * x1 **2 + l * x2 **2 + c * x3 **2) * y3 **2 + (g * x2 **2 + h * x3 **2 + a * x1 **2) * y1 **2 + &
         (g * x1 **2 + l * x3 **2 + b * x2 **2) * y2 **2
C6(2,1) = C6(1,2)

C6(1,3) = 4 * ((e * x3 * z3 + f * x2 * z2) * x1 * z1 + d * x2 * x3 * z2 * z3) + &
         (h * x1 **2 + l * x2 **2 + c * x3 **2) * z3 **2 + (g * x2 **2 + h * x3 **2 + a * x1 **2) * z1 **2 + &
         (g * x1 **2 + l * x3 **2 + b * x2 **2) * z2 **2
C6(3,1) = C6(1,2)

C6(1,4) = 2 * (((x2 * z3 + x3 * z2) * e * x1 + (y2 * z3 + y3 * z2) * d * x2) * x3 + & 
         (x2 * y3 + x3 * y2) * f * x1 * x2) + (h * x1 **2 + l * x2 **2 + c * x3 **2) * y3 * z3 + &
         (g * x2 **2 + h * x3 **2 + a * x1 **2) * y1 * z1 + (g * x1 **2 + l * x3 **2 + b * x2 **2) * y2 * z2
C6(4,1) = C6(1,2)

C6(1,5) = 2 * (((x1 * z3 + x3 * z1) * e * x1 + (y1 * z3 + y3 * z1) * d * x2) * x3 + &
         (x1 * y3 + x3 * y1) * f * x1 * x2) + (h * x1 **2 + l * x2 **2 + c * x3 **2) * x3 * z3 + &
         (g * x2 **2 + h * x3 **2 + a * x1 **2) * x1 * z1 + (g * x1 **2 + l * x3 **2 + b * x2 **2) * x2 * z2
C6(5,1) = C6(1,2)

C6(1,6) = 2 * (((x1 * z2 + x2 * z1) * e * x1 + (y1 * z2 + y2 * z1) * d * x2) * x3 + & 
         (x1 * y2 + x2 * y1) * f * x1 * x2) + (h * x1 **2 + l * x2 **2 + c * x3 **2) * x3 * y3 + & 
         (g * x2 **2 + h * x3 **2 + a * x1 **2) * x1 * y1 + (g * x1 **2 + l * x3 **2 + b * x2 **2) * x2 * y2
C6(6,1) = C6(1,2)

C6(2,2) = 4 * ((e * y3 **2 + f * y2 **2) * y1 **2 + d * y2 **2 * y3 **2) + &
         (h * y1 **2 + l * y2 **2 + c * y3 **2) * y3 **2 + (g * y2 **2 + h * y3 **2 + a * y1 **2) * y1 **2 + &
         (g * y1 **2 + l * y3 **2 + b * y2 **2) * y2 **2

C6(2,3) = 4 * ((e * y3 * z3 + f * y2 * z2) * y1 * z1 + d * y2 * y3 * z2 * z3) + & 
         (h * y1 **2 + l * y2 **2 + c * y3 **2) * z3 **2 + (g * y2 **2 + h * y3 **2 + a * y1 **2) * z1 **2 + &
         (g * y1 **2 + l * y3 **2 + b * y2 **2) * z2 **2
C6(3,2) = C6(2,3)

C6(2,4) = 2 * (((x2 * z3 + x3 * z2) * e * y1 + (y2 * z3 + y3 * z2) * d * y2) * y3 + & 
         (x2 * y3 + x3 * y2) * f * y1 * y2) + (h * y1 **2 + l * y2 **2 + c * y3 **2) * y3 * z3 + &
         (g * y2 **2 + h * y3 **2 + a * y1 **2) * y1 * z1 + (g * y1 **2 + l * y3 **2 + b * y2 **2) * y2 * z2
C6(4,2) = C6(2,4)

C6(2,5) = 2 * (((x1 * z3 + x3 * z1) * e * y1 + (y1 * z3 + y3 * z1) * d * y2) * y3 + (x1 * y3 + x3 * y1) * f * y1 * y2) + &
         (h * y1 **2 + l * y2 **2 + c * y3 **2) * x3 * z3 + (g * y2 **2 + h * y3 **2 + a * y1 **2) * x1 * z1 + & 
         (g * y1 **2 + l * y3 **2 + b * y2 **2) * x2 * z2
C6(5,2) = C6(2,5)
        
C6(2,6) = 2 * (((x1 * z2 + x2 * z1) * e * y1 + (y1 * z2 + y2 * z1) * d * y2) * y3 + (x1 * y2 + x2 * y1) * f * y1 * y2) + &
         (h * y1 **2 + l * y2 **2 + c * y3 **2) * x3 * y3 + (g * y2 **2 + h * y3 **2 + a * y1 **2) * x1 * y1 + & 
         (g * y1 **2 + l * y3 **2 + b * y2 **2) * x2 * y2
C6(6,2) = C6(2,6)
        
C6(3,3) = 4 * ((e * z3 **2 + f * z2 **2) * z1 **2 + d * z2 **2 * z3 **2) + &
         (h * z1 **2 + l * z2 **2 + c * z3 **2) * z3 **2 + (g * z2 **2 + h * z3 **2 + a * z1 **2) * z1 **2 + &
         (g * z1 **2 + l * z3 **2 + b * z2 **2) * z2 **2
        
C6(3,4) = 2 * (((x2 * z3 + x3 * z2) * e * z1 + (y2 * z3 + y3 * z2) * d * z2) * z3 + (x2 * y3 + x3 * y2) * f * z1 * z2) + &
         (h * z1 **2 + l * z2 **2 + c * z3 **2) * y3 * z3 + (g * z2 **2 + h * z3 **2 + a * z1 **2) * y1 * z1 + &
         (g * z1 **2 + l * z3 **2 + b * z2 **2) * y2 * z2
C6(4,3) = C6(3,4)
        
C6(3,5) = 2 * (((x1 * z3 + x3 * z1) * e * z1 + (y1 * z3 + y3 * z1) * d * z2) * z3 + (x1 * y3 + x3 * y1) * f * z1 * z2) + &
         (h * z1 **2 + l * z2 **2 + c * z3 **2) * x3 * z3 + (g * z2 **2 + h * z3 **2 + a * z1 **2) * x1 * z1 + &
         (g * z1 **2 + l * z3 **2 + b * z2 **2) * x2 * z2
C6(5,3) = C6(3,5)
        
C6(3,6) = 2 * (((x1 * z2 + x2 * z1) * e * z1 + (y1 * z2 + y2 * z1) * d * z2) * z3 + (x1 * y2 + x2 * y1) * f * z1 * z2) + &
         (h * z1 **2 + l * z2 **2 + c * z3 **2) * x3 * y3 + (g * z2 **2 + h * z3 **2 + a * z1 **2) * x1 * y1 + &
         (g * z1 **2 + l * z3 **2 + b * z2 **2) * x2 * y2
C6(6,3) = C6(3,6)
        
C6(4,4) = (x2 * z3 + x3 * z2) **2 * e + (y2 * z3 + y3 * z2) **2 * d + (x2 * y3 + x3 * y2) **2 * f + &
         (h * y1 * z1 + l * y2 * z2 + c * y3 * z3) * y3 * z3 + (g * y2 * z2 + h * y3 * z3 + a * y1 * z1) * y1 * z1 + &
         (g * y1 * z1 + l * y3 * z3 + b * y2 * z2) * y2 * z2
        
C6(4,5) = (x1 * z3 + x3 * z1) * (x2 * z3 + x3 * z2) * e + (y1 * z3 + y3 * z1) * (y2 * z3 + y3 * z2) * d + &
         (x1 * y3 + x3 * y1) * (x2 * y3 + x3 * y2) * f + (h * y1 * z1 + l * y2 * z2 + c * y3 * z3) * x3 * z3 + &
         (g * y2 * z2 + h * y3 * z3 + a * y1 * z1) * x1 * z1 + (g * y1 * z1 + l * y3 * z3 + b * y2 * z2) * x2 * z2
C6(5,4) = C6(4,5)
        
C6(4,6) = (x1 * z2 + x2 * z1) * (x2 * z3 + x3 * z2) * e + (y1 * z2 + y2 * z1) * (y2 * z3 + y3 * z2) * d + &
         (x1 * y2 + x2 * y1) * (x2 * y3 + x3 * y2) * f + (h * y1 * z1 + l * y2 * z2 + c * y3 * z3) * x3 * y3 + &
         (g * y2 * z2 + h * y3 * z3 + a * y1 * z1) * x1 * y1 + (g * y1 * z1 + l * y3 * z3 + b * y2 * z2) * x2 * y2
C6(6,4) = C6(4,6)
        
C6(5,5) = (x1 * z3 + x3 * z1) **2 * e + (y1 * z3 + y3 * z1) **2 * d + (x1 * y3 + x3 * y1) **2 * f + &
         (h * x1 * z1 + l * x2 * z2 + c * x3 * z3) * x3 * z3 + (g * x2 * z2 + h * x3 * z3 + a * x1 * z1) * x1 * z1 + &
         (g * x1 * z1 + l * x3 * z3 + b * x2 * z2) * x2 * z2
        
C6(5,6) = (x1 * z2 + x2 * z1) * (x1 * z3 + x3 * z1) * e + (y1 * z2 + y2 * z1) * (y1 * z3 + y3 * z1) * d + &
         (x1 * y2 + x2 * y1) * (x1 * y3 + x3 * y1) * f + (h * x1 * z1 + l * x2 * z2 + c * x3 * z3) * x3 * y3 + &
         (g * x2 * z2 + h * x3 * z3 + a * x1 * z1) * x1 * y1 + (g * x1 * z1 + l * x3 * z3 + b * x2 * z2) * x2 * y2
C6(6,5) = C6(5,6)
        
C6(6,6) = (x1 * z2 + x2 * z1) **2 * e + (y1 * z2 + y2 * z1) **2 * d + (x1 * y2 + x2 * y1) **2 * f + &
         (h * x1 * y1 + l * x2 * y2 + c * x3 * y3) * x3 * y3 + (g * x2 * y2 + h * x3 * y3 + a * x1 * y1) * x1 * y1 + &
         (g * x1 * y1 + l * x3 * y3 + b * x2 * y2) * x2 * y2
    

END SUBROUTINE rotvoigt2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine rotmatrix, build rotation matrix                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rotmatrixZXZ(phi1,theta,phi2,acs)

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
!!! subroutine yin2yang, find coordinates in other spherical grid          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE yin2yang(lr,cr,lr1,cr1)

IMPLICIT NONE

DOUBLE PRECISION :: lr,cr,lr1,cr1,sinl,cosl,sinc,cosc,cosmc

sinl=sin(lr);
cosl=cos(lr);
sinc=sin(cr);
cosc=cos(cr);
cosmc=sinl*sinc;
if(cosmc<-1.0) cosmc=-1.0;
if(cosmc> 1.0) cosmc= 1.0;
! e = Yin; n = Yang
lr1=atan2(cosc,-sinc*cosl);
cr1=acos(cosmc)

END SUBROUTINE yin2yang

