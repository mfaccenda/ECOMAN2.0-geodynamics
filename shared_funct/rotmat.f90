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
!!! subroutine rotvoigt, rotate tnesor in voigt notation (6x6)             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rotvoigt(C0,acsrot,C6)

   USE omp_lib

   IMPLICIT NONE

   INTEGER :: i,j
   DOUBLE PRECISION, DIMENSION (3,3) :: acsrot
   DOUBLE PRECISION, DIMENSION (6,6) :: C0,C6,K6    

   !C0  =  input stiffness tensor in Voigt notation (6x6)
   !C6  = output stiffness tensor in Voigt notation (6x6)
   !acsrot = rotation matrix (3x3) 

   !Rotate stiffness tensor in Voigt notation
   !Build 6x6 rotation matrix (Bond, 1943; Auld, 1973)              
   K6 = 0d0

   !Upper left quadrant
   DO j = 1 , 3; DO i = 1, 3
      K6(i,j) = acsrot(i,j)**2d0
   END DO; END DO
   !Upper rigth quadrant
   K6(1,4) = 2d0*acsrot(1,2)*acsrot(1,3); K6(1,5)= 2d0*acsrot(1,3)*acsrot(1,1); K6(1,6)= 2d0*acsrot(1,1)*acsrot(1,2);
   K6(2,4) = 2d0*acsrot(2,2)*acsrot(2,3); K6(2,5)= 2d0*acsrot(2,3)*acsrot(2,1); K6(2,6)= 2d0*acsrot(2,1)*acsrot(2,2);
   K6(3,4) = 2d0*acsrot(3,2)*acsrot(3,3); K6(3,5)= 2d0*acsrot(3,3)*acsrot(3,1); K6(3,6)= 2d0*acsrot(3,1)*acsrot(3,2);
   !Lower left quadrant
   K6(4,1) = acsrot(2,1)*acsrot(3,1); K6(4,2)= acsrot(2,2)*acsrot(3,2); K6(4,3)= acsrot(2,3)*acsrot(3,3);
   K6(5,1) = acsrot(3,1)*acsrot(1,1); K6(5,2)= acsrot(3,2)*acsrot(1,2); K6(5,3)= acsrot(3,3)*acsrot(1,3);
   K6(6,1) = acsrot(1,1)*acsrot(2,1); K6(6,2)= acsrot(1,2)*acsrot(2,2); K6(6,3)= acsrot(1,3)*acsrot(2,3);
   !Lower right quadrant
   K6(4,4) = acsrot(2,2)*acsrot(3,3)+acsrot(2,3)*acsrot(3,2); K6(4,5)= acsrot(2,3)*acsrot(3,1)+acsrot(2,1)*acsrot(3,3); K6(4,6)= acsrot(2,1)*acsrot(3,2)+acsrot(2,2)*acsrot(3,1);
   K6(5,4) = acsrot(3,2)*acsrot(1,3)+acsrot(3,3)*acsrot(1,2); K6(5,5)= acsrot(3,3)*acsrot(1,1)+acsrot(3,1)*acsrot(1,3); K6(5,6)= acsrot(3,1)*acsrot(1,2)+acsrot(3,2)*acsrot(1,1);
   K6(6,4) = acsrot(1,2)*acsrot(2,3)+acsrot(1,3)*acsrot(2,2); K6(6,5)= acsrot(1,3)*acsrot(2,1)+acsrot(1,1)*acsrot(2,3); K6(6,6)= acsrot(1,1)*acsrot(2,2)+acsrot(1,2)*acsrot(2,1);
 
   C6=MATMUL(MATMUL(K6,C0),TRANSPOSE(K6))
 
   END SUBROUTINE rotvoigt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine rotmatrix, build rotation matrix                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rotmatrixZXZ(phi1,theta,phi2,acsrot)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: phi1,theta,phi2
   DOUBLE PRECISION, DIMENSION(3,3) :: acsrot

   !Transform Euler angles into direction cosine matrix
   !Z-X-Z: Euler angles in radians
   acsrot(1,1)=COS(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*SIN(phi2)
   acsrot(2,1)=COS(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*SIN(phi2)
   acsrot(3,1)=SIN(phi2)*SIN(theta)

   acsrot(1,2)=-SIN(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*COS(phi2)
   acsrot(2,2)=-SIN(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*COS(phi2)
   acsrot(3,2)=COS(phi2)*SIN(theta)

   acsrot(1,3)=SIN(theta)*SIN(phi1)
   acsrot(2,3)=-SIN(theta)*COS(phi1)
   acsrot(3,3)=COS(theta)
            
   RETURN

   END SUBROUTINE rotmatrixZXZ      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine rotmatrix, build rotation matrix                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rotmatrixZYZ(phi1,theta,phi2,acsrot)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: phi1,theta,phi2
   DOUBLE PRECISION, DIMENSION(3,3) :: acsrot

   !Transform Euler angles into direction cosine matrix
   !Z-X-Z: Euler angles in radians
   acsrot(1,1)=COS(phi2)*COS(phi1)*COS(theta)-SIN(phi1)*SIN(phi2)
   acsrot(2,1)=COS(phi1)*SIN(phi2)+COS(theta)*COS(phi2)*SIN(phi1)
   acsrot(3,1)=-COS(phi2)*SIN(theta)

   acsrot(1,2)=-SIN(phi1)*COS(phi2)-COS(theta)*SIN(phi2)*COS(phi1)
   acsrot(2,2)=COS(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*SIN(phi2)
   acsrot(3,2)=SIN(phi2)*SIN(theta)

   acsrot(1,3)=SIN(theta)*COS(phi1)
   acsrot(2,3)=-SIN(theta)*SIN(phi1)
   acsrot(3,3)=COS(theta)
            
   RETURN

   END SUBROUTINE rotmatrixZYZ      

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

   !Reverse rotation
   acsnew = MATMUL(acsold,acsrot)

   END SUBROUTINE rot3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine rot 3D - direction cosine matrix 3D rotation around given axis!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rot3Dspo(acsold,acsnew,a,angle)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: angle
   DOUBLE PRECISION, DIMENSION(3) :: a
   DOUBLE PRECISION, DIMENSION(3,3) :: acsold,acsnew,acsrot

   !!! angle must be in radians

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

   acsnew = MATMUL(acsrot,acsold)

   END SUBROUTINE rot3Dspo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
