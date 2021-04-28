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

 SUBROUTINE rotvoigt(C0,acs,C6)

 USE omp_lib

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
