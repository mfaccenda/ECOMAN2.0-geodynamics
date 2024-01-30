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
!!! subroutine UPPERLEFT, calculation of upper left node                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE upperleft(x10,x20,x30,i1,i2,i3)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: x10,x20,x30

   INTEGER :: i1,i2,i3 

   i1 = 1; i2 = 1; i3 = 1

   DO WHILE (X1(i1+1) .LT. x10 .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

   DO WHILE (X2(i2+1) .LT. x20 .AND. i2 .LT. nx2-1) ; i2 = i2+1 ; ENDDO
      
   DO WHILE (X3(i3+1) .LT. x30 .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO
                 
   RETURN

   END SUBROUTINE upperleft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine VELOCITYCALC, calculation of velocity at a given point      !!!
!!! by interpolation method given in Numerical Recipies, Press et al., p96 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE velocitycalc(X1s,X2s,X3s,U1s,U2s,U3s,i1s,i2s,i3s,yy)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X2s,X3s,X1dum,X2dum,X3dum,shiftX,shiftY,shiftZ
   ! x and y coordinates of the point on the streamline

   DOUBLE PRECISION :: U1s,U2s,U3s
   ! interpolated velocity at the point

   INTEGER :: i1s,i2s,i3s,m1,m2,m3,ix(2),iy(2),iz(2),yy
   ! indices of the UP-LEFT grid point closest to the extrapolation point

   DOUBLE PRECISION :: y1,y2,y3,y4,y5,y6,y7,y8
   ! dummies for interpolation

!	             8----7
!	            /    /|
!	           /    / |
!	          5----6  |
!       +Y	  |  4-|--3
!	|  +Z     | /  | /
!	|  /      |/   |/
!	| /       1----2
!	|/ 
!	----- +X

   IF(basicstag == 1) THEN

   !Vx
   y1 = Ui(yy,1,i1s,i2s,i3s) ; y2 = Ui(yy,1,i1s+1,i2s,i3s)
   y3 = Ui(yy,1,i1s+1,i2s,i3s+1) ; y4 = Ui(yy,1,i1s,i2s,i3s+1)
   y5 = Ui(yy,1,i1s,i2s+1,i3s) ; y6 = Ui(yy,1,i1s+1,i2s+1,i3s)
   y7 = Ui(yy,1,i1s+1,i2s+1,i3s+1) ; y8 = Ui(yy,1,i1s,i2s+1,i3s+1)

   CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,U1s)

   !Vy
   y1 = Ui(yy,2,i1s,i2s,i3s) ; y2 = Ui(yy,2,i1s+1,i2s,i3s)
   y3 = Ui(yy,2,i1s+1,i2s,i3s+1) ; y4 = Ui(yy,2,i1s,i2s,i3s+1)
   y5 = Ui(yy,2,i1s,i2s+1,i3s) ; y6 = Ui(yy,2,i1s+1,i2s+1,i3s)
   y7 = Ui(yy,2,i1s+1,i2s+1,i3s+1) ; y8 = Ui(yy,2,i1s,i2s+1,i3s+1)

   CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,U2s)

   !Vz
   y1 = Ui(yy,3,i1s,i2s,i3s) ; y2 = Ui(yy,3,i1s+1,i2s,i3s)
   y3 = Ui(yy,3,i1s+1,i2s,i3s+1) ; y4 = Ui(yy,3,i1s,i2s,i3s+1)
   y5 = Ui(yy,3,i1s,i2s+1,i3s) ; y6 = Ui(yy,3,i1s+1,i2s+1,i3s)
   y7 = Ui(yy,3,i1s+1,i2s+1,i3s+1) ; y8 = Ui(yy,3,i1s,i2s+1,i3s+1)

   CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,U3s)

   ELSE

!	             8----7
!	            /    /|
!	           /    / |
!	          5----6  |
!       +Y	  |  4-|--3
!	|  +Z     | /  | /
!	|  /      |/   |/
!	| /       1----2
!	|/ 
!	----- +X

   shiftX=(X1(i1s)+X1(i1s+1))/2d0
   shiftY=(X2(i2s)+X2(i2s+1))/2d0
   shiftZ=(X3(i3s)+X3(i3s+1))/2d0
   !Vx
   X1dum = X1s
   X2dum = X2s
   X3dum = X3s
   !Modify indexes for Vx
   m1 = i1s
   if(m1 < 1) m1 = 1 ; if(m1 > nx1-1) m1 = nx1-1
   ix(1) = m1 ; ix(2) = m1 + 1
 
   m2 = i2s
   if(X2s < shiftY) m2 = i2s-1
   if(x2periodic > 0 .and. (m2 < 1 .or. m2 > nx2-2)) then
      if(m2 < 1) X2dum = x2max + X2s
      iy(1) = nx2-1 ; iy(2) = 1
   else
      if(m2 < 1) m2 = 1 ; if(m2 > nx2-2) m2 = nx2-2
      iy(1) = m2 ; iy(2) = m2+1
   end if
   
   m3 = i3s
   if(X3s < shiftZ) m3 = i3s-1
   if(x3periodic > 0 .and. yinyang == 1 .and. (m3 < 1 .or. m3 > nx3-2)) then
      if(m3 < 1) X3dum = x3max + X3s
      iz(1) = nx3-1 ; iz(2) = 1
   else
      if(m3 < 1) m3 = 1 ; if(m3 > nx3-2) m3 = nx3-2
      iz(1) = m3 ; iz(2) = m3+1
   end if

   y1 = Ui(yy,1,ix(1),iy(1),iz(1)) ; y2 = Ui(yy,1,ix(2),iy(1),iz(1))
   y3 = Ui(yy,1,ix(2),iy(1),iz(2)) ; y4 = Ui(yy,1,ix(1),iy(1),iz(2))
   y5 = Ui(yy,1,ix(1),iy(2),iz(1)) ; y6 = Ui(yy,1,ix(2),iy(2),iz(1))
   y7 = Ui(yy,1,ix(2),iy(2),iz(2)) ; y8 = Ui(yy,1,ix(1),iy(2),iz(2))

   CALL interpshift(X1dum,X2dum,X3dum,m1,m2,m3,0d0,shiftY,shiftZ,y1,y2,y3,y4,y5,y6,y7,y8,U1s)

   !Vy
   X1dum = X1s
   X2dum = X2s
   X3dum = X3s
   !Modify indexes for Vy
   m1=i1s
   if(X1s < shiftX) m1 = i1s-1
   if(x1periodic > 0 .and. yinyang == 1 .and. (m1 < 1 .or. m1 > nx1-2)) then
      if(m1 < 1) X1dum = x1max + X1s
      ix(1) = nx1-1 ; ix(2) = 1
   else
      if(m1 < 1) m1=1 ; if(m1 > nx1-2) m1 = nx1-2
      ix(1) = m1 ; ix(2) = m1 + 1
   end if

   
   m2=i2s
   if(m2 < 1) m2=1 ; if(m2 > nx2-1) m2 = nx2-1
   iy(1) = m2 ; iy(2) = m2+1
   
   m3=i3s
   if(X3s < shiftZ) m3 = i3s-1
   if(x3periodic > 0 .and. yinyang == 1 .and. (m3 < 1 .or. m3 > nx3-2)) then
      if(m3 < 1) X3dum = x3max + X3s
      iz(1) = nx3-1 ; iz(2) = 1
   else
      if(m3 < 1) m3=1 ; if(m3 > nx3-2) m3 = nx3-2
      iz(1) = m3 ; iz(2) = m3+1
   end if


   y1 = Ui(yy,2,ix(1),iy(1),iz(1)) ; y2 = Ui(yy,2,ix(2),iy(1),iz(1))
   y3 = Ui(yy,2,ix(2),iy(1),iz(2)) ; y4 = Ui(yy,2,ix(1),iy(1),iz(2))
   y5 = Ui(yy,2,ix(1),iy(2),iz(1)) ; y6 = Ui(yy,2,ix(2),iy(2),iz(1))
   y7 = Ui(yy,2,ix(2),iy(2),iz(2)) ; y8 = Ui(yy,2,ix(1),iy(2),iz(2))

   CALL interpshift(X1dum,X2dum,X3dum,m1,m2,m3,shiftX,0d0,shiftZ,y1,y2,y3,y4,y5,y6,y7,y8,U2s)

   !Vz
   X1dum = X1s
   X2dum = X2s
   X3dum = X3s
   !Modify indexes for Vz
   m1=i1s
   if(X1s < shiftX) m1 = i1s-1
   if(x1periodic > 0 .and. yinyang == 1 .and. (m1 < 1 .or. m1 > nx1-2)) then
      if(m1 < 1) X1dum = x1max + X1s
      ix(1) = nx1-1 ; ix(2) = 1
   else
      if(m1 < 1) m1 = 1 ; if(m1 > nx1-2) m1 = nx1-2
      ix(1) = m1 ; ix(2) = m1 + 1
   end if

   m2=i2s
   if(X2s < shiftY) m2 = i2s-1
   if(x2periodic .ne. 0 .and. (m2 < 1 .or. m2 > nx2-2)) then
      if(m2 < 1) X2dum = x2max + X2s
      iy(1) = nx2-1 ; iy(2) = 1
   else
      if(m2 < 1) m2=1 ; if(m2 > nx2-2) m2 = nx2-2
      iy(1) = m2 ; iy(2) = m2 + 1
   end if
   
   m3=i3s
   if(m3 < 1) m3=1 ; if(m3 > nx3-1) m3 = nx3-1
   iz(1) = m3 ; iz(2) = m3+1

   y1 = Ui(yy,3,ix(1),iy(1),iz(1)) ; y2 = Ui(yy,3,ix(2),iy(1),iz(1))
   y3 = Ui(yy,3,ix(2),iy(1),iz(2)) ; y4 = Ui(yy,3,ix(1),iy(1),iz(2))
   y5 = Ui(yy,3,ix(1),iy(2),iz(1)) ; y6 = Ui(yy,3,ix(2),iy(2),iz(1))
   y7 = Ui(yy,3,ix(2),iy(2),iz(2)) ; y8 = Ui(yy,3,ix(1),iy(2),iz(2))
   
   CALL interpshift(X1dum,X2dum,X3dum,m1,m2,m3,shiftX,shiftY,0d0,y1,y2,y3,y4,y5,y6,y7,y8,U3s)

   END IF

   RETURN

   END SUBROUTINE velocitycalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine interp, bilinear interpolation of a function                !!!
!!! following the method given in Numerical Recipies, Press et al., p96    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,res)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X2s,X3s
   ! x and y coordinates of the point on the streamline

   DOUBLE PRECISION :: res
   ! result of the interpolation at the point

   INTEGER :: i1s,i2s,i3s
   ! indices of the UP-LEFT grid point closest to the considered point

   DOUBLE PRECISION :: y1,y2,y3,y4,y5,y6,y7,y8,xx,yy,zz
   ! dummies for interpolation (numerical recipies, p96)

!	             8----7
!	            /    /|
!	           /    / |
!	          5----6  |
!       +Y	  |  4-|--3
!	|  +Z     | /  | /
!	|  /      |/   |/
!	| /       1----2
!	|/ 
!	----- +X

   xx = (X1s-X1(i1s))/(X1(i1s+1)-X1(i1s))
   yy = (X2s-X2(i2s))/(X2(i2s+1)-X2(i2s))
   zz = (X3s-X3(i3s))/(X3(i3s+1)-X3(i3s))

   res = (1d0-xx)*(1d0-yy)*(1d0-zz)*y1 &
              +xx*(1d0-yy)*(1d0-zz)*y2 &
                    +xx*(1d0-yy)*zz*y3 &
              +(1d0-xx)*(1d0-yy)*zz*y4 &
              +(1d0-xx)*yy*(1d0-zz)*y5 &
                    +xx*yy*(1d0-zz)*y6 &
                          +xx*yy*zz*y7 &
                    +(1d0-xx)*yy*zz*y8

   RETURN

   END SUBROUTINE interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine interp, bilinear interpolation of a function                !!!
!!! following the method given in Numerical Recipies, Press et al., p96    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE interpshift(X1s,X2s,X3s,i1s,i2s,i3s,shiftX,shiftY,shiftZ,y1,y2,y3,y4,y5,y6,y7,y8,res)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X2s,X3s,shiftX,shiftY,shiftZ
!  DOUBLE PRECISION :: x1c,x2c,x3c,x1cc,x2cc,x3cc ! unused
   ! x and y coordinates of the point on the streamline

   DOUBLE PRECISION :: res
   ! result of the interpolation at the point

   INTEGER :: i1s,i2s,i3s
   ! indices of the UP-LEFT grid point closest to the considered point

   DOUBLE PRECISION :: y1,y2,y3,y4,y5,y6,y7,y8,xx,yy,zz
   ! dummies for interpolation (numerical recipies, p96)

!                    8----7
!                   /    /|
!                  /    / |
!                 5----6  |
!       +Y        |  4-|--3
!       |  +Z     | /  | /
!       |  /      |/   |/
!       | /       1----2
!       |/ 
!       ----- +X

   ! initialization to avoid warning during compilation
   xx = 0.d0
   yy = 0.d0
   zz = 0.d0

   if(shiftX .ne. 0) then
      if(x1periodic > 0 .and. yinyang == 1 .and. i1s== nx1-1) then
         xx=(X1s-shiftX)/((X1(i1s+1)-X1(i1s)+X1(2)-X1(1))/2)
      else
         xx=(X1s-shiftX)/((X1(i1s+2)-X1(i1s))/2)
      end if
   end if
   
   if(shiftY .ne. 0) then
      if(x2periodic .ne. 0 .and. i2s== nx2-1) then
         yy=(X2s-shiftY)/((X2(i2s+1)-X2(i2s)+X2(2)-X2(1))/2)
      else
         yy=(X2s-shiftY)/((X2(i2s+2)-X2(i2s))/2)
      end if
   end if
   
   if(shiftZ .ne. 0) then
      if(x3periodic > 0 .and. yinyang == 1 .and. i3s== nx3-1) then
         zz=(X3s-shiftZ)/((X3(i3s+1)-X3(i3s)+X3(2)-X3(1))/2)
      else
         zz=(X3s-shiftZ)/((X3(i3s+2)-X3(i3s))/2)
      end if
   end if
   
   if(shiftX .eq. 0) xx=(X1s-X1(i1s))/(X1(i1s+1)-X1(i1s))
   if(shiftY .eq. 0) yy=(X2s-X2(i2s))/(X2(i2s+1)-X2(i2s))
   if(shiftZ .eq. 0) zz=(X3s-X3(i3s))/(X3(i3s+1)-X3(i3s))

   IF(xx .LT. 0d0) xx = 0d0 ; IF (xx .GT. 1d0) xx = 1d0
   IF(yy .LT. 0d0) yy = 0d0 ; IF (yy .GT. 1d0) yy = 1d0
   IF(zz .LT. 0d0) zz = 0d0 ; IF (zz .GT. 1d0) zz = 1d0

   res = (1d0-xx)*(1d0-yy)*(1d0-zz)*y1 &
              +xx*(1d0-yy)*(1d0-zz)*y2 &
                    +xx*(1d0-yy)*zz*y3 &
              +(1d0-xx)*(1d0-yy)*zz*y4 &
              +(1d0-xx)*yy*(1d0-zz)*y5 &
                    +xx*yy*(1d0-zz)*y6 &
                          +xx*yy*zz*y7 &
                    +(1d0-xx)*yy*zz*y8

   RETURN

   END SUBROUTINE interpshift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine GRADIENTCALC, interpolation of velocity gradient tensor     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE gradientcalc(tid,X1s,X2s,X3s,i1,i2,i3,yy)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X2s,X3s,X1dum,X2dum,X3dum,shiftX,shiftY,shiftZ
   ! x and y coordinates of the point on the streamline

   INTEGER :: i1s,i2s,i3s,i1,i2,i3,ix(2),iy(2),iz(2),yy
   ! indices of the UPPER-LEFT grid point closest to the considered point

   DOUBLE PRECISION :: y1,y2,y3,y4,y5,y6,y7,y8
   ! dummies for interpolation

   DOUBLE PRECISION, DIMENSION(3,3) :: ee
   ! dummy for reference strain rate calculation

   INTEGER :: tid
!  INTEGER :: nrot ! unused
   ! number of rotations for the Jacobi
   DOUBLE PRECISION, DIMENSION(3) :: evals
   DOUBLE PRECISION, DIMENSION(3,3) :: evects
   ! eigen values and vectors in jacobi

   l(tid,:,:) = 0d0 ; e(tid,:,:) = 0d0

!	             8----7
!	            /    /|
!	           /    / |
!	          5----6  |
!       +Y	  |  4-|--3
!	|  +Z     | /  | /
!	|  /      |/   |/
!	| /       1----2
!	|/ 
!	----- +X
  
   IF(basicstag == 1) THEN

   i1s = i1   
   i2s = i2   
   i3s = i3

   IF(i1s < 2) i1s = 2; IF(i1s > nx1-2) i1s = nx1-2   
   IF(i2s < 2) i2s = 2; IF(i2s > nx2-2) i2s = nx2-2   
   IF(i3s < 2) i3s = 2; IF(i3s > nx3-2) i3s = nx3-2   

   !dVx/dx
   y1 = Dij(yy,1,1,i1s,i2s,i3s) ;       y2 = Dij(yy,1,1,i1s+1,i2s,i3s)
   y3 = Dij(yy,1,1,i1s+1,i2s,i3s+1) ;   y4 = Dij(yy,1,1,i1s,i2s,i3s+1)
   y5 = Dij(yy,1,1,i1s,i2s+1,i3s) ;     y6 = Dij(yy,1,1,i1s+1,i2s+1,i3s)
   y7 = Dij(yy,1,1,i1s+1,i2s+1,i3s+1) ; y8 = Dij(yy,1,1,i1s,i2s+1,i3s+1)

   CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,1,1))

   !dVx/dy
   y1 = Dij(yy,1,2,i1s,i2s,i3s) ;       y2 = Dij(yy,1,2,i1s+1,i2s,i3s)
   y3 = Dij(yy,1,2,i1s+1,i2s,i3s+1) ;   y4 = Dij(yy,1,2,i1s,i2s,i3s+1)
   y5 = Dij(yy,1,2,i1s,i2s+1,i3s) ;     y6 = Dij(yy,1,2,i1s+1,i2s+1,i3s)
   y7 = Dij(yy,1,2,i1s+1,i2s+1,i3s+1) ; y8 = Dij(yy,1,2,i1s,i2s+1,i3s+1)

   CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,1,2))

   !dVx/dz
   y1 = Dij(yy,1,3,i1s,i2s,i3s) ;       y2 = Dij(yy,1,3,i1s+1,i2s,i3s)
   y3 = Dij(yy,1,3,i1s+1,i2s,i3s+1) ;   y4 = Dij(yy,1,3,i1s,i2s,i3s+1)
   y5 = Dij(yy,1,3,i1s,i2s+1,i3s) ;     y6 = Dij(yy,1,3,i1s+1,i2s+1,i3s)
   y7 = Dij(yy,1,3,i1s+1,i2s+1,i3s+1) ; y8 = Dij(yy,1,3,i1s,i2s+1,i3s+1)

   CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,1,3))

   !dVy/dx
   y1 = Dij(yy,2,1,i1s,i2s,i3s) ;       y2 = Dij(yy,2,1,i1s+1,i2s,i3s)
   y3 = Dij(yy,2,1,i1s+1,i2s,i3s+1) ;   y4 = Dij(yy,2,1,i1s,i2s,i3s+1)
   y5 = Dij(yy,2,1,i1s,i2s+1,i3s) ;     y6 = Dij(yy,2,1,i1s+1,i2s+1,i3s)
   y7 = Dij(yy,2,1,i1s+1,i2s+1,i3s+1) ; y8 = Dij(yy,2,1,i1s,i2s+1,i3s+1)

   CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,2,1))

   !dVy/dy
   y1 = Dij(yy,2,2,i1s,i2s,i3s) ;       y2 = Dij(yy,2,2,i1s+1,i2s,i3s)
   y3 = Dij(yy,2,2,i1s+1,i2s,i3s+1) ;   y4 = Dij(yy,2,2,i1s,i2s,i3s+1)
   y5 = Dij(yy,2,2,i1s,i2s+1,i3s) ;     y6 = Dij(yy,2,2,i1s+1,i2s+1,i3s)
   y7 = Dij(yy,2,2,i1s+1,i2s+1,i3s+1) ; y8 = Dij(yy,2,2,i1s,i2s+1,i3s+1)

   CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,2,2))

   !dVy/dz
   y1 = Dij(yy,2,3,i1s,i2s,i3s) ;       y2 = Dij(yy,2,3,i1s+1,i2s,i3s)
   y3 = Dij(yy,2,3,i1s+1,i2s,i3s+1) ;   y4 = Dij(yy,2,3,i1s,i2s,i3s+1)
   y5 = Dij(yy,2,3,i1s,i2s+1,i3s) ;     y6 = Dij(yy,2,3,i1s+1,i2s+1,i3s)
   y7 = Dij(yy,2,3,i1s+1,i2s+1,i3s+1) ; y8 = Dij(yy,2,3,i1s,i2s+1,i3s+1)

   CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,2,3))

   !dVz/dx
   y1 = Dij(yy,3,1,i1s,i2s,i3s) ;       y2 = Dij(yy,3,1,i1s+1,i2s,i3s)
   y3 = Dij(yy,3,1,i1s+1,i2s,i3s+1) ;   y4 = Dij(yy,3,1,i1s,i2s,i3s+1)
   y5 = Dij(yy,3,1,i1s,i2s+1,i3s) ;     y6 = Dij(yy,3,1,i1s+1,i2s+1,i3s)
   y7 = Dij(yy,3,1,i1s+1,i2s+1,i3s+1) ; y8 = Dij(yy,3,1,i1s,i2s+1,i3s+1)

   CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,3,1))

   !dVz/dy
   y1 = Dij(yy,3,2,i1s,i2s,i3s) ;       y2 = Dij(yy,3,2,i1s+1,i2s,i3s)
   y3 = Dij(yy,3,2,i1s+1,i2s,i3s+1) ;   y4 = Dij(yy,3,2,i1s,i2s,i3s+1)
   y5 = Dij(yy,3,2,i1s,i2s+1,i3s) ;     y6 = Dij(yy,3,2,i1s+1,i2s+1,i3s)
   y7 = Dij(yy,3,2,i1s+1,i2s+1,i3s+1) ; y8 = Dij(yy,3,2,i1s,i2s+1,i3s+1)

   CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,3,2))

   !dVz/dz
   y1 = Dij(yy,3,3,i1s,i2s,i3s) ;       y2 = Dij(yy,3,3,i1s+1,i2s,i3s)
   y3 = Dij(yy,3,3,i1s+1,i2s,i3s+1) ;   y4 = Dij(yy,3,3,i1s,i2s,i3s+1)
   y5 = Dij(yy,3,3,i1s,i2s+1,i3s) ;     y6 = Dij(yy,3,3,i1s+1,i2s+1,i3s)
   y7 = Dij(yy,3,3,i1s+1,i2s+1,i3s+1) ; y8 = Dij(yy,3,3,i1s,i2s+1,i3s+1)

   CALL interp(X1s,X2s,X3s,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,3,3))

!!! Incompressibility
   !l(tid,3,3) = -l(tid,1,1)-l(tid,2,2)

   ELSE

   !dVx/dx at -stepx1/2,-stepx2/2,-stepx3/2
   X1dum = X1s
   X2dum = X2s
   X3dum = X3s


   !Modify indexes
   i1s = i1
   if(X1s > (X1(i1s)+X1(i1s+1))/2d0) i1s = i1+1
   if(x1periodic .ne. 0 .and. yinyang == 1 .and. (i1s==1 .or. i1s==nx1)) then
      if(i1s < 2) X1dum = x1max + X1s
      ix(1) = nx1; ix(2) = 2;
   else 
      if(i1s < 2) i1s = 2 ; if(i1s > nx1-1) i1s = nx1-1
      ix(1) = i1s ; ix(2) = i1s+1
   end if 
   
   i2s = i2
   if(X2s > (X2(i2s)+X2(i2s+1))/2d0) i2s = i2+1
   if(x2periodic .ne. 0 .and. (i2s==1 .or. i2s==nx2)) then
      if(i2s < 2) X2dum = x2max + X2s
      iy(1) = nx2; iy(2) = 2; 
   else 
      if(i2s < 2) i2s = 2 ; if(i2s > nx2-1) i2s = nx2-1
      iy(1) = i2s ; iy(2) = i2s+1
   end if 
   
   i3s = i3
   if(X3s > (X3(i3s)+X3(i3s+1))/2d0) i3s = i3+1
   if(x3periodic .ne. 0 .and. yinyang == 1 .and. (i3s==1 .or. i3s==nx3)) then
      if(i3s < 2) X3dum = x3max + X3s
      iz(1) = nx3; iz(2) = 2; 
   else 
      if(i3s < 2) i3s = 2 ; if(i3s > nx3-1) i3s = nx3-1
      iz(1) = i3s ; iz(2) = i3s+1
   end if 

   y1 = Dij(yy,1,1,ix(1),iy(1),iz(1)) ; y2 = Dij(yy,1,1,ix(2),iy(1),iz(1))
   y3 = Dij(yy,1,1,ix(2),iy(1),iz(2)) ; y4 = Dij(yy,1,1,ix(1),iy(1),iz(2))
   y5 = Dij(yy,1,1,ix(1),iy(2),iz(1)) ; y6 = Dij(yy,1,1,ix(2),iy(2),iz(1))
   y7 = Dij(yy,1,1,ix(2),iy(2),iz(2)) ; y8 = Dij(yy,1,1,ix(1),iy(2),iz(2))

   CALL interpshift(X1dum,X2dum,X3dum,i1s-1,i2s-1,i3s-1,(X1(i1s)+X1(i1s-1))/2d0,          &
   &                (X2(i2s)+X2(i2s-1))/2d0,(X3(i3s)+X3(i3s-1))/2d0,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,1,1))

   !dVy/dy at -stepx1/2,-stepx2/2,-stepx3/2

   y1 = Dij(yy,2,2,ix(1),iy(1),iz(1)) ; y2 = Dij(yy,2,2,ix(2),iy(1),iz(1))
   y3 = Dij(yy,2,2,ix(2),iy(1),iz(2)) ; y4 = Dij(yy,2,2,ix(1),iy(1),iz(2))
   y5 = Dij(yy,2,2,ix(1),iy(2),iz(1)) ; y6 = Dij(yy,2,2,ix(2),iy(2),iz(1))
   y7 = Dij(yy,2,2,ix(2),iy(2),iz(2)) ; y8 = Dij(yy,2,2,ix(1),iy(2),iz(2))

   CALL interpshift(X1dum,X2dum,X3dum,i1s-1,i2s-1,i3s-1,(X1(i1s)+X1(i1s-1))/2d0,          &
   &                (X2(i2s)+X2(i2s-1))/2d0,(X3(i3s)+X3(i3s-1))/2d0,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,2,2))

   !dVz/dz
   y1 = Dij(yy,3,3,ix(1),iy(1),iz(1)) ; y2 = Dij(yy,3,3,ix(2),iy(1),iz(1))
   y3 = Dij(yy,3,3,ix(2),iy(1),iz(2)) ; y4 = Dij(yy,3,3,ix(1),iy(1),iz(2))
   y5 = Dij(yy,3,3,ix(1),iy(2),iz(1)) ; y6 = Dij(yy,3,3,ix(2),iy(2),iz(1))
   y7 = Dij(yy,3,3,ix(2),iy(2),iz(2)) ; y8 = Dij(yy,3,3,ix(1),iy(2),iz(2))

   CALL interpshift(X1dum,X2dum,X3dum,i1s-1,i2s-1,i3s-1,(X1(i1s)+X1(i1s-1))/2d0,          &
   &                (X2(i2s)+X2(i2s-1))/2d0,(X3(i3s)+X3(i3s-1))/2d0,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,3,3))

!!! Incompressibility
   !l(tid,3,3) = -l(tid,1,1)-l(tid,2,2)

   !dVx/dy at 0,0,stepx3/2
   X1dum = X1s
   X2dum = X2s
   X3dum = X3s

   shiftZ = (X3(i3)+X3(i3+1))/2d0
   !Modify indexes
   i1s = i1
   if(i1s .lt. 1) i1s = 1 ; if(i1s > nx1-1) i1s = nx1-1
   ix(1) = i1s ; ix(2) = i1s+1
   
   i2s = i2
   if(x2periodic .ne. 0) then
      if(i2s < 1) i2s = 1 ; if(i2s > nx2-1) i2s = nx2-1
   else
      if(i2s < 2) i2s = 2 ; if(i2s > nx2-2) i2s = nx2-2
   end if
   iy(1) = i2s ; iy(2) = i2s+1
   
   i3s = i3
   if(X3s < shiftZ) i3s = i3-1
   if(x3periodic .ne. 0 .and. yinyang == 1 .and. (i3s==1 .or. i3s==nx3-1)) then
      if(i3s < 2) X3dum = x3max + X3s
      iz(1) = nx3-1; iz(2) = 1;
   else 
      if(i3s < 1) i3s = 1 ; if(i3s > nx3-2) i3s = nx3-2
      iz(1) = i3s ; iz(2) = i3s+1
   end if

   y1 = Dij(yy,1,2,ix(1),iy(1),iz(1)) ; y2 = Dij(yy,1,2,ix(2),iy(1),iz(1))
   y3 = Dij(yy,1,2,ix(2),iy(1),iz(2)) ; y4 = Dij(yy,1,2,ix(1),iy(1),iz(2))
   y5 = Dij(yy,1,2,ix(1),iy(2),iz(1)) ; y6 = Dij(yy,1,2,ix(2),iy(2),iz(1))
   y7 = Dij(yy,1,2,ix(2),iy(2),iz(2)) ; y8 = Dij(yy,1,2,ix(1),iy(2),iz(2))

   CALL interpshift(X1dum,X2dum,X3dum,i1s,i2s,i3s,0d0,0d0,shiftZ,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,1,2))

   !dVy/dx at 0,0,stepx3/2

   i1s = i1
   if(x1periodic .ne. 0 .and. yinyang == 1) then
      if(i1s < 1) i1s = 1 ; if(i1s > nx1-1) i1s = nx1-1
   else
      if(i1s < 2) i1s = 2 ; if(i1s > nx1-2) i1s = nx1-2
   end if
   ix(1) = i1s ; ix(2) = i1s+1
   
   i2s = i2
   if(i2s < 1) i2s = 1 ; if(i2s > nx2-1) i2s = nx2-1
   iy(1) = i2s ; iy(2) = i2s+1
   
   y1 = Dij(yy,2,1,ix(1),iy(1),iz(1)) ; y2 = Dij(yy,2,1,ix(2),iy(1),iz(1))
   y3 = Dij(yy,2,1,ix(2),iy(1),iz(2)) ; y4 = Dij(yy,2,1,ix(1),iy(1),iz(2))
   y5 = Dij(yy,2,1,ix(1),iy(2),iz(1)) ; y6 = Dij(yy,2,1,ix(2),iy(2),iz(1))
   y7 = Dij(yy,2,1,ix(2),iy(2),iz(2)) ; y8 = Dij(yy,2,1,ix(1),iy(2),iz(2))

   CALL interpshift(X1dum,X2dum,X3dum,i1s,i2s,i3s,0d0,0d0,shiftZ,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,2,1))

   !dVx/dz at 0,stepx2/2,0
   X1dum = X1s
   X2dum = X2s
   X3dum = X3s
   shiftY = (X2(i2)+X2(i2+1))/2d0
   !Modify indexes
   i1s=i1
   if(i1s < 1) i1s=1 ; if(i1s > nx1-1) i1s = nx1-1
   ix(1) = i1s ; ix(2) = i1s+1
   
   i2s = i2
   if(X2s < shiftY) i2s = i2-1
   if(x2periodic .ne. 0 .and. (i2s==1 .or. i2s==nx2-1)) then
      if(i2s < 2) X2dum = x2max + X2s
      iy(1) = nx2-1; iy(2) = 1;
   else 
      if(i2s < 1) i2s = 1 ; if(i2s > nx2-2) i2s = nx2-2
      iy(1) = i2s ; iy(2) = i2s+1
   end if

   i3s = i3
   if(x3periodic .ne. 0 .and. yinyang == 1) then
      if(i3s < 1) i3s = 1 ; if(i3s > nx3-1) i3s = nx3-1
   else
      if(i3s < 2) i3s = 2 ; if(i3s > nx3-2) i3s = nx3-2
   end if
   iz(1) = i3s ; iz(2) = i3s+1

   y1 = Dij(yy,1,3,ix(1),iy(1),iz(1)) ; y2 = Dij(yy,1,3,ix(2),iy(1),iz(1))
   y3 = Dij(yy,1,3,ix(2),iy(1),iz(2)) ; y4 = Dij(yy,1,3,ix(1),iy(1),iz(2))
   y5 = Dij(yy,1,3,ix(1),iy(2),iz(1)) ; y6 = Dij(yy,1,3,ix(2),iy(2),iz(1))
   y7 = Dij(yy,1,3,ix(2),iy(2),iz(2)) ; y8 = Dij(yy,1,3,ix(1),iy(2),iz(2))

   CALL interpshift(X1dum,X2dum,X3dum,i1s,i2s,i3s,0d0,shiftY,0d0,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,1,3))

   !dVz/dx at 0,stepx2/2,0

   i1s = i1
   if(x1periodic .ne. 0 .and. yinyang == 1) then
      if(i1s < 1) i1s = 1 ; if(i1s > nx1-1) i1s = nx1-1
   else
      if(i1s < 2) i1s = 2 ; if(i1s > nx1-2) i1s = nx1-2
   end if
   ix(1) = i1s ; ix(2) = i1s+1
 
   i3s = i3
   if(i3s < 1) i3s = 1 ; if(i3s > nx3-1) i3s = nx3-1
   iz(1) = i3s ; iz(2) = i3s+1

   y1 = Dij(yy,3,1,ix(1),iy(1),iz(1)) ; y2 = Dij(yy,3,1,ix(2),iy(1),iz(1))
   y3 = Dij(yy,3,1,ix(2),iy(1),iz(2)) ; y4 = Dij(yy,3,1,ix(1),iy(1),iz(2))
   y5 = Dij(yy,3,1,ix(1),iy(2),iz(1)) ; y6 = Dij(yy,3,1,ix(2),iy(2),iz(1))
   y7 = Dij(yy,3,1,ix(2),iy(2),iz(2)) ; y8 = Dij(yy,3,1,ix(1),iy(2),iz(2))

   CALL interpshift(X1dum,X2dum,X3dum,i1s,i2s,i3s,0d0,shiftY,0d0,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,3,1))

   !dVy/dz at stepx1/2,0,0
   X1dum = X1s
   X2dum = X2s
   X3dum = X3s
   shiftX = (X1(i1)+X1(i1+1))/2d0
   !Modify indexes
   i1s=i1
   if(X1s < shiftX) i1s=i1-1
   if(x1periodic .ne. 0 .and. yinyang == 1 .and. (i1s==1 .or. i1s==nx1-1)) then
      if(i1s < 2) X1dum = x1max + X1s
      ix(1) = nx1-1; ix(2) = 1;
   else
      if(i1s .lt. 1) i1s=1 ; if(i1s>nx1-2) i1s=nx1-2
      ix(1) = i1s ; ix(2) = i1s+1
   end if
   
   i2s=i2
   if(i2s .lt. 1) i2s=1 ; if(i2s>nx2-1) i2s=nx2-1
   iy(1)=i2s ; iy(2)=i2s+1
   
   i3s=i3
   if(x3periodic .ne. 0 .and. yinyang == 1) then
      if(i3s .lt. 1) i3s=1 ; if(i3s>nx3-1) i3s=nx3-1
   else
      if(i3s .lt. 2) i3s=2 ; if(i3s>nx3-2) i3s=nx3-2
   end if
   iz(1) = i3s ; iz(2) = i3s+1

   y1 = Dij(yy,2,3,ix(1),iy(1),iz(1)) ; y2 = Dij(yy,2,3,ix(2),iy(1),iz(1))
   y3 = Dij(yy,2,3,ix(2),iy(1),iz(2)) ; y4 = Dij(yy,2,3,ix(1),iy(1),iz(2))
   y5 = Dij(yy,2,3,ix(1),iy(2),iz(1)) ; y6 = Dij(yy,2,3,ix(2),iy(2),iz(1))
   y7 = Dij(yy,2,3,ix(2),iy(2),iz(2)) ; y8 = Dij(yy,2,3,ix(1),iy(2),iz(2))

   CALL interpshift(X1dum,X2dum,X3dum,i1s,i2s,i3s,shiftX,0d0,0d0,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,2,3))

   !dVz/dy at stepx1/2,0,0

   i2s=i2
   if(x2periodic .ne. 0) then
      if(i2s .lt. 1) i2s=1 ; if(i2s>nx2-1) i2s=nx2-1
   else
      if(i2s .lt. 2) i2s=2 ; if(i2s>nx2-2) i2s=nx2-2
   end if
   iy(1)=i2s ; iy(2)=i2s+1

   i3s=i3
   if(i3s .lt. 1) i3s=1 ; if(i3s>nx3-1) i3s=nx3-1
   iz(1) = i3s ; iz(2) = i3s+1

   y1 = Dij(yy,3,2,ix(1),iy(1),iz(1)) ; y2 = Dij(yy,3,2,ix(2),iy(1),iz(1))
   y3 = Dij(yy,3,2,ix(2),iy(1),iz(2)) ; y4 = Dij(yy,3,2,ix(1),iy(1),iz(2))
   y5 = Dij(yy,3,2,ix(1),iy(2),iz(1)) ; y6 = Dij(yy,3,2,ix(2),iy(2),iz(1))
   y7 = Dij(yy,3,2,ix(2),iy(2),iz(2)) ; y8 = Dij(yy,3,2,ix(1),iy(2),iz(2))

   CALL interpshift(X1dum,X2dum,X3dum,i1s,i2s,i3s,shiftX,0d0,0d0,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,3,2))

   END IF

!!! strain rate tensor
   e(tid,1,1) = l(tid,1,1) ; e(tid,3,3) = l(tid,3,3) ; e(tid,2,2) = l(tid,2,2)
   e(tid,3,1) = (l(tid,1,3)+l(tid,3,1))/2d0 ; e(tid,1,3) = e(tid,3,1)
   e(tid,1,2) = (l(tid,2,1)+l(tid,1,2))/2d0 ; e(tid,2,1) = e(tid,1,2)
   e(tid,2,3) = (l(tid,3,2)+l(tid,2,3))/2d0 ; e(tid,3,2) = e(tid,2,3)

!! reference strain rate
   ee = e(tid,:,:)
   if(isnan(ee(1,1))) print *,i1,i2,i3,tid,ee,l(tid,:,:),Dij(yy,:,:,i1,i2,i3)

   CALL DSYEVQ3(ee,evects,evals)
   epsnot(tid) = MAXVAL(ABS(evals))

   RETURN
   
   END SUBROUTINE gradientcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine PATHLINE - Calculation of tracers pathlines                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE advection(m)

   USE comvar 
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: m,i1,i2,i3,yy
   ! indices of the closest upper-left grid point

   DOUBLE PRECISION :: X1i,X2i,X3i,U1i,U2i,U3i,lr,cr
   ! intermediate velocity components on the streamline

   DOUBLE PRECISION :: kx1,kx2,kx3,kx4,dx
   DOUBLE PRECISION :: ky1,ky2,ky3,ky4,dy
   DOUBLE PRECISION :: kz1,kz2,kz3,kz4,dz
   ! RGK intermediate for streamlines calculation 

   ! value of the dimensionless strain for which any former
   ! LPO has been erased. From Kaminski and Ribe 2002, 10
   ! is a safe value. Should be tested for a given flow though.


!!! X PERIODIC MUST BE TAKEN INTO ACCOUNT

   X1i = mx1(m) ; X2i = mx2(m) ; X3i = mx3(m)

   yy = mYY(m)
   CALL upperleft(X1i,X2i,X3i,i1,i2,i3)
   CALL velocitycalc(X1i,X2i,X3i,U1i,U2i,U3i,i1,i2,i3,yy)

   kx1 = U1i*dt ; ky1 = U2i*dt ; kz1 = U3i*dt

   X1i = mx1(m) + 0.5d0*kx1
   X2i = mx2(m) + 0.5d0*ky1
   X3i = mx3(m) + 0.5d0*kz1

!!! Check if marker position out of domain 
   IF(yinyang == 1) THEN
      
   !!! Check for periodicity along x1
   IF(x1periodic > 0) THEN
      IF(cartspher == 2) THEN
         !Check if colat is < 0°
         IF(X1i < x1min) THEN
            X1i = - X1i
            X2i = X2i + pi
         END IF
         !Check if colat is > 90°
         IF(X1i > x1max) THEN
            X1i = 2*pi - X1i
            X2i = X2i + pi
         END IF
      ELSE
         IF(X1i < x1min) X1i = x1max + X1i - x1min
         IF(X1i > x1max) X1i = x1min - X1i - x1max
      END IF
   ELSE
      !IF(X1i < x1min) X1i = x1min ; IF (X1i > x1max) X1i = x1max
   END IF
   !!! Check for periodicity along x2
   IF(x2periodic > 0) THEN
      IF(X2i < x2min) X2i = x2max + X2i - x2min
      IF(X2i > x2max) X2i = x2min + X2i - x2max
   ELSE
      !IF(X2i < x2min) X2i = x2min ; IF (X2i > x2max) X2i = x2max
   END IF
   !!! Check for periodicity along x3
   IF(x3periodic > 0) THEN
      IF(X3i < x3min) X3i = x3max + X3i - x3min
      IF(X3i > x3max) X3i = x3min + X3i - x3max
   ELSE
      !IF(X3i < x3min) X3i = x3min ; IF (X3i > x3max) X3i = x3max
   END IF

   END IF

!!! find the UPPER-LEFT grid point closest to the point of calculation
      
   CALL upperleft(X1i,X2i,X3i,i1,i2,i3)
   CALL velocitycalc(X1i,X2i,X3i,U1i,U2i,U3i,i1,i2,i3,yy)

   kx2 = U1i*dt ; ky2 = U2i*dt ; kz2 = U3i*dt

   X1i = mx1(m) + 0.5d0*kx2
   X2i = mx2(m) + 0.5d0*ky2
   X3i = mx3(m) + 0.5d0*kz2

!!! Check if marker position out of domain 
   IF(yinyang == 1) THEN
      
   !!! Check for periodicity along x1
   IF(x1periodic > 0) THEN
      IF(cartspher == 2) THEN
         !Check if colat is < 0°
         IF(X1i < x1min) THEN
            X1i = - X1i
            X2i = X2i + pi
         END IF
         !Check if colat is > 90°
         IF(X1i > x1max) THEN
            X1i = 2*pi - X1i
            X2i = X2i + pi
         END IF
      ELSE
         IF(X1i < x1min) X1i = x1max + X1i - x1min
         IF(X1i > x1max) X1i = x1min - X1i - x1max
      END IF
   ELSE
      !IF(X1i < x1min) X1i = x1min ; IF (X1i > x1max) X1i = x1max
   END IF
   !!! Check for periodicity along x2
   IF(x2periodic > 0) THEN
      IF(X2i < x2min) X2i = x2max + X2i - x2min
      IF(X2i > x2max) X2i = x2min + X2i - x2max
   ELSE
      !IF(X2i < x2min) X2i = x2min ; IF (X2i > x2max) X2i = x2max
   END IF
   !!! Check for periodicity along x3
   IF(x3periodic > 0) THEN
      IF(X3i < x3min) X3i = x3max + X3i - x3min
      IF(X3i > x3max) X3i = x3min + X3i - x3max
   ELSE
      !IF(X3i < x3min) X3i = x3min ; IF (X3i > x3max) X3i = x3max
   END IF

   END IF

!!! find the UPPER-LEFT grid point closest to the point of calculation
      
   CALL upperleft(X1i,X2i,X3i,i1,i2,i3)
   CALL velocitycalc(X1i,X2i,X3i,U1i,U2i,U3i,i1,i2,i3,yy)

   kx3 = U1i*dt ; ky3 = U2i*dt ; kz3 = U3i*dt

   X1i = mx1(m) + kx3
   X2i = mx2(m) + ky3
   X3i = mx3(m) + kz3

!!! Check if marker position out of domain 
   IF(yinyang == 1) THEN
      
   !!! Check for periodicity along x1
   IF(x1periodic > 0) THEN
      IF(cartspher == 2) THEN
         !Check if colat is < 0°
         IF(X1i < x1min) THEN
            X1i = - X1i
            X2i = X2i + pi
         END IF
         !Check if colat is > 90°
         IF(X1i > x1max) THEN
            X1i = 2*pi - X1i
            X2i = X2i + pi
         END IF
      ELSE
         IF(X1i < x1min) X1i = x1max + X1i - x1min
         IF(X1i > x1max) X1i = x1min - X1i - x1max
      END IF
   ELSE
      !IF(X1i < x1min) X1i = x1min ; IF (X1i > x1max) X1i = x1max
   END IF
   !!! Check for periodicity along x2
   IF(x2periodic > 0) THEN
      IF(X2i < x2min) X2i = x2max + X2i - x2min
      IF(X2i > x2max) X2i = x2min + X2i - x2max
   ELSE
      !IF(X2i < x2min) X2i = x2min ; IF (X2i > x2max) X2i = x2max
   END IF
   !!! Check for periodicity along x3
   IF(x3periodic > 0) THEN
      IF(X3i < x3min) X3i = x3max + X3i - x3min
      IF(X3i > x3max) X3i = x3min + X3i - x3max
   ELSE
      !IF(X3i < x3min) X3i = x3min ; IF (X3i > x3max) X3i = x3max
   END IF

   END IF

!!! find the UPPER-LEFT grid point closest to the point of calculation
      
   CALL upperleft(X1i,X2i,X3i,i1,i2,i3)
   CALL velocitycalc(X1i,X2i,X3i,U1i,U2i,U3i,i1,i2,i3,yy)

   kx4 = U1i*dt ; ky4 = U2i*dt ; kz4 = U3i*dt

   dx = (kx1/2d0+kx2+kx3+kx4/2d0)/3d0
   dy = (ky1/2d0+ky2+ky3+ky4/2d0)/3d0
   dz = (kz1/2d0+kz2+kz3+kz4/2d0)/3d0

   mx1(m) = mx1(m) + dx
   mx2(m) = mx2(m) + dy
   mx3(m) = mx3(m) + dz

   !Change rocktype to compute (< 10 ) or not ( > 10) the LPO and FSE for aggregates that, respectively, enter or leave the domain
   IF(rocktype(m) < 10) THEN

      !!! Check if particle moved to other grid
      IF(yinyang == 2) THEN

         IF(mx1(m) < X1(2) .OR. mx1(m) > X1(nx1-1) .OR. mx3(m) < X3(2) .OR. mx3(m) > X3(nx3-1)) THEN
            CALL yin2yang(mx1(m),mx3(m),lr,cr)
            mx1(m) = lr 
            mx3(m) = cr 
            yy = mYY(m)
            IF(yy == 1) THEN
               mYY(m) = 2
            ELSE
               mYY(m) = 1
            END IF
         END IF

      ELSE

      !!! Check for periodicity along x1
      IF(x1periodic > 0) THEN
         IF(cartspher == 2) THEN
            !Check if colat is < 0°
            IF(mx1(m) < x1min) THEN
               mx1(m) = - mx1(m)
               mx2(m) = mx2(m) + pi
            END IF
            !Check if colat is > 90°
            IF(mx1(m) > x1max) THEN
               mx1(m) = 2*pi - mx1(m)
               mx2(m) = mx2(m) + pi
            END IF
         ELSE
            IF(mx1(m) < x1min) mx1(m) = x1max + mx1(m) - x1min
            IF(mx1(m) > x1max) mx1(m) = x1min + mx1(m) - x1max
         END IF
      ELSE
         IF(mx1(m) < x1min .OR. mx1(m) > x1max) rocktype(m) = rocktype(m) + 10
      END IF

      !!! Check for periodicity along x2
      IF(x2periodic > 0) THEN
         IF(mx2(m) < x2min) mx2(m) = x2max + mx2(m) - x2min
         IF(mx2(m) > x2max) mx2(m) = x2min + mx2(m) - x2max
      ELSE
         IF(mx2(m) < x2min .OR. mx2(m) > x2max) rocktype(m) = rocktype(m) + 10
      END IF

      !!! Check for periodicity along x3
      IF(x3periodic > 0) THEN
         IF(mx3(m) < x3min) mx3(m) = x3max + mx3(m) - x3min
         IF(mx3(m) > x3max) mx3(m) = x3min + mx3(m) - x3max
      ELSE
         IF(mx3(m) < x3min .OR. mx3(m) > x3max) rocktype(m) = rocktype(m) + 10
      END IF

      END IF

   ELSE IF(rocktype(m) > 10 .AND. rocktype(m) < 100 .AND. mx1(m) > x1min .AND.            &
   &       mx1(m) < x1max .AND. mx2(m) > x2min .AND. mx2(m) < x2max .AND. mx3(m) > x3min .AND. mx3(m) < x3max) THEN 

      rocktype(m) = rocktype(m) - 10

   END IF

   RETURN

   END SUBROUTINE advection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE disldiff(m,i1,i2,i3,fractdisl)

   USE comvar
   USE omp_lib

   IMPLICIT NONE
  
   DOUBLE PRECISION fractdisl,y1,y2,y3,y4,y5,y6,y7,y8 
   INTEGER m,i1,i2,i3,yy
!  INTEGER n1,n2 ! unused

   yy = mYY(m)

   !!! Temperature
   y1 = Fd(yy,i1,i2,i3) ;       y2 = Fd(yy,i1+1,i2,i3)
   y3 = Fd(yy,i1+1,i2,i3+1) ;   y4 = Fd(yy,i1,i2,i3+1)
   y5 = Fd(yy,i1,i2+1,i3) ;     y6 = Fd(yy,i1+1,i2+1,i3)
   y7 = Fd(yy,i1+1,i2+1,i3+1) ; y8 = Fd(yy,i1,i2+1,i3+1)

   CALL interp(mx1(m),mx2(m),mx3(m),i1,i2,i3,y1,y2,y3,y4,y5,y6,y7,y8,fractdisl)
 
   !Scaling factor for deformation accommodated by anisotropic phase
   fractdisl = fractdisl * fractdislrock(rocktype(m))

   END SUBROUTINE disldiff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
