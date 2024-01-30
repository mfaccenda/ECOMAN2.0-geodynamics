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
!!! subroutine UPPERLEFT2D, calculation of upper left node for 2D model    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE upperleft2D(x10,x20,i1,i2)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: x10,x20

   INTEGER :: i1,i2 

   i1 = 1; i2 = 1

   DO WHILE (X1(i1+1) .LT. x10 .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

   DO WHILE (X2(i2+1) .LT. x20 .AND. i2 .LT. nx2-1) ; i2 = i2+1 ; ENDDO
      
   RETURN

   END SUBROUTINE upperleft2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine VELOCITYCALC2D, calculation of velocity at a given points   !!!
!!! by interpolation method given in Numerical Recipies, Press et al., p96 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE velocitycalc2D(X1s,X2s,U1s,U2s,i1s,i2s)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X2s
   ! x and y coordinates of the point on the streamline

   DOUBLE PRECISION :: U1s,U2s
   ! interpolated velocity at the point

   INTEGER :: i1s,i2s,m1,m2
   ! indices of the UP-LEFT grid point closest to the extrapolation point

   DOUBLE PRECISION :: y1,y2,y3,y4
   ! dummies for interpolation

!                          
!                          
!                 2----4   
!       +Y        |    |   
!       |         |    |  
!       |         |    | 
!       |         1----3
!       |  
!       ----- +X

   IF(basicstag == 1) THEN
   !Vx
   y1 = Ui(1,1,i1s,i2s  ,1) ; y3 = Ui(1,1,i1s+1,i2s  ,1)
   y2 = Ui(1,1,i1s,i2s+1,1) ; y4 = Ui(1,1,i1s+1,i2s+1,1)

   CALL interp2D(X1s,X2s,i1s,i2s,y1,y2,y3,y4,U1s)

   !Vy
   y1 = Ui(1,2,i1s,i2s  ,1) ; y3 = Ui(1,2,i1s+1,i2s  ,1)
   y2 = Ui(1,2,i1s,i2s+1,1) ; y4 = Ui(1,2,i1s+1,i2s+1,1)

   CALL interp2D(X1s,X2s,i1s,i2s,y1,y2,y3,y4,U2s)

   ELSE

   !Vx

   !Modify indexes for Vx
   m1=i1s
   if(m1 .lt. 1) m1=1 ; if(m1 .gt. nx1-1) m1 = nx1-1
   m2=i2s
   if(X2s .lt. (X2(m2)+X2(m2+1))/2d0) m2 = i2s-1
   if(m2 .lt. 1) m2=1 ; if(m2 .gt. nx2-2) m2 = nx2-2

   y1 = Ui(1,1,m1,m2  ,1) ; y3 = Ui(1,1,m1+1,m2  ,1)
   y2 = Ui(1,1,m1,m2+1,1) ; y4 = Ui(1,1,m1+1,m2+1,1)

   CALL interpshift2D(X1s,X2s,m1,m2,0d0,(X2(m2)+X2(m2+1))/2d0,y1,y2,y3,y4,U1s)

   !Vy (should modify when periodic)

   !Modify indexes for Vy
   m1=i1s
   if(X1s .lt. (X1(m1)+X1(m1+1))/2d0) m1 = i1s-1
   if(m1 .lt. 1) m1=1 ; if(m1 .gt. nx1-2) m1 = nx1-2
   m2=i2s
   if(m2 .lt. 1) m2=1 ; if(m2 .gt. nx2-1) m2 = nx2-1

   y1 = Ui(1,2,m1,m2  ,1) ; y3 = Ui(1,2,m1+1,m2  ,1)
   y2 = Ui(1,2,m1,m2+1,1) ; y4 = Ui(1,2,m1+1,m2+1,1)

   CALL interpshift2D(X1s,X2s,m1,m2,(X1(m1)+X1(m1+1))/2d0,0d0,y1,y2,y3,y4,U2s)

   END IF

   RETURN

   END SUBROUTINE velocitycalc2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine interp, bilinear interpolation of a function                !!!
!!! following the method given in Numerical Recipies, Press et al., p96    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE interp2D(X1s,X2s,i1s,i2s,y1,y2,y3,y4,res)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X2s
   ! x and y coordinates of the point on the streamline

   DOUBLE PRECISION :: res
   ! result of the interpolation at the point

   INTEGER :: i1s,i2s
   ! indices of the UP-LEFT grid point closest to the considered point

   DOUBLE PRECISION :: y1,y2,y3,y4,xx,yy
   ! dummies for interpolation (numerical recipies, p96)

!                          
!                          
!                 2----4   
!       +Y        |    |   
!       |         |    |  
!       |         |    | 
!       |         1----3
!       |  
!       ----- +X

   xx = (X1s-X1(i1s))/(X1(i1s+1)-X1(i1s))
   yy = (X2s-X2(i2s))/(X2(i2s+1)-X2(i2s))

   res = (1d0-xx)*(1d0-yy)*y1 &
              +(1d0-xx)*yy*y2 &
              +xx*(1d0-yy)*y3 &
                    +xx*yy*y4

   RETURN

   END SUBROUTINE interp2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine interp, bilinear interpolation of a function                !!!
!!! following the method given in Numerical Recipies, Press et al., p96    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE interpshift2D(X1s,X2s,i1s,i2s,shiftX,shiftY,y1,y2,y3,y4,res)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X2s,shiftX,shiftY
!  DOUBLE PRECISION :: x1c,x2c,x1cc,x2ccc ! unused
   ! x and y coordinates of the point on the streamline

   DOUBLE PRECISION :: res
   ! result of the interpolation at the point

   INTEGER :: i1s,i2s
   ! indices of the UP-LEFT grid point closest to the considered point

   DOUBLE PRECISION :: y1,y2,y3,y4,xx,yy
   ! dummies for interpolation (numerical recipies, p96)

!                          
!                 2----4   
!       +Y        |    |   
!       |         |    |  
!       |         |    | 
!       |         1----3
!       |  
!       ----- +X
!                          

   ! initialization to avoid warning during compilation
   xx = 0.d0
   yy = 0.d0

   if(shiftX .ne. 0) then
      if(x1periodic .ne. 0 .and. i1s== nx1-1) then
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
   
   if(shiftX .eq. 0) xx=(X1s-X1(i1s))/(X1(i1s+1)-X1(i1s))
   if(shiftY .eq. 0) yy=(X2s-X2(i2s))/(X2(i2s+1)-X2(i2s))

   IF(xx .LT. 0d0) xx = 0d0 ; IF (xx .GT. 1d0) xx = 1d0
   IF(yy .LT. 0d0) yy = 0d0 ; IF (yy .GT. 1d0) yy = 1d0

   res = (1d0-xx)*(1d0-yy)*y1 &
              +(1d0-xx)*yy*y2 &  
              +xx*(1d0-yy)*y3 &
                    +xx*yy*y4  

   RETURN

   END SUBROUTINE interpshift2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine GRADIENTCALC, interpolation of velocity gradient tensor     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE gradientcalc2D(tid,X1s,X2s,i1,i2)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X2s
   ! x and y coordinates of the point on the streamline

   INTEGER :: i1s,i2s,i1,i2
!  INTEGER :: m  ! unused
   ! indices of the UPPER-LEFT grid point closest to the considered point

   DOUBLE PRECISION :: y1,y2,y3,y4
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

!                          
!                 2----4   
!       +Y        |    |   
!       |         |    |  
!       |         |    | 
!       |         1----3
!       |  
!       ----- +X
!                          

   IF(basicstag == 1) THEN
   
   !dVx/dx
   y1 = Dij(1,1,1,i1,i2  ,1) ;     y3 = Dij(1,1,1,i1+1,i2  ,1)
   y2 = Dij(1,1,1,i1,i2+1,1) ;     y4 = Dij(1,1,1,i1+1,i2+1,1)

   CALL interp2D(X1s,X2s,i1,i2,y1,y2,y3,y4,l(tid,1,1))

   !dVx/dy
   y1 = Dij(1,1,2,i1,i2  ,1) ;     y3 = Dij(1,1,2,i1+1,i2  ,1)
   y2 = Dij(1,1,2,i1,i2+1,1) ;     y4 = Dij(1,1,2,i1+1,i2+1,1)

   CALL interp2D(X1s,X2s,i1,i2,y1,y2,y3,y4,l(tid,1,2))

   !dVy/dx
   y1 = Dij(1,2,1,i1,i2  ,1) ;     y3 = Dij(1,2,1,i1+1,i2  ,1)
   y2 = Dij(1,2,1,i1,i2+1,1) ;     y4 = Dij(1,2,1,i1+1,i2+1,1)

   CALL interp2D(X1s,X2s,i1,i2,y1,y2,y3,y4,l(tid,2,1))

   !dVy/dy
   y1 = Dij(1,2,2,i1,i2  ,1) ;     y3 = Dij(1,2,2,i1+1,i2  ,1)
   y2 = Dij(1,2,2,i1,i2+1,1) ;     y4 = Dij(1,2,2,i1+1,i2+1,1)

   CALL interp2D(X1s,X2s,i1,i2,y1,y2,y3,y4,l(tid,2,2))

   ELSE

   !dVx/dx at -stepx1/2,-stepx2/2

   !Modify indexes
   i1s=i1
   if(X1s .gt. (X1(i1s)+X1(i1s+1))/2d0) i1s=i1+1
   if(i1s .lt. 2) i1s=2 ; if(i1s>nx1-1) i1s=nx1-1
   i2s=i2
   if(X2s .gt. (X2(i2s)+X2(i2s+1))/2d0) i2s=i2+1
   if(i2s .lt. 2) i2s=2 ; if(i2s>nx2-1) i2s=nx2-1

   y1 = Dij(1,1,1,i1s,i2s  ,1) ;     y3 = Dij(1,1,1,i1s+1,i2s  ,1)
   y2 = Dij(1,1,1,i1s,i2s+1,1) ;     y4 = Dij(1,1,1,i1s+1,i2s+1,1)

   CALL interpshift2D(X1s,X2s,i1s-1,i2s-1,(X1(i1s)+X1(i1s-1))/2d0,(X2(i2s)+X2(i2s-1))/2d0,y1,y2,y3,y4,l(tid,1,1))

   !dVy/dy at -stepx1/2,-stepx2/2

   y1 = Dij(1,2,2,i1s,i2s  ,1) ;     y3 = Dij(1,2,2,i1s+1,i2s  ,1)
   y2 = Dij(1,2,2,i1s,i2s+1,1) ;     y4 = Dij(1,2,2,i1s+1,i2s+1,1)

   CALL interpshift2D(X1s,X2s,i1s-1,i2s-1,(X1(i1s)+X1(i1s-1))/2d0,(X2(i2s)+X2(i2s-1))/2d0,y1,y2,y3,y4,l(tid,2,2))

   !dVx/dy at 0,0

   !Modify indexes
   i1s=i1
   if(i1s .lt. 1) i1s=1 ; if(i1s>nx1-1) i1s=nx1-1
   i2s=i2
   if(i2s .lt. 2) i2s=2 ; if(i2s>nx2-2) i2s=nx2-2

   y1 = Dij(1,1,2,i1s,i2s  ,1) ;     y3 = Dij(1,1,2,i1s+1,i2s  ,1)
   y2 = Dij(1,1,2,i1s,i2s+1,1) ;     y4 = Dij(1,1,2,i1s+1,i2s+1,1)

   CALL interpshift2D(X1s,X2s,i1s-1,i2s-1,0d0,0d0,y1,y2,y3,y4,l(tid,1,2))

   !dVy/dx at 0,0

   i1s=i1
   if(i1s .lt. 2) i1s=2 ; if(i1s>nx1-2) i1s=nx1-2
   i2s=i2
   if(i2s .lt. 1) i2s=1 ; if(i2s>nx2-1) i2s=nx2-1

   y1 = Dij(1,2,1,i1s,i2s  ,1) ;     y3 = Dij(1,2,1,i1s+1,i2s  ,1)
   y2 = Dij(1,2,1,i1s,i2s+1,1) ;     y4 = Dij(1,2,1,i1s+1,i2s+1,1)

   CALL interpshift2D(X1s,X2s,i1s-1,i2s-1,0d0,0d0,y1,y2,y3,y4,l(tid,2,1))

   END IF

!!! strain rate tensor
   e(tid,1,1) = l(tid,1,1) ; e(tid,2,2) = l(tid,2,2)
   e(tid,1,2) = (l(tid,2,1)+l(tid,1,2))/2d0 ; e(tid,2,1) = e(tid,1,2)

!! reference strain rate
   ee = e(tid,:,:)
   CALL DSYEVQ3(ee,evects,evals)
   epsnot(tid) = MAXVAL(ABS(evals))

   RETURN
   
   END SUBROUTINE gradientcalc2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine PATHLINE - Calculation of tracers pathlines                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE advection2D(m)

   USE comvar 
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: m,i1,i2
   ! indices of the closest upper-left grid point

   DOUBLE PRECISION :: X1i,X2i,U1i,U2i
   ! intermediate velocity components on the streamline

   DOUBLE PRECISION :: kx1,kx2,kx3,kx4,dx
   DOUBLE PRECISION :: ky1,ky2,ky3,ky4,dy
   ! RGK intermediate for streamlines calculation 

   ! value of the dimensionless strain for which any former
   ! LPO has been erased. From Kaminski and Ribe 2002, 10
   ! is a safe value. Should be tested for a given flow though.


!!! X PERIODIC MUST BE TAKEN INTO ACCOUNT

   X1i = mx1(m) ; X2i = mx2(m) 

   CALL upperleft2D(X1i,X2i,i1,i2)
   CALL velocitycalc2D(X1i,X2i,U1i,U2i,i1,i2)

   kx1 = U1i*dt ; ky1 = U2i*dt 

   X1i = mx1(m) + 0.5d0*kx1
   X2i = mx2(m) + 0.5d0*ky1

!!! Check if marker position out of domain 
   !!! Check for periodicity along x1
   IF(x1periodic > 0) THEN
      IF(X1i < x1min) X1i = x1max + X1i - x1min
      IF(X1i > x1max) X1i = x1min + X1i - x1max
   ELSE
      IF(X1i < x1min) X1i = x1min ; IF (X1i > x1max) X1i = x1max
   END IF
   !!! Check for periodicity along x2
   IF(x2periodic > 0) THEN
      IF(X2i < x2min) X2i = x2max + X2i - x2min
      IF(X2i > x2max) X2i = x2min - X2i - x2max
   ELSE
      IF(X2i < x2min) X2i = x2min ; IF (X2i > x2max) X2i = x2max
   END IF

!!! find the UPPER-LEFT grid point closest to the point of calculation
      
   CALL upperleft2D(X1i,X2i,i1,i2)
   CALL velocitycalc2D(X1i,X2i,U1i,U2i,i1,i2)

   kx2 = U1i*dt ; ky2 = U2i*dt

   X1i = mx1(m) + 0.5d0*kx2
   X2i = mx2(m) + 0.5d0*ky2

!!! Check if marker position out of domain 
   !!! Check for periodicity along x1
   IF(x1periodic > 0) THEN
      IF(X1i < x1min) X1i = x1max + X1i - x1min
      IF(X1i > x1max) X1i = x1min + X1i - x1max
   ELSE
      IF(X1i < x1min) X1i = x1min ; IF (X1i > x1max) X1i = x1max
   END IF
   !!! Check for periodicity along x2
   IF(x2periodic > 0) THEN
      IF(X2i < x2min) X2i = x2max + X2i - x2min
      IF(X2i > x2max) X2i = x2min - X2i - x2max
   ELSE
      IF(X2i < x2min) X2i = x2min ; IF (X2i > x2max) X2i = x2max
   END IF

!!! find the UPPER-LEFT grid point closest to the point of calculation
      
   CALL upperleft2D(X1i,X2i,i1,i2)
   CALL velocitycalc2D(X1i,X2i,U1i,U2i,i1,i2)

   kx3 = U1i*dt ; ky3 = U2i*dt 

   X1i = mx1(m) + kx3
   X2i = mx2(m) + ky3

!!! Check if marker position out of domain 
   !!! Check for periodicity along x1
   IF(x1periodic > 0) THEN
      IF(X1i < x1min) X1i = x1max + X1i - x1min
      IF(X1i > x1max) X1i = x1min + X1i - x1max
   ELSE
      IF(X1i < x1min) X1i = x1min ; IF (X1i > x1max) X1i = x1max
   END IF
   !!! Check for periodicity along x2
   IF(x2periodic > 0) THEN
      IF(X2i < x2min) X2i = x2max + X2i - x2min
      IF(X2i > x2max) X2i = x2min - X2i - x2max
   ELSE
      IF(X2i < x2min) X2i = x2min ; IF (X2i > x2max) X2i = x2max
   END IF

!!! find the UPPER-LEFT grid point closest to the point of calculation
      
   CALL upperleft2D(X1i,X2i,i1,i2)
   CALL velocitycalc2D(X1i,X2i,U1i,U2i,i1,i2)

   kx4 = U1i*dt ; ky4 = U2i*dt

   dx = (kx1/2d0+kx2+kx3+kx4/2d0)/3d0
   dy = (ky1/2d0+ky2+ky3+ky4/2d0)/3d0

   !Orthogonal motion of external aggregates
   IF(rocktype(m) < 10) THEN
      IF(mx1(m)+dx < x1min .OR. mx1(m)+dx > x1max) dy = 0
      IF(mx2(m)+dy < x2min .OR. mx2(m)+dy > x2max) dx = 0
   END IF
   IF(rocktype(m) > 10) THEN
      IF(mx1(m) < x1min .OR. mx1(m) > x1max) dy = 0
      IF(mx2(m) < x2min .OR. mx2(m) > x2max) dx = 0
   END IF

   mx1(m) = mx1(m) + dx
   mx2(m) = mx2(m) + dy

   !Change rocktype to compute (< 10 ) or not ( > 10) the LPO and FSE for aggregates that, respectively, enter or leave the domain
   IF(rocktype(m) < 10) THEN

      !!! Check for periodicity along x1
      IF(x1periodic > 0) THEN
         IF(mx1(m) < x1min) mx1(m) = x1max + mx1(m) - x1min
         IF(mx1(m) > x1max) mx1(m) = x1min + mx1(m) - x1max
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

   ELSE IF(rocktype(m) > 10 .AND. rocktype(m) < 100 .AND. mx1(m) > x1min .AND.            &
   &       mx1(m) < x1max .AND. mx2(m) > x2min .AND. mx2(m) < x2max) THEN 

      rocktype(m) = rocktype(m) - 10 

   END IF

   RETURN

   END SUBROUTINE advection2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE disldiff2D(m,i1,i2,fractdisl)

   USE comvar
   USE omp_lib

   IMPLICIT NONE
  
   DOUBLE PRECISION fractdisl,y1,y2,y3,y4
!  DOUBLE PRECISION y5,y6,y7,y8 ! unused
   INTEGER m,i1,i2
!  INTEGER n1,n2,i3 ! unused

   
   !!! Fraction of dislocation creep
   y1 = Fd(1,i1  ,i2,1) ;     y2 = Fd(1,i1  ,i2+1,1)
   y3 = Fd(1,i1+1,i2,1) ;     y4 = Fd(1,i1+1,i2+1,1)

   CALL interp2D(mx1(m),mx2(m),i1,i2,y1,y2,y3,y4,fractdisl)
 
   !Scaling factor for deformation accommodated by anisotropic phase
   fractdisl = fractdisl * fractdislrock(rocktype(m))

   END SUBROUTINE disldiff2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rhopt(m,mtk,mpgpa)

   USE comvar
   USE omp_lib

   IMPLICIT NONE
  
   DOUBLE PRECISION mtk,mpa,mpgpa,mpb,ee,n,R0,R1,R2,R3,y1,y2,y3,y4,y5,y6,y7,y8 
   INTEGER m,n1,n2,i1,i2,i3,yy

!                          
!                          
!                 2----4   
!       +Y        |    |   
!       |         |    |  
!       |         |    | 
!       |         1----3
!       |  
!       ----- +X
 
   ! 2D model
   IF(dimensions == 2) THEN

   CALL upperleft2D(mx1(m),mx2(m),i1,i2)
  
   !!! Temperature
   y1 = Tk(1,i1,i2,1) ;       y2 = Tk(1,i1,i2+1,1)
   y3 = Tk(1,i1+1,i2,1) ;     y4 = Tk(1,i1+1,i2+1,1)

   CALL interp2D(mx1(m),mx2(m),i1,i2,y1,y2,y3,y4,mtk)

   !!! Pressure    
   y1 = Pa(1,i1,i2,1) ;       y2 = Pa(1,i1,i2+1,1)
   y3 = Pa(1,i1+1,i2,1) ;     y4 = Pa(1,i1+1,i2+1,1)

   CALL interp2D(mx1(m),mx2(m),i1,i2,y1,y2,y3,y4,mpa)

   ! 3D model
   ELSE

   yy = mYY(m)

   CALL upperleft(mx1(m),mx2(m),mx3(m),i1,i2,i3)

   !!! Temperature
   y1 = Tk(yy,i1,i2,i3) ;       y2 = Tk(yy,i1+1,i2,i3)
   y3 = Tk(yy,i1+1,i2,i3+1) ;   y4 = Tk(yy,i1,i2,i3+1)
   y5 = Tk(yy,i1,i2+1,i3) ;     y6 = Tk(yy,i1+1,i2+1,i3)
   y7 = Tk(yy,i1+1,i2+1,i3+1) ; y8 = Tk(yy,i1,i2+1,i3+1)

   CALL interp(mx1(m),mx2(m),mx3(m),i1,i2,i3,y1,y2,y3,y4,y5,y6,y7,y8,mtk)

   !!! Pressure    
   y1 = Pa(yy,i1,i2,i3) ;       y2 = Pa(yy,i1+1,i2,i3)
   y3 = Pa(yy,i1+1,i2,i3+1) ;   y4 = Pa(yy,i1,i2,i3+1)
   y5 = Pa(yy,i1,i2+1,i3) ;     y6 = Pa(yy,i1+1,i2+1,i3)
   y7 = Pa(yy,i1+1,i2+1,i3+1) ; y8 = Pa(yy,i1,i2+1,i3+1)

   CALL interp(mx1(m),mx2(m),mx3(m),i1,i2,i3,y1,y2,y3,y4,y5,y6,y7,y8,mpa)

   END IF

   IF(mtk<300d0) mtk=300d0
   IF(mpa<1d5) mpa=1d5

   ! Transform pressure in GPa
   mpgpa = mpa / 1d9

   !!! Density

   ! Transform pressure from Pascal to bar
   mpb = mpa / 1d5

   ! ABCD-4Cell Number 
   ee = (mtk-tkmin)/tkstp
   IF(ee < 0) ee=0
   IF(ee > REAL(tknum)) ee=REAL(tknum)
   n=(mpb-pbmin)/pbstp 
   IF(n < 0) n=0
   IF(n > REAL(pbnum)) n=REAL(pbnum)
   n1=FLOOR(ee) + 1
   IF(n1 < 1) n1=1
   IF(n1 > tknum-1) n1=tknum-1
   n2=FLOOR(n)+1
   IF(n2 < 1) n2=1
   IF(n2 > pbnum-1) n2=pbnum-1
   ! Calc normalized distances 
   ee=(ee-REAL(n1-1))
   n=(n-REAL(n2-1))
   ! Ro values
   ! 0 2
   ! 1 3 
   R0=td_rho(n1  ,n2  )
   R1=td_rho(n1  ,n2+1)
   R2=td_rho(n1+1,n2  )
   R3=td_rho(n1+1,n2+1)
   
   rho(m)=((R0*(1.0-n)+R1*n)*(1.0-ee)+(R2*(1.0-n)+R3*n)*ee)

   END SUBROUTINE rhopt           

