
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine for elastic decomposition
!  into transverse isotropy
!
!  DECMOD   :   module of global parameters
!  DECSYM   :   decomposition into transverse
!               isotropy symmetry
!
!  Library of independent subroutines
!  for elasticity and tensor
!
!  FULLSYM6 :   upper right 6x6 matrix into
!               full symmetric 6x6 matrix
!  NDELTA   :   Kronecker function
!  TENS4    :   transform 6x6 matrix
!               into 4th order tensor
!  JACOBI   :   get eigenvectors and eigenvalues
!               of a real symmetric matrix
!  EIGSRT   :   order eigenvectors with
!               increasing eigenvalues
!  MAT6     :   transform 4th order tensor
!               into 6x6 matrix
!  PERMUT   :   permutation of index 123
!  ROT4     :   rotation of 4th order tensor 
!
!  Library of subroutines
!  for symmetry decomposition
!  (elastic coefficients in GPa)
!
!  SCCA       : form symmetry cartesian system
!  V21D       : form 21D vector from 6x6 matrix (independent)
!  PROJECTI   : transverse isotropy projector   (independent)
!
!****************************************************************
!

      module DECMOD

      !double precision, dimension(3,3) :: SCC
      !double precision, dimension(6,6) :: CE1
      !double precision, dimension(3,3,3,3) :: EL1
      !double precision, dimension(21) :: XEC
      !double precision :: XN,ANIS

      end module DECMOD

!
!****************************************************************
!

      subroutine DECSYM(Sav1,rt,PERC,TIAXIS,m,ANIS,Dhexa,Dtetra,Dortho,Dmono,Dtri)

      use DECMOD

      implicit none

      double precision, dimension(6,6) :: Sav1
      double precision :: PERC,DC5,PI
      double precision, dimension(3) :: TIAXIS
      double precision, dimension(3,3) :: SCC
      double precision, dimension(6,6) :: CE1
      double precision, dimension(3,3,3,3) :: EL1
      double precision, dimension(21) :: XEC
      double precision, dimension(5) :: DEV
      double precision :: Diso,Dhexa,Dtetra,Dortho,Dmono,Dtri
      double precision :: XN,ANIS
      INTEGER :: m,rt

      XN = 0d0; ANIS=0D0
      Diso=0D0 ; Dhexa = 0d0; Dtetra = 0d0; Dortho = 0d0; Dmono = 0d0; Dtri = 0d0
      PI=acos(-1D0)
      CE1=Sav1
      EL1=0D0
      call FULLSYM6(CE1)
      call TENS4(CE1,EL1)
      call V21D(CE1,XEC)
      XN=sqrt(dot_product(XEC,XEC))

      call SCCA(SCC,CE1,EL1,XEC,XN,ANIS)
      
      call PROJECTI(XEC,DC5)

      PERC=(ANIS-DC5)/XN*100D0

      TIAXIS=SCC(3,:)
      TIAXIS=TIAXIS/sqrt(sum(TIAXIS*TIAXIS))
      !Choosing Z unit vector
      !INCLTI=asin(TIAXIS(3))
      
      !Reset negative hexagonal anisotropy
      IF(PERC < 0d0) then
         PERC = 0d0
         TIAXIS = 0d0
      END IF

      !print *,m,rt,PERC,TIAXIS(1),TIAXIS(2),TIAXIS(3)

      !!! Elastic tensor decomposition
      call DECSYM2(XEC,DEV) 
      ANIS  = ANIS/XN*100d0
      DEV   = DEV/XN*100d0

      !if(ABS(DEV(5)-ANIS)>1d0 .AND. PERC>0)then
      !print *,m,100d0-ANIS,100d0-DEV(5),ANIS-DC5/XN*100D0
      !read(*,*)
      !end if

      Diso  = 100.0d0 - ANIS
      !The isotropic component of the tensor computed in SCCA (from dilational
      !and Voigt tensors) differ fom that computed in PROJECISO
      !Dhexa = DEV(5)-DEV(4)
      IF(PERC > 0d0) THEN
         Dhexa = ANIS  -DEV(4)
         Dtetra= DEV(4)-DEV(3)
         Dortho= DEV(3)-DEV(2)
         Dmono = DEV(2)-DEV(1)
         Dtri  = DEV(1)
      END IF

      Dhexa = Dhexa /ANIS*100d0
      Dtetra= Dtetra/ANIS*100d0
      Dortho= Dortho/ANIS*100d0
      Dmono = Dmono /ANIS*100d0
      Dtri  = Dtri  /ANIS*100d0

      return

      end  subroutine DECSYM

!
!****************************************************************
!

      subroutine FULLSYM6(C)

      implicit none

      double precision, dimension(6,6) :: C

      C(3,2)=C(2,3)
      C(3,1)=C(1,3)
      C(2,1)=C(1,2)

      C(4,1)=C(1,4)
      C(5,1)=C(1,5)
      C(6,1)=C(1,6)
      C(4,2)=C(2,4)
      C(5,2)=C(2,5)
      C(6,2)=C(2,6)
      C(4,3)=C(3,4)
      C(5,3)=C(3,5)
      C(6,3)=C(3,6)

      C(6,5)=C(5,6)
      C(6,4)=C(4,6)
      C(5,4)=C(4,5)

      return

      end subroutine FULLSYM6

!
!****************************************************************
!

      subroutine TENS4(C,C4)

      implicit none

      integer :: i,j,k,l
      integer :: p,q
      integer :: NDELTA
      double precision, dimension(6,6) :: C
      double precision, dimension(3,3,3,3) :: C4

      C4=0D0

      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3

         p=NDELTA(i,j)*i+(1-NDELTA(i,j))*(9-i-j)
         q=NDELTA(k,l)*k+(1-NDELTA(k,l))*(9-k-l)
         C4(i,j,k,l)=C(p,q)

      end do
      end do
      end do
      end do

      end subroutine TENS4

!
!****************************************************************
!

      function NDELTA(i,j)

      implicit none

      integer :: i,j
      integer :: NDELTA

      NDELTA=0
      if (i==j) NDELTA=1

      end function NDELTA

!
!****************************************************************
!

      subroutine EIGSRT(d,v,n,np)

      ! Order eigenvalues and eigenvectors
      ! 1 : max
      ! 2 : mid
      ! 3 : min

      implicit none

      integer :: np,n
      integer :: i,j,k
      double precision, dimension(np) :: d
      double precision, dimension(np,np) :: v
      double precision :: p

      do i=1,n-1
         k=i
         p=d(i)
         do j=i+1,n
            if (d(j)>=p) then
               k=j
               p=d(j)
            end if
         end do
         if (k/=i) then
            d(k)=d(i)
            d(i)=p
            do j=1,n
               p=v(j,i)
               v(j,i)=v(j,k)
               v(j,k)=p
            end do
         end if
      end do

      return

      end subroutine EIGSRT

!
!****************************************************************
!

      subroutine MAT6(C4,C)

      implicit none

      integer :: i
      double precision, dimension(6,6) :: C
      double precision, dimension(3,3,3,3) :: C4

      C = 0D0

      do i=1,3
         C(i,i)=C4(i,i,i,i)
      end do
      do i=2,3
         C(1,i)=(C4(1,1,i,i)+C4(i,i,1,1))/2D0
         C(i,1)=C(1,i)
      end do
      C(2,3)=(C4(2,2,3,3)+C4(3,3,2,2))/2D0
      C(3,2)=C(2,3)

      do i=1,3
         C(i,4)=(C4(i,i,2,3)+C4(i,i,3,2)+                       &
                 C4(2,3,i,i)+C4(3,2,i,i))/4D0
         C(4,i)=C(i,4)
      end do
      do i=1,3
         C(i,5)=(C4(i,i,1,3)+C4(i,i,3,1)+                       &
                 C4(1,3,i,i)+C4(3,1,i,i))/4D0
         C(5,i)=C(i,5)
      end do
      do i=1,3
         C(i,6)=(C4(i,i,1,2)+C4(i,i,2,1)+                       &
                 C4(1,2,i,i)+C4(2,1,i,i))/4D0
         C(6,i)=C(i,6)
      end do

      C(4,4)=(C4(2,3,2,3)+C4(2,3,3,2)+                          &
              C4(3,2,2,3)+C4(3,2,3,2))/4D0
      C(5,5)=(C4(1,3,1,3)+C4(1,3,3,1)+                          &
              C4(3,1,1,3)+C4(3,1,3,1))/4D0
      C(6,6)=(C4(2,1,2,1)+C4(2,1,1,2)+                          &
              C4(1,2,2,1)+C4(1,2,1,2))/4D0
      C(4,5)=(C4(2,3,1,3)+C4(2,3,3,1)+                          &
              C4(3,2,1,3)+C4(3,2,3,1)+                          &
              C4(1,3,2,3)+C4(1,3,3,2)+                          &
              C4(3,1,2,3)+C4(3,1,3,2))/8D0

      C(5,4)=C(4,5)
      C(4,6)=(C4(2,3,1,2)+C4(2,3,2,1)+                          &
              C4(3,2,1,2)+C4(3,2,2,1)+                          &
              C4(1,2,2,3)+C4(1,2,3,2)+                          &
              C4(2,1,2,3)+C4(2,1,3,2))/8D0
      C(6,4)=C(4,6)
      C(5,6)=(C4(1,3,1,2)+C4(1,3,2,1)+                          &
              C4(3,1,1,2)+C4(3,1,2,1)+                          &
              C4(1,2,1,3)+C4(1,2,3,1)+                          &
              C4(2,1,1,3)+C4(2,1,3,1))/8D0
      C(6,5)=C(5,6)

      return

      end subroutine MAT6

!
!****************************************************************
!

      subroutine PERMUT(INDEX,PERM)

      implicit none

      integer :: INDEX
      integer, dimension(3) :: PERM

      if (INDEX==1) then
         PERM(1)=1
         PERM(2)=2
         PERM(3)=3
      end if
      if (INDEX==2) then
         PERM(1)=2
         PERM(2)=3
         PERM(3)=1
      end if
      if (INDEX==3) then 
         PERM(1)=3
         PERM(2)=1
         PERM(3)=2
      endif

      return

      end subroutine PERMUT

!
!****************************************************************
!

      subroutine ROT4(C4,R,C4C)

      implicit none

      integer :: i1,i2,i3,i4,j1,j2,j3,j4
      double precision, dimension(3,3,3,3) :: C4,C4C
      double precision, dimension(3,3) :: R

      C4C = 0D0

      do i1=1,3
      do i2=1,3
      do i3=1,3
      do i4=1,3

         do j1=1,3
         do j2=1,3
         do j3=1,3
         do j4=1,3

            C4C(i1,i2,i3,i4) = C4C(i1,i2,i3,i4) +               &
            R(i1,j1)*R(i2,j2)*R(i3,j3)*R(i4,j4)*C4(j1,j2,j3,j4)

         end do
         end do
         end do
         end do

      end do
      end do
      end do
      end do

      return

      end subroutine ROT4

!
!****************************************************************
!

      subroutine SCCA(SCC,CE1,EL1,XEC,XN,ANIS)

      use DECMOD

      implicit none

      integer :: i,NROT,i1,i2,NDVC
      integer, dimension(3) :: IHS
      double precision, dimension(3) :: EGDI,EGVO
      double precision, dimension(3,3) :: DI,VO,VECDI,VECVO
      double precision, dimension(6,6) :: CEC
      double precision, dimension(3,3,3,3) :: ELC
      double precision, dimension(21) :: XH,XD
      double precision :: SDV,ADV,ADVC,SCN,DEV,K,G
      double precision, dimension(3,3) :: SCC
      double precision, dimension(6,6) :: CE1
      double precision, dimension(3,3,3,3) :: EL1
      double precision, dimension(21) :: XEC
      double precision :: XN,ANIS

      DI=0D0
      VO=0D0
      K=0D0
      G=0D0

      do i=1,3
         DI(1,1)=CE1(1,i)+DI(1,1)
         DI(2,2)=CE1(2,i)+DI(2,2)
         DI(3,3)=CE1(3,i)+DI(3,3)
         DI(2,1)=CE1(6,i)+DI(2,1)
         DI(3,1)=CE1(5,i)+DI(3,1)
         DI(3,2)=CE1(4,i)+DI(3,2)
      end do
      DI(1,2)=DI(2,1)
      DI(1,3)=DI(3,1)
      DI(2,3)=DI(3,2)

      VO(1,1)=CE1(1,1)+CE1(6,6)+CE1(5,5)
      VO(2,2)=CE1(6,6)+CE1(2,2)+CE1(4,4)
      VO(3,3)=CE1(5,5)+CE1(4,4)+CE1(3,3)
      VO(2,1)=CE1(1,6)+CE1(2,6)+CE1(4,5)
      VO(1,2)=VO(2,1)
      VO(3,1)=CE1(1,5)+CE1(3,5)+CE1(4,6)
      VO(1,3)=VO(3,1)
      VO(3,2)=CE1(2,4)+CE1(3,4)+CE1(5,6)
      VO(2,3)=VO(3,2)

      do i=1,3
         K=K+DI(i,i)
         G=G+VO(i,i)
      end do
      K=K/9D0
      G=G/10D0-3D0*K/10D0

      ! Anisotropy

      ANIS=0D0
      XH=0D0
      XD=0D0
      XH(1) = K + 4D0*G/3D0
      XH(2) = K + 4D0*G/3D0
      XH(3) = K + 4D0*G/3D0
      XH(4) = sqrt(2D0)*(K-2D0*G/3D0)
      XH(5) = sqrt(2D0)*(K-2D0*G/3D0)
      XH(6) = sqrt(2D0)*(K-2D0*G/3D0)
      XH(7) = 2D0*G
      XH(8) = 2D0*G
      XH(9) = 2D0*G

      XD=XEC-XH
      ANIS=sqrt(dot_product(XD,XD))

      ! Dil. and Voigt axes

      CALL DSYEVQ3(DI,VECDI,EGDI)
      call EIGSRT(EGDI,VECDI,3,3)
      CALL DSYEVQ3(VO,VECVO,EGVO)
      call EIGSRT(EGVO,VECVO,3,3)

      ! Search for SCCA directions

      do i1=1,3
         NDVC=0
         ADVC=10D0
         SCN=0D0
         do i2=1,3
            SDV=dot_product(VECDI(:,i1),VECVO(:,i2))
            if (abs(SDV)>=1D0) SDV=sign(1D0,SDV)
            ADV=acos(SDV)
            if (SDV<0D0) ADV=acos(-1D0)-ADV
            if (ADV<ADVC) then
               NDVC=sign(1D0,SDV)*i2
               ADVC=ADV
            endif
         end do

         VECDI(:,i1)=(VECDI(:,i1)+NDVC*VECVO(:,abs(NDVC)))/2D0
         SCN=sqrt(VECDI(1,i1)**2D0+VECDI(2,i1)**2D0+            &
                  VECDI(3,i1)**2D0)
         VECDI(:,i1)=VECDI(:,i1)/SCN
      end do

      ! Higher symmetry axis

      SCC = transpose(VECDI)
      ELC = 0D0
      SDV = XN
      do i=1,3
         call PERMUT(i,IHS)
         do i1=1,3
            VECDI(i1,:)=SCC(IHS(i1),:)
         end do
         call ROT4(EL1,VECDI,ELC)
         call MAT6(ELC,CEC)
         call V21D(CEC,XEC)
         call PROJECTI(XEC,DEV)
         if (DEV<SDV) then
             SDV=DEV
             NDVC=i
         endif
      end do

      VECDI=SCC
      call PERMUT(NDVC,IHS)
      do i1=1,3
         SCC(i1,:)=VECDI(IHS(i1),:)
      end do

      ! Rotate in SCCA

      call ROT4(EL1,SCC,ELC)
      call MAT6(ELC,CEC)
      call V21D(CEC,XEC)

      return

      end subroutine SCCA

!
!****************************************************************
!

      subroutine V21D(C,X)

      implicit none

      double precision, dimension(6,6) :: C
      double precision, dimension(21) :: X

      X  = 0D0
      X(1)  = C(1,1)
      X(2)  = C(2,2)
      X(3)  = C(3,3)
      X(4)  = sqrt(2D0)*C(2,3)
      X(5)  = sqrt(2D0)*C(1,3)
      X(6)  = sqrt(2D0)*C(1,2)
      X(7)  = 2D0*C(4,4)
      X(8)  = 2D0*C(5,5)
      X(9)  = 2D0*C(6,6)
      X(10) = 2D0*C(1,4)
      X(11) = 2D0*C(2,5)
      X(12) = 2D0*C(3,6)
      X(13) = 2D0*C(3,4)
      X(14) = 2D0*C(1,5)
      X(15) = 2D0*C(2,6)
      X(16) = 2D0*C(2,4)
      X(17) = 2D0*C(3,5)
      X(18) = 2D0*C(1,6)
      X(19) = 2D0*sqrt(2D0)*C(5,6)
      X(20) = 2D0*sqrt(2D0)*C(4,6)
      X(21) = 2D0*sqrt(2D0)*C(4,5)

      return

      end subroutine V21D

!
!***************************************************************
!

      subroutine PROJECTI(X,DEV)

      implicit none

      double precision :: DEV
      double precision, dimension(21) :: X,XH,XD

      XH=0D0
      XD=0D0
      DEV=0D0

      XH(1)=3D0/8D0*(X(1)+X(2))+X(6)/4D0/sqrt(2D0)+X(9)/4D0
      XH(2)=XH(1)
      XH(3)=X(3)
      XH(4)=(X(4)+X(5))/2D0
      XH(5)=XH(4)
      XH(6)=(X(1)+X(2))/4D0/sqrt(2D0)+3D0/4D0*X(6)              &
            -X(9)/2D0/sqrt(2D0)
      XH(7)=(X(7)+X(8))/2D0
      XH(8)=XH(7)
      XH(9)=(X(1)+X(2))/4D0-X(6)/2D0/sqrt(2D0)+X(9)/2D0

      XD=X-XH
      DEV=sqrt(dot_product(XD,XD))

      return

      end subroutine PROJECTI
!
!****************************************************************
!

      subroutine DECSYM2(XEC,DEV)!Dhexa,Dtetra,Dortho,Dmono,Dtri)

      !use comvar

      implicit none

      double precision,dimension(5) :: DEV
      double precision, dimension(21) :: XEC

      DEV = 0D0
      call PROJECMONO(XEC,DEV(1))
      call PROJECORTHO(XEC,DEV(2))
      call PROJECTETRA(XEC,DEV(3))
      call PROJECHEXA(XEC,DEV(4))
      call PROJECISO(XEC,DEV(5))

      
!!! Save anisotropy components
      !WRITE(11,"(101(1PE12.4))"),Dtri,Dmono,Dortho,Dtetra,Dhexa      

      return

      end  subroutine DECSYM2

!
!***************************************************************
!

      subroutine PROJECISO(X,DEV)

      implicit none

      double precision :: DEV
      double precision, dimension(21) :: X,XH,XD

      XH=0D0
      XD=0D0
      DEV=0D0

      XH(1)=3D0/15D0*(X(1)+X(2)+X(3))+sqrt(2D0)/15D0*(X(4)+X(5)+X(6))+2D0/15D0*(X(7)+X(8)+X(9))
      XH(2)=XH(1)
      XH(3)=XH(1)
      XH(4)=sqrt(2D0)/15D0*(X(1)+X(2)+X(3))+4D0/15D0*(X(4)+X(5)+X(6))-sqrt(2D0)/15D0*(X(7)+X(8)+X(9))
      !XH(4)=XH(4)/sqrt(2d0)
      XH(5)=XH(4)
      XH(6)=XH(4)
      XH(7)=2D0/15D0*(X(1)+X(2)+X(3))-sqrt(2D0)/15D0*(X(4)+X(5)+X(6))+1D0/5D0*(X(7)+X(8)+X(9))
      !XH(7)=XH(7)/2d0
      XH(8)=XH(7)
      XH(9)=XH(7)

!      XH(4)=XH(4)/sqrt(2d0)
!      XH(5)=XH(4)
!      XH(6)=XH(4)
!      XH(7)=XH(7)/2d0
!      XH(8)=XH(7)
!      XH(9)=XH(7)

!print *, XH
      XD=X-XH
      DEV=sqrt(dot_product(XD,XD))
!      print *,'XN,ANIS(%)=',sqrt(dot_product(X,X)),DEV/sqrt(dot_product(X,X))*100d0
!      print *,'Norm Isotropy=',sqrt(dot_product(XH,XH))
!      print *,'Norm Anisotropy=',sqrt(dot_product(XD,XD))
      
      return

      end subroutine PROJECISO
!
!***************************************************************
!

      subroutine PROJECHEXA(X,DEV)

      implicit none

      double precision :: DEV
      double precision, dimension(21) :: X,XH,XD

      XH=0D0
      XD=0D0
      DEV=0D0

      XH(1)=3D0/8D0*(X(1)+X(2))+X(6)/4D0/sqrt(2D0)+X(9)/4D0
      XH(2)=XH(1)
      XH(3)=X(3)
      XH(4)=(X(4)+X(5))/2D0
      XH(5)=XH(4)
      XH(6)=(X(1)+X(2))/4D0/sqrt(2D0)+3D0/4D0*X(6)              &
            -X(9)/2D0/sqrt(2D0)
      XH(7)=(X(7)+X(8))/2D0
      XH(8)=XH(7)
      XH(9)=(X(1)+X(2))/4D0-X(6)/2D0/sqrt(2D0)+X(9)/2D0

!print *, XH
      XD=X-XH
      DEV=sqrt(dot_product(XD,XD))
!      print *,DEV,DEV/sqrt(dot_product(X,X))
!      print *,'Norm Hexagonal=',sqrt(dot_product(XH,XH))
!      print *,'Norm Anisotropy=',sqrt(dot_product(XD,XD))
      
      return

      end subroutine PROJECHEXA
      
!
!***************************************************************
!

      subroutine PROJECTETRA(X,DEV)

      implicit none

      double precision :: DEV
      double precision, dimension(21) :: X,XH,XD

      XH=0D0
      XD=0D0
      DEV=0D0

      XH(1)=5D-1*(X(1)+X(2))
      XH(2)=XH(1)
      XH(3)=X(3)
      XH(4)=5D-1*(X(4)+X(5))
      XH(5)=XH(4)
      XH(6)=X(6)
      XH(7)=5D-1*(X(7)+X(8))
      XH(8)=XH(7)
      XH(9)=X(9)

!print *, 'Tetragonal=',XH
      XD=X-XH
      DEV=sqrt(dot_product(XD,XD))
!      print *,DEV,DEV/sqrt(dot_product(X,X))
!      print *,'Norm Tetragonal=',sqrt(dot_product(XH,XH))
!      print *,'Norm Anisotropy=',sqrt(dot_product(XD,XD))
      
      return

      end subroutine PROJECTETRA
!
!***************************************************************
!

      subroutine PROJECORTHO(X,DEV)

      implicit none

      double precision :: DEV
      double precision, dimension(21) :: X,XH,XD

      XH=0D0
      XD=0D0
      DEV=0D0

      XH(1)=X(1)
      XH(2)=X(2)
      XH(3)=X(3)
      XH(4)=X(4)
      XH(5)=X(5)
      XH(6)=X(6)
      XH(7)=X(7)
      XH(8)=X(8)
      XH(9)=X(9)

!print *, 'Orthogonal=',XH
      XD=X-XH
      DEV=sqrt(dot_product(XD,XD))
!      print *,DEV,DEV/sqrt(dot_product(X,X))
!      print *,'Norm Orthogonal=',sqrt(dot_product(XH,XH))
!      print *,'Norm Anisotropy=',sqrt(dot_product(XD,XD))
      
      return

      end subroutine PROJECORTHO
!
!***************************************************************
!

      subroutine PROJECMONO(X,DEV)

      implicit none

      double precision :: DEV
      double precision, dimension(21) :: X,XH,XD

      XH=0D0
      XD=0D0
      DEV=0D0

      XH(1)=X(1)
      XH(2)=X(2)
      XH(3)=X(3)
      XH(4)=X(4)
      XH(5)=X(5)
      XH(6)=X(6)
      XH(7)=X(7)
      XH(8)=X(8)
      XH(9)=X(9)
      XH(12)=X(12)
      XH(15)=X(15)
      XH(18)=X(18)
      XH(21)=X(21)

!print *, 'Monoclinic=',XH
      XD=X-XH
      DEV=sqrt(dot_product(XD,XD))
      !print *,DEV,DEV/sqrt(dot_product(X,X))
      !print *,'Norm Monoclinic=',sqrt(dot_product(XH,XH))
      !print *,'Norm Anisotropy=',sqrt(dot_product(XD,XD))
      
      return

      end subroutine PROJECMONO
