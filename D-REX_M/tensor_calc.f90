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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Calculates elastic tensor cav_{ijkl} for an aggregate using Voigt-Reuss averaging  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE tensorscalc(m,mtk,mpgpa,Sav1)

   USE comvar

   !M3E!!!!!!!!!!!!!!
   use class_DistGrid
   !M3E!!!!!!!!!!!!!!

   IMPLICIT NONE

   INTEGER :: i,j,k,ll,nu,m
   DOUBLE PRECISION :: mtk,mpgpa,volfract
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C0,Cav,CMP,CMP0,CMP2
   DOUBLE PRECISION, DIMENSION(6,6) :: Sav1,CMP6,Voigt,Reuss
   DOUBLE PRECISION, DIMENSION(nboxnum) :: f1_ol,f1_opx
  
   volfract = 1d0/REAL(size3**3)

   !Run SBFTEX
   !IF(sbfmod > 0 .AND. rocktype(m) == 1 ) CALL SBFTEX(m,f1_ol,f1_opx)

10 C0 = 0d0 ; Cav = 0d0 ; CMP = 0d0 ; CMP0 = 0d0 ; CMP2 = 0d0 
   Sav1 = 0d0 ; CMP6 = 0d0 ; Voigt = 0d0 ; Reuss = 0d0 ; 

   !!! Major phase

   !!! converting (6,6) to (3,3,3,3) matrix
   CALL S0_2_C0(1,m,mtk,mpgpa,C0,CMP0)

   !SBFTEX
   IF(sbfmod > 0) THEN

      IF(rocktype(m) == 1) THEN
         DO nu = 1 , nboxnum
            IF(fractvoigt > 0d0) CALL Caverage(f1_ol(nu)*Xol(rocktype(m)),acs0sbf(nu,:,:),C0,Cav)
            IF(fractvoigt < 1d0) CALL Caverage(f1_ol(nu)*Xol(rocktype(m)),acs0sbf(nu,:,:),CMP0,CMP) 
         END DO
      ELSE
         IF(ptmod > 0) THEN
            CALL isotensor(m,mtk,mpgpa,1,Sav1)
            RETURN 
         ELSE
            CALL c2sav(C0,Sav1)
            CALL isotensor(m,mtk,mpgpa,0,Sav1)
            Voigt = Sav1*Xol(rocktype(m))
         END IF
      END IF

   !D-REX_M
   ELSE

      IF(rocktype(m) == 3) THEN

         IF(ptmod > 0) THEN 
            CALL isotensor(m,mtk,mpgpa,1,Sav1)
            RETURN
         ELSE 
            CALL c2sav(C0,Sav1)
            CALL isotensor(m,mtk,mpgpa,0,Sav1)
            Voigt = Sav1*Xol(rocktype(m))
         END IF

      ELSE

         DO nu = 1 , size
            IF(fractvoigt > 0d0) CALL Caverage(odf(m,nu)*Xol(rocktype(m)),acs(nu,:,:,m),C0,Cav)
            IF(fractvoigt < 1d0) CALL Caverage(odf(m,nu)*Xol(rocktype(m)),acs(nu,:,:,m),CMP0,CMP)
         END DO

      END IF

   END IF

   !!! Minor phase
   IF(Xol(rocktype(m)) < 1d0) THEN

   !!! converting (6,6) to (3,3,3,3) matrix
   CALL S0_2_C0(2,m,mtk,mpgpa,C0,CMP0)
   
   IF(sbfmod > 0) THEN

      IF(rocktype(m) == 1 ) THEN
         DO nu = 1 , nboxnum
            IF(fractvoigt > 0d0) CALL Caverage(f1_opx(nu)*(1.0-Xol(rocktype(m))),acs0sbf(nu,:,:),C0,Cav)
            IF(fractvoigt < 1d0) CALL Caverage(f1_opx(nu)*(1.0-Xol(rocktype(m))),acs0sbf(nu,:,:),CMP0,CMP)
         END DO
      ELSE
         CALL c2sav(C0,Sav1)
         CALL isotensor(m,mtk,mpgpa,0,Sav1)
         Sav1 = Voigt + Sav1*(1.0-Xol(rocktype(m)))
         RETURN
      END IF

   ELSE

      IF(rocktype(m) == 1) THEN

         DO nu = 1 , size
            IF(fractvoigt > 0d0) CALL Caverage(odf_ens(m,nu)*(1.0-Xol(rocktype(m))),acs_ens(nu,:,:,m),C0,Cav)
            IF(fractvoigt < 1d0) CALL Caverage(odf_ens(m,nu)*(1.0-Xol(rocktype(m))),acs_ens(nu,:,:,m),CMP0,CMP)
         END DO

      ELSE

         CALL c2sav(C0,Sav1)
         CALL isotensor(m,mtk,mpgpa,0,Sav1)
         !For isotropic tensors, stiffness and compliance tensors are the same
         DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
            Cav(i,j,k,ll) = Cav(i,j,k,ll) + Sav1(ijkl(i,j),ijkl(k,ll))*(1.0-Xol(rocktype(m)))
            CMP(i,j,k,ll) = CMP(i,j,k,ll) + Sav1(ijkl(i,j),ijkl(k,ll))*(1.0-Xol(rocktype(m)))
         END DO ; END DO ; END DO ; END DO
      
      END IF

   END IF

   END IF

   !Reset fabric if bulk moduli are NaN
   IF(isnan(Cav(1,1,1,1)) .OR. (Cav(1,1,1,1)>400d0 .AND. rocktype(m)==1)) THEN
      IF(sbfmod == 0 ) THEN
         print *,"Cav Nan at",m,rocktype(m),mx1(m),mx2(m),mx3(m),mtk,mpgpa,odf(m,1),acs(1,1,1,m),acs_ens(1,1,1,m),Cav(1,1,1,1)
         acs(:,:,:,m) = acs0
         acs_ens(:,:,:,m) = acs0
         odf(m,:) = volfract           
         odf_ens(m,:) = odf(m,:)
      ELSE
         !IF(rocktype(m)==1) THEN
            print *,"Cav Nan at",m,rocktype(m),mx1(m),mx2(m),mx3(m),mtk,mpgpa,f1_ol(1),f1_opx(1),Cav(1,1,1,1)
            f1_ol = 1d0/REAL(nboxnum)
            f1_opx= 1d0/REAL(nboxnum)
         !END IF
      END IF
      goto 10
   END IF

   !!! Average stiffness matrix
   CALL c2sav(Cav,Voigt)

   !!! Inverse of CMP
   IF(fractvoigt < 1d0) THEN
   
      DO i = 1 , 6 ; DO j = 1 , 6
         CMP6(i,j) = CMP(l1(i),l2(i),l1(j),l2(j))*mandel_scale(i,j)
      END DO ; END DO

      CALL inverse(CMP6,Reuss,6)

      DO i = 1 , 6 ; DO j = 1 , 6
         Reuss(i,j) = Reuss(i,j)/mandel_scale(i,j)
      END DO ; END DO

      DO i = 1 , 6 ; DO j = 1 , 6
         Sav1(i,j) = fractvoigt*Voigt(i,j) + (1d0-fractvoigt)*Reuss(i,j)
      END DO ; END DO

   ELSE

      Sav1 = Voigt  

   END IF

   !!! Correct isotropic component of aggregate elastic tensor with that from thermodynamic database
   IF(ptmod > 0) CALL isocorrection(m,mtk,mpgpa,Sav1)

   RETURN

   END SUBROUTINE tensorscalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE Caverage(w,acss,C0,Cav)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i,j,k,ll,p,q,r,ss,phase,m,nu
   DOUBLE PRECISION :: w
   DOUBLE PRECISION, DIMENSION (3,3) :: acss
   DOUBLE PRECISION, DIMENSION (3,3,3,3) :: C0,Cav,Cav2

   !!! f_sbf : weigth of crystal in the ODF from sbf (Ribe et al., 2019, GJI)

   !!! Cav2 : oriented with the external axis
   !!! acs(k,j) = cosine of the angle between the kth crystallographic axes and
   !the jth external axes 
   Cav2 = 0d0

   DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
      DO p = 1 , 3 ; DO q = 1 , 3 ; DO r = 1 , 3 ; DO ss = 1 , 3
         Cav2(i,j,k,ll) = Cav2(i,j,k,ll) + &
         acss(p,i)*acss(q,j)*acss(r,k)*acss(ss,ll)*C0(p,q,r,ss)
      END DO ; END DO ; END DO ; END DO
      Cav(i,j,k,ll) = Cav(i,j,k,ll) + Cav2(i,j,k,ll)*w
   END DO ; END DO ; END DO ; END DO

   END SUBROUTINE Caverage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE S0_2_C0(phase,m,mtk,mpgpa,C0,CMP0)

   USE comvar

   IMPLICIT NONE

   INTEGER :: phase,mtype,i,j,k,ll,m
   DOUBLE PRECISION :: mpgpa,mtk,DP,DTk
   DOUBLE PRECISION, DIMENSION (6,6) :: S0dummy,S0pt,dSdt,dSdp,CMP,CMPpt
   DOUBLE PRECISION, DIMENSION (3,3,3,3) :: C0,CMP0

   S0pt = 0d0 ; dSdt = 0d0 ; dSdp = 0d0 ; CMP = 0d0 ; CMPpt = 0d0 ; C0 = 0d0 ; CMP0 = 0d0
 
   mtype = single_crystal_elastic_db(rocktype(m),phase)
   
!!! Room P-T 
   IF(ptmod == 0 .OR. mtype == 8 .OR. mtype >= 17) THEN

      CALL sav2c(S0(mtype,:,:),C0)
  
      IF(fractvoigt < 1d0) THEN
         S0dummy=S0(mtype,:,:)
         DO i = 1, 6 ; DO j = 1 , 6
            S0dummy(i,j) = S0dummy(i,j)*mandel_scale(i,j)
         END DO ; END DO

         CALL inverse(S0dummy,CMP,6)
        
         DO i = 1, 6 ; DO j = 1 , 6
            CMP(i,j) = CMP(i,j)/mandel_scale(i,j)
         END DO ; END DO

         DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
            CMP0(i,j,k,ll) = CMP(ijkl(i,j),ijkl(k,ll))
         END DO ; END DO ; END DO ; END DO
 
      END IF
 
!!! Correction of elastic moduli for pressure and temperature derivatives
   ELSE
     
      DP=mpgpa-1d-4
      DTk=mtk-298d0

      !!! Interpolate pressure and temperature derivatives between 0 and 100 GPa
      IF(mtype == 7 .and. mpgpa <= 100d0) THEN
         DO i = 1, 6 ; DO j = 1 , 6
            dSdp(i,j) = dS0dp5(i,j) - (dS0dp5(i,j)-dS0dp(mtype,i,j))/100d0*DP
            dSdt(i,j) = dS0dt5(i,j) - (dS0dt5(i,j)-dS0dt(mtype,i,j))/100d0*DP
         END DO ; END DO
      ELSE 
         DO i = 1, 6 ; DO j = 1 , 6
            dSdp(i,j) = dS0dp(mtype,i,j) 
            dSdt(i,j) = dS0dt(mtype,i,j) 
         END DO ; END DO
      END IF
    
      IF(mtype == 15) THEN
         DTk=mtk-1500
         DP=mpgpa-34
      END IF
      IF(mtype == 16) THEN
         DTk=mtk-1500
         DP=mpgpa-48
      END IF
      IF(mtype == 11) THEN
         DTk=mtk - 1500d0
         DP=mpgpa - 99.95d0
      END IF

      IF(DP < 0d0) DP = 0d0!; END IF
      IF(DTk< 0d0) DTk= 0d0!; END IF

      DO i = 1, 6 ; DO j = 1 , 6
         S0pt(i,j) = S0(mtype,i,j) - dSdt(i,j)*DTk + &
               dSdp(i,j)*DP + 0.5*dS0dp2(mtype,i,j)*(DP**2) + dS0dpt(mtype,i,j)*DP*DTk
      END DO ; END DO

      CALL sav2c(S0pt,C0)
   
      IF(fractvoigt < 1d0) THEN
         DO i = 1, 6 ; DO j = 1 , 6
            S0pt(i,j) = S0pt(i,j)*mandel_scale(i,j)
         END DO ; END DO

         CALL inverse(S0pt,CMPpt,6)

         DO i = 1, 6 ; DO j = 1 , 6
            CMPpt(i,j) = CMPpt(i,j)/mandel_scale(i,j)
         END DO ; END DO

         DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
            CMP0(i,j,k,ll) = CMPpt(ijkl(i,j),ijkl(k,ll))
         END DO ; END DO ; END DO ; END DO
 
      END IF

   END IF

   END SUBROUTINE S0_2_C0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE isocorrection(m,mtk,mpgpa,Sav1)

   IMPLICIT NONE
   
   INTEGER :: m
   DOUBLE PRECISION :: mtk,mpgpa,Mcur,Gcur,Lcur,Mnew,Gnew,Lnew 
   DOUBLE PRECISION, DIMENSION(21) :: XE
   DOUBLE PRECISION, DIMENSION(6,6) :: Sav1

   !Find current isotropic elastic moduli
   XE = 0d0 
   CALL V21DD(Sav1,XE)
   Mcur = 3D0/15D0*(XE(1)+XE(2)+XE(3))+sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+2D0/15D0*(XE(7)+XE(8)+XE(9))
   Gcur = (2D0/15D0*(XE(1)+XE(2)+XE(3))-sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+1D0/5D0*(XE(7)+XE(8)+XE(9)))/2
   !Lambda
   Lcur = Mcur - 2*Gcur

   !Find new isotropic elastic moduli
   CALL isotropicmoduli(m,mtk,mpgpa,Mnew,Gnew)
   !Lambda
   Lnew = Mnew - 2*Gnew

   Sav1(1,1) = Sav1(1,1) + (Mnew - Mcur)  
   Sav1(2,2) = Sav1(2,2) + (Mnew - Mcur)  
   Sav1(3,3) = Sav1(3,3) + (Mnew - Mcur)  
   Sav1(1,2) = Sav1(1,2) + (Lnew - Lcur)  
   Sav1(2,1) = Sav1(1,2)
   Sav1(1,3) = Sav1(1,3) + (Lnew - Lcur)  
   Sav1(3,1) = Sav1(1,3)
   Sav1(2,3) = Sav1(2,3) + (Lnew - Lcur)  
   Sav1(3,2) = Sav1(2,3)
   Sav1(4,4) = Sav1(4,4) + (Gnew - Gcur)  
   Sav1(5,5) = Sav1(5,5) + (Gnew - Gcur)  
   Sav1(6,6) = Sav1(6,6) + (Gnew - Gcur)  

   END SUBROUTINE isocorrection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE isotensor(m,mtk,mpgpa,minrock,Sav1)

   USE comvar

   IMPLICIT NONE
   
   INTEGER :: m,minrock
   DOUBLE PRECISION :: mtk,mpgpa,Mnew,Gnew,Lnew 
   DOUBLE PRECISION, DIMENSION(21) :: XE
   DOUBLE PRECISION, DIMENSION(6,6) :: Sav1

   !Find new isotropic elastic moduli
   IF(minrock > 0) THEN
      CALL isotropicmoduli(m,mtk,mpgpa,Mnew,Gnew)
   ELSE
      XE = 0d0
      CALL V21DD(Sav1,XE)
      Mnew = 3D0/15D0*(XE(1)+XE(2)+XE(3))+sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+2D0/15D0*(XE(7)+XE(8)+XE(9))
      Gnew = (2D0/15D0*(XE(1)+XE(2)+XE(3))-sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+1D0/5D0*(XE(7)+XE(8)+XE(9)))/2   
   END IF

   !Lambda
   Lnew = Mnew - 2*Gnew

   Sav1(1,1) = Mnew   
   Sav1(2,2) = Sav1(1,1)
   Sav1(3,3) = Sav1(1,1)
   Sav1(1,2) = Lnew
   Sav1(2,1) = Sav1(1,2)
   Sav1(1,3) = Sav1(1,2)
   Sav1(3,1) = Sav1(1,2)
   Sav1(2,3) = Sav1(1,2)
   Sav1(3,2) = Sav1(1,2)
   Sav1(4,4) = Gnew
   Sav1(5,5) = Sav1(4,4)
   Sav1(6,6) = Sav1(4,4)

   END SUBROUTINE isotensor    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE V21DD(C,X)

   IMPLICIT NONE

   DOUBLE PRECISION, DIMENSION(6,6) :: C
   DOUBLE PRECISION, DIMENSION(21) :: X

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

   RETURN

   END SUBROUTINE V21DD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE isotropicmoduli(m,mtk,mpgpa,K,G) 

   USE comvar
   USE omp_lib

   IMPLICIT NONE
  
   DOUBLE PRECISION mtk,mpa,mpgpa,mpb,ee,n,Vp,Vs,K,K0,K1,K2,K3,G,G0,G1,G2,G3
   INTEGER m,n1,n2

!                          
!                          
!                 2----4   
!       +Y        |    |   
!       |         |    |  
!       |         |    | 
!       |         1----3
!       |  
!       ----- +X
 
   mpa = mpgpa*1d9

   IF(mtk<300d0) mtk=300d0
   IF(mpa<1d5) mpa=1d5

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
   ! Vs,Vp values
   ! 0 2
   ! 1 3 
   G0=td_vs(n1  ,n2  )
   G1=td_vs(n1  ,n2+1)
   G2=td_vs(n1+1,n2  )
   G3=td_vs(n1+1,n2+1)
   K0=td_vp(n1  ,n2  )
   K1=td_vp(n1  ,n2+1)
   K2=td_vp(n1+1,n2  )
   K3=td_vp(n1+1,n2+1)
   
   !Vs (km/s)
   Vs = ((G0*(1.0-n)+G1*n)*(1.0-ee)+(G2*(1.0-n)+G3*n)*ee)
   !Swave modulus (GPa)
   G = rho(m)*(Vs)**2/1d3

   !Vp (km/s)
   Vp = ((K0*(1.0-n)+K1*n)*(1.0-ee)+(K2*(1.0-n)+K3*n)*ee)
   !Pwave modulus (GPa)
   K = rho(m)*(Vp)**2/1d3 

   RETURN

   END SUBROUTINE isotropicmoduli 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE sav2c(Sav1,C0)

   USE comvar
   
   IMPLICIT NONE
  
   INTEGER :: i,j,k,ll
   DOUBLE PRECISION, DIMENSION(6,6) :: Sav1
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C0
  
   DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
      C0(i,j,k,ll) = Sav1(ijkl(i,j),ijkl(k,ll))
   END DO ; END DO ; END DO ; END DO

   END SUBROUTINE sav2c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE c2sav(C0,Sav1)

   USE comvar
   
   IMPLICIT NONE
  
   INTEGER :: i,j
   DOUBLE PRECISION, DIMENSION(6,6) :: Sav1
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C0
  
   DO i = 1 , 6 ; DO j = 1 , 6
      Sav1(i,j) = C0(l1(i),l2(i),l1(j),l2(j))
   END DO ; END DO

   END SUBROUTINE c2sav
