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
!!! Load/Save dataset, format double precision                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE loadsave_double(data_attr,rank,loc_id,dims,memtype,buf1,path,mode)

   USE comvar
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

   USE comvar
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine voigt - Calculates elastic tensor cav_{ijkl} for an olivine !!!
!!! aggregate using Voigt averaging                                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE stifftenz(m)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i,j,k,ll,n,nu,p,q,r,ss,m
   DOUBLE PRECISION :: Mcur,Gcur,Lcur,Mnew,Gnew,Lnew 
   DOUBLE PRECISION, DIMENSION(21) :: XE
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C0,Cav,CMP,CMP0
   DOUBLE PRECISION, DIMENSION(6,6) :: Sav1,CMP6
  
   C0 = 0d0 ; Cav = 0d0 ; CMP = 0d0 ; CMP0 = 0d0    
   Sav1 = 0d0 ; CMP6 = 0d0 ; Voigt = 0d0 ; Reuss = 0d0 ; Mixed = 0d0 

   !!! Major phase

   !!! converting (6,6) to (3,3,3,3) matrix
   CALL S0_2_C0(1,m,mtk(m),mpgpa(m),C0,CMP0)

   DO nu = 1 , size
      
      !IF(fractvoigt > 0d0) THEN 
         !!! VOIGT
         !!! C0: single crystal stiffness tensor
         !!! Cav: aggregate stiffness tensor
         CALL Caverage(1,m,nu,C0,Cav) 
      !END IF

      !IF(fractvoigt < 1d0) THEN
         !!! REUSS
         !!! CMP0: single crystal inverse stiffness tensor
         !!! CMP: aggregate inverse stiffness tensor
         CALL Caverage(1,m,nu,CMP0,CMP)
      !END IF 
   
   END DO

   !!! Minor phase
   IF(Xol(rocktype(m)) < 1d0) THEN

   !!! converting (6,6) to (3,3,3,3) matrix
   CALL S0_2_C0(2,m,mtk(m),mpgpa(m),C0,CMP0)
   
      DO nu = 1 , size

         !IF(fractvoigt > 0d0) THEN 
            !!! VOIGT
            !!! C0: single crystal stiffness tensor
            !!! Cav: aggregate stiffness tensor
            CALL Caverage(2,m,nu,C0,Cav) 
         !END IF
         
         !IF(fractvoigt < 1d0) THEN 
            !!! REUSS
            !!! CMP0: single crystal inverse stiffness tensor
            !!! CMP: aggregate inverse stiffness tensor
            CALL Caverage(2,m,nu,CMP0,CMP) 
         !END IF
   
      END DO

   END IF

   !!! Average stiffness matrix
   DO i = 1 , 6 ; DO j = 1 , 6
      Voigt(i,j) = Cav(l1(i),l2(i),l1(j),l2(j))
   END DO ; END DO

   !!! Inverse of CMP
   !IF(fractvoigt < 1d0) THEN
   
      DO i = 1 , 6 ; DO j = 1 , 6
         CMP6(i,j) = CMP(l1(i),l2(i),l1(j),l2(j))*mandel_scale(i,j)
      END DO ; END DO

      CALL inverse(CMP6,Reuss,6)

   !ELSE

   !END IF

!!! Correct isotropic component of single crystal elastic tensor with that from thermodynamic database
   IF(ptmod > 0) THEN

      !Find current isotropic elastic moduli
      Sav1 = Voigt  
      XE = 0d0 
      CALL V21D(Sav1,XE)
      Mcur = 3D0/15D0*(XE(1)+XE(2)+XE(3))+sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+2D0/15D0*(XE(7)+XE(8)+XE(9))
      Gcur = (2D0/15D0*(XE(1)+XE(2)+XE(3))-sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+1D0/5D0*(XE(7)+XE(8)+XE(9)))/2
      !Lambda
      Lcur = Mcur - 2*Gcur

      !Find new isotropic elastic moduli
      CALL isotropicmoduli(m,mtk(m),mpgpa(m),Mnew,Gnew)
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
      Voigt = Sav1

      Sav1 = Reuss  
      XE = 0d0 
      CALL V21D(Sav1,XE)
      Mcur = 3D0/15D0*(XE(1)+XE(2)+XE(3))+sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+2D0/15D0*(XE(7)+XE(8)+XE(9))
      Gcur = (2D0/15D0*(XE(1)+XE(2)+XE(3))-sqrt(2D0)/15D0*(XE(4)+XE(5)+XE(6))+1D0/5D0*(XE(7)+XE(8)+XE(9)))/2
      !Lambda
      Lcur = Mcur - 2*Gcur

      !Find new isotropic elastic moduli
      CALL isotropicmoduli(m,mtk(m),mpgpa(m),Mnew,Gnew)
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
      Reuss = Sav1

   END IF 
!!! End correct isotropic component of single crystal elastic tensor with that from thermodynamic database


   DO i = 1 , 6 ; DO j = 1 , 6
      Mixed(i,j) = fractvoigt*Voigt(i,j) + (1d0-fractvoigt)*Reuss(i,j)
   END DO ; END DO

   RETURN

   END SUBROUTINE stifftenz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE Caverage(phase,m,nu,C0,Cav)

   USE comvar

   IMPLICIT NONE
   
   INTEGER :: i,j,k,ll,p,q,r,ss,phase,m,nu
   DOUBLE PRECISION, DIMENSION (3,3,3,3) :: C0,Cav,Cav2

   !!! Cav2 : oriented with the external axis
   !!! acs(k,j) = cosine of the angle between the kth crystallographic axes and the jth external axes 
   Cav2 = 0d0

   IF(phase == 1) THEN

      IF(rocktype(m) == 3) THEN
   
         ! Isotropic LPO for Ringwoodite
         DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
            DO p = 1 , 3 ; DO q = 1 , 3 ; DO r = 1 , 3 ; DO ss = 1 , 3
               Cav2(i,j,k,ll) = Cav2(i,j,k,ll) + &
               acs0(nu,p,i)*acs0(nu,q,j)*acs0(nu,r,k)*acs0(nu,ss,ll)*C0(p,q,r,ss)
            END DO ; END DO ; END DO ; END DO
            Cav(i,j,k,ll) = Cav(i,j,k,ll) + Cav2(i,j,k,ll)*odf(m,nu)*Xol(rocktype(m))
         END DO ; END DO ; END DO ; END DO
   
      ELSE

         DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
            DO p = 1 , 3 ; DO q = 1 , 3 ; DO r = 1 , 3 ; DO ss = 1 , 3
               Cav2(i,j,k,ll) = Cav2(i,j,k,ll) + &
               acs(nu,p,i,m)*acs(nu,q,j,m)*acs(nu,r,k,m)*acs(nu,ss,ll,m)*C0(p,q,r,ss)
            END DO ; END DO ; END DO ; END DO
            Cav(i,j,k,ll) = Cav(i,j,k,ll) + Cav2(i,j,k,ll)*odf(m,nu)*Xol(rocktype(m))
         END DO ; END DO ; END DO ; END DO

      END IF

   ELSE

      IF(rocktype(m) == 1) THEN

         DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
            DO p = 1 , 3 ; DO q = 1 , 3 ; DO r = 1 , 3 ; DO ss = 1 , 3
               Cav2(i,j,k,ll) = Cav2(i,j,k,ll) + &
               acs_ens(nu,p,i,m)*acs_ens(nu,q,j,m)*acs_ens(nu,r,k,m)*acs_ens(nu,ss,ll,m)*C0(p,q,r,ss)
            END DO ; END DO ; END DO ; END DO
            Cav(i,j,k,ll) = Cav(i,j,k,ll) + Cav2(i,j,k,ll)*odf_ens(m,nu)*(1d0-Xol(rocktype(m)))
         END DO ; END DO ; END DO ; END DO

      ELSE
  
         ! Isotropic LPO for Garnet,Rw and MgO
         DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
            DO p = 1 , 3 ; DO q = 1 , 3 ; DO r = 1 , 3 ; DO ss = 1 , 3
               Cav2(i,j,k,ll) = Cav2(i,j,k,ll) + &
               acs0(nu,p,i)*acs0(nu,q,j)*acs0(nu,r,k)*acs0(nu,ss,ll)*C0(p,q,r,ss)
            END DO ; END DO ; END DO ; END DO
            Cav(i,j,k,ll) = Cav(i,j,k,ll) + Cav2(i,j,k,ll)*odf_ens(m,nu)*(1d0-Xol(rocktype(m)))
         END DO ; END DO ; END DO ; END DO

      END IF

   END IF

   END SUBROUTINE Caverage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE S0_2_C0(phase,m,mtk0,mpgpa0,C0,CMP0)

   USE comvar

   IMPLICIT NONE

   INTEGER :: phase,mtype,i,j,k,ll,m
   DOUBLE PRECISION :: mpgpa0,mtk0,DP,DTk
   DOUBLE PRECISION, DIMENSION (6,6) :: S0dummy,S0pt,dSdt,dSdp,CMP,CMPpt
   DOUBLE PRECISION, DIMENSION (3,3,3,3) :: C0,CMP0

   S0pt = 0d0 ; dSdt = 0d0 ; dSdp = 0d0 ; CMP = 0d0 ; CMPpt = 0d0 ; C0 = 0d0 ; CMP0 = 0d0
 
   mtype = single_crystal_elastic_db(rocktype(m),phase)
   
!!! Room P-T 
   IF(ptmod == 0 .OR. mtype == 8 .OR. mtype >= 17) THEN

      DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
         C0(i,j,k,ll) = S0(mtype,ijkl(i,j),ijkl(k,ll))
      END DO ; END DO ; END DO ; END DO
  
      !IF(fractvoigt < 1d0) THEN
         S0dummy=S0(mtype,:,:)
         CALL inverse(S0dummy,CMP,6)
        
         DO i = 1, 6 ; DO j = 1 , 6
            CMP(i,j) = CMP(i,j)/mandel_scale(i,j)
         END DO ; END DO

         DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
            CMP0(i,j,k,ll) = CMP(ijkl(i,j),ijkl(k,ll))
         END DO ; END DO ; END DO ; END DO
 
      !END IF
 
!!! Correction of elastic moduli for pressure and temperature derivatives
   ELSE
     
      DP=mpgpa0-1d-4
      DTk=mtk0-298d0

      !!! Interpolate pressure and temperature derivatives between 0 and 100 GPa
      IF(mtype == 7 .and. mpgpa0 <= 100d0) THEN
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
         DTk=mtk0-1500
         DP=mpgpa0-34
      END IF
      IF(mtype == 16) THEN
         DTk=mtk0-1500
         DP=mpgpa0-48
      END IF

      DO i = 1, 6 ; DO j = 1 , 6
         S0pt(i,j) = S0(mtype,i,j) - dSdt(i,j)*DTk + &
               dSdp(i,j)*DP + 0.5*dS0dp2(mtype,i,j)*(DP**2) + dS0dpt(mtype,i,j)*DP*DTk
      END DO ; END DO

      DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
         C0(i,j,k,ll) = S0pt(ijkl(i,j),ijkl(k,ll))
      END DO ; END DO ; END DO ; END DO
   
      !IF(fractvoigt < 1d0) THEN
         CALL inverse(S0pt,CMPpt,6)

         DO i = 1, 6 ; DO j = 1 , 6
            CMPpt(i,j) = CMPpt(i,j)/mandel_scale(i,j)
         END DO ; END DO

         DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
            CMP0(i,j,k,ll) = CMPpt(ijkl(i,j),ijkl(k,ll))
         END DO ; END DO ; END DO ; END DO
 
      !END IF

   END IF

   END SUBROUTINE S0_2_C0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE isotropicmoduli(m,mtk0,mpgpa0,K,G)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION mtk0,mpa,mpgpa0,mpb,ee,n,Vp,Vs
   DOUBLE PRECISION K,K0,K1,K2,K3,G,G0,G1,G2,G3,R0,R1,R2,R3
   INTEGER m,n1,n2

   mpa = mpgpa0 * 1d9
   !Transform pressure from Pascal to bar
   mpb = mpa / 1d5

   ! ABCD-4Cell Number 
   ee = (mtk0-tkmin)/tkstp
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
   R0=td_rho(n1  ,n2  )
   R1=td_rho(n1  ,n2+1)
   R2=td_rho(n1+1,n2  )
   R3=td_rho(n1+1,n2+1)

   rho(1)=((R0*(1.0-n)+R1*n)*(1.0-ee)+(R2*(1.0-n)+R3*n)*ee)

   !Vs (km/s)
   Vs = ((G0*(1.0-n)+G1*n)*(1.0-ee)+(G2*(1.0-n)+G3*n)*ee)
   !Swave modulus (GPa)
   G = rho(1)*(Vs)**2/1d3

   !Vp (km/s)
   Vp = ((K0*(1.0-n)+K1*n)*(1.0-ee)+(K2*(1.0-n)+K3*n)*ee)
   !Pwave modulus (GPa)
   K = rho(1)*(Vp)**2/1d3

   RETURN

   END SUBROUTINE isotropicmoduli

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine output, saving fse, hexagonal symmetry axis                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE output

   USE comvar

   IMPLICIT NONE

   INTEGER :: t
   DOUBLE PRECISION :: Ax,Ay,Az,Bx,By,Bz
   !Buffers
   
   DOUBLE PRECISION, DIMENSION(3,3) :: LSij
   ! left-strech tensor for FSE calculation

   INTEGER :: nrot
   ! number of rotations for the Jacobi

   DOUBLE PRECISION, DIMENSION(3) :: evals
   DOUBLE PRECISION, DIMENSION(3,3) :: evects
   ! eigen values and vectors in jacobi

!!! Left-stretch tensor for FSE calculation

   CALL eigen(1,Fij(:,:,1),evects,evals,1)

!!! Pick up the orientation of the long axis of the FSE 

   IF (evals(1) .EQ. MAXVAL(evals)) THEN
      phi_fse(:) = evects(:,1)
      ELSE IF(evals(2) .EQ. MAXVAL(evals)) THEN
         phi_fse(:) = evects(:,2)
         ELSE
         phi_fse(:) = evects(:,3)
   END IF

!!! natural strain = ln(a/c) where a is the long axis = maxval(evals)**0.5

   ln_fse = 0.5*LOG(MAXVAL(evals)/MINVAL(evals))

!!! Save only if enough strain rate (otherwise anisotropy data may be wrong)
!!! ln_fse = 0.5 it is enough to overprint any preexisting fabric (Ribe, 1992, JGR; Becker, 2006, EPSL) 

   print *,'Max strain:',max_strain,timesum,ln_fse,perc_a      
   
   IF(ln_fse .GE. 0.0) THEN

      !anis_x
      Bx = perc_a*phi_a(1)
      !anis_y
      By = perc_a*phi_a(2)
      !anis_z
      Bz = perc_a*phi_a(3)
   
      WRITE(10,"(101(1PE12.4))"), max_strain,timesum,ln_fse,perc_a      

   END IF

   END SUBROUTINE output 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

      double precision, dimension(3,3) :: SCC
      double precision, dimension(6,6) :: CE1
      double precision, dimension(3,3,3,3) :: EL1
      double precision, dimension(21) :: XEC
      double precision :: XN,ANIS

      end module DECMOD

!
!****************************************************************
!

      subroutine DECSYM(CED,PERC,TIAXIS)

      use DECMOD

      implicit none

      double precision, dimension(6,6) :: CED
      double precision :: PERC,INCLTI,DC5,PI
      double precision, dimension(3) :: TIAXIS

      PI = 3.141592653589793238462643383279
      CE1=CED
      EL1=0D0
      call FULLSYM6(CE1)
      call TENS4(CE1,EL1)
      call V21D(CE1,XEC)
      !XN=sqrt(dot_product(XEC,XEC))

      call SCCA
      call PROJECTI(XEC,DC5)

      !!! After rotation of the elastic tensor toward TI symmetry axis
      XN=sqrt(dot_product(XEC,XEC))
      !!! Elastic tensor decomposition
      call DECSYM2

      PERC=(ANIS-DC5)/XN*100D0
!print *,'Anis, DC5=',ANIS/XN*100D0,DC5/XN*100D0
      TIAXIS=SCC(3,:)
      TIAXIS=TIAXIS/sqrt(sum(TIAXIS*TIAXIS))
      !INCLTI=asin(TIAXIS(3))

      !Reset negative hexagonal anisotropy
      IF(PERC < 0d0) then
         PERC = 0d0
         TIAXIS = 0d0
      END IF
     
      print *,'Transverse anisotropy:',PERC,TIAXIS(1),TIAXIS(2),TIAXIS(3)
      
      return

      end  subroutine DECSYM

!
!****************************************************************
!

      subroutine DECSYM2

      use DECMOD

      implicit none

      double precision :: Diso,Dhexa,Dtetra,Dortho,Dmono,Dtri
      double precision,dimension(5) :: DEV

      Diso=0d0
      Dhexa=0d0
      Dtetra=0d0
      Dortho=0d0
      Dmono=0d0
      Dtri=0d0
      DEV=0d0

      call PROJECMONO(XEC,DEV(1))
      call PROJECORTHO(XEC,DEV(2))
      call PROJECTETRA(XEC,DEV(3))
      call PROJECHEXA(XEC,DEV(4))
      call PROJECISO(XEC,DEV(5))

      DEV=DEV/XN*100d0

!print *,'DEV=',DEV
!print *,'XN=',XN    
      
      Dtri=DEV(1)
      Dmono=DEV(2)-DEV(1)
      Dortho=DEV(3)-DEV(2)
      Dtetra=DEV(4)-DEV(3)
      Dhexa=DEV(5)-DEV(4)

print *,'Anis components (%)=',Dtri,Dmono,Dortho,Dtetra,Dhexa
      
!!! Save anisotropy components
      !WRITE(11,"(101(1PE12.4))"),Dtri,Dmono,Dortho,Dtetra,Dhexa      

      return

      end  subroutine DECSYM2

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

      subroutine JACOBI(a,n,np,d,v,nrot)

      ! Jacobi algorithm for real symmetric matrix
      ! Gives eigenvalues and orthonormalized eigenvectors
      ! Half of the input matrix is crushed

      implicit none

      integer :: n,np,nrot
      integer, parameter :: NMAX=500,IDP=kind(1D0)

      double precision, dimension(np,np) :: a,v
      double precision, dimension(np) :: d
      double precision, dimension(NMAX) :: b,z

      integer :: i,ip,iq,j
      double precision :: c,g,h,s,sm,t,tau,theta,tresh

      do ip=1,n
         do iq=1,n
            v(ip,iq)=0D0
         enddo
         v(ip,ip)=1D0
      enddo
      do ip=1,n
         b(ip)=a(ip,ip)
         d(ip)=b(ip)
         z(ip)=0D0
      enddo
      nrot=0
      do i=1,50
         sm=0D0
         do ip=1,n-1
         do iq=ip+1,n
            sm=sm+abs(a(ip,iq))
         enddo
         enddo
         if (sm==0D0) return
         if (i<4) then
            tresh=0.2D0*sm/real(n,IDP)**2D0
         else
            tresh=0D0
         endif
         do ip=1,n-1
         do iq=ip+1,n
            g=100D0*abs(a(ip,iq))
            if ((i>4).and.(abs(d(ip))+                          &
            g==abs(d(ip))).and.(abs(d(iq))+g==abs(d(iq)))) then
               a(ip,iq)=0D0
            else if (abs(a(ip,iq))>tresh) then
               h=d(iq)-d(ip)
               if (abs(h)+g==abs(h)) then
                  t=a(ip,iq)/h
               else
                  theta=0.5D0*h/a(ip,iq)
                  t=1D0/(abs(theta)+sqrt(1D0+theta**2D0))
                  if (theta<0D0) t=-t
               endif
               c=1D0/sqrt(1D0+t**2D0)
               s=t*c
               tau=s/(1D0+c)
               h=t*a(ip,iq)
               z(ip)=z(ip)-h
               z(iq)=z(iq)+h
               d(ip)=d(ip)-h
               d(iq)=d(iq)+h
               a(ip,iq)=0D0
               do j=1,ip-1
                  g=a(j,ip)
                  h=a(j,iq)
                  a(j,ip)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
               enddo
               do j=ip+1,iq-1
                  g=a(ip,j)
                  h=a(j,iq)
                  a(ip,j)=g-s*(h+g*tau)
                  a(j,iq)=h+s*(g-h*tau)
               enddo
               do j=iq+1,n
                  g=a(ip,j)
                  h=a(iq,j)
                  a(ip,j)=g-s*(h+g*tau)
                  a(iq,j)=h+s*(g-h*tau)
               enddo
               do j=1,n
                  g=v(j,ip)
                  h=v(j,iq)
                  v(j,ip)=g-s*(h+g*tau)
                  v(j,iq)=h+s*(g-h*tau)
               enddo
               nrot=nrot+1
            endif
         enddo
         enddo
         do ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)
            z(ip)=0D0
         enddo
      enddo
      write(6,'(''Too many iterations in JACOBI'')')

      return

      end subroutine JACOBI

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

      subroutine SCCA

      use DECMOD

      implicit none

      integer :: i,NROT,i1,i2,NDVC
      integer, dimension(3) :: IHS
      double precision, dimension(3) :: EGDI,EGVO
      double precision, dimension(3,3) :: DI,VO,VECDI,VECVO
      double precision, dimension(6,6) :: CEC
      double precision, dimension(3,3,3,3) :: ELC
      double precision, dimension(21) :: XH,XD
      double precision :: SDV,ADV,ADVC,SCN,DEV,DEV2,K,G

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

      !XD=XEC-XH
      !ANIS=sqrt(dot_product(XD,XD))
      !print *,'Norm Isotropy=',sqrt(dot_product(XH,XH))
      !print *,'Norm Anisotropy=',ANIS/XN*100d0

      ! Dil. and Voigt axes

      call JACOBI(DI,3,3,EGDI,VECDI,NROT)
      call EIGSRT(EGDI,VECDI,3,3)
      call JACOBI(VO,3,3,EGVO,VECVO,NROT)
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

      XD=XEC-XH
      ANIS=sqrt(dot_product(XD,XD))
      
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

      double precision :: DEV,DEV2
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
      DEV2=sqrt(dot_product(XH,XH))
!print *,DEV/sqrt(dot_product(X,X)),DEV2/sqrt(dot_product(X,X))
      
      return

      end subroutine PROJECTI
      
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
      XH(5)=XH(4)
      XH(6)=XH(4)
      XH(7)=2D0/15D0*(X(1)+X(2)+X(3))-sqrt(2D0)/15D0*(X(4)+X(5)+X(6))+1D0/5D0*(X(7)+X(8)+X(9))
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

      double precision :: DEV,DEV2
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
 !     print *,DEV,DEV/sqrt(dot_product(X,X))
 !     print *,'Norm Monoclinic=',sqrt(dot_product(XH,XH))
 !     print *,'Norm Anisotropy=',sqrt(dot_product(XD,XD))
      
      return

      end subroutine PROJECMONO
