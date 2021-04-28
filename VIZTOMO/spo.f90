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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Read input file inputspo.dat and allocate matrices    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE initspo

   USE comvar
   USE omp_lib
   
   INTEGER :: ii,jj,kk,ll,j1,j2,rt,nt

   OPEN(15,file='spo_input.dat')

!!! Read input values for SPO model       

   !Grain-scale / rock-scale SPO
   call read_par_int(15,spograinmod)
   call read_par_int(15,sporockmod)
   call read_par_double(15,volfractrock(1))
   call read_par_double(15,volfractrock(2))
   call read_par_double(15,volfractrock(3))
   call read_par_double(15,volfractrock(4))
   call read_par_double(15,volfractrock(5))

   !Porosity from geodynamic model
   call read_par_int(15,meltspomod)
   call read_par_char(15,meltfilename)

   !Constant Cij
   call read_par_double6(15,Cback(1,1:6))
   call read_par_double6(15,Cback(2,1:6))
   call read_par_double6(15,Cback(3,1:6))
   call read_par_double6(15,Cback(4,1:6))
   call read_par_double6(15,Cback(5,1:6))
   call read_par_double6(15,Cback(6,1:6))
   call read_par_double(15,ro_back)
   call read_par_double6(15,Cinc(1,1:6))
   call read_par_double6(15,Cinc(2,1:6))
   call read_par_double6(15,Cinc(3,1:6))
   call read_par_double6(15,Cinc(4,1:6))
   call read_par_double6(15,Cinc(5,1:6))
   call read_par_double6(15,Cinc(6,1:6))
   call read_par_double(15,ro_incl)

   !DEM parameter
   call read_par_double(15,esa(1))
   call read_par_double(15,esa(2))
   call read_par_double(15,esa(3))
   call read_par_double(15,Vmax)   
   call read_par_double(15,Vstp)   

   !Angle between SPO fabric and max FSE or max compression
   call read_par_double(15,phi_spo)

   if(spomod > 0) then

      write(*,*)
      write(*,"(a)"),' SPO MODELLING IS ACTIVE'
      write(*,*)
      if(spomod == 1) then
         write(*,"(a,i1,a)"),' spomod :          ',spomod,' --> STILWE'
         write(*,"(a,i1)")  ,' spograinmod :     ',spograinmod
         write(*,"(a,i1)")  ,' sporockmod :      ',sporockmod
      end if
      if(spomod >  1) then
         write(*,"(a,i1,a)"),' spomod :          ',spomod,' --> DEM'
         write(*,"(a,i1)")  ,' meltspomod :      ',meltspomod
      end if
      write(*,"(a,f6.2)")   ,' phi_spo :  ',phi_spo
      write(*,*)

      !Convert degrees to radians
      phi_spo = phi_spo*deg2rad

      !!! STILWE models
      if(spomod == 1 .and. ptmod == 0) then

         write(*,"(a)"),' STILWE with constant Cijs and Vmax from inputspo.dat'

      else if(spomod == 1 .and. ptmod > 0 .and. spograinmod > 0 .and. sporockmod == 0) then

         write(*,"(a)"),' Grain-scale SPO. LPO fabrics are ignored'
         write(*,"(a)"),' STILWE with isotropic moduli and density from databases f(P,T)'
         write(*,*)
         write(*,"(a)")                         ,' Rock types           :   Dunite Harzb. Pyrol.  MORB  Pyrox.'
         write(*,"(a,f7.2,f7.2,f7.2,f7.2,f7.2)"),' Rock vol. fract. (%) : ',volfractrock

         if(sum(volfractrock) /= 100.0) then
            write(*,*)
            write(*,"(a,f8.2)"),' Sum of rock volume fractions should be 100 (%), but is instead',sum(volfractrock)
            write(*,*)
            stop
         end if

         if(volfractrock(1) > 0) then
            print *,' No grain-scale SPO for dunite as it is made by olivine only' 
            stop
         end if

      else if(spomod == 1 .and. ptmod > 0 .and. spograinmod == 0 .and. sporockmod > 0) then

         write(*,"(a)"),' Rock-scale SPO. LPO fabrics are ignored'
         write(*,"(a)"),' STILWE with isotropic moduli and density from databases f(P,T)'
         write(*,*)
         write(*,"(a)")                         ,' Rock types           :   Dunite Harzb. Pyrol.  MORB  Pyrox.'
         write(*,"(a,f7.2,f7.2,f7.2,f7.2,f7.2)"),' Rock vol. fract. (%) : ',volfractrock

         if(sum(volfractrock) /= 100.0) then
            write(*,*)
            write(*,"(a,f8.2)"),' Sum of rock volume fractions should be 100 (%), but is instead',sum(volfractrock)
            write(*,*)
            stop
         end if

      !!! DEM models
      else if(spomod == 2 .and. meltspomod == 0) then

         write(*,"(a)"),' DEM with constant Cij, ellipsoid shape and Vmax from inputspo.dat'
         write(*,*)
         write(*,'(a,2f10.3)') ' Vamx, Vstp: ',Vmax,Vstp
         write(*,*)
         write(*,'(a,3f10.3)') ' FSE semiaxes max, med, min: ',esa   

      else if(spomod == 2 .and. meltspomod >  0) then

         write(*,"(a)"),' DEM with constant Cij, ellipsoid shape and melt content from geodynamic model'
         write(*,*)
         write(*,'(a,2f10.3)') ' Vamx, Vstp: ',Vmax,Vstp
         write(*,*)
         write(*,'(a,3f10.3)') ' FSE semiaxes max, med, min: ',esa   

      else if(spomod == 3 .and. meltspomod == 0) then

         write(*,"(a)"),' DEM with Cij from LPO for matrix, constant Cij for inclusions, ellipsoid shape and Vmax from inputspo.dat'
         write(*,*)
         write(*,'(a,2f10.3)') ' Vmax, Vstp: ',Vmax,Vstp
         write(*,*)
         write(*,'(a,3f10.3)') ' FSE semiaxes max, med, min: ',esa   

      else if(spomod == 3 .and. meltspomod >  0) then

         write(*,"(a)"),' DEM with Cij from LPO for matrix, constant Cij for inclusions, ellipsoid shape and melt content from geodynamic model'
         write(*,*)
         write(*,'(a,2f10.3)') ' Vamx, Vstp: ',Vmax,Vstp
         write(*,*)
         write(*,'(a,3f10.3)') ' FSE semiaxes max, med, min: ',esa   

      else

         write(*,"(a)"),'NO allowed SPO model' 
         stop

      end if

      volfractrock = volfractrock / 100.0

      !write(*,*)
      !write(*,"(a)"),'--------------------------------------------------------'

   end if
 
   Is = 0d0 ; kron = 0d0

   !Kronecker delta
   kron(1,1)=1;
   kron(2,2)=1;
   kron(3,3)=1;
   
   !Symmetricic symmetric fourth-rank unit tensor
   DO ii = 1,3; DO jj = 1,3; DO kk = 1,3; DO ll=1,3
      Is(ii,jj,kk,ll) = Is(ii,jj,kk,ll) + 0.5*(kron(ii,kk)*kron(jj,ll)+kron(ii,ll)*kron(jj,kk));
   END DO; END DO; END DO; END DO

   !Make spherical grid
   pi = 3.141592653589793238462643383279
   stp=pi/100
   thetanum = INT(pi/stp)+1
   ALLOCATE(thetaspo(thetanum))
   DO ii = 1, thetanum
      thetaspo(ii) = stp*(ii)
   END DO

   phinum = INT(2*pi/stp)+1
   ALLOCATE(phispo(phinum))
   DO ii = 1, phinum
      phispo(ii) = stp*(ii)
   END DO

   IF(spomod == 1 .AND. ptmod > 0.AND. (spograinmod > 0 .OR. sporockmod > 0)) THEN

   !P-T grid
   tnum = 85
   pnum=1401
   ptnum = pnum*tnum

   !!! Grain-scale SPO  
   IF(spograinmod > 0) THEN

      numminer(:) = 21
      numminer(2) = 20 !Hartzburgite 
      numminermax = MAXVAL(numminer)
      ALLOCATE(rhominer(rocknum,tnum,pnum,numminermax),vsminer(rocknum,tnum,pnum,numminermax),&
           vpminer(rocknum,tnum,pnum,numminermax),volfractminer(rocknum,tnum,pnum,numminermax))

      !!! Load databases for minerals made with HeFESTo (Stixrude & Lithgow-Bertelloni, 2005, 2011, GJI) 
      CALL loadHefestograin

   END IF
   
   !!! Rock-scale SPO  
   IF(sporockmod > 0) THEN

      ALLOCATE(rhorock(rocknum,tnum,pnum),vsrock(rocknum,tnum,pnum),vprock(rocknum,tnum,pnum),srock(rocknum,tnum,pnum))

      !!! Load databases for minerals made with HeFESTo (Stixrude & Lithgow-Bertelloni, 2005, 2011, GJI) 
      CALL loadHefestorock

   END IF

   END IF
   
   END SUBROUTINE initspo  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  DEM model !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DEM(m,Vmaxm,Cdem)              

   USE comvar   

   IMPLICIT NONE

   INTEGER :: m,ii,jj,kk,ll,cyc,ncyc
   DOUBLE PRECISION :: vol_fract,Vmin,Vmaxm
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Vol
   DOUBLE PRECISION, DIMENSION(1000,6,6) :: Cdem
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: c0dem,Cback3,Cinc3,Ctemp,dC1,dC2,dC3,dC4

   Vmin = 0d0
   ncyc=INT((Vmaxm-Vmin)/Vstp)+1
   ALLOCATE(Vol(ncyc))
   DO ii = 1,ncyc
      Vol(ii) = Vmin + Vstp*(ii-1)
   END DO
       
   IF(spomod == 3) Cback = Sav(:,:,m) !Set matrix Cij

   !Convert from Voigt notation to 4th order tensor
   DO ii = 1,3; DO jj = 1,3; DO kk = 1,3; DO ll=1,3
      Cback3(ii,jj,kk,ll) = Cback(ijkl(ii,jj),ijkl(kk,ll))
      Cinc3(ii,jj,kk,ll)  = Cinc(ijkl(ii,jj),ijkl(kk,ll))
   END DO; END DO; END DO; END DO

   DO cyc=1,ncyc

      IF(Vol(cyc) .EQ. 0d0) THEN

         c0dem = Cback3

      ELSE IF(Vol(cyc) .EQ. 1d0) THEN

         c0dem = Cinc3

      ELSE

         vol_fract = log(1d0+Vstp/(1d0-Vol(cyc))) 

         !Calculate symmetrical tensor Green function
         !4th order Runge-Kutta integration scheme    
         Ctemp = c0dem 
         CALL greenfunction(Ctemp,Cinc3,vol_fract,dC1)
         Ctemp = c0dem + 0.5d0*dC1
         CALL greenfunction(Ctemp,Cinc3,vol_fract,dC2)
         Ctemp = c0dem + 0.5d0*dC2
         CALL greenfunction(Ctemp,Cinc3,vol_fract,dC3)
         Ctemp = c0dem + dC3
         CALL greenfunction(Ctemp,Cinc3,vol_fract,dC4)
         
         c0dem = c0dem + (dC1 + 2d0*dC2 + 2d0*dC3 + dC4)/6d0
 
      END IF
        
      !Convert to Voigt notation
      DO ii = 1, 6 ; DO jj = 1 , 6
         Cdem(cyc,ii,jj) = c0dem(l1(ii),l2(ii),l1(jj),l2(jj))
      END DO ; END DO

   END DO 

   DEALLOCATE(Vol)
   
   END SUBROUTINE DEM    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Double index 4th order tensor contraction                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE tensorcontraction(A,B,C)

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

   SUBROUTINE greenfunction(c0dem,Cinc3,vol_fract,dc0dem)

   USE comvar    

   IMPLICIT NONE

   INTEGER :: ii,jj,kk,ll
   DOUBLE PRECISION :: vol_fract             
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: c0dem,Cinc3,dc0dem,Ai,Ainv,dC,Gs,G,SdC
   DOUBLE PRECISION, DIMENSION(6,6) :: Ai_mandel,Ainv_mandel 
  
   !Calculate interaction tensor
   CALL interactiontensor(c0dem,G)

   ! Symmetric interaction tensor
   Gs = 0d0
   DO ii = 1,3; DO jj = 1,3; DO kk = 1,3; DO ll=1,3
      !Mainprice
      !Gs(ii,kk,jj,ll) = 0.5d0*(G(ii,kk,jj,ll)+G(jj,kk,ii,ll))
      !Hornby
      Gs(ii,jj,kk,ll) = 0.5d0*(G(ii,kk,jj,ll)+G(jj,kk,ii,ll))
   END DO; END DO; END DO; END DO
   
   !Difference in viscous tensor between inclusion and matrix 
   dC = Cinc3 - c0dem

   !Eshelby tensor = S = Jd*T*Cm = Ts*Cm 
   CALL tensorcontraction(Gs,dC,SdC)
   
   !Add symmetric identity 4th-order tensor
   Ai =  Is + SdC

   !Convert Ai to Mandel notation
   DO ii = 1, 6 ; DO jj = 1 , 6
      Ai_mandel(ii,jj) = Ai(l1(ii),l2(ii),l1(jj),l2(jj))*mandel_scale(ii,jj)
   END DO ; END DO
   
   !Invert Ai
   CALL inverse(Ai_mandel,Ainv_mandel,6)

   !Convert from Mandel notation to 4th order tensor
   DO ii = 1 , 3 ; DO jj = 1 , 3 ; DO kk = 1 , 3 ; DO ll = 1 , 3
      Ainv(ii,jj,kk,ll) = Ainv_mandel(ijkl(ii,jj),ijkl(kk,ll))/mandel_scale(ijkl(ii,jj),ijkl(kk,ll))
   END DO ; END DO ; END DO ; END DO
         
   !New composite stiffness tensor
   CALL tensorcontraction(dC,Ainv,dc0dem)
   
   !Multiply by volume fraction increment
   dc0dem = vol_fract*dc0dem

   END SUBROUTINE greenfunction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE interactiontensor(c0dem,G)
   
   USE comvar    

   IMPLICIT NONE

   INTEGER :: ii,jj,kk,ll,pp,qq,rr,ss,tt,uu
   DOUBLE PRECISION :: SurfArea,x(3)  
   DOUBLE PRECISION, DIMENSION(3,3) :: Kinv,K    
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: c0dem,G
  
   G = 0d0
   SurfArea = 0d0
   DO pp=1,thetanum; DO qq=1,phinum
      
      !Directions
      x(1)=sin(thetaspo(pp))*cos(phispo(qq))/esa(1)
      x(2)=sin(thetaspo(pp))*sin(phispo(qq))/esa(2)
      x(3)=cos(thetaspo(pp))/esa(3)
      
      !Christoffel stiffness tensor K
      K = 0d0
      DO rr = 1,3; DO ss = 1,3; DO tt = 1,3; DO uu=1,3
         K(rr,tt)=K(rr,tt)+c0dem(rr,ss,tt,uu)*x(ss)*x(uu)
      END DO; END DO; END DO; END DO
       
      !Find inverse of K
      CALL inverse(K,Kinv,3)

      DO ii = 1,3; DO kk = 1,3; DO jj = 1,3; DO ll=1,3
         !Hornby, 1994, Geophysics
         G(ii,jj,kk,ll) = G(ii,jj,kk,ll) + Kinv(ii,jj)*x(kk)*x(ll)*sin(thetaspo(pp))!*stp*stp
         !Mainprice, 2015, Treat. Geophys.
         !G(ii,kk,jj,ll) = G(ii,kk,jj,ll) + Kinv(ii,jj)*x(kk)*x(ll)*sin(thetaspo(pp))!*stp*stp
      END DO; END DO; END DO; END DO
      SurfArea = SurfArea + sin(thetaspo(pp))!*stp*stp
   END DO; END DO

   G = G/SurfArea

   END SUBROUTINE interactiontensor
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  STILWE    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE stilwe(m,C1)              

   USE comvar

   IMPLICIT NONE

   INTEGER :: m,n1,n2
   DOUBLE PRECISION :: tmin,tstp,pmin,pstp,ee,nnn
   DOUBLE PRECISION :: A,B,C,F,Lh,N,R,S,T,RRhom
   DOUBLE PRECISION :: A0,A1,A2,A3       
   DOUBLE PRECISION :: B0,B1,B2,B3       
   DOUBLE PRECISION :: CC0,CC1,CC2,CC3       
   DOUBLE PRECISION :: F0,F1,F2,F3       
   DOUBLE PRECISION :: LL0,LL1,LL2,LL3       
   DOUBLE PRECISION :: NN0,NN1,NN2,NN3       
   DOUBLE PRECISION :: RR0,RR1,RR2,RR3       
   DOUBLE PRECISION :: K1,K2,mu1,mu2,lambda1,lambda2
   DOUBLE PRECISION :: theta1,theta2,fract1,fract2
   DOUBLE PRECISION, DIMENSION(6,6) :: C1

   !!! Smoothed, transversely isotropic, long-wavelength equivalent medium
   !!! (Backus, 1962, JGR)

   IF(ptmod > 0 .AND. (spograinmod > 0 .OR. sporockmod > 0)) THEN

   tmin = 300 ; tstp = 50
   pmin = 0 ; pstp = 0.1

   !Find nearest nodes in the P,T grid
   ! ABCD-4Cell Number 
   ee = (mtk(m)-tmin)/tstp
   IF(ee < 0) ee=0
   IF(ee > REAL(tnum)) ee=REAL(tnum)
   nnn=(mpgpa(m)-pmin)/pstp
   IF(nnn < 0) nnn=0
   IF(nnn > REAL(pnum)) nnn=REAL(pnum)
   n1=FLOOR(ee) + 1
   IF(n1 < 1) n1=1
   IF(n1 > tnum-1) n1=tnum-1
   n2=FLOOR(nnn)+1
   IF(n2 < 1) n2=1
   IF(n2 > pnum-1) n2=pnum-1
   ! Calc normalized distances 
   ee=(ee-REAL(n1-1))
   nnn=(nnn-REAL(n2-1))

   ! Ro values
   ! 0 2
   ! 1 3 
   A0=AA(n1  ,n2  )
   A1=AA(n1  ,n2+1)
   A2=AA(n1+1,n2  )
   A3=AA(n1+1,n2+1)

   A = ((A0*(1.0-nnn)+A1*nnn)*(1.0-ee)+(A2*(1.0-nnn)+A3*nnn)*ee) 

   B0=BB(n1  ,n2  )
   B1=BB(n1  ,n2+1)
   B2=BB(n1+1,n2  )
   B3=BB(n1+1,n2+1)

   B = ((B0*(1.0-nnn)+B1*nnn)*(1.0-ee)+(B2*(1.0-nnn)+B3*nnn)*ee) 

   CC0=CC(n1  ,n2  )
   CC1=CC(n1  ,n2+1)
   CC2=CC(n1+1,n2  )
   CC3=CC(n1+1,n2+1)

   C = ((CC0*(1.0-nnn)+CC1*nnn)*(1.0-ee)+(CC2*(1.0-nnn)+CC3*nnn)*ee) 

   F0=FF(n1  ,n2  )
   F1=FF(n1  ,n2+1)
   F2=FF(n1+1,n2  )
   F3=FF(n1+1,n2+1)

   F = ((F0*(1.0-nnn)+F1*nnn)*(1.0-ee)+(F2*(1.0-nnn)+F3*nnn)*ee) 

   LL0=LLL(n1  ,n2  )
   LL1=LLL(n1  ,n2+1)
   LL2=LLL(n1+1,n2  )
   LL3=LLL(n1+1,n2+1)

   Lh = ((LL0*(1.0-nnn)+LL1*nnn)*(1.0-ee)+(LL2*(1.0-nnn)+LL3*nnn)*ee) 

   NN0=NN(n1  ,n2  )
   NN1=NN(n1  ,n2+1)
   NN2=NN(n1+1,n2  )
   NN3=NN(n1+1,n2+1)

   N = ((NN0*(1.0-nnn)+NN1*nnn)*(1.0-ee)+(NN2*(1.0-nnn)+NN3*nnn)*ee) 

   RR0=RRho(n1  ,n2  )
   RR1=RRho(n1  ,n2+1)
   RR2=RRho(n1+1,n2  )
   RR3=RRho(n1+1,n2+1)

   RRhom = ((RR0*(1.0-nnn)+RR1*nnn)*(1.0-ee)+(RR2*(1.0-nnn)+RR3*nnn)*ee) 

   ELSE

   !!! Phase 1 
   C1 = Cback
   mu1 = (C1(1,1)+C1(2,2)+C1(3,3)-(C1(1,2)+C1(1,3)+C1(2,3))+3*(C1(4,4)+C1(5,5)+C1(6,6)))/15
   K1  = (C1(1,1)+C1(2,2)+C1(3,3)+2*(C1(1,2)+C1(1,3)+C1(2,3)))/9
   lambda1= K1-2*mu1/3
   
   !!! Phase 2
   C1 = Cinc 
   mu2 = (C1(1,1)+C1(2,2)+C1(3,3)-(C1(1,2)+C1(1,3)+C1(2,3))+3*(C1(4,4)+C1(5,5)+C1(6,6)))/15
   K2  = (C1(1,1)+C1(2,2)+C1(3,3)+2*(C1(1,2)+C1(1,3)+C1(2,3)))/9
   lambda2= K2-2*mu2/3

   theta1=mu1/(lambda1+2d0*mu1)
   theta2=mu2/(lambda2+2d0*mu2)
   fract1 = 1d0 - Vmax
   fract2 = Vmax

   R=fract1*theta1/mu1+fract2*theta2/mu2
   S=fract1*theta1*mu1+fract2*theta2*mu2
   T=fract1*theta1+fract2*theta2

   N=fract1*mu1+fract2*mu2 !Vsh^2*rho

   Lh=(fract1/mu1+fract2/mu2) !Vsv^2*rho
   Lh=1/Lh

   IF(mu1==0 .OR. mu2==0) THEN
      IF(mu1==0) print *,'WARNING: shear modulus of medium 1 is 0, then L = 0 and radial anisotropy is infinite'
      IF(mu2==0) print *,'WARNING: shear modulus of medium 2 is 0, then L = 0 and radial anisotropy is infinite'
      stop
   END IF

   C=1d0/R !Vph^2*rho
   F=C*(1d0-2d0*T)
   B=2d0*N - 4d0*S + C*(1d0-2d0*T)**2d0
   A=B+2d0*N !Vpv^2*rho

   RRhom = ro_back*fract1 + ro_incl*fract2

   END IF
   
   !New elastic tensor
   C1(:,:) = 0d0
   C1(1,1) = A ; C1(1,2) = B ; C1(1,3) = F
   C1(2,1) = B ; C1(2,2) = A ; C1(2,3) = F
   C1(3,1) = F ; C1(3,2) = F ; C1(3,3) = C
   C1(4,4) = Lh; C1(5,5) = Lh; C1(6,6) = N

   !New density
   rho(m) = RRhom

   END SUBROUTINE stilwe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Load/Save dataset, format double precision                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE loadHefestograin

   USE comvar
   USE hdf5

   IMPLICIT NONE

   INTEGER :: i,j,k,rt,gi,dum_int(2)
   DOUBLE PRECISION :: dum(3),A,B,C,F,Lh,N,S,R,T,RRhom
   DOUBLE PRECISION :: volfractmin,rhomin,vsmin,vpmin,mu,theta
   DOUBLE PRECISION, DIMENSION(ptnum,21) :: Rhodum,Vpdum,Vsdum,Vdum
   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl,memtype !Handles
   INTEGER     ::   error  ! Error flag

   !Read mineral density, vs, vp and volume fraction         
   dum_int(1)=ptnum
   dum_int(2)=21

   DO rt = 1, rocknum

   IF(volfractrock(rt) /= 0) THEN

   !Dunite
   IF(rt == 1) THEN

      !Dunite   
      CALL H5open_f (error)
      CALL H5Fopen_f("../DATABASES/MMA-EoS/dunite.h5", H5F_ACC_RDONLY_F, file_id, error)
      CALL H5Gopen_f(file_id, "/minerals", group_id, error)

      !Read mineral density, vs, vp and volume fraction         
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Rhodum(:,:),'Rho',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vpdum(:,:),'Vp',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vsdum(:,:),'Vs',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vdum(:,:),'Volume',0)

      CALL H5Gclose_f(group_id, error)
      CALL H5Fclose_f(file_id, error)
      CALL H5close_f(error)
      !Close FORTRAN interface.
  
      DO j = 1 , tnum ; DO i = 1 , pnum
         gi = j + (i-1)*tnum
         vsminer(rt,j,i,:) = Vsdum(gi,:)
         vpminer(rt,j,i,:) = Vpdum(gi,:)
         rhominer(rt,j,i,:) = Rhodum(gi,:)
         volfractminer(rt,j,i,:) = Vdum(gi,:)
         DO k = 1, numminer(rt)
            IF (isnan(vsminer(rt,j,i,k))) vsminer(rt,j,i,k)=0d0
            IF (isnan(vpminer(rt,j,i,k))) vpminer(rt,j,i,k)=0d0
            IF (isnan(rhominer(rt,j,i,k))) rhominer(rt,j,i,k)=0d0
            IF (isnan(volfractminer(rt,j,i,k))) volfractminer(rt,j,i,k)=0d0
         END DO
      END DO ; END DO

   END IF
   
   !Harzburgite (Xu et al., 2008, EPSL)
   IF(rt == 2) THEN

      !Read mineral density, vs, vp and volume fraction         
      CALL H5open_f (error)
      CALL H5Fopen_f("../DATABASES/MMA-EoS/hartzburgite.h5", H5F_ACC_RDONLY_F, file_id, error)
      CALL H5Gopen_f(file_id, "/minerals", group_id, error)

      !Read mineral density, vs, vp and volume fraction         
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Rhodum(:,:),'Rho',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vpdum(:,:),'Vp',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vsdum(:,:),'Vs',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vdum(:,:),'Volume',0)

      CALL H5Gclose_f(group_id, error)
      CALL H5Fclose_f(file_id, error)
      CALL H5close_f(error)
      !Close FORTRAN interface.

      DO j = 1 , tnum ; DO i = 1 , pnum
         gi = j + (i-1)*tnum
         vsminer(rt,j,i,:) = Vsdum(gi,:)
         vpminer(rt,j,i,:) = Vpdum(gi,:)
         rhominer(rt,j,i,:) = Rhodum(gi,:)
         volfractminer(rt,j,i,:) = Vdum(gi,:)
         DO k = 1, numminer(rt)
            IF (isnan(vsminer(rt,j,i,k))) vsminer(rt,j,i,k)=0d0
            IF (isnan(vpminer(rt,j,i,k))) vpminer(rt,j,i,k)=0d0
            IF (isnan(rhominer(rt,j,i,k))) rhominer(rt,j,i,k)=0d0
            IF (isnan(volfractminer(rt,j,i,k))) volfractminer(rt,j,i,k)=0d0
         END DO
      END DO ; END DO

   END IF
   
   !Pyrolite (Xu et al., 2008, EPSL)
   IF(rt == 3) THEN

      !Read mineral density, vs, vp and volume fraction         
      CALL H5open_f (error)
      CALL H5Fopen_f("../DATABASES/MMA-EoS/pyrolite.h5", H5F_ACC_RDONLY_F, file_id, error)
      CALL H5Gopen_f(file_id, "/minerals", group_id, error)

      !Read mineral density, vs, vp and volume fraction         
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Rhodum(:,:),'Rho',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vpdum(:,:),'Vp',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vsdum(:,:),'Vs',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vdum(:,:),'Volume',0)

      CALL H5Gclose_f(group_id, error)
      CALL H5Fclose_f(file_id, error)
      CALL H5close_f(error)
      !Close FORTRAN interface.

      DO j = 1 , tnum ; DO i = 1 , pnum
         gi = j + (i-1)*tnum
         vsminer(rt,j,i,:) = Vsdum(gi,:)
         vpminer(rt,j,i,:) = Vpdum(gi,:)
         rhominer(rt,j,i,:) = Rhodum(gi,:)
         volfractminer(rt,j,i,:) = Vdum(gi,:)
         DO k = 1, numminer(rt)
            IF (isnan(vsminer(rt,j,i,k))) vsminer(rt,j,i,k)=0d0
            IF (isnan(vpminer(rt,j,i,k))) vpminer(rt,j,i,k)=0d0
            IF (isnan(rhominer(rt,j,i,k))) rhominer(rt,j,i,k)=0d0
            IF (isnan(volfractminer(rt,j,i,k))) volfractminer(rt,j,i,k)=0d0
         END DO
      END DO ; END DO

   END IF
   
   !Basalt (MORB) (Xu et al., 2008, EPSL))
   IF(rt == 4) THEN

      !Read mineral density, vs, vp and volume fraction         
      CALL H5open_f (error)
      CALL H5Fopen_f("../DATABASES/MMA-EoS/morb.h5", H5F_ACC_RDONLY_F, file_id, error)
      CALL H5Gopen_f(file_id, "/minerals", group_id, error)

      !Read mineral density, vs, vp and volume fraction         
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Rhodum(:,:),'Rho',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vpdum(:,:),'Vp',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vsdum(:,:),'Vs',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vdum(:,:),'Volume',0)

      CALL H5Gclose_f(group_id, error)
      CALL H5Fclose_f(file_id, error)
      CALL H5close_f(error)
      !Close FORTRAN interface.

      DO j = 1 , tnum ; DO i = 1 , pnum
         gi = j + (i-1)*tnum
         vsminer(rt,j,i,:) = Vsdum(gi,:)
         vpminer(rt,j,i,:) = Vpdum(gi,:)
         rhominer(rt,j,i,:) = Rhodum(gi,:)
         volfractminer(rt,j,i,:) = Vdum(gi,:)
         DO k = 1, numminer(rt)
            IF (isnan(vsminer(rt,j,i,k))) vsminer(rt,j,i,k)=0d0
            IF (isnan(vpminer(rt,j,i,k))) vpminer(rt,j,i,k)=0d0
            IF (isnan(rhominer(rt,j,i,k))) rhominer(rt,j,i,k)=0d0
            IF (isnan(volfractminer(rt,j,i,k))) volfractminer(rt,j,i,k)=0d0
         END DO
      END DO ; END DO

   END IF
   
   !Pyroxenite (Hirschmann, 2006, Geology)
   IF(rt == 5) THEN

      !Read mineral density, vs, vp and volume fraction         
      CALL H5open_f (error)
      CALL H5Fopen_f("../DATABASES/MMA-EoS/pyroxenite.h5", H5F_ACC_RDONLY_F, file_id, error)
      CALL H5Gopen_f(file_id, "/minerals", group_id, error)

      !Read mineral density, vs, vp and volume fraction         
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Rhodum(:,:),'Rho',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vpdum(:,:),'Vp',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vsdum(:,:),'Vs',0)
      CALL loadsave_double(0,2,group_id,dum_int,H5T_NATIVE_DOUBLE,Vdum(:,:),'Volume',0)

      CALL H5Gclose_f(group_id, error)
      CALL H5Fclose_f(file_id, error)
      CALL H5close_f(error)
      !Close FORTRAN interface.

      DO j = 1 , tnum ; DO i = 1 , pnum
         gi = j + (i-1)*tnum
         vsminer(rt,j,i,:) = Vsdum(gi,:)
         vpminer(rt,j,i,:) = Vpdum(gi,:)
         rhominer(rt,j,i,:) = Rhodum(gi,:)
         volfractminer(rt,j,i,:) = Vdum(gi,:)
         DO k = 1, numminer(rt)
            IF (isnan(vsminer(rt,j,i,k))) vsminer(rt,j,i,k)=0d0
            IF (isnan(vpminer(rt,j,i,k))) vpminer(rt,j,i,k)=0d0
            IF (isnan(rhominer(rt,j,i,k))) rhominer(rt,j,i,k)=0d0
            IF (isnan(volfractminer(rt,j,i,k))) volfractminer(rt,j,i,k)=0d0
         END DO
      END DO ; END DO

   END IF
   
   END IF

   END DO

   !Compute STILWE f(P,T)
   ALLOCATE(AA(tnum,pnum),BB(tnum,pnum),CC(tnum,pnum),FF(tnum,pnum),LLL(tnum,pnum),NN(tnum,pnum),RRho(tnum,pnum))

   DO j = 1 , tnum ; DO i = 1 , pnum

      N = 0d0 ; Lh = 0d0 ; R = 0d0 ; S = 0d0 ; T = 0d0

      DO rt = 1, rocknum

         IF(volfractrock(rt) /= 0) THEN

         DO k = 1, numminer(rt)
         IF(volfractminer(rt,j,i,k) .GT. 0d0) THEN
            !Mineral volume fraction
            volfractmin = volfractminer(rt,j,i,k)*volfractrock(rt)
            !Mineral density
            rhomin = rhominer(rt,j,i,k)
            !Mineral Vs
            vsmin = vsminer(rt,j,i,k) 
            !Mineral Vp
            vpmin = vpminer(rt,j,i,k)
            !Shear modulus
            mu = rhomin*vsmin**2d0 
            !Vs^2/Vp^2
            theta = (vsmin/vpmin)**2d0
            !Love parameters
            N = N + volfractmin*mu
            Lh = Lh + volfractmin/mu
            R = R + volfractmin*theta/mu
            S = S + volfractmin*theta*mu
            T = T + volfractmin*theta
            RRhom = RRhom + volfractmin*rhomin
         END IF
         END DO
         
         END IF

      END DO

      Lh=1d0/Lh;
      C=1d0/R !Vph^2*rho
      F=C*(1d0-2d0*T)
      B=2d0*N - 4d0*S + C*(1d0-2d0*T)**2d0
      A=B+2d0*N !Vpv^2*rho

      !Polarization anisotropy
      AA(j,i) = A
      BB(j,i) = B
      CC(j,i) = C
      FF(j,i) = F
      LLL(j,i) = Lh
      NN(j,i) = N
      RRho(j,i) = RRhom

   END DO ; END DO

   END SUBROUTINE loadHefestograin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Load/Save dataset, format double precision                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE loadHefestorock

   USE comvar
   USE hdf5

   IMPLICIT NONE

   INTEGER :: i,j,rt,gi
   DOUBLE PRECISION :: dum(3),A,B,C,F,Lh,N,S,R,T,RRhom
   DOUBLE PRECISION :: volfractmin,rhomin,vsmin,vpmin,mu,theta
   DOUBLE PRECISION, DIMENSION(5,ptnum) :: Rhodum,Vpdum,Vsdum
   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl,memtype !Handles
   INTEGER     ::   error  ! Error flag

   !Read mineral density, vs, vp

   !Dunite   
   CALL H5open_f (error)
   CALL H5Fopen_f("../DATABASES/MMA-EoS/dunite.h5", H5F_ACC_RDONLY_F, file_id, error)
   CALL H5Gopen_f(file_id, "/rock", group_id, error)

   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Rhodum(1,:),'Rho',0)
   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vpdum(1,:),'Vp',0)
   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vsdum(1,:),'Vs',0)

   CALL H5Gclose_f(group_id, error)
   CALL H5Fclose_f(file_id, error)
   CALL H5close_f(error)
   !Close FORTRAN interface.
  
   !Hartzburgite
   CALL H5open_f (error)
   CALL H5Fopen_f("../DATABASES/MMA-EoS/hartzburgite.h5", H5F_ACC_RDONLY_F, file_id, error)
   CALL H5Gopen_f(file_id, "/rock", group_id, error)

   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Rhodum(2,:),'Rho',0)
   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vpdum(2,:),'Vp',0)
   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vsdum(2,:),'Vs',0)

   CALL H5Gclose_f(group_id, error)
   CALL H5Fclose_f(file_id, error)
   CALL H5close_f(error)
   !Close FORTRAN interface.
  
   !Pyrolite
   CALL H5open_f (error)
   CALL H5Fopen_f("../DATABASES/MMA-EoS/pyrolite.h5", H5F_ACC_RDONLY_F, file_id, error)
   CALL H5Gopen_f(file_id, "/rock", group_id, error)

   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Rhodum(3,:),'Rho',0)
   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vpdum(3,:),'Vp',0)
   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vsdum(3,:),'Vs',0)

   CALL H5Gclose_f(group_id, error)
   CALL H5Fclose_f(file_id, error)
   CALL H5close_f(error)
   !Close FORTRAN interface.

   !Basalt
   CALL H5open_f (error)
   CALL H5Fopen_f("../DATABASES/MMA-EoS/morb.h5", H5F_ACC_RDONLY_F, file_id, error)
   CALL H5Gopen_f(file_id, "/rock", group_id, error)

   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Rhodum(4,:),'Rho',0)
   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vpdum(4,:),'Vp',0)
   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vsdum(4,:),'Vs',0)

   CALL H5Gclose_f(group_id, error)
   CALL H5Fclose_f(file_id, error)
   CALL H5close_f(error)
   !Close FORTRAN interface.

   !Pyroxenite
   CALL H5open_f (error)
   CALL H5Fopen_f("../DATABASES/MMA-EoS/pyroxenite.h5", H5F_ACC_RDONLY_F, file_id, error)
   CALL H5Gopen_f(file_id, "/rock", group_id, error)

   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Rhodum(5,:),'Rho',0)
   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vpdum(5,:),'Vp',0)
   CALL loadsave_double(0,1,group_id,ptnum,H5T_NATIVE_DOUBLE,Vsdum(5,:),'Vs',0)

   CALL H5Gclose_f(group_id, error)
   CALL H5Fclose_f(file_id, error)
   CALL H5close_f(error)
   !Close FORTRAN interface.

   !Compute STILWE f(P,T)
   IF(spograinmod == 0) ALLOCATE(AA(tnum,pnum),BB(tnum,pnum),CC(tnum,pnum),FF(tnum,pnum),LLL(tnum,pnum),NN(tnum,pnum),RRho(tnum,pnum))

   DO j = 1 , tnum ; DO i = 1 , pnum

      N = 0d0 ; Lh = 0d0 ; R = 0d0 ; S = 0d0 ; T = 0d0 ; RRhom = 0d0

      DO rt = 1, rocknum

         gi = j + (i-1)*tnum
         rhorock(rt,j,i) = Rhodum(rt,gi)
         vsrock(rt,j,i)  = Vpdum(rt,gi) 
         vprock(rt,j,i)  = Vsdum(rt,gi)

         IF(volfractrock(rt) > 0d0) THEN
            !Mineral density
            rhomin = rhorock(rt,j,i)
            !Mineral Vs
            vsmin = vsrock(rt,j,i)
            !Mineral Vp
            vpmin = vprock(rt,j,i)
            !Shear modulus
            mu = rhomin*vsmin**2d0
            !Vs^2/Vp^2
            theta = (vsmin/vpmin)**2d0
            !Love parameters
            N = N + volfractrock(rt)*mu
            Lh = Lh + volfractrock(rt)/mu
            R = R + volfractrock(rt)*theta/mu
            S = S + volfractrock(rt)*theta*mu
            T = T + volfractrock(rt)*theta
            RRhom = RRhom + volfractrock(rt)*rhomin
         END IF

      END DO

      Lh=1d0/Lh;
      C=1d0/R !Vph^2*rho
      F=C*(1d0-2d0*T)
      B=2d0*N - 4d0*S + C*(1d0-2d0*T)**2d0
      A=B+2d0*N !Vpv^2*rho

      !Polarization anisotropy
      AA(j,i) = A
      BB(j,i) = B
      CC(j,i) = C
      FF(j,i) = F
      LLL(j,i) = Lh
      NN(j,i) = N
      RRho(j,i) = RRhom

   END DO ; END DO

   END SUBROUTINE loadHefestorock

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
