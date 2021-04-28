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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine lagrangian_output, saving fse, hexagonal symmetry axis, etc. !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE lagrangian_output(dt_str4,cijkl_dir,output_dir)

   USE comvar
   USE omp_lib
   USE hdf5   

   IMPLICIT NONE

   INTEGER :: m,t,tid,nt,i,i1,i2,i3,j,j1,j2,k,ll,dum_int(3),ln_fse_num,nx(2),ti(1)
   !Buffers
   INTEGER, DIMENSION(:), ALLOCATABLE :: idx  
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: mxyz,vpmax,dvsmax,LS
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi_fse_max0,phi_fse_min0
   ! orientation of the long axis of the FSE 

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi_a,phi_a0
   ! average orientation of a-axis 

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: perc_a,perc_a0
   ! percentage of S wave anisotropy
   
   DOUBLE PRECISION :: vp_anis,dvs_anis
   DOUBLE PRECISION, DIMENSION(3) :: evals,vpvect,dvsvect,acsavg,c2
   DOUBLE PRECISION, DIMENSION(3,3) :: evects,acsV,cV

   DOUBLE PRECISION :: phi1,theta,phi2,phi,dm,Vmaxm,a0,w3,e3,vortnum
   DOUBLE PRECISION :: pphi1,ttheta,pphi2,lr,cr
   DOUBLE PRECISION, DIMENSION(3,3) :: acs,acs1,Rotm1,lrot
   ! eulerian angles

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phix0,phiy0,phiz0,phiz1
   ! spherical to cartesian coordinates arrays
   
   DOUBLE PRECISION, DIMENSION(6,6) :: Cstilwe
   DOUBLE PRECISION, DIMENSION(1000,6,6) :: Cdem

   CHARACTER (500) :: filename,filenamexmf,filenamecijkl,command
   CHARACTER (len=*) :: cijkl_dir,output_dir
   CHARACTER (len(trim(cijkl_dir))) :: str1
   CHARACTER (len(trim(output_dir))) :: str2
   CHARACTER (4) :: dt_str4
   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl,memtype ! Handles
   INTEGER     ::   error  ! Error flag
   DOUBLE PRECISION :: dum_db(5)
   
   str1=trim(cijkl_dir)
   str2=trim(output_dir)

   filename = str2//'lagrangian'//dt_str4//'.h5'

!!! Select aggregates with sufficient deformation
   ln_fse_num = 0
   DO m = 1,marknum
      IF(ln_fse(m)>=ln_fse_min .AND. Sav(1,1,m) > 0) THEN 
      IF(mx1(m) >= n1first .AND. mx1(m) <= n1last) THEN
      IF(mx2(m) >= n2first .AND. mx2(m) <= n2last) THEN
      IF(mx3(m) >= n3first .AND. mx3(m) <= n3last) THEN
      IF(uppermantlemod==0 .OR. (uppermantlemod>0 .AND. rocktype(m)==1)) THEN
         ln_fse_num = ln_fse_num + 1
      END IF;END IF;END IF;END IF;END IF
   END DO

   ALLOCATE(idx(ln_fse_num))

   ln_fse_num = 0
   DO m = 1,marknum
      IF(ln_fse(m)>=ln_fse_min .AND. Sav(1,1,m) > 0) THEN
      IF(mx1(m) >= n1first .AND. mx1(m) <= n1last) THEN
      IF(mx2(m) >= n2first .AND. mx2(m) <= n2last) THEN
      IF(mx3(m) >= n3first .AND. mx3(m) <= n3last) THEN
      IF(uppermantlemod==0 .OR. (uppermantlemod>0 .AND. rocktype(m)==1)) THEN
         ln_fse_num = ln_fse_num + 1
         idx(ln_fse_num) = m
      END IF;END IF;END IF;END IF;END IF
   END DO

   IF(uppermantlemod == 0) write(*,'(a,i0)'),' Number of aggregates sufficiently deformed (ln_fse > ln_fse_min) = ',ln_fse_num
   IF(uppermantlemod > 0)  write(*,'(a,i0)'),' Number of UPPER MANTLE aggregates sufficiently deformed (ln_fse > ln_fse_min) = ',ln_fse_num
   print *

   IF(ln_fse_num==0) THEN
      print *,''
      print *,'No sufficiently deformed aggregate'
      !DEALLOCATE(phi_fse_min,phi_fse_max,ln_fse,idx)
      RETURN
   END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Save min axis of FSE !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(fseminmod .NE. 0) THEN
  
   ALLOCATE(phi_fse_min0(3,ln_fse_num))

   DO m = 1,ln_fse_num
      
      i=idx(m)
      
      !fse_min_x
      phi_fse_min0(1,m) = ln_fse(i)*Rotm(1,1,i)
      !fse_min_y
      phi_fse_min0(2,m) = ln_fse(i)*Rotm(2,1,i)
      !fse_min_z
      phi_fse_min0(3,m) = ln_fse(i)*Rotm(3,1,i)

    END DO

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Save max axis of FSE !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(fsemaxmod .NE. 0) THEN
   
   ALLOCATE(phi_fse_max0(3,ln_fse_num))

   DO m = 1,ln_fse_num
      
      i=idx(m)
      
      !fse_min_x
      phi_fse_max0(1,m) = ln_fse(i)*Rotm(1,3,i)
      !fse_min_y
      phi_fse_max0(2,m) = ln_fse(i)*Rotm(2,3,i)     
      !fse_min_z
      phi_fse_max0(3,m) = ln_fse(i)*Rotm(3,3,i)      

    END DO

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  SPO calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(spomod > 0) THEN
  
   !!! Read input file infos for DEM model and allocate TD database matrices
   CALL initspo

   !$omp parallel    
   nt = OMP_GET_NUM_THREADS()
   !$omp end parallel

   print *,' '
   write(*,'(a,i0)'),' Number of threads = ',nt
   print *,' '

   !!! Reset all elastic tensors to check only the effect of SPO 
   !!! for grain-scale or rock-scale layering
   IF(spomod == 1 .and. ptmod > 0 .and. (spograinmod > 0 .and. sporockmod > 0)) Sav = 0d0

   !Load velocity field and compute Dij 
   IF(spomod > 1) THEN

      CALL load(dt_str4,cijkl_dir)

      !Read porosity from geodynamic model
      IF(meltspomod > 0) CALL readmeltcontent(cijkl_dir)    

      ALLOCATE(l(nt,3,3),e(nt,3,3),epsnot(nt))

   END IF

   !!! SPO from constant tensors defined in inputspo.dat
   IF(spomod == 1 .AND. ptmod == 0) CALL stilwe(1,Cstilwe) 
   IF(spomod == 2) CALL DEM(1,Vmax,Cdem) 

   !$omp parallel & 
   !$omp shared(idx,Sav,meltmark,Rotm,rho,mx1,mx2,mx3,mYY,Ui,Ux,Dij,e,l,epsnot) &
   !$omp private(i,j,m,tid,i1,i2,i3,ti,dm,Vmaxm,phi1,theta,phi2,acs,pphi1,ttheta,pphi2,acs1,Rotm1,evects,evals,c2,lrot,w3,a0,Cdem) &    
   !$omp firstprivate(ln_fse_num,spomod,ptmod,meltspomod,spograinmod,sporockmod,dimensions,Vstp,Cstilwe,phi_spo,ro_back,ro_incl)
   !$omp do schedule(guided,8)
   DO m = 1,ln_fse_num
      
      i = idx(m)  

      !Assign constant tensors defined in inputspo.dat
      !STILWE (Backus, 1962, JGR)
      IF(spomod == 1) THEN

         IF(ptmod == 0) THEN
        
            Sav(:,:,i) = Cstilwe 

         !Use elastic tensor f(P,T) taken from HeFESTo
         !STILWE (Backus, 1962, JGR)
         ELSEIF(ptmod > 0 .AND. (spograinmod > 0 .OR. sporockmod > 0)) THEN
   
            CALL stilwe(i,Sav(:,:,i)) 
 
         END IF

         !Save FSE semiaxes orientation
         !1st column: orientaton mininum axis
         !2nd column: orientaton medium  axis
         !3rd column: orientaton maximum axis
         Rotm1(:,:) = Rotm(:,:,i)
         IF(phi_spo /=0) CALL rot3Dspo(Rotm(:,:,i),Rotm1(:,:),Rotm(:,2,i),phi_spo)
          
         !Rotate tensor normal to min axis of FSE
         phi1 = atan2(Rotm1(1,1),-Rotm1(2,1)) !Angle of FSE max semiaxis with x1 axis 
         theta = acos(Rotm1(3,1)) !Angle of FSE min semiaxis with x3 axis
         phi2  = atan2(Rotm1(3,3),Rotm1(3,2))
         CALL rotmatrixZXZ(phi1,theta,phi2,acs) !Find rot. matrix
         CALL rotvoigt(Sav(:,:,i),acs,Sav(:,:,i)) !Rotate

      END IF

      !DEM (Mainprice, 1997, EPSL)
      IF(spomod == 2 .OR. (spomod == 3 .AND. (meltspomod == 0 .OR. (meltspomod > 0 .AND. meltmark(i) > Vstp)))) THEN

         tid = OMP_GET_THREAD_NUM() + 1

         !Calculate strain rate tensor
         IF(dimensions == 2) THEN
            CALL upperleft2Dn(mx1(i),mx2(i),i1,i2)
            CALL gradientcalc2D(tid,mx1(i),mx2(i),i1,i2,mYY(i))
         END IF
         IF(dimensions == 3) THEN
            CALL upperleft3Dn(mx1(i),mx2(i),mx3(i),i1,i2,i3)
            CALL gradientcalc(tid,mx1(i),mx2(i),mx3(i),i1,i2,i3,mYY(i))
         END IF
         CALL DSYEVQ3(e(tid,:,:)/epsnot(tid),evects,evals)

         !Order from smallest to largest semiaxis
         DO j = 1,3
            ti = MINLOC(evals) ;acs(:,j) = evects(:,ti(1)) ; c2(j) = evals(ti(1)) ; evals(ti(1))= 1d60
         END DO

         !Save semiaxes orientation
         !1st column: orientaton compressive  axis
         !2nd column: orientaton intermediate axis
         !3rd column: orientaton extensional  axis
         Rotm1 = acs

         !Rotate tensor such that compressive (min) and extensional (max)
         !strain rate axes coincide with the inclusion x1 and x3 axes
         phi1  = atan2(Rotm1(1,3),-Rotm1(2,3)) !Angle of compressive strain rate axis with x1 axis 
         theta =  acos(Rotm1(3,3)) !Angle of extensional strain rate axis with x3 axis
         phi2  = atan2(Rotm1(3,1),Rotm1(3,2))

         !Additional rotation from sigma 1
         IF(phi_spo /=0) THEN

            !Rotate velocity gradient tensor parallel to sigma 2
            pphi1  = atan2(Rotm1(1,2),-Rotm1(2,2))
            ttheta =  acos(Rotm1(3,2))
            pphi2  = 0d0
            CALL rotmatrixZXZ(pphi1,ttheta,pphi2,acs)
            !L' = R'*L*R (inverse rotation)
            lrot = MATMUL(TRANSPOSE(acs),l(tid,:,:))
            lrot = MATMUL(lrot,acs)
            !Vorticity around Z', rigth-hand rule
            w3 = 0.5*(lrot(2,1) - lrot(1,2))
            !e3 = 0.5*(lrot(2,1) + lrot(1,2))
            !vortnum = ABS(w3/e3)
            a0 = phi_spo
            IF(w3 < 0) a0 = -phi_spo !Opposite rotation when vorticity is positive

         END IF

         !First rotate backward the matrix fabric as inclusions are assumed to
         !be initially coaxial to the X-Y-Z axes
         CALL rotmatrixZXZ(phi1,theta,phi2,acs)
         CALL rotvoigt(Sav(:,:,i),TRANSPOSE(acs),Sav(:,:,i)) !Reverse rotation with transpose rot. matrix

         !Now rotate around Y axis if phi_spo different from 0
         IF(phi_spo /= 0) THEN
            CALL rotmatrixZYZ(0d0,a0,0d0,acs1)
            CALL rotvoigt(Sav(:,:,i),TRANSPOSE(acs1),Sav(:,:,i)) !Reverse rotation with transpose rot. matrix
         END IF

         !Assign constant tensors defined in inputspo.dat
         IF(spomod == 2) THEN

         !Find DEM elastic tensor according to porosity of the geodynamic model
         IF(meltspomod > 0) THEN
            IF(meltmark(i) > 0) THEN
               j = FLOOR(meltmark(i)/Vstp) + 1
               dm = (meltmark(i) - Vstp*(j-1))/Vstp
               Sav(:,:,i) = Cdem(j,:,:)*(1-dm) + Cdem(j+1,:,:)*dm 
               rho(i) = ro_back*(1.0-meltmark(i)) + ro_incl*meltmark(i)
            END IF
         !Find DEM elastic tensor with constant porosity = Vmax
         ELSE
            j = INT(Vmax/Vstp) + 1
            Sav(:,:,i) = Cdem(j,:,:)
            rho(i) = ro_back*(1.0-Vmax) + ro_incl*Vmax
         END IF

         !Use elastic tensor from DREX_S for the matrix
         ELSEIF(spomod == 3) THEN

            !Find DEM elastic tensor according to porosity of the geodynamic model
            IF(meltspomod > 0) THEN
               IF(meltmark(i) >= Vstp) THEN
               j = FLOOR(meltmark(i)/Vstp) + 1
               Vmaxm = j*Vstp !Set porosity
               CALL DEM(i,Vmaxm,Cdem)
               dm = (meltmark(i) - Vstp*(j-1))/Vstp
               Sav(:,:,i) = Cdem(j,:,:)*(1.0-dm) + Cdem(j+1,:,:)*dm 
               rho(i) = rho(i)*(1.0-meltmark(i)) + ro_incl*meltmark(i)
               END IF
            !Find DEM elastic tensor with constant porosity = Vmax
            ELSE
               CALL DEM(i,Vmax,Cdem)
               j = INT(Vmax/Vstp) + 1
               Sav(:,:,i) = Cdem(j,:,:)
               rho(i) = rho(i)*(1.0-Vmax) + ro_incl*Vmax
            END IF

         END IF


         IF(phi_spo /= 0) THEN
            CALL rotvoigt(Sav(:,:,i),acs1,Sav(:,:,i))
         END IF

         !Rotate tensor toward Sigma 1
         CALL rotvoigt(Sav(:,:,i),acs,Sav(:,:,i))

      END IF
   
   END DO
   !$omp end do
   !$omp end parallel


   IF(spomod > 1 .AND. meltspomod > 0) DEALLOCATE(meltmark)

   !!! Write infos in hdf5 format for SKS splitting modeling
   command = 'mkdir '//str2//'SPO/'
   CALL SYSTEM(command)
   filename = str2//'SPO/Cijkl'//dt_str4//'.h5'

   CALL H5open_f (error)

   ! Create a new file using default properties.

   CALL H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

   !Lagrangian particles (crystal aggregates)
   CALL H5Gcreate_f(file_id, "/Particles", group_id, error)
   dum_int(1)=marknum    
   CALL loadsave_integer(1,1,group_id,1,H5T_NATIVE_INTEGER,dum_int(1:1),'marknum',1)

   !Rocktype
   CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mx1,'mx1',1)
   CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mx2,'mx2',1)
   CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mx2,'mx3',1)
   if(yinyang == 2) CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_INTEGER,mYY,'mYY',1)
   CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,rho,'rho',1)

   ALLOCATE(X1save(21,marknum))
   X1save  = 0d0 

   DO m = 1, marknum    
      !Calculate T,P

      !Save elastic tensor: only 21 elastic moduli
      X1save(1,m)  = Sav(1,1,m)
      X1save(2,m)  = Sav(2,2,m)
      X1save(3,m)  = Sav(3,3,m)
      X1save(4,m)  = Sav(2,3,m)
      X1save(5,m)  = Sav(1,3,m)
      X1save(6,m)  = Sav(1,2,m)
      X1save(7,m)  = Sav(4,4,m)
      X1save(8,m)  = Sav(5,5,m)
      X1save(9,m)  = Sav(6,6,m)
      X1save(10,m) = Sav(1,4,m)
      X1save(11,m) = Sav(2,5,m)
      X1save(12,m) = Sav(3,6,m)
      X1save(13,m) = Sav(3,4,m)
      X1save(14,m) = Sav(1,5,m)
      X1save(15,m) = Sav(2,6,m)
      X1save(16,m) = Sav(2,4,m)
      X1save(17,m) = Sav(3,5,m)
      X1save(18,m) = Sav(1,6,m)
      X1save(19,m) = Sav(5,6,m)
      X1save(20,m) = Sav(4,6,m)
      X1save(21,m) = Sav(4,5,m)
   END DO

   !Sav(GPa)
   nx(1)=21
   nx(2)=marknum
   CALL loadsave_double(0,2,group_id,nx(1:2),H5T_NATIVE_DOUBLE,X1save,'Sav',1)

   DEALLOCATE(X1save,X1,X2,X3)

   CALL H5Gclose_f(group_id, error)
   !Terminate access to the file.

   CALL H5Fclose_f(file_id, error)

   CALL H5close_f(error)
   !Close FORTRAN interface.

   filename = str2//'lagrangian'//dt_str4//'.h5'

   IF(fseminmod == 0 .AND. fsemaxmod == 0 .AND. TIaxismod == 0 .AND. vpmaxmod == 0 .AND. dvsmaxmod == 0) RETURN

   IF(spomod==1) write(*,'(a)'),' LPO fabric reset for all aggregates, replaced with SPO fabric (STILWE) for sufficiently deformed aggregates'
   IF(spomod==2) write(*,'(a)'),' LPO fabric reset for all aggregates, replaced with SPO fabric (DEM) for sufficiently deformed aggregates'
   write(*,*)
   
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Save symmetry axis of transverse isotropy sub-tensor !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(TIaxismod .NE. 0) THEN
   
   ALLOCATE(phi_a(ln_fse_num,3),perc_a(ln_fse_num),phi_a0(3,ln_fse_num))
   
   phi_a = 0d0 ; perc_a = 0d0 ; phi_a0 = 0d0

   !$omp parallel & 
   !$omp shared(rocktype,ln_fse,Sav,phi_a,perc_a,phi_a0,idx) &
   !$omp private(m,i) &    
   !$omp firstprivate(ln_fse_num)
   !$omp do schedule(guided,8)
   DO m = 1,ln_fse_num
      i=idx(m)

      !!! Percentage of anisotropy and orientation of axis of hexagonal symmetry
      !No hexagonal anisotropy for lower transition zones
      IF(rocktype(i) .NE. 3 ) THEN
         CALL DECSYM(Sav(:,:,i),rocktype(i),perc_a(m),phi_a(m,:),i)
         !anis_x
         phi_a0(1,m) = perc_a(m)*phi_a(m,1)
         !anis_y
         phi_a0(2,m) = perc_a(m)*phi_a(m,2)
         !anis_z
         phi_a0(3,m) = perc_a(m)*phi_a(m,3)
      END IF

   END DO

   !$omp end do
   !$omp end parallel

END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Save directions of max Vp and max dVs, and Gitudes of anisotropy !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(vpmaxmod .NE. 0 .or. dvsmaxmod .NE. 0) THEN  
   
   ALLOCATE(vpmax(3,ln_fse_num),dvsmax(3,ln_fse_num))

   degstp=15 !!! step in degrees
   nxy=360/degstp
   nz=90/degstp + 1

   ALLOCATE(phix0(nxy),phiy0(nxy),phiz0(nz),phiz1(nz))
   ALLOCATE(phix(nz,nxy),phiy(nz,nxy),phiz(nz,nxy))

   DO i = 1 , nz
      phi = degstp*(i-1)*deg2rad
      phiz0(i) = cos(phi);
      phiz1(i) = sin(phi);
   END DO

   DO i = 1 , nxy
      phi = degstp*(i-1)*deg2rad
      phix0(i) = cos(phi);
      phiy0(i) = sin(phi);
   END DO

   DO i = 1 , nz
      DO j = 1 , nxy
      phix(i,j) = phix0(j) * phiz1(i)
      phiy(i,j) = phiy0(j) * phiz1(i)
      phiz(i,j) = phiz0(i)
      END DO
   END DO

   DEALLOCATE(phix0,phiy0,phiz0,phiz1)
!!! Save only if enough strain rate (otherwise anisotropy data may be wrong)
!!! ln_fse = 0.5 it is enough to overprint any preexisting fabric (Ribe, 1992, JGR; Becker, 2006, EPSL) 
   DO m = 1,ln_fse_num
      
      i=idx(m)
      
      !!! Find min/max Vp, dVs direction 
      IF(rocktype(i) .NE. 3 ) THEN

      CALL MAX_VP_DVS(i,vp_anis,vpvect,dvs_anis,dvsvect)
         
      IF(vpmaxmod .NE. 0) THEN
         vpmax(1,m) = vp_anis * vpvect(1) 
         vpmax(2,m) = vp_anis * vpvect(2) 
         vpmax(3,m) = vp_anis * vpvect(3)
      END IF

      IF(dvsmaxmod .NE. 0) THEN
         dvsmax(1,m) = dvs_anis * dvsvect(1) 
         dvsmax(2,m) = dvs_anis * dvsvect(2) 
         dvsmax(3,m) = dvs_anis * dvsvect(3)
      END IF
         
      END IF
         
   END DO

END IF

   !!! Write infos in hdf5 format
   !Initialize FORTRAN interface.

   CALL H5open_f (error)

   ! Create a new file using default properties.

   CALL H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

   ALLOCATE(mxyz(3,ln_fse_num))
  
   DO m = 1,ln_fse_num
      i=idx(m)
      IF(cartspher==1) THEN
         mxyz(1,m)=mx1(i)
         mxyz(2,m)=mx2(i)   
         mxyz(3,m)=mx3(i) 
      ELSE IF(cartspher==2) THEN
         IF(yinyang == 2 .AND. mYY(i) == 2) THEN
            CALL yin2yang(mx1(i),mx3(i),lr,cr)
            mxyz(1,m)=mx2(i)*sin(cr)*cos(lr)
            mxyz(2,m)=mx2(i)*sin(cr)*sin(lr)
            mxyz(3,m)=mx2(i)*cos(cr)  
         ELSE
            mxyz(1,m)=mx2(i)*sin(mx3(i))*cos(mx1(i))
            mxyz(2,m)=mx2(i)*sin(mx3(i))*sin(mx1(i))
            mxyz(3,m)=mx2(i)*cos(mx3(i))  
         END IF    
      END IF
   END DO

   dum_int(1)=3
   dum_int(2)=ln_fse_num
   CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,mxyz,'xyz',1)

   !Save rocktype
   IF(rocktypemod) THEN
      print *,'Save rocktype'
      print * 
      ALLOCATE(rt1(ln_fse_num))
      DO m = 1,ln_fse_num
         i=idx(m)
         rt1(m)  = rocktype(i) 
      END DO
      CALL loadsave_integer(0,1,file_id,ln_fse_num,H5T_NATIVE_INTEGER,rt1,'rocktype',1)
      DEALLOCATE(rt1)
   END IF
   DEALLOCATE(idx,mxyz)

   IF(fseminmod .NE. 0) THEN
      print *,'Save FSE short axis'
      print * 
      CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,phi_fse_min0,'fse_min',1)
      DEALLOCATE(phi_fse_min0)
   END IF
   IF(fsemaxmod .NE. 0) THEN
      print *,'Save FSE long axis'
      print * 
      CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,phi_fse_max0,'fse_max',1)
      DEALLOCATE(phi_fse_max0)
   END IF
   IF(TIaxismod .NE. 0) THEN 
      print *,'Save TI symmetry axis'
      print * 
      CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,phi_a0,'TIaxis',1)
      DEALLOCATE(phi_a,perc_a,phi_a0)
   END IF
   IF(vpmaxmod .NE. 0)  THEN
      print *,'Save Vp max direction'
      print * 
      CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,vpmax,'vpmax',1)
      DEALLOCATE(vpmax)
      DEALLOCATE(phix,phiy,phiz)
   END IF
   IF(dvsmaxmod .NE. 0) THEN
      print *,'Save dVs max direction'
      print * 
      CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,dvsmax,'dvsmax',1)
      DEALLOCATE(dvsmax)
      IF(vpmaxmod .EQ. 0) DEALLOCATE(phix,phiy,phiz)
   END IF

   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   !Close FORTRAN interface.
   CALL H5close_f(error)
  
   !Make xmf file for visualization in Paraview
   filenamexmf = str2//'XDMF.lagrangian'//dt_str4//'.xmf'

   OPEN(19,file=filenamexmf,status='replace')
   WRITE(19,'(a)') '<?xml version="1.0" ?>'
   WRITE(19,'(a)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') ' <Domain>'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '   <Grid Name="materialSwarm" >'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a,1f10.8,a)') '      <Time Value="',timesum,'" />'
   WRITE(19,*) ' '
   WRITE(19,'(a,i0,a)') '          <Topology Type="POLYVERTEX" NodesPerElement="',ln_fse_num,'"> </Topology>'
   WRITE(19,'(a)') '          <Geometry Type="XYZ">'
   WRITE(19,'(a,i0,a,a,a)') '             <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',ln_fse_num,' 3">lagrangian',dt_str4,'.h5:/xyz </DataItem>'
   WRITE(19,'(a)') '          </Geometry>'
   WRITE(19,'(a)') ' '
   
   IF(rocktypemod .NE. 0) THEN
      WRITE(19,'(a)') '          <Attribute Type="Scalar" Center="Node" Name="rocktype">' 
      WRITE(19,'(a,i0,a,a,a)') '             <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',ln_fse_num,' 1">lagrangian',dt_str4,'.h5:/rocktype</DataItem>'
      WRITE(19,'(a)') '          </Attribute>' 
   END IF

   IF(fseminmod .NE. 0) THEN
      WRITE(19,'(a)') '          <Attribute Type="Vector" Center="Node" Name="fsemin">' 
      WRITE(19,'(a,i0,a,a,a)') '             <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',ln_fse_num,' 3">lagrangian',dt_str4,'.h5:/fse_min</DataItem>'
      WRITE(19,'(a)') '          </Attribute>' 
   END IF
   
   IF(fsemaxmod .NE. 0) THEN
      WRITE(19,'(a)') '          <Attribute Type="Vector" Center="Node" Name="fsemax">' 
      WRITE(19,'(a,i0,a,a,a)') '             <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',ln_fse_num,' 3">lagrangian',dt_str4,'.h5:/fse_max</DataItem>'
      WRITE(19,'(a)') '          </Attribute>' 
   END IF
   
   IF(TIaxismod .NE. 0) THEN
      WRITE(19,'(a)') '          <Attribute Type="Vector" Center="Node" Name="TIaxis">' 
      WRITE(19,'(a,i0,a,a,a)') '             <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',ln_fse_num,' 3">lagrangian',dt_str4,'.h5:/TIaxis</DataItem>'
      WRITE(19,'(a)') '          </Attribute>' 
   END IF
   
   IF(vpmaxmod .NE. 0) THEN
      WRITE(19,'(a)') '          <Attribute Type="Vector" Center="Node" Name="vpmax">' 
      WRITE(19,'(a,i0,a,a,a)') '             <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',ln_fse_num,' 3">lagrangian',dt_str4,'.h5:/vpmax</DataItem>'
      WRITE(19,'(a)') '          </Attribute>' 
   END IF
   
   IF(dvsmaxmod .NE. 0) THEN
      WRITE(19,'(a)') '          <Attribute Type="Vector" Center="Node" Name="dvsmax">' 
      WRITE(19,'(a,i0,a,a,a)') '             <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',ln_fse_num,' 3">lagrangian',dt_str4,'.h5:/dvsmax</DataItem>'
      WRITE(19,'(a)') '          </Attribute>' 
   END IF
   
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '   </Grid>'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') ' </Domain>'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '</Xdmf>'

   CLOSE(19)
    
   IF(fse3Dmod > 0) THEN

   print *,'Save 3D FSE'
   print * 

   !$omp parallel &
   !$omp shared(rocktype,Fij,Rotm,amin,amed,amax) &
   !$omp private(m,cV,acsV) &    
   !$omp firstprivate(marknum)
   !$omp do schedule(guided,8)
   DO m = 1,marknum

      IF(rocktype(m) < 10) THEN
      IF(mx1(m) >= n1first .AND. mx1(m) <= n1last) THEN
      IF(mx2(m) >= n2first .AND. mx2(m) <= n2last) THEN
      IF(mx3(m) >= n3first .AND. mx3(m) <= n3last) THEN

      cV(1,1)= amax(m) ; cV(2,2) = amed(m) ; cV(3,3) = amin(m)
      
      acsV(:,1) = Rotm(:,3,m)
      acsV(:,2) = Rotm(:,2,m)
      acsV(:,3) = Rotm(:,1,m)

      !Given a matrix A, the eigenvalues and eigenvectors are found
      !such that V*acsV = acsV*cV, and since acsV' = inv(acsV), V = acsV * cV * acsV'
      !Calculate the left stretch tensor V
      Fij(:,:,m) = MATMUL(MATMUL(acsV , cV),TRANSPOSE(acsV))

      END IF;END IF;END IF;END IF

   END DO
!$omp end do
!$omp end parallel

   !!! Write infos in hdf5 format
   filename = str2//'3Dfse'//dt_str4//'.h5'

   !Initialize FORTRAN interface.

   CALL H5open_f (error)

   ! Create a new file using default properties.

   CALL H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

   ALLOCATE(mxyz(3,marknum),LS(9,marknum),rt1(marknum))

   i = 0
   DO m = 1,marknum

      IF(rocktype(m) < 10) THEN
      IF(mx1(m) >= n1first .AND. mx1(m) <= n1last) THEN
      IF(mx2(m) >= n2first .AND. mx2(m) <= n2last) THEN
      IF(mx3(m) >= n3first .AND. mx3(m) <= n3last) THEN

      i = i + 1
      IF(cartspher==1) THEN
         mxyz(1,i)=mx1(m)
         mxyz(2,i)=mx2(m)
         mxyz(3,i)=mx3(m) 
      ELSE IF(cartspher==2) THEN
         IF(yinyang == 2 .AND. mYY(m) == 2) THEN
            CALL yin2yang(mx1(m),mx3(m),lr,cr)
            mxyz(1,i)=mx2(m)*sin(cr)*cos(lr)
            mxyz(2,i)=mx2(m)*sin(cr)*sin(lr)
            mxyz(3,i)=mx2(m)*cos(cr)  
         ELSE
            mxyz(1,i)=mx2(m)*sin(mx3(m))*cos(mx1(m))
            mxyz(2,i)=mx2(m)*sin(mx3(m))*sin(mx1(m))
            mxyz(3,i)=mx2(m)*cos(mx3(m))      
         END IF
      END IF
      LS(1,i) = Fij(1,1,m) ; LS(2,i) = Fij(1,2,m) ; LS(3,i) = Fij(1,3,m) 
      LS(4,i) = Fij(2,1,m) ; LS(5,i) = Fij(2,2,m) ; LS(6,i) = Fij(2,3,m) 
      LS(7,i) = Fij(3,1,m) ; LS(8,i) = Fij(3,2,m) ; LS(9,i) = Fij(3,3,m) 
      rt1(i)  = rocktype(m) 

      END IF;END IF;END IF;END IF

   END DO

   !Save coordinatesh tensor
   dum_int(1)=3
   dum_int(2)=marknum
   CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,mxyz,'xyz',1)
   DEALLOCATE(mxyz)

   !Save left stretch tensor
   dum_int(1)=9
   dum_int(2)=marknum
   CALL loadsave_double(0,2,file_id,dum_int(1:2),H5T_NATIVE_DOUBLE,LS,'V',1)
   DEALLOCATE(LS)

   !Save rocktype
   CALL loadsave_integer(0,1,file_id,marknum,H5T_NATIVE_INTEGER,rt1,'rocktype',1)
   DEALLOCATE(rt1)
   
   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   !Close FORTRAN interface.
   CALL H5close_f(error)
  
   !Make xmf file for visualization in Paraview
   filenamexmf = str2//'XDMF.3Dfse'//dt_str4//'.xmf'

   OPEN(19,file=filenamexmf,status='replace')
   WRITE(19,'(a)') '<?xml version="1.0" ?>'
   WRITE(19,'(a)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') ' <Domain>'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '   <Grid Name="materialSwarm" >'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a,1f10.8,a)') '      <Time Value="',timesum,'" />'
   WRITE(19,*) ' '
   WRITE(19,'(a,i0,a)') '          <Topology Type="POLYVERTEX" NodesPerElement="',marknum,'"> </Topology>'
   WRITE(19,'(a)') '          <Geometry Type="XYZ">'
   WRITE(19,'(a,i0,a,a,a)') '             <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',marknum,' 3">3Dfse',dt_str4,'.h5:/xyz </DataItem>'
   WRITE(19,'(a)') '          </Geometry>'
   WRITE(19,'(a)') ' '
   

   WRITE(19,'(a)') '          <Attribute Type="Tensor" Center="Node" Name="Left stretch tensor">'
   WRITE(19,'(a,i0,a,a,a)') '             <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',marknum,' 9">3Dfse',dt_str4,'.h5:/V</DataItem>'
   WRITE(19,'(a)') '          </Attribute>'

   WRITE(19,'(a)') '          <Attribute Type="Scalar" Center="Node" Name="rocktype">' 
   WRITE(19,'(a,i0,a,a,a)') '             <DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="',marknum,' 1">3Dfse',dt_str4,'.h5:/rocktype</DataItem>'
   WRITE(19,'(a)') '          </Attribute>' 

   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '   </Grid>'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') ' </Domain>'
   WRITE(19,'(a)') ' '
   WRITE(19,'(a)') '</Xdmf>'

   CLOSE(19)
    
   END IF
   
   END SUBROUTINE lagrangian_output 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine yin2yang, find coordinates in other spherical grid          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE yin2yang(lr,cr,lr1,cr1)

   USE comvar

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Load/Save dataset, format double precision                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE readmeltcontent(cijkl_dir)

   USE comvar
   USE hdf5

   IMPLICIT NONE

   CHARACTER(500) :: filename
   CHARACTER(4) :: dt_str4
   CHARACTER (len=*) :: cijkl_dir
   INTEGER :: i,j,k,xnum,ynum,znum,m,m1,dum_int(3)
   DOUBLE PRECISION :: xstp,ystp,zstp,dx,dy,dz
   DOUBLE PRECISION ,DIMENSION(:), ALLOCATABLE :: dum_melt
   DOUBLE PRECISION ,DIMENSION(:,:), ALLOCATABLE :: Xsave
   DOUBLE PRECISION ,DIMENSION(:,:,:), ALLOCATABLE :: Melt
   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl,memtype !Handles
   INTEGER     ::   error  ! Error flag


   filename=trim(cijkl_dir)//trim(meltfilename)

   write(*,*)
   write(*,'(a)'),' READ MELT CONTENT INFOS FROM FILE:'
   write(*,*)
   write(*,'(a,a)'),' ',trim(filename)
   write(*,*)

   !Initialize FORTRAN interface.

   CALL H5open_f (error)

   ! Create a new file using default properties.

   CALL H5Fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

   !Load input parameters
   CALL loadsave_integer(0,1,file_id,3,H5T_NATIVE_INTEGER,dum_int,'Gridnum',0)
   xnum=dum_int(1)
   ynum=dum_int(2)
   znum=dum_int(3)
   ALLOCATE(Melt(xnum,ynum,znum),dum_melt(xnum*ynum*znum),meltmark(marknum))
   CALL loadsave_double(0,1,file_id,xnum*ynum*znum,H5T_NATIVE_DOUBLE,dum_melt,'Field',0)

   IF(dimensions == 2) THEN

   DO i=1, xnum ; DO j=1, ynum
      m = j + (i-1)*ynum
      Melt(i,j,1) = dum_melt(m)
   END DO ; END DO

   !Interpolate melt content to aggregates
   xstp = x10size/(xnum-1)
   ystp = x20size/(ynum-1)
   DO m=1, marknum
      i = FLOOR(mx1(m)/xstp) + 1
      j = FLOOR(mx2(m)/ystp) + 1
      IF(cartspher == 2) j = FLOOR((6371d3-mx2(m))/ystp) + 1
      IF(i<1) i=1 ; IF(i> xnum-1) i = xnum-1
      IF(j<1) j=1 ; IF(j> ynum-1) j = ynum-1
      dx = (mx1(m) - xstp*(i-1))/xstp
      dy = (mx2(m) - ystp*(j-1))/ystp
      IF(cartspher == 2) dy = ((6371d3-mx2(m)) - ystp*(j-1))/ystp
      IF(dx<0) dx=0 ; IF(dx> 1d0) dx=1d0
      IF(dy<0) dy=0 ; IF(dy> 1d0) dy=1d0
      meltmark(m) = Melt(i,j,1)*(1d0-dx)*(1d0-dy) + Melt(i+1,j,1)*dx*(1d0-dy) + &
                 Melt(i,j+1,1)*(1d0-dx)*dy + Melt(i+1,j+1,1)*dx*dy 
!if(meltmark(m)>0) then
!print *,m,mx1(m),mx2(m),xstp,ystp,i,j,dx,dy,meltmark(m)
!read(*,*)
!end if
   END DO
   
   ELSE

   DO i=1, xnum ; DO j=1, ynum ; DO k=1, znum
      m = j + (i-1)*ynum + (k-1)*ynum*xnum
      Melt(i,j,k) = dum_melt(m)
   END DO ; END DO ; END DO

   !Interpolate melt content to aggregates
   xstp = x10size/(xnum-1)
   ystp = x20size/(ynum-1)
   zstp = x30size/(znum-1)
   DO m=1, marknum
      i = FLOOR(mx1(m)/xstp) + 1
      j = FLOOR(mx2(m)/ystp) + 1
      k = FLOOR(mx3(m)/zstp) + 1
      IF(cartspher == 2) j = FLOOR((6371d3-mx2(m))/ystp) + 1
      IF(i<1) i=1 ; IF(i> xnum-1) i = xnum-1
      IF(j<1) j=1 ; IF(j> ynum-1) j = ynum-1
      IF(k<1) k=1 ; IF(k> znum-1) k = znum-1
      dx = (mx1(m) - xstp*(i-1))/xstp
      dy = (mx2(m) - ystp*(j-1))/ystp
      dz = (mx3(m) - zstp*(k-1))/zstp
      IF(cartspher == 2) dy = ((6371d3-mx2(m)) - ystp*(j-1))/ystp
      IF(dx<0) dx=0 ; IF(dx> 1d0) dx=1d0
      IF(dy<0) dy=0 ; IF(dy> 1d0) dy=1d0
      IF(dz<0) dz=0 ; IF(dz> 1d0) dz=1d0
      meltmark(m) = Melt(i,j,k)*(1d0-dx)*(1d0-dy)*(1d0-dz) + Melt(i+1,j,k)*dx*(1d0-dy)*(1d0-dz) + &
                 Melt(i,j+1,k)*(1d0-dx)*dy*(1d0-dz) + Melt(i+1,j+1,k)*dx*dy*(1d0-dz) + & 
                 Melt(i,j,k)*(1d0-dx)*(1d0-dy)*dz + Melt(i+1,j,k)*dx*(1d0-dy)*dz + &
                 Melt(i,j+1,k)*(1d0-dx)*dy*dz + Melt(i+1,j+1,k)*dx*dy*dz 
   END DO

   END IF

   DEALLOCATE(dum_melt,Melt)

   END SUBROUTINE readmeltcontent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE MAX_VP_DVS(m,vp_anis,vpvect,vs_anis,dvsvect)

   USE comvar
   USE omp_lib

   IMPLICIT NONE
  
   DOUBLE PRECISION c(3,3,3,3),n(3),mtens(3,3),evals(3),evects(3,3),vpvect(3),dvsvect(3)
   DOUBLE PRECISION vpmax,vpmin,vp_anis,dvsmax,vs_anis,vp,vsmax,vsmin,dvs
   INTEGER i,j,k,ll,a,b,p,nrot,m
    
   !Convert 6x6 to 3x3x3x3
   CALL TENS4(Sav(:,:,m),c)
   
   vpmax = -1d0 ; vpmin = 1d30 ; dvsmax = -1d0
   vpvect = 0d0 ; dvsvect = 0d0
 
   DO a = 1 , nz
   DO b = 1 , nxy
       
   !Direction of p wave incidence
   n(1)=phix(a,b);
   n(2)=phiy(a,b);
   n(3)=phiz(a,b);

   !Christoffel tensor
   mtens = 0d0
   DO i = 1 , 3
      DO j = 1 , 3
         DO k = 1 , 3
            DO ll = 1 , 3
               mtens(i,k) = mtens(i,k) + c(i,j,k,ll)*n(j)*n(ll)
            END DO
         END DO
      END DO
   END DO

   !Find 3 wave moduli
   CALL DSYEVQ3(mtens,evects,evals)

   !Vp has largest eigenvalue 
   vp = 0.d0 ; p = 1
   DO i = 1 , 3
      IF (evals(i) .GT. vp) THEN
         vp = evals(i)
         p = i
      END IF
   END DO
   evals(p) = 0d0

   IF(vp .GT. vpmax) THEN
      vpmax = vp
      vpvect = n !evects(:,p)
   END IF

   IF(vp .LT. vpmin) THEN
      vpmin = vp
   END IF

   !Vs1
   vsmax = 0.d0 ; p = 1
   DO i = 1 , 3
      IF (evals(i) .GT. vsmax) THEN
         vsmax = evals(i)
         p = i
      END IF
   END DO
   evals(p) = 0d0

   !Vs2
   vsmin = 0.d0
   DO i = 1 , 3
      IF (evals(i) .GT. vsmin) THEN
         vsmin = evals(i)
      END IF
   END DO
   
   vsmax = sqrt(vsmax/rho(m)*1e+9) 
   vsmin = sqrt(vsmin/rho(m)*1e+9)
   dvs = vsmax - vsmin

   IF(dvs .GT. dvsmax) THEN
      dvsmax = dvs
      dvsvect = n !evects(:,p)
      vs_anis = dvs/(vsmax + vsmin)*200
   END IF

   END DO ; END DO

   !Convert into velocity km/s
   vpmax = sqrt(vpmax/rho(m)*1e+9)
   vpmin = sqrt(vpmin/rho(m)*1e+9)
   vp_anis = (vpmax - vpmin)/(vpmax + vpmin)*200

   END SUBROUTINE MAX_VP_DVS
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

      subroutine DECSYM(Sav1,rt,PERC,TIAXIS,m)

      use DECMOD

      implicit none

      double precision, dimension(6,6) :: Sav1
      double precision :: PERC,DC5,PI
      double precision, dimension(3) :: TIAXIS
      double precision, dimension(3,3) :: SCC
      double precision, dimension(6,6) :: CE1
      double precision, dimension(3,3,3,3) :: EL1
      double precision, dimension(21) :: XEC
      double precision :: XN,ANIS
      INTEGER :: m,rt

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
