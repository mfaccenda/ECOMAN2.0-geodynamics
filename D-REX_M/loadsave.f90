 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !!
 !!    Copyright (c) 2018-2023, Università di Padova, Manuele Faccenda
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
!!! Loading hdf5 files                                                     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE load(t,input_dir)

   USE comvar
   use omp_lib
   USE hdf5   

   IMPLICIT NONE

   !M3E!!!!!!!!!!!!
   include 'mpif.h'
   !M3E!!!!!!!!!!!!

   INTEGER :: i1,i2,i3,i4,t,yy
!  INTEGER :: n1,n2,n3,nx(3) ! unused
   CHARACTER (4) :: dt_str4
   CHARACTER (500) :: filename,str
   CHARACTER (len=*) :: input_dir
   CHARACTER (len(trim(input_dir))) :: str1  
   DOUBLE PRECISION, DIMENSION(2,2) :: Adummy2,dV2
   DOUBLE PRECISION, DIMENSION(3,3) :: Adummy,dV
   DOUBLE PRECISION :: datadb(2),lr,rr,cr,x2y2,x2z2,xx,yyy,zz,Uxx(3)
!  DOUBLE PRECISION :: dum,sii,sigmin,sigmax ! unused
   DOUBLE PRECISION ,DIMENSION(:), ALLOCATABLE :: V1,V2,V3,Fd0,Tk0,Pa0
   INTEGER(HID_T)  :: file_id !Handles
!  INTEGER(HID_T)  :: group_id, dataspace_id, dataset_id, attr_id,dcpl,memtype ! unused
   INTEGER     ::   error  ! Error flag

   INTEGER, DIMENSION(1:1) :: dims1D ! to avoid warning (rank mismatch) during compilation

   !M3E!!!!!!!!!!!!!!!!!!!!!
   integer :: errMPI,rankMPI
   !M3E!!!!!!!!!!!!!!!!!!!!!

   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Get MPI info
   call MPI_COMM_RANK(MPI_COMM_WORLD,rankMPI,errMPI)
   rankMPI = rankMPI + 1
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (rankMPI .eq. 1 ) then
      OPEN(19,file='cycle.txt',status='replace')
      write(19,"(a)")  'VTP file number'
      write(19,'(i4)') t
      write(19,"(a)")  'Cycle'
      CLOSE(19)
   endif

   str1=trim(input_dir)
   write(dt_str4,'(i4.4)') t
   str='vtp'//dt_str4//'.h5'

   if (rankMPI .eq. 1 ) then
      write(*,"(a)") '--------------------------------------------------------'
      write(*,*)
      write(*,'(a,a)') ' LOAD FILE ',trim(str)
      write(*,*)
   endif

   filename = str1//str

   !Load input files in hdf5 format

   !Initialize FORTRAN interface.

   CALL H5open_f (error)

   ! Open a new file using default properties.

   CALL H5Fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

   !Load input parameters
   dims1D = 2
   CALL loadsave_double(1,1,file_id,dims1D,H5T_NATIVE_DOUBLE,datadb,'Time',0)
   dt0=datadb(1)
   timesum=datadb(2)

   IF(t == 1) timesum0 = 0
   IF(Tstep > 1) dt0=timesum - timesum0
  
   timesum0 = timesum
 
   !Load velocity field
   yy = yinyang
   ALLOCATE(V1(nodenum*yy),V2(nodenum*yy))
   dims1D = nodenum*yy
   CALL loadsave_double(0,1,file_id,dims1D,H5T_NATIVE_DOUBLE,V1,'Nodes/V1',0)!Vx/VPhi
   CALL loadsave_double(0,1,file_id,dims1D,H5T_NATIVE_DOUBLE,V2,'Nodes/V2',0)!Vy/Vradial
   IF(dimensions == 3) THEN
      ALLOCATE(V3(nodenum*yy))
      CALL loadsave_double(0,1,file_id,dims1D,H5T_NATIVE_DOUBLE,V3,'Nodes/V3',0)!Vz/Vcolat 
   END IF

   !Load temperature and pressure
   IF(ptmod > 0) THEN 
      ALLOCATE(Tk0(nodenum*yy),Pa0(nodenum*yy)) 
      CALL loadsave_double(0,1,file_id,dims1D,H5T_NATIVE_DOUBLE,Tk0,'Nodes/Tk',0)
      CALL loadsave_double(0,1,file_id,dims1D,H5T_NATIVE_DOUBLE,Pa0,'Nodes/P',0)
   END IF

   !Load fraction of deformation accommodated by dislocation creep
   IF(fractdislmod > 0) THEN 
      ALLOCATE(Fd0(nodenum*yy))
      CALL loadsave_double(0,1,file_id,dims1D,H5T_NATIVE_DOUBLE,Fd0,'Nodes/Fd',0)
   END IF

   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   !Close FORTRAN interface.
   CALL H5close_f(error)

   Ui = 0d0 ; Ux = 0d0

   DO yy = 1, yinyang

   !$omp parallel do & 
   !$omp schedule(dynamic) &
   !$omp shared(V1,V2,V3,Ui,Tk0,Tk,Pa0,Pa,Fd0,Fd) &
   !$omp private(i1,i2,i3,i4) &    
   !$omp firstprivate(nx1,nx2,nx3,nodenum,yy,dimensions,ptmod,fractdislmod)
   DO i1 = 1, nx1
      DO i2 = 1, nx2
         DO i3 = 1, nx3
         
         !Global index
         i4 = nodenum*(yy-1) + i2 + (i1-1)*nx2 + (i3-1)*nx1*nx2

         Ui(yy,1,i1,i2,i3)=V1(i4)
         Ui(yy,2,i1,i2,i3)=V2(i4)
         IF(dimensions == 3) Ui(yy,3,i1,i2,i3)=V3(i4)

         IF(ptmod > 0) THEN
            if(Tk0(i4) < 0.0) Tk0(i4) = 0.0
            Tk(yy,i1,i2,i3) = Tk0(i4)
            if(Pa0(i4) < 0.0) Pa0(i4) = 0.0
            Pa(yy,i1,i2,i3) = Pa0(i4)
         END IF

         IF(fractdislmod > 0) THEN
            Fd(yy,i1,i2,i3) = Fd0(i4)
            IF(Fd(yy,i1,i2,i3)<0.0) write(*,"(a,i4,i4,i4,f6.3)") 'Fraction of disl. creep < 0 at node ',i1,i2,i3,Fd(yy,i1,i2,i3)
            IF(Fd(yy,i1,i2,i3)>1.0) write(*,"(a,i4,i4,i4,f6.3)") 'Fraction of disl. creep > 1 at node ',i1,i2,i3,Fd(yy,i1,i2,i3)
         END IF

         END DO
      END DO
   END DO
   !$omp end parallel do

   END DO

   DEALLOCATE(V1,V2)
   IF(dimensions == 3) DEALLOCATE(V3)
   IF(ptmod > 0) DEALLOCATE(Tk0,Pa0)
   IF(fractdislmod > 0) DEALLOCATE(Fd0)

   IF(cartspher == 1) Ux = Ui
   
!!! 2D model
   IF(dimensions == 2) THEN

!!! Convert velocities (m/s) from polar to cartesian coordinate system
   IF(cartspher == 2) THEN

   DO i1 = 1, nx1   
      DO i2 = 1, nx2  

         Adummy2 = 0d0
         
         !Matrix for converting velocities (m/s) from polar to cartesian
         !coordinate system. 
         !|Vx| = | -sin(Phi)  cos(Phi)| |VPhi|
         !|Vy| = |  cos(Phi)  sin(Phi)| |VR  |
         Adummy2(1,1) =-sin(X1(i1))
         Adummy2(1,2) = cos(X1(i1))

         Adummy2(2,1) = cos(X1(i1))
         Adummy2(2,2) = sin(X1(i1))
                
         Ux(1,1:2,i1,i2,1) = MATMUL(Adummy2,Ui(1,1:2,i1,i2,1))

      END DO
   END DO

   END IF

!!! Velocity gradient tensor
   Dij = 0d0

   !$omp parallel do & 
   !$omp schedule(dynamic) &
   !$omp shared(V1,V2,Ui,Dij) &
   !$omp private(i1,i2,dV2,Adummy2) &    
   !$omp firstprivate(nx1,nx2,X1,X2,x1periodic,x2periodic,basicstag,cartspher,sbfmod)
   DO i1 = 1, nx1
      DO i2 = 1, nx2
             
            dV2 = 0d0

            !Staggered grid  
            IF(basicstag == 2) THEN

               IF(cartspher == 2) THEN; print *,'Staggered grid not available yet for spherical coordinates'; stop; END IF

               !dVi/dx
               IF (i1 .LT. nx1) THEN
                  !At (+xstp/2,+ystp/2)
                  IF (i2 < nx2) Dij(1,1,1,i1+1,i2+1,1) = (Ui(1,1,i1+1,i2,1)-Ui(1,1,i1,i2,1))/(X1(i1+1)-X1(i1))
                  !At (0,0)
                  IF (i1 > 1) Dij(1,2,1,i1,i2,1) = (Ui(1,2,i1,i2,1)-Ui(1,2,i1-1,i2,1))/(X1(i1+1)-X1(i1-1))*2d0
                  IF (x1periodic > 0 .AND. i1 == 1) THEN
                     Dij(1,2,1,i1,i2,1) = (Ui(1,2,i1,i2,1)-Ui(1,2,nx1-1,i2,1))/(X1(i1+1)-X1(i1)+X1(nx1)-X1(nx1-1))*2d0
                  END IF
               END IF

               !dVi/dy
               IF (i2 .LT. nx2) THEN
                  !At (0,0)
                  IF (i2 > 1) Dij(1,1,2,i1,i2,1) = (Ui(1,1,i1,i2,1)-Ui(1,1,i1,i2-1,1))/(X2(i2+1)-X2(i2-1))*2d0
                  !At (+xstp/2,+ystp/2)
                  IF (i1 < nx1) Dij(1,2,2,i1+1,i2+1,1) = (Ui(1,2,i1,i2+1,1)-Ui(1,2,i1,i2,1))/(X2(i2+1)-X2(i2))
               END IF

            !Non-staggered grid
            ELSE

                !dVx in polar coordinates
                !dVx/dPhi
                IF(i1 > 1 .AND. i1 < nx1) THEN
                   dV2(1,1) = ((Ux(1,1,i1+1,i2,1)-Ux(1,1,i1  ,i2,1))/(X1(i1+1)-X1(i1  )) +&
                   &           (Ux(1,1,i1  ,i2,1)-Ux(1,1,i1-1,i2,1))/(X1(i1  )-X1(i1-1)))/2
                ENDIF
                IF(x1periodic > 0 .AND. (i1 == 1 .OR. i1 == nx1)) THEN
                   dV2(1,1) = ((Ux(1,1,2,i2,1)-Ux(1,1,1,i2,1))/(X1(2)-X1(1)) +            &
                   &           (Ux(1,1,nx1,i2,1)-Ux(1,1,nx1-1,i2,1))/(X1(nx1)-X1(nx1-1)))/2
                ENDIF
                !dVx/dR    
                IF (i2 > 1 .AND. i2 < nx2) THEN
                   dV2(1,2) = ((Ux(1,1,i1,i2+1,1)-Ux(1,1,i1,i2  ,1))/(X2(i2+1)-X2(i2  )) +&
                   &           (Ux(1,1,i1,i2  ,1)-Ux(1,1,i1,i2-1,1))/(X2(i2  )-X2(i2-1)))/2
                ENDIF
                IF(x2periodic > 0 .AND. (i2 == 1 .OR. i2 == nx2)) THEN
                   dV2(1,2) = ((Ux(1,1,i1,2,1)-Ux(1,1,i1,1,1))/(X2(2)-X2(1)) +            &
                   &           (Ux(1,1,i1,nx2,1)-Ux(1,1,i1,nx2-1,1))/(X2(nx2)-X2(nx2-1)))/2
                ENDIF
                !dVy in polar coordinates
                !dVy/dPhi  
                IF(i1 > 1 .AND. i1 < nx1) THEN
                   dV2(2,1) = ((Ux(1,2,i1+1,i2,1)-Ux(1,2,i1  ,i2,1))/(X1(i1+1)-X1(i1  )) +&
                   &           (Ux(1,2,i1  ,i2,1)-Ux(1,2,i1-1,i2,1))/(X1(i1  )-X1(i1-1)))/2
                ENDIF
                IF(x1periodic > 0 .AND. (i1 == 1 .OR. i1 == nx1)) THEN
                   dV2(2,1) = ((Ux(1,2,2,i2,1)-Ux(1,2,1,i2,1))/(X1(2)-X1(1)) +            &
                   &           (Ux(1,2,nx1,i2,1)-Ux(1,2,nx1-1,i2,1))/(X1(nx1)-X1(nx1-1)))/2
                ENDIF
                !dVy/dR    
                IF (i2 > 1 .AND. i2 < nx2) THEN
                   dV2(2,2) = ((Ux(1,2,i1,i2+1,1)-Ux(1,2,i1,i2  ,1))/(X2(i2+1)-X2(i2  )) +&
                   &           (Ux(1,2,i1,i2  ,1)-Ux(1,2,i1,i2-1,1))/(X2(i2  )-X2(i2-1)))/2
                ENDIF
                IF(x2periodic > 0 .AND. (i2 == 1 .OR. i2 == nx2)) THEN
                   dV2(2,2) = ((Ux(1,2,i1,2,1)-Ux(1,2,i1,1,1))/(X2(2)-X2(1)) +            &
                   &           (Ux(1,2,i1,nx2,1)-Ux(1,2,i1,nx2-1,1))/(X2(nx2)-X2(nx2-1)))/2
                ENDIF
                !Incompressibility
                IF(sbfmod > 0) dV2(2,2) = -dV2(1,1)
                IF(cartspher == 1) Dij(1,1:2,1:2,i1,i2,1) = dV2
      
                !Now in cartesian coordinates
                IF(cartspher == 2) THEN

                   Adummy2 = 0d0

                   !Matrix for converting derivatives from polar to cartesian
                   !coordinate system. For example, for Vx:
                   !dVx/dx = dVx/dPhi*dPhi/dx + dVx/dR*dR/dx
                   !dVx/dy = dVx/dPhi*dPhi/dy + dVx/dR*dR/dy
                   !|dPhi/dx dPhi/dy | =  |-sin(Phi)/R  cos(Phi)/R|
                   !|  dR/dx   dR/dy | =  | cos(Phi)    sin(Phi)  |
                   Adummy2(1,1) =-sin(X1(i1))/X2(i2) 
                   Adummy2(1,2) = cos(X1(i1))/X2(i2) 

                   Adummy2(2,1) = cos(X1(i1))
                   Adummy2(2,2) = sin(X1(i1))

                   Dij(1,1:2,1:2,i1,i2,1) = MATMUL(dV2,Adummy2)
            
                END IF

            END IF

      END DO
   END DO
   !$omp end parallel do

!!! Convert Vlong from m/s to rad/s for advection
   IF(cartspher == 2) THEN

      DO i1 = 1, nx1   
         DO i2 = 1, nx2    
            !VPhi(rad/s) = VPhi(m/s)/R
            IF(Ui(1,1,i1,i2,1) .NE. 0) Ui(1,1,i1,i2,1)=Ui(1,1,i1,i2,1)/X2(i2)               
         END DO
      END DO

   END IF

   dt = 1d30
   DO i1 = 1, nx1
      DO i2 = 1, nx2
         if(Ui(1,1,i1,i2,1)/= 0 .AND. dt > ABS(x1stp/2/Ui(1,1,i1,i2,1))) dt = ABS(x1stp/2/Ui(1,1,i1,i2,1)) 
         if(Ui(1,2,i1,i2,1)/= 0 .AND. dt > ABS(x2stp/2/Ui(1,2,i1,i2,1))) dt = ABS(x2stp/2/Ui(1,2,i1,i2,1)) 
      END DO
   END DO

   END IF

!!! END 2D model

!!! 3D model
   IF(dimensions == 3) THEN

   IF(cartspher == 2) THEN

   DO yy = 1, yinyang

   !$omp parallel do & 
   !$omp schedule(dynamic) &
   !$omp shared(Ui,Ux) &
   !$omp private(i1,i2,i3,lr,cr,Adummy) &    
   !$omp firstprivate(nx1,nx2,nx3,X1,X3,yy)
   DO i1 = 1, nx1
      DO i2 = 1, nx2
         DO i3 = 1, nx3

         !Matrix for converting velocities (m/s) from spherical to cartesian
         !coordinate system. 
         !|Vx| |-sin(Phi) sin(Theta)cos(Phi) cos(Theta)cos(Phi)| | Vphi |
         !|Vy|=| cos(Phi) sin(Theta)sin(Phi) cos(Theta)sin(Phi)|*|  VR  |
         !|Vz| |        0         cos(Theta)        -sin(Theta)| |Vtheta|
         IF(yy == 1) THEN
            lr = X1(i1)
            cr = X3(i3)
         ELSE
            CALL yin2yang(X1(i1),X3(i3),lr,cr)
         END IF

         Adummy(1,1) =-sin(lr)    
         Adummy(2,1) = cos(lr)    
         Adummy(3,1) = 0d0

         Adummy(1,2) = sin(cr)*cos(lr)
         Adummy(2,2) = sin(cr)*sin(lr)
         Adummy(3,2) = cos(cr)

         Adummy(1,3) = cos(cr)*cos(lr)
         Adummy(2,3) = cos(cr)*sin(lr)
         Adummy(3,3) =-sin(cr)

         Ux(yy,:,i1,i2,i3) = MATMUL(Adummy,Ui(yy,:,i1,i2,i3))

         END DO
      END DO
   END DO
   !$omp end parallel do

   END DO

   END IF

!!! Velocity gradient tensor
   Dij = 0d0

   DO yy = 1, yinyang

   !$omp parallel do & 
   !$omp schedule(dynamic) &
   !$omp shared(V1,V2,V3,Ui,Ux,Dij) &
   !$omp private(i1,i2,i3,dV,Adummy,Uxx,lr,cr,rr,xx,yyy,zz,x2z2) &    
   !$omp firstprivate(nx1,nx2,nx3,X1,X2,X3,x1periodic,x2periodic,x3periodic,yy,x1stp,x3stp,yinyang,basicstag,cartspher,sbfmod)
   DO i1 = 1, nx1
      DO i2 = 1, nx2
         DO i3 = 1, nx3
             
            dV = 0     
       
            !Staggered grid  
            IF(basicstag == 2) THEN

            IF(cartspher == 2) THEN; print *,'Staggered grid not available yet for spherical coordinates'; stop; END IF

            !At (+xstp/2,+ystp/2,+zstp/2)
            IF (i1 < nx1 .AND. i2 < nx2 .AND. i3 < nx3) THEN 
               !dVx/dx
               Dij(yy,1,1,i1+1,i2+1,i3+1) = (Ui(yy,1,i1+1,i2,i3)-Ui(yy,1,i1,i2,i3))/(X1(i1+1)-X1(i1))
               !dVy/dy
               Dij(yy,2,2,i1+1,i2+1,i3+1) = (Ui(yy,2,i1,i2+1,i3)-Ui(yy,2,i1,i2,i3))/(X2(i2+1)-X2(i2))
               !dVz/dz
               Dij(yy,3,3,i1+1,i2+1,i3+1) = (Ui(yy,3,i1,i2,i3+1)-Ui(yy,3,i1,i2,i3))/(X3(i3+1)-X3(i3))
               !!! Incompressibility           
               !Dij(3,3,i1+1,i2+1,i3+1) = -Dij(1,1,i1+1,i2+1,i3+1)-Dij(2,2,i1+1,i2+1,i3+1)
            END IF

            !dVy/dx
            !At (0,0,+zstp/2)
            IF (x1periodic > 0 .and. yinyang == 1 .and. i1 == 1 .and. i3 < nx3) THEN
                Dij(yy,2,1,i1,i2,i3) = (Ui(yy,2,1,i2,i3)-Ui(yy,2,nx1-1,i2,i3))/(X1(i1+1)-X1(i1)+X1(nx1)-X1(nx1-1))*2d0
                Dij(yy,2,1,nx1,i2,i3) = Dij(yy,2,1,i1,i2,i3)
            ELSE 
                IF (i1 > 1 .AND. i1 < nx1 .AND. i3 < nx3) THEN
                   Dij(yy,2,1,i1,i2,i3) = (Ui(yy,2,i1,i2,i3)-Ui(yy,2,i1-1,i2,i3))/(X1(i1+1)-X1(i1-1))*2d0
                END IF
            END IF
            !dVz/dx
            !At (0,+ystp/2,0)
            IF (x1periodic > 0 .and. yinyang == 1 .and. i1 == 1 .and. i2 < nx2) THEN
                Dij(yy,3,1,i1,i2,i3) = (Ui(yy,3,1,i2,i3)-Ui(yy,3,nx1-1,i2,i3))/(X1(i1+1)-X1(i1)+X1(nx1)-X1(nx1-1))*2d0
                Dij(yy,3,1,nx1,i2,i3) = Dij(yy,3,1,i1,i2,i3)
            ELSE 
                IF (i1 > 1 .AND. i1 < nx1 .AND. i2 < nx2) THEN
                   Dij(yy,3,1,i1,i2,i3) = (Ui(yy,3,i1,i2,i3)-Ui(yy,3,i1-1,i2,i3))/(X1(i1+1)-X1(i1-1))*2d0
                END IF
            END IF

            !dVx/dy
            !At (0,0,+zstp/2)
            IF (x2periodic > 0 .and. i2 == 1 .and. i3 < nx3) THEN
                Dij(yy,1,2,i1,i2,i3) = (Ui(yy,1,i1,1,i3)-Ui(yy,1,i1,nx2-1,i3))/(X2(i2+1)-X2(i2)+X2(nx2)-X2(nx2-1))*2d0
                Dij(yy,1,2,i1,nx2,i3) = Dij(yy,1,2,i1,i2,i3)
            ELSE 
                IF (i2 > 1 .AND. i2 < nx2 .AND. i3 < nx3) THEN
                   Dij(yy,1,2,i1,i2,i3) = (Ui(yy,1,i1,i2,i3)-Ui(yy,1,i1,i2-1,i3))/(X2(i2+1)-X2(i2-1))*2d0
                END IF
            END IF
    
            !dVz/dy
            !At (+xstp/2,0,0)
            IF (x2periodic > 0 .and. i2 == 1 .and. i1 < nx1) THEN
                Dij(yy,3,2,i1,i2,i3) = (Ui(yy,3,i1,1,i3)-Ui(yy,3,i1,nx2-1,i3))/(X2(i2+1)-X2(i2)+X2(nx2)-X2(nx2-1))*2d0
                Dij(yy,3,2,i1,nx2,i3) = Dij(yy,3,2,i1,i2,i3)
            ELSE 
                IF (i2 > 1 .AND. i2 < nx2 .AND. i1 < nx1) THEN
                   Dij(yy,3,2,i1,i2,i3) = (Ui(yy,3,i1,i2,i3)-Ui(yy,3,i1,i2-1,i3))/(X2(i2+1)-X2(i2-1))*2d0
                END IF
            END IF

            !dVx/dz
            !At (0,+ystp/2,0)
            IF (x3periodic > 0 .and. yinyang == 1 .and. i3 == 1 .and. i2 < nx2) THEN
                Dij(yy,1,3,i1,i2,i3) = (Ui(yy,1,i1,i2,1)-Ui(yy,1,i1,i2,nx3-1))/(X3(i3+1)-X3(i3)+X3(nx3)-X3(nx3-1))*2d0
                Dij(yy,1,3,i1,i2,nx3) = Dij(yy,1,3,i1,i2,i3)
            ELSE 
                IF (i3 > 1 .AND. i3 < nx3 .AND. i2 < nx2) THEN
                   Dij(yy,1,3,i1,i2,i3) = (Ui(yy,1,i1,i2,i3)-Ui(yy,1,i1,i2,i3-1))/(X3(i3+1)-X3(i3-1))*2d0
                END IF
            END IF

            !dVy/dz
            !At (+xstp/2,0,0)
            IF (x3periodic > 0 .and. yinyang == 1 .and. i3 == 1 .and. i1 < nx1) THEN
                Dij(yy,2,3,i1,i2,i3) = (Ui(yy,2,i1,i2,1)-Ui(yy,2,i1,i2,nx3-1))/(X3(i3+1)-X3(i3)+X3(nx3)-X3(nx3-1))*2d0
                Dij(yy,2,3,i1,i2,nx3) = Dij(yy,2,3,i1,i2,i3)
            ELSE 
                IF (i3 > 1 .AND. i3 < nx3 .AND. i1 .LT. nx1) THEN
                   Dij(yy,2,3,i1,i2,i3) = (Ui(yy,2,i1,i2,i3)-Ui(yy,2,i1,i2,i3-1))/(X3(i3+1)-X3(i3-1))*2d0
                END IF
            END IF


            !Non-staggered grid
            ELSE

                !dVi/dPhi   
                IF(i1 > 1 .AND. i1 < nx1) THEN
                   !dVx/dPhi  
                   dV(1,1) = ((Ux(yy,1,i1+1,i2,i3)-Ux(yy,1,i1,i2,i3))/(X1(i1+1)-X1(i1)) + &
                              (Ux(yy,1,i1,i2,i3)-Ux(yy,1,i1-1,i2,i3))/(X1(i1)-X1(i1-1)))/2
                   !dVy/dPhi  
                   dV(2,1) = ((Ux(yy,2,i1+1,i2,i3)-Ux(yy,2,i1,i2,i3))/(X1(i1+1)-X1(i1)) + &
                              (Ux(yy,2,i1,i2,i3)-Ux(yy,2,i1-1,i2,i3))/(X1(i1)-X1(i1-1)))/2
                   !dVz/dPhi  
                   dV(3,1) = ((Ux(yy,3,i1+1,i2,i3)-Ux(yy,3,i1,i2,i3))/(X1(i1+1)-X1(i1)) + &
                              (Ux(yy,3,i1,i2,i3)-Ux(yy,3,i1-1,i2,i3))/(X1(i1)-X1(i1-1)))/2
                END IF
                !Interpolate from other grid
                IF(i1 == 1 .AND. yinyang == 2) THEN

                   CALL interpUx(X1(i1)-x1stp,i2,X3(i3),yy,Uxx)

                   !dVx/dPhi  
                   dV(1,1) = ((Ux(yy,1,i1+1,i2,i3)-Ux(yy,1,i1,i2,i3))/(X1(i1+1)-X1(i1)) + &
                              (Ux(yy,1,i1,i2,i3)-Uxx(1))/x1stp)/2
                   !dVy/dPhi  
                   dV(2,1) = ((Ux(yy,2,i1+1,i2,i3)-Ux(yy,2,i1,i2,i3))/(X1(i1+1)-X1(i1)) + &
                              (Ux(yy,2,i1,i2,i3)-Uxx(2))/x1stp)/2
                   !dVz/dPhi  
                   dV(3,1) = ((Ux(yy,3,i1+1,i2,i3)-Ux(yy,3,i1,i2,i3))/(X1(i1+1)-X1(i1)) + &
                              (Ux(yy,3,i1,i2,i3)-Uxx(3))/x1stp)/2
                END IF
                !Interpolate from other grid
                IF(i1 == nx1 .AND. yinyang == 2) THEN

                   CALL interpUx(X1(i1)+x1stp,i2,X3(i3),yy,Uxx)

                   !dVx/dPhi  
                   dV(1,1) = ((Uxx(1)-Ux(yy,1,i1,i2,i3))/x1stp + &
                              (Ux(yy,1,i1,i2,i3)-Ux(yy,1,i1-1,i2,i3))/(X1(i1)-X1(i1-1)))/2
                   !dVy/dPhi  
                   dV(2,1) = ((Uxx(2)-Ux(yy,2,i1,i2,i3))/x1stp + &
                              (Ux(yy,2,i1,i2,i3)-Ux(yy,2,i1-1,i2,i3))/(X1(i1)-X1(i1-1)))/2
                   !dVz/dPhi  
                   dV(3,1) = ((Uxx(3)-Ux(yy,3,i1,i2,i3))/x1stp + &
                              (Ux(yy,3,i1,i2,i3)-Ux(yy,3,i1-1,i2,i3))/(X1(i1)-X1(i1-1)))/2
                END IF
                IF(cartspher == 1 .AND. x1periodic > 0 .AND. (i1 == 1 .OR. i1 == nx1)) THEN
                   !dVx/dPhi  
                   dV(1,1) = ((Ux(yy,1,2,i2,i3)-Ux(yy,1,1,i2,i3))/(X1(2)-X1(1)) + &
                              (Ux(yy,1,nx1,i2,i3)-Ux(yy,1,nx1-1,i2,i3))/(X1(nx1)-X1(nx1-1)))/2
                   !dVy/dPhi  
                   dV(2,1) = ((Ux(yy,2,2,i2,i3)-Ux(yy,2,1,i2,i3))/(X1(2)-X1(1)) + &
                              (Ux(yy,2,nx1,i2,i3)-Ux(yy,2,nx1-1,i2,i3))/(X1(nx1)-X1(nx1-1)))/2
                   !dVz/dPhi  
                   dV(3,1) = ((Ux(yy,3,2,i2,i3)-Ux(yy,3,1,i2,i3))/(X1(2)-X1(1)) + &
                              (Ux(yy,3,nx1,i2,i3)-Ux(yy,3,nx1-1,i2,i3))/(X1(nx1)-X1(nx1-1)))/2
                END IF

                !dVi/dR   
                IF (i2 > 1 .AND. i2 < nx2) THEN
                   !dVx/dR    
                   dV(1,2) = ((Ux(yy,1,i1,i2+1,i3)-Ux(yy,1,i1,i2,i3))/(X2(i2+1)-X2(i2)) + &
                              (Ux(yy,1,i1,i2,i3)-Ux(yy,1,i1,i2-1,i3))/(X2(i2)-X2(i2-1)))/2
                   !dVy/dR    
                   dV(2,2) = ((Ux(yy,2,i1,i2+1,i3)-Ux(yy,2,i1,i2,i3))/(X2(i2+1)-X2(i2)) + &
                              (Ux(yy,2,i1,i2,i3)-Ux(yy,2,i1,i2-1,i3))/(X2(i2)-X2(i2-1)))/2
                   !dVz/dR    
                   dV(3,2) = ((Ux(yy,3,i1,i2+1,i3)-Ux(yy,3,i1,i2,i3))/(X2(i2+1)-X2(i2)) + &
                              (Ux(yy,3,i1,i2,i3)-Ux(yy,3,i1,i2-1,i3))/(X2(i2)-X2(i2-1)))/2
                END IF
                IF (x2periodic > 0 .AND. (i2 == 1 .OR. i2 == nx2)) THEN 
                   !dVx/dR    
                   dV(1,2) = ((Ux(yy,1,i1,2,i3)-Ux(yy,1,i1,1,i3))/(X2(2)-X2(1)) + &
                              (Ux(yy,1,i1,nx2,i3)-Ux(yy,1,i1,nx2-1,i3))/(X2(nx2)-X2(nx2-1)))/2
                   !dVy/dR    
                   dV(2,2) = ((Ux(yy,2,i1,2,i3)-Ux(yy,2,i1,1,i3))/(X2(2)-X2(1)) + &
                              (Ux(yy,2,i1,nx2,i3)-Ux(yy,2,i1,nx2-1,i3))/(X2(nx2)-X2(nx2-1)))/2
                   !dVz/dR    
                   dV(3,2) = ((Ux(yy,3,i1,2,i3)-Ux(yy,3,i1,1,i3))/(X2(2)-X2(1)) + &
                              (Ux(yy,3,i1,nx2,i3)-Ux(yy,3,i1,nx2-1,i3))/(X2(nx2)-X2(nx2-1)))/2
                END IF

                !dVi/dTheta
                IF (i3 > 1 .AND. i3 < nx3) THEN
                   !dVx/dTheta
                   dV(1,3) = ((Ux(yy,1,i1,i2,i3+1)-Ux(yy,1,i1,i2,i3))/(X3(i3+1)-X3(i3)) + &
                              (Ux(yy,1,i1,i2,i3)-Ux(yy,1,i1,i2,i3-1))/(X3(i3)-X3(i3-1)))/2
                   !dVy/dTheta
                   dV(2,3) = ((Ux(yy,2,i1,i2,i3+1)-Ux(yy,2,i1,i2,i3))/(X3(i3+1)-X3(i3)) + &
                              (Ux(yy,2,i1,i2,i3)-Ux(yy,2,i1,i2,i3-1))/(X3(i3)-X3(i3-1)))/2
                   !dVz/dTheta
                   dV(3,3) = ((Ux(yy,3,i1,i2,i3+1)-Ux(yy,3,i1,i2,i3))/(X3(i3+1)-X3(i3)) + &
                              (Ux(yy,3,i1,i2,i3)-Ux(yy,3,i1,i2,i3-1))/(X3(i3)-X3(i3-1)))/2
                END IF
                !Interpolate from other grid
                IF(i3 == 1 .AND. yinyang == 2) THEN

                   CALL interpUx(X1(i1),i2,X3(i3)-x3stp,yy,Uxx)

                   !dVx/dPhi  
                   dV(1,3) = ((Ux(yy,1,i1,i2,i3+1)-Ux(yy,1,i1,i2,i3))/(X3(i3+1)-X3(i3)) + &
                              (Ux(yy,1,i1,i2,i3)-Uxx(1))/x3stp)/2
                   !dVy/dPhi  
                   dV(2,3) = ((Ux(yy,2,i1,i2,i3+1)-Ux(yy,2,i1,i2,i3))/(X3(i3+1)-X3(i3)) + &
                              (Ux(yy,2,i1,i2,i3)-Uxx(2))/x3stp)/2
                   !dVz/dPhi  
                   dV(3,3) = ((Ux(yy,3,i1,i2,i3+1)-Ux(yy,3,i1,i2,i3))/(X3(i3+1)-X3(i3)) + &
                              (Ux(yy,3,i1,i2,i3)-Uxx(3))/x3stp)/2
                END IF
                !Interpolate from other grid
                IF(i3 == nx3 .AND. yinyang == 2) THEN

                   CALL interpUx(X1(i1),i2,X3(i3)+x3stp,yy,Uxx)

                   !dVx/dTheta
                   dV(1,3) = ((Uxx(1)-Ux(yy,1,i1,i2,i3))/x3stp + &
                              (Ux(yy,1,i1,i2,i3)-Ux(yy,1,i1,i2,i3-1))/(X3(i3)-X3(i3-1)))/2
                   !dVy/dTheta
                   dV(2,3) = ((Uxx(2)-Ux(yy,2,i1,i2,i3))/x3stp + &
                              (Ux(yy,2,i1,i2,i3)-Ux(yy,2,i1,i2,i3-1))/(X3(i3)-X3(i3-1)))/2
                   !dVz/dTheta
                   dV(3,3) = ((Uxx(3)-Ux(yy,3,i1,i2,i3))/x3stp + &
                              (Ux(yy,3,i1,i2,i3)-Ux(yy,3,i1,i2,i3-1))/(X3(i3)-X3(i3-1)))/2
                END IF
                IF(cartspher == 1 .AND. x3periodic > 0 .AND. (i3 == 1 .OR. i3 == nx3)) THEN
                   !dVx/dTheta
                   dV(1,3) = ((Ux(yy,1,i1,i2,2)-Ux(yy,1,i1,i2,1))/(X3(2)-X3(1)) + &
                              (Ux(yy,1,i1,i2,nx3)-Ux(yy,1,i1,i2,nx3-1))/(X3(nx3)-X3(nx3-1)))/2
                   !dVy/dTheta
                   dV(2,3) = ((Ux(yy,2,i1,i2,2)-Ux(yy,2,i1,i2,1))/(X3(2)-X3(1)) + &
                              (Ux(yy,2,i1,i2,nx3)-Ux(yy,2,i1,i2,nx3-1))/(X3(nx3)-X3(nx3-1)))/2
                   !dVz/dTheta
                   dV(3,3) = ((Ux(yy,3,i1,i2,2)-Ux(yy,3,i1,i2,1))/(X3(2)-X3(1)) + &
                              (Ux(yy,3,i1,i2,nx3)-Ux(yy,3,i1,i2,nx3-1))/(X3(nx3)-X3(nx3-1)))/2
                END IF

                IF(cartspher == 1) Dij(yy,:,:,i1,i2,i3) = dV
      
                !Now in cartesian coordinates
                IF(cartspher == 2) THEN

                   !Matrix for converting derivatives from spherical to cartesian
                   !coordinate system. For example, for Vx:
                   !dVx/dx = dVx/dPhi*dPhi/dx + dVx/dR*dR/dx + dVx/dTheta*dTheta/dx
                   !dVx/dy = dVx/dPhi*dPhi/dy + dVx/dR*dR/dy + dVx/dTheta*dTheta/dy
                   !dVx/dz = dVx/dPhi*dPhi/dz + dVx/dR*dR/dz + dVx/dTheta*dTheta/dz
                   !|  dPhi/dx   dPhi/dy   dPhi/dz| = |-sin(Phi)/sin(Theta)/R  cos(Phi)/sin(Theta)/R             0|
                   !|    dR/dx     dR/dy     dR/dz|   |   sin(Theta)cos(Phi)     sin(Theta)sin(Phi)    cos(Theta)|
                   !|dTheta/dx dTheta/dy dTheta/dz|   | cos(Theta)cos(Phi)/R   cos(Theta)sin(Phi)/R -sin(Theta)/R|
                   IF(yy == 1) THEN
                   
                   !Matrix in spherical coordinates
                   Adummy(1,1) =-sin(X1(i1))/sin(X3(i3))/X2(i2) 
                   Adummy(1,2) = cos(X1(i1))/sin(X3(i3))/X2(i2) 
                   Adummy(1,3) = 0d0                 

                   Adummy(2,1) = sin(X3(i3))*cos(X1(i1))
                   Adummy(2,2) = sin(X3(i3))*sin(X1(i1))
                   Adummy(2,3) = cos(X3(i3))

                   Adummy(3,1) = cos(X3(i3))*cos(X1(i1))/X2(i2) 
                   Adummy(3,2) = cos(X3(i3))*sin(X1(i1))/X2(i2) 
                   Adummy(3,3) = -sin(X3(i3))/X2(i2) 

                   !Matrix in cartesian coordinates
                   if(0==1) then
                   rr = X2(i2)
                   xx = rr*sin(X3(i3))*cos(X1(i1))
                   yyy= rr*sin(X3(i3))*sin(X1(i1))
                   zz = rr*cos(X3(i3))
                   x2y2 = (xx**2.0 + yyy**2.0)**0.5
                   !dPhi_n/di_e
                   Adummy(1,1) =-yyy/(x2y2**2.0)
                   Adummy(1,2) = xx/(x2y2**2.0)
                   Adummy(1,3) = 0d0

                   !dR_n/di_e
                   Adummy(2,1) = xx/rr
                   Adummy(2,2) = yyy/rr
                   Adummy(2,3) = zz/rr

                   !dTheta_n/di_e
                   Adummy(3,1) = xx*zz/((rr**2.0)*x2y2)
                   Adummy(3,2) = zz*yyy/((rr**2.0)*x2y2)
                   Adummy(3,3) =-x2y2/(rr**2.0)                    
                   end if

                   !Yang grid
                   ELSE
                
                   !Matrix for converting derivatives from spherical to cartesian
                   !coordinate system. For example, for Vx:
                   !dVx_e/dx_e = dVx_e/dPhi_n*dPhi_n/dx_e + dVx_e/dR_n*dR_n/dx_e + dVx_e/dTheta_n*dTheta_n/dx_e
                   !dVx_e/dy_e = dVx_e/dPhi_n*dPhi_n/dy_e + dVx_e/dR_n*dR_n/dy_e + dVx_e/dTheta_n*dTheta_n/dy_e
                   !dVx_e/dz_e = dVx_e/dPhi_n*dPhi_n/dz_e + dVx_e/dR_n*dR_n/dz_e + dVx_e/dTheta_n*dTheta_n/dz_e

                   !Matrix in cartesian coordinates
                   CALL yin2yang(X1(i1),X3(i3),lr,cr)
                   rr = X2(i2)
                   xx = rr*sin(cr)*cos(lr)
                   yyy= rr*sin(cr)*sin(lr)
                   zz = rr*cos(cr)
                   x2z2 = (xx**2.0 + zz**2.0)**0.5
                   !dPhi_n/di_e
                   Adummy(1,1) =-zz/(x2z2**2.0)
                   Adummy(1,2) = 0d0
                   Adummy(1,3) = xx/(x2z2**2.0)

                   !dR_n/di_e
                   !Adummy(2,1) = sin(cr)*cos(lr)
                   !Adummy(2,2) = sin(cr)*sin(lr)
                   !Adummy(2,3) = cos(cr)
                   Adummy(2,1) = xx/rr
                   Adummy(2,2) = yyy/rr
                   Adummy(2,3) = zz/rr

                   !dTheta_n/di_e
                   Adummy(3,1) = xx*yyy/((rr**2.0)*x2z2)
                   Adummy(3,2) =-x2z2/(rr**2.0)                    
                   Adummy(3,3) = zz*yyy/((rr**2.0)*x2z2)

                   dV(:,1) = -dV(:,1)

                   END IF

                   Dij(yy,:,:,i1,i2,i3) = MATMUL(dV,Adummy)
            
                   !Incompressibility
                   IF(sbfmod > 0) Dij(yy,3,3,i1,i2,i3) = -Dij(yy,1,1,i1,i2,i3)-Dij(yy,2,2,i1,i2,i3)

                END IF

            END IF

         END DO
      END DO
   END DO
   !$omp end parallel do

   END DO

!!! Convert Vcolat and Vphi from m/s to rad/s for advection
   IF(cartspher == 2) THEN

   DO yy = 1, yinyang

   !$omp parallel do & 
   !$omp schedule(dynamic) &
   !$omp shared(Ui) &
   !$omp private(i1,i2,i3,lr,cr) &    
   !$omp firstprivate(nx1,nx2,nx3,X1,X2,X3,yy,pi)
   DO i1 = 1, nx1   
      DO i2 = 1, nx2    
         DO i3 = 1, nx3    
            IF(yy == 1) THEN
               !Vcolat(rad/s) = Vcolat(m/s)/R
               IF(Ui(yy,3,i1,i2,i3) /= 0) Ui(yy,3,i1,i2,i3)=Ui(yy,3,i1,i2,i3)/X2(i2)               
               !Vphi(rad/s) = Vphi(m/s)/R/sin(theta)
               IF(Ui(yy,1,i1,i2,i3) /= 0 .AND. X3(i3) /= 0 .AND. X3(i3) /= pi) THEN
                  Ui(yy,1,i1,i2,i3)=Ui(yy,1,i1,i2,i3)/X2(i2)/sin(X3(i3))
               ENDIF
            ELSE
               CALL yin2yang(X1(i1),X3(i3),lr,cr)
               !Vcolat(rad/s) = Vcolat(m/s)/R
               IF(Ui(yy,3,i1,i2,i3) /= 0) Ui(yy,3,i1,i2,i3)=Ui(yy,3,i1,i2,i3)/X2(i2)               
               !Vphi(rad/s) = Vphi(m/s)/R/sin(theta)
               IF(Ui(yy,1,i1,i2,i3) /= 0) Ui(yy,1,i1,i2,i3)=Ui(yy,1,i1,i2,i3)/X2(i2)/sin(cr)        
            END IF
         END DO
      END DO
   END DO
   !$omp end parallel do

   END DO

   END IF

!!! Check advection timestep 
   dt = 1d30
   
   DO yy = 1, yinyang

   DO i1 = 1, nx1
      DO i2 = 1, nx2
         DO i3 = 1, nx3
            if(Ui(yy,1,i1,i2,i3)/= 0 .AND. dt > ABS(x1stp/2/Ui(yy,1,i1,i2,i3))) dt = ABS(x1stp/2/Ui(yy,1,i1,i2,i3)) 
            if(Ui(yy,2,i1,i2,i3)/= 0 .AND. dt > ABS(x2stp/2/Ui(yy,2,i1,i2,i3))) dt = ABS(x2stp/2/Ui(yy,2,i1,i2,i3)) 
            if(Ui(yy,3,i1,i2,i3)/= 0 .AND. dt > ABS(x3stp/2/Ui(yy,3,i1,i2,i3))) dt = ABS(x3stp/2/Ui(yy,3,i1,i2,i3)) 
         END DO
      END DO
   END DO

   END DO

   END IF

!!! END 3D model

   IF(dt >= dt0) THEN
      dt = dt0
      numdt = 1
   ELSE
      numdt = FLOOR(dt0/dt) + 1
      IF(dt*numdt == dt0) numdt = numdt - 1
   END IF

   IF(Tinit == Tend) THEN
      numdt = 1
      if ( rankMPI .eq. 1 ) then
         write(*,'(a,1es13.6)') ' TOTAL TIMESTEP     = ',dt0
         write(*,*)
         write(*,'(a,1es13.6,a,i13)') ' ADVECTION TIMESTEP = ',dt,', NUMBER OF CYCLES = ',NINT(timemax/dt)
         write(*,*)
      endif
   ELSE
      if ( rankMPI .eq. 1 ) then
         write(*,'(a,1es13.6)') ' TOTAL TIMESTEP     = ',dt0
         write(*,*)
         write(*,'(a,1es13.6,a,i13)') ' ADVECTION TIMESTEP = ',dt,', NUMBER OF CYCLES = ',numdt
         write(*,*)
      endif
   END IF

   END SUBROUTINE load

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
!!! subroutine yin2yang, find coordinates in other spherical grid          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE interpUx(lr,i2s,cr,yy,Uxx)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i1s,i2s,i3s,i20s,yy,yyy
   DOUBLE PRECISION :: lr,cr,lr1,cr1,Uxx(3)
   DOUBLE PRECISION :: y1,y2,y3,y4,y5,y6,y7,y8

   CALL yin2yang(lr,cr,lr1,cr1)
   CALL upperleft(lr1,X2(i2s),cr1,i1s,i20s,i3s)
   
   IF(yy == 1) yyy = 2
   IF(yy == 2) yyy = 1

   !Vx
   y1 = Ux(yyy,1,i1s,i2s,i3s) ;       y2 = Ux(yyy,1,i1s+1,i2s,i3s)
   y3 = Ux(yyy,1,i1s+1,i2s,i3s+1) ;   y4 = Ux(yyy,1,i1s,i2s,i3s+1)
   y5 = Ux(yyy,1,i1s,i2s+1,i3s) ;     y6 = Ux(yyy,1,i1s+1,i2s+1,i3s)
   y7 = Ux(yyy,1,i1s+1,i2s+1,i3s+1) ; y8 = Ux(yyy,1,i1s,i2s+1,i3s+1)

   CALL interp(lr1,X2(i2s),cr1,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,Uxx(1))

   !Vy
   y1 = Ux(yyy,2,i1s,i2s,i3s) ;       y2 = Ux(yyy,2,i1s+1,i2s,i3s)
   y3 = Ux(yyy,2,i1s+1,i2s,i3s+1) ;   y4 = Ux(yyy,2,i1s,i2s,i3s+1)
   y5 = Ux(yyy,2,i1s,i2s+1,i3s) ;     y6 = Ux(yyy,2,i1s+1,i2s+1,i3s)
   y7 = Ux(yyy,2,i1s+1,i2s+1,i3s+1) ; y8 = Ux(yyy,2,i1s,i2s+1,i3s+1)

   CALL interp(lr1,X2(i2s),cr1,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,Uxx(2))

   !Vz
   y1 = Ux(yyy,3,i1s,i2s,i3s) ;       y2 = Ux(yyy,3,i1s+1,i2s,i3s)
   y3 = Ux(yyy,3,i1s+1,i2s,i3s+1) ;   y4 = Ux(yyy,3,i1s,i2s,i3s+1)
   y5 = Ux(yyy,3,i1s,i2s+1,i3s) ;     y6 = Ux(yyy,3,i1s+1,i2s+1,i3s)
   y7 = Ux(yyy,3,i1s+1,i2s+1,i3s+1) ; y8 = Ux(yyy,3,i1s,i2s+1,i3s+1)

   CALL interp(lr1,X2(i2s),cr1,i1s,i2s,i3s,y1,y2,y3,y4,y5,y6,y7,y8,Uxx(3))

   END SUBROUTINE interpUx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine anisotropy, calculte hexagonal symmetry axis and Cijkl tensor!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine stifftenzsave_mpi(t,output_dir,ncyc2)

   use comvar
   use omp_lib
   use hdf5    

   !M3E!!!!!!!!!!!!!!
   use class_DistGrid
   !M3E!!!!!!!!!!!!!!

   implicit none

   ! Input variables
   integer,           intent(in) :: t
   character (len=*), intent(in) :: output_dir
   integer,           intent(in) :: ncyc2

   ! Local variables
   double precision                              :: wtime
   character (500)                               :: filename
   character (len(trim(output_dir)))             :: str1
   character (4)                                 :: dt_str4
   integer                                       :: m,j,m1
   double precision, dimension(:),   allocatable :: mtk,mpgpa
   double precision, dimension(:,:), allocatable :: Xsave
   integer                                       :: dum_int(3)
   double precision                              :: dum_db(3)
   integer                        :: rankMPI,sizeMPI,errMPI
   integer                        :: marknum00
   integer                        :: mrkdist,mrkdump1,mrkdump2
   integer                        :: del1,del2
   integer,           allocatable :: alldel(:,:),off1(:),off2(:)
   integer                        :: allmarknumsav
   integer                        :: errH5
   integer(hid_t)                 :: plistH5,fileH5,groupH5
   integer                        :: dims(3),offs(3)
   integer,           allocatable :: rkscr1(:),rkscr2(:)
   integer,           allocatable :: scrI1(:)
   real(kind=double), allocatable :: scrR11(:),scrR12(:),scrR13(:)
   real(kind=double), allocatable :: scrR2(:,:)
   real(kind=double), allocatable :: scrR3(:,:,:)

   ! Get MPI info
   call MPI_COMM_RANK(MPI_COMM_WORLD,rankMPI,errMPI)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeMPI,errMPI)
   rankMPI = rankMPI + 1

   ! Calculate Anisotropy 
   if (rankMPI .eq. 1) then
      write(*,*)
      write(*,"(a)") '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
      write(*,"(a)") '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
      write(*,*)
      write(*,'(a,i0,a)') ' CALCULATE ELASTIC TENSOR AND SAVE OUTPUT AFTER ',ncyc2,' CYCLES'
   endif

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime = MPI_Wtime() 

   IF(fsemod == 0 .OR. sbfmod > 0) THEN

      ALLOCATE(Sav(6,6,marknum))
      IF(ptmod > 0) THEN
         ALLOCATE(mtk(marknum),mpgpa(marknum))
         mtk = 0d0 ; mpgpa = 0d0
      ENDIF

      !$omp parallel do & 
      !$omp schedule(dynamic) &
      !$omp shared(mx1,mx2,mx3,mYY,rocktype,Xol,odf,odf_ens,acs,acs_ens,Fij,Tk,Pa,Fd,td_rho,rho) &
      !$omp shared(mtk,mpgpa,Sav,S0,dS0dp,dS0dp2,dS0dt,dS0dpt,dS0dp5,dS0dt5,ijkl) &
      !$omp private(m) &    
      !$omp firstprivate(marknum,ptmod,uppermantlemod,sbfmod)
      DO m = 1 , marknum

         IF(rocktype(m) > 0 .AND. rocktype(m) < 10) then

            !!! P-T conditions and rho
            IF(ptmod > 0) CALL rhopt(m,mtk(m),mpgpa(m))

            ! Cijkl tensor (using Voigt-Reuss average)
            CALL tensorscalc(m,mtk(m),mpgpa(m),Sav(:,:,m))

         ELSE

            do j = 1,6
               Sav(1,j,m) = 0.d0
               Sav(2,j,m) = 0.d0
               Sav(3,j,m) = 0.d0
               Sav(4,j,m) = 0.d0
               Sav(5,j,m) = 0.d0
               Sav(6,j,m) = 0.d0
            end do

         END IF

      END DO
      !$omp end parallel do

      ! Lagrangian particles (crystal aggregates)

      ALLOCATE(Xsave(21,marknum))
      Xsave  = 0d0

      DO m = 1, marknum
         ! Save elastic tensor: only 21 elastic moduli
         Xsave(1,m)  = Sav(1,1,m)
         Xsave(2,m)  = Sav(2,2,m)
         Xsave(3,m)  = Sav(3,3,m)
         Xsave(4,m)  = Sav(2,3,m)
         Xsave(5,m)  = Sav(1,3,m)
         Xsave(6,m)  = Sav(1,2,m)
         Xsave(7,m)  = Sav(4,4,m)
         Xsave(8,m)  = Sav(5,5,m)
         Xsave(9,m)  = Sav(6,6,m)
         Xsave(10,m) = Sav(1,4,m)
         Xsave(11,m) = Sav(2,5,m)
         Xsave(12,m) = Sav(3,6,m)
         Xsave(13,m) = Sav(3,4,m)
         Xsave(14,m) = Sav(1,5,m)
         Xsave(15,m) = Sav(2,6,m)
         Xsave(16,m) = Sav(2,4,m)
         Xsave(17,m) = Sav(3,5,m)
         Xsave(18,m) = Sav(1,6,m)
         Xsave(19,m) = Sav(5,6,m)
         Xsave(20,m) = Sav(4,6,m)
         Xsave(21,m) = Sav(4,5,m)
      END DO

   END IF

   ! Get number of marknum without yinyang
   marknum00 = marknum
   if ( yinyang .eq. 2 ) marknum00 = marknum / 2

   ! Get number of distributed markers
   mrkdist = ptnmrk(rankMPI+1) - ptnmrk(rankMPI)

   ! Rocktype communications
   allocate(rkscr1(mrkdist))
   call DistMrkValI(1,1,rocktype(1:marknum00),rkscr1)
   if ( yinyang .eq. 2 ) then
      allocate(rkscr2(mrkdist))
      call DistMrkValI(1,1,rocktype(marknum00+1:marknum),rkscr2)
   endif 

   ! Compute deleted markers (i.e., with rocktype >= 100)
   mrkdump1 = 0
   do m = 1,mrkdist
      if( rkscr1(m) .lt. 100 ) then
         mrkdump1 = mrkdump1 + 1
      end if
   end do
   del1 = mrkdist - mrkdump1

   mrkdump2 = 0
   del2 = 0
   if ( yinyang .eq. 2 ) then
      do m = 1,mrkdist
         if( rkscr2(m) .lt. 100 ) then
            mrkdump2 = mrkdump2 + 1
         end if
      end do
      del2 = mrkdist - mrkdump2
   endif

   ! Communicate deleted markers
   allocate(alldel(2,sizeMPI))
   dum_int(1) = del1
   dum_int(2) = del2
   call MPI_Allgather(dum_int,2,MPI_INTEGER,alldel,2,MPI_INTEGER,MPI_COMM_WORLD,errMPI)

   ! Compute printing offsets

   ! Allocate offset scratches
   allocate(off1(sizeMPI+1),off2(sizeMPI+2))

   ! Initialization
   do m = 1,sizeMPI+1
      off1(m) = ptnmrk(m)
      off2(m) = ptnmrk(m) + ptnmrk(sizeMPI+1) - 1
   end do

   ! Remove deleted markers
   del1 = 0
   do m = 2,sizeMPI+1
      del1 = del1 + alldel(1,m-1)
      off1(m) = off1(m) - del1
   end do
   del2 = 0
   do m = 2,sizeMPI+1
      del2 = del2 + alldel(2,m-1)
      off2(m) = off2(m) - del1 - del2 
   end do

   allmarknumsav = off1(sizeMPI+1) - 1
   if ( yinyang .eq. 2 ) allmarknumsav = off2(sizeMPI+1) - 1 

   if ( rankMPI .eq. 1 ) then
      write(*,*)
      write(*,'(a,i0)') ' NUMBER OF AGGREGATES:       ',allmarknumsav+del1+del2
      write(*,'(a,i0)') ' NUMBER OF AGGREGATES SAVED: ',allmarknumsav
   endif

   ! Delete communication scratch
   deallocate(alldel)

   ! Get filename
   str1=trim(output_dir)
   write(dt_str4,'(i4.4)') t  
   filename=str1//'Cijkl'//dt_str4//'.h5'

   ! Write infos in hdf5 format
   if (rankMPI .eq. 1 ) then
      write(*,*)
      write(*,"(a)",ADVANCE='no') ' SAVE DATA TO ',trim(filename)
      write(*,*)
   endif

   ! Initialize FORTRAN interface
   call h5open_f(errH5);

   ! Setup file access property list with parallel I/O access
   call h5pcreate_f(H5P_FILE_ACCESS_F,plistH5,errH5);
   call h5pset_fapl_mpio_f(plistH5,MPI_COMM_WORLD,MPI_INFO_NULL,errH5);

   ! Create the file collectively
   call h5fcreate_f(filename,H5F_ACC_TRUNC_F,fileH5,errH5,access_prp=plistH5);
   call h5pclose_f(plistH5,errH5);

   ! Run input parameters
   dims(1) = 3
   dum_int(1) = dimensions
   dum_int(2) = cartspher
   dum_int(3) = basicstag
   call SaveAttributeMPI_integer(fileH5,'CoordinateSystem',1,dims,dum_int)
   dum_int(1) = x1periodic
   dum_int(2) = x2periodic
   dum_int(3) = x3periodic
   call SaveAttributeMPI_integer(fileH5,'Periodic',1,dims,dum_int)
   dims(1) = 2
   dum_db(1) = dt
   dum_db(2) = timesum

   call SaveAttributeMPI_double(fileH5,'Time',1,dims,dum_db)

   ! Eulerian Nodes
   call H5gcreate_f(fileH5,"/Nodes",groupH5,errH5)

   ! Dimensions
   dims(1) = 3
   dum_int(1) = nx1
   dum_int(2) = nx2
   dum_int(3) = nx3
   call SaveAttributeMPI_integer(groupH5,'nxyz',1,dims,dum_int)

   ! Coordinates
   dims(1) = nx1
   offs(1) = 0
   call SaveDatasetMPI_double(.true.,groupH5,'gx',1,dims,nx1,offs,X1)
   dims(1) = nx2
   call SaveDatasetMPI_double(.true.,groupH5,'gy',1,dims,nx2,offs,X2)
   dims(1) = nx3
   call SaveDatasetMPI_double(.true.,groupH5,'gz',1,dims,nx3,offs,X3)

   ! Close eulerian group
   call H5gclose_f(groupH5,errH5)

   ! Lagrangian nodes
   call H5Gcreate_f(fileH5,"/Particles",groupH5,errH5)

   ! Dimensions
   dims(1) = 1
   offs(1) = 0
   dum_int(1) = allmarknumsav 
   call SaveAttributeMPI_integer(groupH5,'marknum',1,dims,dum_int)
   dims(1) = 3
   dum_db(1) = CompDom%mx1stp
   dum_db(2) = CompDom%mx2stp
   dum_db(3) = CompDom%mx3stp
   call SaveDatasetMPI_double(.true.,groupH5,'mxstp',1,dims,3,offs,dum_db)
   dum_db(1) = CompDom%mx1min
   dum_db(2) = CompDom%mx2min
   dum_db(3) = CompDom%mx3min
   call SaveDatasetMPI_double(.true.,groupH5,'mxmin',1,dims,3,offs,dum_db)
   dum_db(1) = CompDom%mx1max
   dum_db(2) = CompDom%mx2max
   dum_db(3) = CompDom%mx3max
   call SaveDatasetMPI_double(.true.,groupH5,'mxmax',1,dims,3,offs,dum_db)

   ! Rocktype
   allocate(scrI1(mrkdist))
   m1 = 0
   do m = 1,mrkdist
      if(rkscr1(m) < 100) then
         m1 = m1 + 1
         scrI1(m1) = rkscr1(m) 
      endif 
   end do
   dims(1) = off1(rankMPI+1) - off1(rankMPI)
   offs(1) = off1(rankMPI) - 1
   call SaveDatasetMPI_integer(.true.,groupH5,'rocktype',1,dims,allmarknumsav,offs,scrI1)

   if(yinyang == 2) then

      m1 = 0
      do m = 1,mrkdist
         if(rkscr2(m) < 100) then
            m1 = m1 + 1
            scrI1(m1) = rkscr2(m) 
         endif 
      end do
      dims(1) = off2(rankMPI+1) - off2(rankMPI)
      offs(1) = off2(rankMPI) - 1
      call SaveDatasetMPI_integer(.false.,groupH5,'rocktype',1,dims,allmarknumsav,offs,scrI1)

   endif

   deallocate(scrI1)

   ! yinyang
   if(yinyang == 2) then

      allocate(scrI1(mrkdist))

      call DistMrkValI(1,1,mYY(1:marknum00),scrI1)
      m1 = 0
      do m = 1,mrkdist
         if(rkscr1(m) < 100) then
            m1 = m1 + 1
            scrI1(m1) = scrI1(m)
         endif 
      end do
      dims(1) = off1(rankMPI+1) - off1(rankMPI)
      offs(1) = off1(rankMPI) - 1
      call SaveDatasetMPI_integer(.true.,groupH5,'mYY',1,dims,allmarknumsav,offs,scrI1)

      call DistMrkValI(1,1,mYY(marknum00+1:marknum),scrI1)
      m1 = 0
      do m = 1,mrkdist
         if(rkscr2(m) < 100) then
            m1 = m1 + 1
            scrI1(m1) = scrI1(m)
         endif 
      end do
      dims(1) = off2(rankMPI+1) - off2(rankMPI)
      offs(1) = off2(rankMPI) - 1
      call SaveDatasetMPI_integer(.false.,groupH5,'mYY',1,dims,allmarknumsav,offs,scrI1)

      deallocate(scrI1)

   endif

   ! Coordinates
   allocate(scrR11(mrkdist),scrR12(mrkdist),scrR13(mrkdist))
   call DistMrkValR(1,1,mx1(1:marknum00),scrR11)
   call DistMrkValR(1,1,mx2(1:marknum00),scrR12)
   call DistMrkValR(1,1,mx3(1:marknum00),scrR13)
   m1 = 0
   do m = 1,mrkdist
      if(rkscr1(m) < 100) then
         m1 = m1 + 1
         scrR11(m1) = scrR11(m)
         scrR12(m1) = scrR12(m)
         scrR13(m1) = scrR13(m)
      endif 
   end do
   dims(1) = off1(rankMPI+1) - off1(rankMPI)
   offs(1) = off1(rankMPI) - 1
   call SaveDatasetMPI_double(.true.,groupH5,'mx1',1,dims,allmarknumsav,offs,scrR11)
   call SaveDatasetMPI_double(.true.,groupH5,'mx2',1,dims,allmarknumsav,offs,scrR12)
   call SaveDatasetMPI_double(.true.,groupH5,'mx3',1,dims,allmarknumsav,offs,scrR13)
 
   if(yinyang == 2) then

      call DistMrkValR(1,1,mx1(marknum00+1:marknum),scrR11)
      call DistMrkValR(1,1,mx2(marknum00+1:marknum),scrR12)
      call DistMrkValR(1,1,mx3(marknum00+1:marknum),scrR13)
      m1 = 0
      do m = 1,mrkdist
         if(rkscr2(m) < 100) then
            m1 = m1 + 1
            scrR11(m1) = scrR11(m)
            scrR12(m1) = scrR12(m)
            scrR13(m1) = scrR13(m)
         endif 
      end do
      dims(1) = off2(rankMPI+1) - off2(rankMPI)
      offs(1) = off2(rankMPI) - 1
      call SaveDatasetMPI_double(.false.,groupH5,'mx1',1,dims,allmarknumsav,offs,scrR11)
      call SaveDatasetMPI_double(.false.,groupH5,'mx2',1,dims,allmarknumsav,offs,scrR12)
      call SaveDatasetMPI_double(.false.,groupH5,'mx3',1,dims,allmarknumsav,offs,scrR13)

   endif

   deallocate(scrR11,scrR12,scrR13)

   ! FSE
   allocate(scrR3(3,3,mrkdist))
   call DistMrkValR(3,3,Fij(1:3,1:3,1:marknum00),scrR3)
   m1 = 0
   do m = 1,mrkdist
      if(rkscr1(m) < 100) then
         m1 = m1 + 1
         do j = 1,3
            scrR3(1,j,m1) = scrR3(1,j,m)
            scrR3(2,j,m1) = scrR3(2,j,m)
            scrR3(3,j,m1) = scrR3(3,j,m)
         end do
      endif 
   end do
   dims(1) = 3 
   dims(2) = 3 
   dims(3) = off1(rankMPI+1) - off1(rankMPI)
   offs(1) = 0 
   offs(2) = 0 
   offs(3) = off1(rankMPI) - 1
   call SaveDatasetMPI_double(.true.,groupH5,'fse',3,dims,allmarknumsav,offs,scrR3)
 
   if(yinyang == 2) then

      call DistMrkValR(3,3,Fij(1:3,1:3,marknum00+1:marknum),scrR3)
      m1 = 0
      do m = 1,mrkdist
         if(rkscr2(m) < 100) then
            m1 = m1 + 1
            do j = 1,3
               scrR3(1,j,m1) = scrR3(1,j,m)
               scrR3(2,j,m1) = scrR3(2,j,m)
               scrR3(3,j,m1) = scrR3(3,j,m)
            end do
         endif 
      end do
      dims(1) = 3
      dims(2) = 3
      dims(3) = off2(rankMPI+1) - off2(rankMPI)
      offs(1) = 0
      offs(2) = 0
      offs(3) = off2(rankMPI) - 1
      call SaveDatasetMPI_double(.false.,groupH5,'fse',3,dims,allmarknumsav,offs,scrR3)

   endif

   deallocate(scrR3)

   ! Rho / T / P / Sav
   IF(fsemod == 0 .OR. sbfmod > 0) THEN

      ! Rho (kg/m^3)
      allocate(scrR11(mrkdist))
      call DistMrkValR(1,1,rho(1:marknum00),scrR11)
      m1 = 0
      do m = 1,mrkdist
         if(rkscr1(m) < 100) then
            m1 = m1 + 1
            scrR11(m1) = scrR11(m)
         endif 
      end do
      dims(1) = off1(rankMPI+1) - off1(rankMPI)
      offs(1) = off1(rankMPI) - 1
      call SaveDatasetMPI_double(.true.,groupH5,'rho',1,dims,allmarknumsav,offs,scrR11)
      if(yinyang == 2) then
         call DistMrkValR(1,1,rho(marknum00+1:marknum),scrR11)
         m1 = 0
         do m = 1,mrkdist
            if(rkscr2(m) < 100) then
               m1 = m1 + 1
               scrR11(m1) = scrR11(m)
            endif 
         end do
         dims(1) = off2(rankMPI+1) - off2(rankMPI)
         offs(1) = off2(rankMPI) - 1
         call SaveDatasetMPI_double(.false.,groupH5,'rho',1,dims,allmarknumsav,offs,scrR11)
      endif
      deallocate(scrR11)

      if(ptmod > 0) then

         ! T(K) / P(GPa)
         allocate(scrR11(mrkdist),scrR12(mrkdist))

         call DistMrkValR(1,1,mtk(1:marknum00),  scrR11)
         call DistMrkValR(1,1,mpgpa(1:marknum00),scrR12)
         m1 = 0
         do m = 1,mrkdist
            if(rkscr1(m) < 100) then
               m1 = m1 + 1
               scrR11(m1) = scrR11(m)
               scrR12(m1) = scrR12(m)
            endif 
         end do
         dims(1) = off1(rankMPI+1) - off1(rankMPI)
         offs(1) = off1(rankMPI) - 1
         call SaveDatasetMPI_double(.true.,groupH5,'mtk',  1,dims,allmarknumsav,offs,scrR11)
         call SaveDatasetMPI_double(.true.,groupH5,'mpgpa',1,dims,allmarknumsav,offs,scrR12)
  
         if(yinyang == 2) then

            call DistMrkValR(1,1,mtk(marknum00+1:marknum),  scrR11)
            call DistMrkValR(1,1,mpgpa(marknum00+1:marknum),scrR12)
            m1 = 0
            do m = 1,mrkdist
               if(rkscr2(m) < 100) then
                  m1 = m1 + 1
                  scrR11(m1) = scrR11(m)
                  scrR12(m1) = scrR12(m)
               endif 
            end do
            dims(1) = off2(rankMPI+1) - off2(rankMPI)
            offs(1) = off2(rankMPI) - 1
            call SaveDatasetMPI_double(.false.,groupH5,'mtk',  1,dims,allmarknumsav,offs,scrR11)
            call SaveDatasetMPI_double(.false.,groupH5,'mpgpa',1,dims,allmarknumsav,offs,scrR12)

         endif
      
         deallocate(scrR11,scrR12)

      endif

      ! Sav(GPa)
      allocate(scrR2(21,mrkdist))

      call DistMrkValR(1,21,Xsave(1:21,1:marknum00),scrR2)
      m1 = 0
      do m = 1,mrkdist
         if(rkscr1(m) < 100) then
            m1 = m1 + 1
            do j = 1,21
               scrR2(j,m1) = scrR2(j,m)
            end do
         endif 
      end do
      dims(1) = 21 
      dims(2) = off1(rankMPI+1) - off1(rankMPI)
      offs(1) = 0
      offs(2) = off1(rankMPI) - 1
      call SaveDatasetMPI_double(.true.,groupH5,'Sav',2,dims,allmarknumsav,offs,scrR2)
  
      if(yinyang == 2) then

         call DistMrkValR(1,21,Xsave(1:21,marknum00+1:marknum),scrR2)
         m1 = 0
         do m = 1,mrkdist
            if(rkscr2(m) < 100) then
               m1 = m1 + 1
               do j = 1,21
                  scrR2(j,m1) = scrR2(j,m)
               end do
            endif 
         end do
         dims(1) = 21 
         dims(2) = off2(rankMPI+1) - off2(rankMPI)
         offs(1) = 0 
         offs(2) = off2(rankMPI) - 1
         call SaveDatasetMPI_double(.false.,groupH5,'Sav',2,dims,allmarknumsav,offs,scrR2)

      endif

      deallocate(scrR2)

      IF(ptmod > 0) DEALLOCATE(mtk,mpgpa)
      DEALLOCATE(Sav,Xsave)

   end if

   ! Close lagrangian group
   call h5gclose_f(groupH5,errH5);

   ! Close file
   call h5fclose_f(fileH5,errH5);

   ! Close FORTRAN interface
   call h5close_f(errH5);

   ! Remove communications scratches (rocktype)
   deallocate(rkscr1)
   if ( yinyang .eq. 2 ) deallocate(rkscr2)

   ! Remove offset scratches
   deallocate(off1,off2)

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime = MPI_Wtime() - wtime

   if (rankMPI .eq. 1) then
      write(*,*)
      write(*,'(a,1f10.2,a)') ' CALCULATE ELASTIC TENSOR AND SAVE OUTPUT OK! (',wtime,' sec)'
      write(*,*)
      write(*,"(a)") '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
      write(*,"(a)") '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
      write(*,*)
   endif

   end subroutine stifftenzsave_mpi

!M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
   SUBROUTINE stifftenzsave(t,output_dir,ncyc2)

   USE comvar
   USE omp_lib
   USE hdf5    

   IMPLICIT NONE

   INTEGER :: m,m1,t,nx(3),ncyc2,dum_int(3)
!  INTEGER :: i,j ! unused
   DOUBLE PRECISION :: wtime,dum_db(5)
!  DOUBLE PRECISION :: mtkm,mpgpam ! unused
   DOUBLE PRECISION ,DIMENSION(:), ALLOCATABLE :: mtk,mpgpa
!  DOUBLE PRECISION ,DIMENSION(:), ALLOCATABLE :: rt1 ! unused
   DOUBLE PRECISION ,DIMENSION(:,:), ALLOCATABLE :: Xsave
   CHARACTER (500) :: filename
   CHARACTER (len=*) :: output_dir
   CHARACTER (len(trim(output_dir))) :: str1
   CHARACTER (4) :: dt_str4
   INTEGER(HID_T)  :: file_id, group_id ! Handles
!  INTEGER(HID_T)  :: dataspace_id, dataset_id, attr_id, dcpl,memtype ! unused
   INTEGER     ::   error  ! Error flag

   ! to avoid warning (rank mismatch) during compilation
   INTEGER, DIMENSION(1:1) :: dims1D

   str1=trim(output_dir)

   write(dt_str4,'(i4.4)') t  
   filename=str1//'Cijkl'//dt_str4//'.h5'

!!! Calculate Anisotropy 
   
   write(*,*)
   write(*,"(a)") '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
   write(*,"(a)") '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
   write(*,*)
   write(*,'(a,i0,a)') ' CALCULATE ELASTIC TENSOR AND SAVE OUTPUT AFTER ',ncyc2,' CYCLES'
   
   wtime = OMP_GET_WTIME()

   IF(fsemod == 0 .OR. sbfmod > 0) THEN

      ALLOCATE(Sav(6,6,marknum))
      IF(ptmod > 0) THEN
         ALLOCATE(mtk(marknum),mpgpa(marknum))
         mtk = 0d0 ; mpgpa = 0d0
      ENDIF

      !$omp parallel do & 
      !$omp schedule(guided,8) &
      !$omp shared(mx1,mx2,mx3,mYY,rocktype,Xol,odf,odf_ens,acs,acs_ens,Tk,Pa,Fd,td_rho,rho) &
      !$omp shared(mtk,mpgpa,Sav,S0,dS0dp,dS0dp2,dS0dt,dS0dpt,dS0dp5,dS0dt5,ijkl) &
      !$omp private(m) &    
      !$omp firstprivate(size,marknum,ptmod,uppermantlemod)
      DO m = 1 , marknum
 
         IF(rocktype(m) > 0 .AND. rocktype(m) < 10) then

            !!! P-T conditions and rho
            IF(ptmod > 0) CALL rhopt(m,mtk(m),mpgpa(m))

            !!! Cijkl tensor (using Voigt-Reuss average)
            CALL tensorscalc(m,mtk(m),mpgpa(m),Sav(:,:,m))

         END IF

      END DO
      !$omp end parallel do

       ALLOCATE(Xsave(21,marknum))
       Xsave  = 0d0

       DO m = 1, marknum
          !Save elastic tensor: only 21 elastic moduli
          Xsave(1,m)  = Sav(1,1,m)
          Xsave(2,m)  = Sav(2,2,m)
          Xsave(3,m)  = Sav(3,3,m)
          Xsave(4,m)  = Sav(2,3,m)
          Xsave(5,m)  = Sav(1,3,m)
          Xsave(6,m)  = Sav(1,2,m)
          Xsave(7,m)  = Sav(4,4,m)
          Xsave(8,m)  = Sav(5,5,m)
          Xsave(9,m)  = Sav(6,6,m)
          Xsave(10,m) = Sav(1,4,m)
          Xsave(11,m) = Sav(2,5,m)
          Xsave(12,m) = Sav(3,6,m)
          Xsave(13,m) = Sav(3,4,m)
          Xsave(14,m) = Sav(1,5,m)
          Xsave(15,m) = Sav(2,6,m)
          Xsave(16,m) = Sav(2,4,m)
          Xsave(17,m) = Sav(3,5,m)
          Xsave(18,m) = Sav(1,6,m)
          Xsave(19,m) = Sav(5,6,m)
          Xsave(20,m) = Sav(4,6,m)
          Xsave(21,m) = Sav(4,5,m)
       END DO

    END IF

    !!! Write infos in hdf5 format
    write(*,*)
    write(*,"(a)",ADVANCE='no') ' SAVE DATA TO ',trim(filename)
    write(*,*)

    !Initialize FORTRAN interface.

    CALL H5open_f (error)

    ! Create a new file using default properties.

    CALL H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

    !Run input parameters
    dum_int(1)=dimensions
    dum_int(2)=cartspher
    dum_int(3)=basicstag
    dims1D = 3
    CALL loadsave_integer(1,1,file_id,dims1D,H5T_NATIVE_INTEGER,dum_int,'CoordinateSystem',1)
    dum_int(1)=x1periodic
    dum_int(2)=x2periodic
    dum_int(3)=x3periodic
    CALL loadsave_integer(1,1,file_id,dims1D,H5T_NATIVE_INTEGER,dum_int,'Periodic',1)
    dum_db(1)=dt
    dum_db(2)=timesum
    dims1D = 2
    CALL loadsave_double(1,1,file_id,dims1D,H5T_NATIVE_DOUBLE,dum_db,'Time',1)

    !Eulerian Nodes
    CALL H5Gcreate_f(file_id, "/Nodes", group_id, error)
    !Create Attribute
    dum_int(1)=nx1
    dum_int(2)=nx2
    dum_int(3)=nx3
    dims1D = 3
    CALL loadsave_integer(1,1,group_id,dims1D,H5T_NATIVE_INTEGER,dum_int,'nxyz',1)
    dims1D = nx1
    CALL loadsave_double(0,1,group_id,dims1D,H5T_NATIVE_DOUBLE,X1,'gx',1)
    dims1D = nx2
    CALL loadsave_double(0,1,group_id,dims1D,H5T_NATIVE_DOUBLE,X2,'gy',1)
    dims1D = nx3
    CALL loadsave_double(0,1,group_id,dims1D,H5T_NATIVE_DOUBLE,X3,'gz',1)

    CALL H5Gclose_f(group_id, error)

    !Lagrangian particles (crystal aggregates)

    !Replace deleted markers (i.e., with rocktype > 100)
    m1 = 0
    DO m = 1, marknum 
       IF(rocktype(m) < 100) THEN
          m1 = m1 + 1
          rocktype(m1) = rocktype(m)
          mx1(m1) = mx1(m)
          mx2(m1) = mx2(m)
          mx3(m1) = mx3(m)
          mYY(m1) = mYY(m)
          Fij(:,:,m1) = Fij(:,:,m)
          IF(fsemod == 0) THEN
             rho(m1) = rho(m)
             IF(ptmod > 0 ) THEN
                mtk(m1) = mtk(m)
                mpgpa(m1) = mpgpa(m)
             END IF
             Xsave(:,m1) = Xsave(:,m)
          END IF
       END IF
    END DO

    marknum = m1

    write(*,*)
    write(*,'(a,i0)') ' NUMBER OF AGGREGATES SAVED: ',marknum

    CALL H5Gcreate_f(file_id, "/Particles", group_id, error)
    dum_int(1)=marknum
    dims1D = 1
    CALL loadsave_integer(1,1,group_id,dims1D,H5T_NATIVE_INTEGER,dum_int(1:1),'marknum',1)
    dum_db(1)=mx1stp
    dum_db(2)=mx2stp
    dum_db(3)=mx3stp
    dims1D = 3
    CALL loadsave_double(0,1,group_id,dims1D,H5T_NATIVE_DOUBLE,dum_db(1:3),'mxstp',1)
    dum_db(1)=mx1min
    dum_db(2)=mx2min
    dum_db(3)=mx3min
    CALL loadsave_double(0,1,group_id,dims1D,H5T_NATIVE_DOUBLE,dum_db(1:3),'mxmin',1)
    dum_db(1)=mx1max
    dum_db(2)=mx2max
    dum_db(3)=mx3max
    CALL loadsave_double(0,1,group_id,dims1D,H5T_NATIVE_DOUBLE,dum_db(1:3),'mxmax',1)

    !Rocktype
    dims1D = marknum 
    CALL loadsave_integer(0,1,group_id,dims1D,H5T_NATIVE_INTEGER,rocktype(1:marknum),'rocktype',1)
    IF(yinyang == 2) CALL loadsave_integer(0,1,group_id,dims1D,H5T_NATIVE_INTEGER,mYY(1:marknum),'mYY',1)
    CALL loadsave_double(0,1,group_id,dims1D,H5T_NATIVE_DOUBLE,mx1(1:marknum),'mx1',1)
    CALL loadsave_double(0,1,group_id,dims1D,H5T_NATIVE_DOUBLE,mx2(1:marknum),'mx2',1)
    CALL loadsave_double(0,1,group_id,dims1D,H5T_NATIVE_DOUBLE,mx3(1:marknum),'mx3',1)

    !FSE
    nx(1)=3
    nx(2)=3
    nx(3)=marknum
    CALL loadsave_double(0,3,group_id,nx,H5T_NATIVE_DOUBLE,Fij(:,:,1:marknum),'fse',1)

    IF(fsemod == 0) THEN

    !Rho (kg/m^3)
    CALL loadsave_double(0,1,group_id,dims1D,H5T_NATIVE_DOUBLE,rho(1:marknum),'rho',1)
    !T(K)
    IF(ptmod > 0) CALL loadsave_double(0,1,group_id,dims1D,H5T_NATIVE_DOUBLE,mtk(1:marknum),'mtk',1)
    !P(GPa)
    IF(ptmod > 0) CALL loadsave_double(0,1,group_id,dims1D,H5T_NATIVE_DOUBLE,mpgpa(1:marknum),'mpgpa',1)
    !Sav(GPa)
    nx(1)=21
    nx(2)=marknum
    CALL loadsave_double(0,2,group_id,nx(1:2),H5T_NATIVE_DOUBLE,Xsave(:,1:marknum),'Sav',1)

    DEALLOCATE(Sav,Xsave)
    IF(ptmod > 0) DEALLOCATE(mtk,mpgpa)

    END IF

    CALL H5Gclose_f(group_id, error)
    !Terminate access to the file.

    CALL H5Fclose_f(file_id, error)

    CALL H5close_f(error)
    !Close FORTRAN interface.

   write(*,*)
   write(*,'(a,1f6.2,a)') ' CALCULATE ELASTIC TENSOR AND SAVE OUTPUT OK! (',OMP_GET_WTIME() - wtime,' sec)'
   write(*,*)
   write(*,"(a)") '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
   write(*,"(a)") '+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
   write(*,*)

   END SUBROUTINE stifftenzsave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Load/Save dataset, format double precision                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine SaveAttributeMPI_double(fileH5,dsetname,rank,dims,buffer)

   use class_precision
   use hdf5

   implicit none

   ! Input variables
   integer(hid_t),    intent(in) :: fileH5
   character(len=*),  intent(in) :: dsetname
   integer,           intent(in) :: rank
   integer,           intent(in) :: dims(rank)
   real(kind=double), intent(in) :: buffer(*)

   ! Local variables
   integer           :: i
   integer           :: errH5
   integer(hsize_t)  :: dimsH5(rank)
   integer(hid_t)    :: filespaceH5,dsetH5

   ! Retrieve total dimensions
   do i = 1,rank
      dimsH5(i) = dims(i)
   end do

   ! Create the data space for the dataset
   call h5screate_simple_f(rank,dimsH5,filespaceH5,errH5)

   ! Create the dataset with default properties
   call h5acreate_f(fileH5,dsetname,H5T_NATIVE_DOUBLE,filespaceH5,dsetH5,errH5)
   call h5sclose_f(filespaceH5,errH5)

   ! Write the dataset collectively
   call h5awrite_f(dsetH5,H5T_NATIVE_DOUBLE,buffer,dimsH5,errH5)

   ! Close resources
   call h5aclose_f(dsetH5,errH5)

   end subroutine SaveAttributeMPI_double

!M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine SaveAttributeMPI_integer(fileH5,dsetname,rank,dims,buffer)

   use hdf5

   implicit none

   ! Input variables
   integer(hid_t),   intent(in) :: fileH5
   character(len=*), intent(in) :: dsetname
   integer,          intent(in) :: rank
   integer,          intent(in) :: dims(rank)
   integer,          intent(in) :: buffer(*)

   ! Local variables
   integer           :: i
   integer           :: errH5
   integer(hsize_t)  :: dimsH5(rank)
   integer(hid_t)    :: filespaceH5,dsetH5

   ! Retrieve total dimensions
   do i = 1,rank
      dimsH5(i) = dims(i)
   end do

   ! Create the data space for the dataset
   call h5screate_simple_f(rank,dimsH5,filespaceH5,errH5)

   ! Create the dataset with default properties
   call h5acreate_f(fileH5,dsetname,H5T_NATIVE_INTEGER,filespaceH5,dsetH5,errH5)
   call h5sclose_f(filespaceH5,errH5)

   ! Write the dataset collectively
   call h5awrite_f(dsetH5,H5T_NATIVE_INTEGER,buffer,dimsH5,errH5)

   ! Close resources
   call h5aclose_f(dsetH5,errH5)

   end subroutine SaveAttributeMPI_integer

!M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine SaveDatasetMPI_double(create,fileH5,dsetname,rank,dims,totdim,offset,buffer)

   use class_precision
   use hdf5

   implicit none

   ! Input variables
   logical,           intent(in) :: create
   integer(hid_t),    intent(in) :: fileH5
   character(len=*),  intent(in) :: dsetname
   integer,           intent(in) :: rank
   integer,           intent(in) :: dims(rank)
   integer,           intent(in) :: totdim
   integer,           intent(in) :: offset(rank)
   real(kind=double), intent(in) :: buffer(*)

   ! Local variables
   integer           :: i
   integer           :: errH5
   integer(hsize_t)  :: dimsH5(rank)
   integer(hid_t)    :: filespaceH5,dsetH5,memspaceH5,plistH5
   integer(hssize_t) :: offsetH5(rank)

   ! Retrieve total dimensions
   do i = 1,rank-1
      dimsH5(i) = dims(i)
   end do
   dimsH5(rank) = totdim

   if (create) then

      ! Create the data space for the dataset
      call h5screate_simple_f(rank,dimsH5,filespaceH5,errH5)

      ! Create the dataset with default properties
      call h5dcreate_f(fileH5,dsetname,H5T_NATIVE_DOUBLE,filespaceH5,dsetH5,errH5)
      call h5sclose_f(filespaceH5,errH5)

   else

      ! Open dataset with default properties
      call h5dopen_f(fileH5,dsetname,dsetH5,errH5)

   endif

   ! Retrieve local dimensions
   dimsH5(rank) = dims(rank) 

   ! Create memspace
   CALL h5screate_simple_f(rank,dimsH5,memspaceH5,errH5)

   ! Retrieve local offset
   do i = 1,rank
      offsetH5(i) = offset(i)
   end do

   ! Select hyperslab in the file
   call h5dget_space_f(dsetH5,filespaceH5,errH5)
   call h5sselect_hyperslab_f(filespaceH5,H5S_SELECT_SET_F,offsetH5,dimsH5,errH5)

   ! Create property list for collective dataset write
   call h5pcreate_f(H5P_DATASET_XFER_F,plistH5,errH5) 
   call h5pset_dxpl_mpio_f(plistH5,H5FD_MPIO_COLLECTIVE_F,errH5)
 
   ! Retrieve total dimensions
   dimsH5(rank) = totdim 
  
   ! Write the dataset collectively
   call h5dwrite_f(dsetH5,H5T_NATIVE_DOUBLE,buffer,dimsH5,errH5,                          &
   &               file_space_id=filespaceH5,mem_space_id=memspaceH5,xfer_prp=plistH5)

   ! Close resources
   call h5sclose_f(filespaceH5,errH5)
   call h5sclose_f(memspaceH5,errH5)
   call h5dclose_f(dsetH5,errH5)
   call h5pclose_f(plistH5,errH5)

   end subroutine SaveDatasetMPI_double

!M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine SaveDatasetMPI_integer(create,fileH5,dsetname,rank,dims,totdim,offset,buffer)

   use hdf5

   implicit none

   ! Input variables
   logical,          intent(in) :: create
   integer(hid_t),   intent(in) :: fileH5
   character(len=*), intent(in) :: dsetname
   integer,          intent(in) :: rank
   integer,          intent(in) :: dims(rank)
   integer,          intent(in) :: totdim
   integer,          intent(in) :: offset(rank)
   integer,          intent(in) :: buffer(*)

   ! Local variables
   integer           :: i
   integer           :: errH5
   integer(hsize_t)  :: dimsH5(rank)
   integer(hid_t)    :: filespaceH5,dsetH5,memspaceH5,plistH5
   integer(hssize_t) :: offsetH5(rank)

   ! Retrieve total dimensions
   do i = 1,rank-1
      dimsH5(i) = dims(i)
   end do
   dimsH5(rank) = totdim

   if (create) then

      ! Create the data space for the dataset
      call h5screate_simple_f(rank,dimsH5,filespaceH5,errH5)

      ! Create the dataset with default properties
      call h5dcreate_f(fileH5,dsetname,H5T_NATIVE_INTEGER,filespaceH5,dsetH5,errH5)
      call h5sclose_f(filespaceH5,errH5)

   else

      ! Open dataset with default properties
      call h5dopen_f(fileH5,dsetname,dsetH5,errH5)

   endif

   ! Retrieve local dimensions
   dimsH5(rank) = dims(rank) 

   ! Create memspace
   CALL h5screate_simple_f(rank,dimsH5,memspaceH5,errH5)

   ! Retrieve local offset
   do i = 1,rank
      offsetH5(i) = offset(i)
   end do

   ! Select hyperslab in the file
   call h5dget_space_f(dsetH5,filespaceH5,errH5)
   call h5sselect_hyperslab_f(filespaceH5,H5S_SELECT_SET_F,offsetH5,dimsH5,errH5)

   ! Create property list for collective dataset write
   call h5pcreate_f(H5P_DATASET_XFER_F,plistH5,errH5) 
   call h5pset_dxpl_mpio_f(plistH5,H5FD_MPIO_COLLECTIVE_F,errH5)
 
   ! Retrieve total dimensions
   dimsH5(rank) = totdim 
  
   ! Write the dataset collectively
   call h5dwrite_f(dsetH5,H5T_NATIVE_INTEGER,buffer,dimsH5,errH5,                         &
   &               file_space_id=filespaceH5,mem_space_id=memspaceH5,xfer_prp=plistH5)

   ! Close resources
   call h5sclose_f(filespaceH5,errH5)
   call h5sclose_f(memspaceH5,errH5)
   call h5dclose_f(dsetH5,errH5)
   call h5pclose_f(plistH5,errH5)

   end subroutine SaveDatasetMPI_integer

!M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE loadsave_double(data_attr,rank,loc_id,dims,memtype,buf1,path,mode)

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

