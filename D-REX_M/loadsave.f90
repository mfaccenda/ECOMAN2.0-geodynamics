 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !! ---------------------------------------------------------------------------
 !!
 !!    Copyright (c) 2018-2020, Università di Padova, Manuele Faccenda
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
   USE hdf5   

   IMPLICIT NONE
 
   INTEGER :: i1,i2,i3,i4,t,n1,n2,n3,nx(3),yy
   CHARACTER (4) :: dt_str4
   CHARACTER (500) :: filename,str
   CHARACTER (len=*) :: input_dir
   CHARACTER (len(trim(input_dir))) :: str1  
   DOUBLE PRECISION, DIMENSION(2,2) :: Adummy2,dV2
   DOUBLE PRECISION, DIMENSION(3,3) :: Adummy,dV
   DOUBLE PRECISION :: dum,sii,sigmin,sigmax,datadb(2),lr,rr,cr,x2y2,x2z2,xx,yyy,zz,Uxx(3)
   DOUBLE PRECISION ,DIMENSION(:), ALLOCATABLE :: V1,V2,V3,Fd0,Tk0,Pa0
   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id,dcpl,memtype !Handles
   INTEGER     ::   error  ! Error flag

   OPEN(19,file='cycle.txt',status='replace')

   write(19,"(a)"), 'VTP file number'
   write(19,'(i4)') t
   write(19,"(a)"), 'Cycle'
   CLOSE(19)

   str1=trim(input_dir)
   write(dt_str4,'(i4.4)') t
   str='vtp'//dt_str4//'.h5'
   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)
   write(*,'(a,a)'),' LOAD FILE ',trim(str)
   write(*,*)

   filename = str1//str

   !Load input files in hdf5 format

   !Initialize FORTRAN interface.

   CALL H5open_f (error)

   ! Open a new file using default properties.

   CALL H5Fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

   !Load input parameters
   CALL loadsave_double(1,1,file_id,2,H5T_NATIVE_DOUBLE,datadb,'Time',0)
   dt0=datadb(1)
   timesum=datadb(2)

   IF(t == 1) timesum0 = 0
   IF(Tstep > 1) dt0=timesum - timesum0
  
   timesum0 = timesum
 
   !Load velocity field
   yy = yinyang
   ALLOCATE(V1(nodenum*yy),V2(nodenum*yy),V3(nodenum*yy))
   CALL loadsave_double(0,1,file_id,nodenum*yy,H5T_NATIVE_DOUBLE,V1,'Nodes/V1',0)!Vx/VPhi
   CALL loadsave_double(0,1,file_id,nodenum*yy,H5T_NATIVE_DOUBLE,V2,'Nodes/V2',0)!Vy/Vradial
   IF(dimensions == 3) CALL loadsave_double(0,1,file_id,nodenum*yy,H5T_NATIVE_DOUBLE,V3,'Nodes/V3',0)!Vz/Vcolat 

   !Load temperature and pressure
   IF(ptmod > 0) THEN 
      ALLOCATE(Tk0(nodenum*yy),Pa0(nodenum*yy)) 
      CALL loadsave_double(0,1,file_id,nodenum*yy,H5T_NATIVE_DOUBLE,Tk0,'Nodes/Tk',0)
      CALL loadsave_double(0,1,file_id,nodenum*yy,H5T_NATIVE_DOUBLE,Pa0,'Nodes/P',0)
   END IF

   !Load fraction of deformation accommodated by dislocation creep
   IF(fractdislmod > 0) THEN 
      ALLOCATE(Fd0(nodenum*yy))
      CALL loadsave_double(0,1,file_id,nodenum*yy,H5T_NATIVE_DOUBLE,Fd0,'Nodes/Fd',0)
   END IF

   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   !Close FORTRAN interface.
   CALL H5close_f(error)

   Ui = 0d0 ; Ux = 0d0

   DO yy = 1, yinyang

   DO i1 = 1, nx1
      DO i2 = 1, nx2
         DO i3 = 1, nx3
         
         !Global index
         i4 = nodenum*(yy-1) + i2 + (i1-1)*nx2 + (i3-1)*nx1*nx2

         Ui(yy,1,i1,i2,i3)=V1(i4)
         Ui(yy,2,i1,i2,i3)=V2(i4)
         IF(dimensions == 3) Ui(yy,3,i1,i2,i3)=V3(i4)

         IF(ptmod > 0) THEN
            Tk(yy,i1,i2,i3) = Tk0(i4)
            Pa(yy,i1,i2,i3) = Pa0(i4)
         END IF

         IF(fractdislmod > 0) THEN
            Fd(yy,i1,i2,i3) = Fd0(i4)
            IF(Fd(yy,i1,i2,i3)<0.0) write(*,"(a,i4,i4,i4,f6.3)"),'Fraction of disl. creep < 0 at node ',i1,i2,i3,Fd(yy,i1,i2,i3)
            IF(Fd(yy,i1,i2,i3)>1.0) write(*,"(a,i4,i4,i4,f6.3)"),'Fraction of disl. creep > 1 at node ',i1,i2,i3,Fd(yy,i1,i2,i3)
         END IF

         END DO
      END DO
   END DO

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
                  IF (x1periodic > 0 .AND. i1 == 1) Dij(1,2,1,i1,i2,1) = (Ui(1,2,i1,i2,1)-Ui(1,2,nx1-1,i2,1))/(X1(i1+1)-X1(i1)+X1(nx1)-X1(nx1-1))*2d0
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
                IF(i1 > 1 .AND. i1 < nx1) dV2(1,1) = ((Ux(1,1,i1+1,i2,1)-Ux(1,1,i1,i2,1))/(X1(i1+1)-X1(i1)) + (Ux(1,1,i1,i2,1)-Ux(1,1,i1-1,i2,1))/(X1(i1)-X1(i1-1)))/2

                IF(x1periodic > 0 .AND. (i1 == 1 .OR. i1 == nx1)) dV2(1,1) = ((Ux(1,1,2,i2,1)-Ux(1,1,1,i2,1))/(X1(2)-X1(1)) + (Ux(1,1,nx1,i2,1)-Ux(1,1,nx1-1,i2,1))/(X1(nx1)-X1(nx1-1)))/2

                !dVx/dR    
                IF (i2 > 1 .AND. i2 < nx2) dV2(1,2) = ((Ux(1,1,i1,i2+1,1)-Ux(1,1,i1,i2,1))/(X2(i2+1)-X2(i2)) + (Ux(1,1,i1,i2,1)-Ux(1,1,i1,i2-1,1))/(X2(i2)-X2(i2-1)))/2

                IF(x2periodic > 0 .AND. (i2 == 1 .OR. i2 == nx2)) dV2(1,2) = ((Ux(1,1,i1,2,1)-Ux(1,1,i1,1,1))/(X2(2)-X2(1)) + (Ux(1,1,i1,nx2,1)-Ux(1,1,i1,nx2-1,1))/(X2(nx2)-X2(nx2-1)))/2
                
                !dVy in polar coordinates
                !dVy/dPhi  
                IF(i1 > 1 .AND. i1 < nx1) dV2(2,1) = ((Ux(1,2,i1+1,i2,1)-Ux(1,2,i1,i2,1))/(X1(i1+1)-X1(i1)) + (Ux(1,2,i1,i2,1)-Ux(1,2,i1-1,i2,1))/(X1(i1)-X1(i1-1)))/2

                IF(x1periodic > 0 .AND. (i1 == 1 .OR. i1 == nx1)) dV2(2,1) = ((Ux(1,2,2,i2,1)-Ux(1,2,1,i2,1))/(X1(2)-X1(1)) + (Ux(1,2,nx1,i2,1)-Ux(1,2,nx1-1,i2,1))/(X1(nx1)-X1(nx1-1)))/2

                !dVy/dR    
                IF (i2 > 1 .AND. i2 < nx2) dV2(2,2) = ((Ux(1,2,i1,i2+1,1)-Ux(1,2,i1,i2,1))/(X2(i2+1)-X2(i2)) + (Ux(1,2,i1,i2,1)-Ux(1,2,i1,i2-1,1))/(X2(i2)-X2(i2-1)))/2

                IF(x2periodic > 0 .AND. (i2 == 1 .OR. i2 == nx2)) dV2(2,2) = ((Ux(1,2,i1,2,1)-Ux(1,2,i1,1,1))/(X2(2)-X2(1)) + (Ux(1,2,i1,nx2,1)-Ux(1,2,i1,nx2-1,1))/(X2(nx2)-X2(nx2-1)))/2

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

   END DO

   END IF

!!! Velocity gradient tensor
   Dij = 0d0

   DO yy = 1, yinyang

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
                IF (i1 > 1 .AND. i1 < nx1 .AND. i3 < nx3) Dij(yy,2,1,i1,i2,i3) = (Ui(yy,2,i1,i2,i3)-Ui(yy,2,i1-1,i2,i3))/(X1(i1+1)-X1(i1-1))*2d0
            END IF
            !dVz/dx
            !At (0,+ystp/2,0)
            IF (x1periodic > 0 .and. yinyang == 1 .and. i1 == 1 .and. i2 < nx2) THEN
                Dij(yy,3,1,i1,i2,i3) = (Ui(yy,3,1,i2,i3)-Ui(yy,3,nx1-1,i2,i3))/(X1(i1+1)-X1(i1)+X1(nx1)-X1(nx1-1))*2d0
                Dij(yy,3,1,nx1,i2,i3) = Dij(yy,3,1,i1,i2,i3)
            ELSE 
                IF (i1 > 1 .AND. i1 < nx1 .AND. i2 < nx2) Dij(yy,3,1,i1,i2,i3) = (Ui(yy,3,i1,i2,i3)-Ui(yy,3,i1-1,i2,i3))/(X1(i1+1)-X1(i1-1))*2d0
            END IF

            !dVx/dy
            !At (0,0,+zstp/2)
            IF (x2periodic > 0 .and. i2 == 1 .and. i3 < nx3) THEN
                Dij(yy,1,2,i1,i2,i3) = (Ui(yy,1,i1,1,i3)-Ui(yy,1,i1,nx2-1,i3))/(X2(i2+1)-X2(i2)+X2(nx2)-X2(nx2-1))*2d0
                Dij(yy,1,2,i1,nx2,i3) = Dij(yy,1,2,i1,i2,i3)
            ELSE 
                IF (i2 > 1 .AND. i2 < nx2 .AND. i3 < nx3) Dij(yy,1,2,i1,i2,i3) = (Ui(yy,1,i1,i2,i3)-Ui(yy,1,i1,i2-1,i3))/(X2(i2+1)-X2(i2-1))*2d0
            END IF
    
            !dVz/dy
            !At (+xstp/2,0,0)
            IF (x2periodic > 0 .and. i2 == 1 .and. i1 < nx1) THEN
                Dij(yy,3,2,i1,i2,i3) = (Ui(yy,3,i1,1,i3)-Ui(yy,3,i1,nx2-1,i3))/(X2(i2+1)-X2(i2)+X2(nx2)-X2(nx2-1))*2d0
                Dij(yy,3,2,i1,nx2,i3) = Dij(yy,3,2,i1,i2,i3)
            ELSE 
                IF (i2 > 1 .AND. i2 < nx2 .AND. i1 < nx1) Dij(yy,3,2,i1,i2,i3) = (Ui(yy,3,i1,i2,i3)-Ui(yy,3,i1,i2-1,i3))/(X2(i2+1)-X2(i2-1))*2d0
            END IF

            !dVx/dz
            !At (0,+ystp/2,0)
            IF (x3periodic > 0 .and. yinyang == 1 .and. i3 == 1 .and. i2 < nx2) THEN
                Dij(yy,1,3,i1,i2,i3) = (Ui(yy,1,i1,i2,1)-Ui(yy,1,i1,i2,nx3-1))/(X3(i3+1)-X3(i3)+X3(nx3)-X3(nx3-1))*2d0
                Dij(yy,1,3,i1,i2,nx3) = Dij(yy,1,3,i1,i2,i3)
            ELSE 
                IF (i3 > 1 .AND. i3 < nx3 .AND. i2 < nx2) Dij(yy,1,3,i1,i2,i3) = (Ui(yy,1,i1,i2,i3)-Ui(yy,1,i1,i2,i3-1))/(X3(i3+1)-X3(i3-1))*2d0
            END IF

            !dVy/dz
            !At (+xstp/2,0,0)
            IF (x3periodic > 0 .and. yinyang == 1 .and. i3 == 1 .and. i1 < nx1) THEN
                Dij(yy,2,3,i1,i2,i3) = (Ui(yy,2,i1,i2,1)-Ui(yy,2,i1,i2,nx3-1))/(X3(i3+1)-X3(i3)+X3(nx3)-X3(nx3-1))*2d0
                Dij(yy,2,3,i1,i2,nx3) = Dij(yy,2,3,i1,i2,i3)
            ELSE 
                IF (i3 > 1 .AND. i3 < nx3 .AND. i1 .LT. nx1) Dij(yy,2,3,i1,i2,i3) = (Ui(yy,2,i1,i2,i3)-Ui(yy,2,i1,i2,i3-1))/(X3(i3+1)-X3(i3-1))*2d0
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
                IF (x2periodic .AND. (i2 == 1 .OR. i2 == nx2)) THEN 
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
                   !Dij(yy,3,3,i1,i2,i3) = -Dij(yy,1,1,i1,i2,i3)-Dij(yy,2,2,i1,i2,i3)

                END IF

            END IF

         END DO
      END DO
   END DO

   END DO

!!! Convert Vcolat and Vphi from m/s to rad/s for advection
   IF(cartspher == 2) THEN

   DO yy = 1, yinyang

   DO i1 = 1, nx1   
      DO i2 = 1, nx2    
         DO i3 = 1, nx3    
            IF(yy == 1) THEN
               !Vcolat(rad/s) = Vcolat(m/s)/R
               IF(Ui(yy,3,i1,i2,i3) /= 0) Ui(yy,3,i1,i2,i3)=Ui(yy,3,i1,i2,i3)/X2(i2)               
               !Vphi(rad/s) = Vphi(m/s)/R/sin(theta)
               IF(Ui(yy,1,i1,i2,i3) /= 0 .AND. X3(i3) /= 0 .AND. X3(i3) /= pi) Ui(yy,1,i1,i2,i3)=Ui(yy,1,i1,i2,i3)/X2(i2)/sin(X3(i3))    
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
      write(*,'(a,1es13.6)'),' TOTAL TIMESTEP     = ',dt0
      write(*,*)
      write(*,'(a,1es13.6,a,i)'),' ADVECTION TIMESTEP = ',dt,', NUMBER OF CYCLES = ',NINT(timemax/dt)
      write(*,*)
   ELSE
      write(*,'(a,1es13.6)'),' TOTAL TIMESTEP     = ',dt0
      write(*,*)
      write(*,'(a,1es13.6,a,i)'),' ADVECTION TIMESTEP = ',dt,', NUMBER OF CYCLES = ',numdt
      write(*,*)
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

   SUBROUTINE stifftenzsave(t,output_dir,ncyc2)

   USE comvar
   USE omp_lib
   USE hdf5    

   IMPLICIT NONE

   INTEGER :: i,j,m,m1,t,nx(3),ncyc2,dum_int(3)
   DOUBLE PRECISION :: wtime,mtkm,mpgpam,dum_db(5)
   DOUBLE PRECISION ,DIMENSION(:), ALLOCATABLE :: mtk,mpgpa,rt1
   DOUBLE PRECISION ,DIMENSION(:,:), ALLOCATABLE :: Xsave
   CHARACTER (500) :: filename
   CHARACTER (len=*) :: output_dir
   CHARACTER (len(trim(output_dir))) :: str1
   CHARACTER (4) :: dt_str4
   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id, dcpl,memtype ! Handles
   INTEGER     ::   error  ! Error flag


   str1=trim(output_dir)

   write(dt_str4,'(i4.4)') t  
   filename=str1//'Cijkl'//dt_str4//'.h5'

!!! Calculate Anisotropy 
   
   write(*,*)
   write(*,"(a)"),'+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
   write(*,"(a)"),'+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
   write(*,*)
   write(*,'(a,i0,a)'),' CALCULATE ELASTIC TENSOR AND SAVE OUTPUT AFTER ',ncyc2,' CYCLES'
   
   wtime = OMP_GET_WTIME()

   IF(fsemod == 0) THEN

   ALLOCATE(Sav(6,6,marknum))

!$omp parallel & 
!$omp shared(mx1,mx2,mx3,mYY,rocktype,Xol,odf,odf_ens,acs,acs_ens,Tk,Pa,Fd,td_rho,rho,Sav,S0,dS0dp,dS0dp2,dS0dt,dS0dpt,dS0dp5,dS0dt5,ijkl) &
!$omp private(m) &    
!$omp firstprivate(size,marknum,ptmod,uppermantlemod)
!$omp do schedule(guided,8)
   DO m = 1 , marknum

      IF(rocktype(m) > 0 .AND. rocktype(m) < 10) then

!!! Cijkl tensor (using Voigt average)
      CALL stifftenz(m,Sav(:,:,m))

      END IF

   END DO

!$omp end do
!$omp end parallel

    END IF

    !!! Write infos in hdf5 format
    write(*,*)
    write(*,"(a)",ADVANCE='no'),' SAVE DATA TO ',trim(filename)
    write(*,*)

    !Initialize FORTRAN interface.

    CALL H5open_f (error)

    ! Create a new file using default properties.

    CALL H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

    !Run input parameters
    dum_int(1)=dimensions
    dum_int(2)=cartspher
    dum_int(3)=basicstag
    CALL loadsave_integer(1,1,file_id,3,H5T_NATIVE_INTEGER,dum_int,'CoordinateSystem',1)
    dum_int(1)=x1periodic
    dum_int(2)=x2periodic
    dum_int(3)=x3periodic
    CALL loadsave_integer(1,1,file_id,3,H5T_NATIVE_INTEGER,dum_int,'Periodic',1)
    dum_db(1)=dt
    dum_db(2)=timesum
    CALL loadsave_double(1,1,file_id,2,H5T_NATIVE_DOUBLE,dum_db,'Time',1)

    !Eulerian Nodes
    CALL H5Gcreate_f(file_id, "/Nodes", group_id, error)
    !Create Attribute
    dum_int(1)=nx1
    dum_int(2)=nx2
    dum_int(3)=nx3
    CALL loadsave_integer(1,1,group_id,3,H5T_NATIVE_INTEGER,dum_int,'nxyz',1)
    CALL loadsave_double(0,1,group_id,nx1,H5T_NATIVE_DOUBLE,X1,'gx',1)
    CALL loadsave_double(0,1,group_id,nx2,H5T_NATIVE_DOUBLE,X2,'gy',1)
    CALL loadsave_double(0,1,group_id,nx3,H5T_NATIVE_DOUBLE,X3,'gz',1)

    CALL H5Gclose_f(group_id, error)

    !Lagrangian particles (crystal aggregates)

    IF(fsemod == 0) THEN

    ALLOCATE(Xsave(21,marknum))
    IF(ptmod > 0) ALLOCATE(mtk(marknum),mpgpa(marknum))
    Xsave  = 0d0 ; mtk = 0d0 ; mpgpa = 0d0

    DO m = 1, marknum
       !Calculate T,P

       IF(ptmod > 0) CALL rhopt(m,mtk(m),mpgpa(m))

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
    write(*,'(a,i0)'),' NUMBER OF AGGREGATES SAVED: ',marknum

    CALL H5Gcreate_f(file_id, "/Particles", group_id, error)
    dum_int(1)=marknum
    CALL loadsave_integer(1,1,group_id,1,H5T_NATIVE_INTEGER,dum_int(1:1),'marknum',1)
    dum_db(1)=mx1stp
    dum_db(2)=mx2stp
    dum_db(3)=mx3stp
    CALL loadsave_double(0,1,group_id,3,H5T_NATIVE_DOUBLE,dum_db(1:3),'mxstp',1)
    dum_db(1)=mx1min
    dum_db(2)=mx2min
    dum_db(3)=mx3min
    CALL loadsave_double(0,1,group_id,3,H5T_NATIVE_DOUBLE,dum_db(1:3),'mxmin',1)
    dum_db(1)=mx1max
    dum_db(2)=mx2max
    dum_db(3)=mx3max
    CALL loadsave_double(0,1,group_id,3,H5T_NATIVE_DOUBLE,dum_db(1:3),'mxmax',1)

    !Rocktype
    CALL loadsave_integer(0,1,group_id,marknum,H5T_NATIVE_INTEGER,rocktype(1:marknum),'rocktype',1)
    IF(yinyang == 2) CALL loadsave_integer(0,1,group_id,marknum,H5T_NATIVE_INTEGER,mYY(1:marknum),'mYY',1)
    CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mx1(1:marknum),'mx1',1)
    CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mx2(1:marknum),'mx2',1)
    CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mx3(1:marknum),'mx3',1)

    !FSE
    nx(1)=3
    nx(2)=3
    nx(3)=marknum 
    CALL loadsave_double(0,3,group_id,nx,H5T_NATIVE_DOUBLE,Fij(:,:,1:marknum),'fse',1)

    IF(fsemod == 0) THEN

    !Rho (kg/m^3)
    CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,rho(1:marknum),'rho',1)
    !T(K)
    IF(ptmod > 0) CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mtk(1:marknum),'mtk',1)
    !P(GPa)
    IF(ptmod > 0) CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mpgpa(1:marknum),'mpgpa',1)
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
   write(*,'(a,1f6.2,a)'),' CALCULATE ELASTIC TENSOR AND SAVE OUTPUT OK! (',OMP_GET_WTIME() - wtime,' sec)'
   write(*,*)
   write(*,"(a)"),'+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
   write(*,"(a)"),'+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+'
   write(*,*)

   END SUBROUTINE stifftenzsave

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

   SUBROUTINE stifftenz(m,Sav1)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i,j,k,ll,n,nu,p,q,r,ss,m
   DOUBLE PRECISION :: mtk,mpgpa,Mcur,Gcur,Lcur,Mnew,Gnew,Lnew 
   DOUBLE PRECISION, DIMENSION(21) :: XE
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C0,Cav,Cav2,CMP,CMP0,CMP2
   DOUBLE PRECISION, DIMENSION(6,6) :: Sav1,CMP6,Voigt,Reuss
  
10 C0 = 0d0 ; Cav = 0d0 ; CMP = 0d0 ; CMP0 = 0d0 ; CMP2 = 0d0 
   Sav1 = 0d0 ; CMP6 = 0d0 ; Voigt = 0d0 ; Reuss = 0d0 ; 

   !!! P-T conditions and rho
   mtk = 0d0 ; mpgpa = 0d0
   IF(ptmod > 0) CALL rhopt(m,mtk,mpgpa)

   !!! Major phase

   !!! converting (6,6) to (3,3,3,3) matrix
   CALL S0_2_C0(1,m,mtk,mpgpa,C0,CMP0)

   DO nu = 1 , size
      
      IF(fractvoigt > 0d0) THEN 
         !!! VOIGT
         !!! C0: single crystal stiffness tensor
         !!! Cav: aggregate stiffness tensor
         CALL Caverage(1,m,nu,C0,Cav) 
      END IF

      IF(fractvoigt < 1d0) THEN
         !!! REUSS
         !!! CMP0: single crystal inverse stiffness tensor
         !!! CMP: aggregate inverse stiffness tensor
         CALL Caverage(1,m,nu,CMP0,CMP)
      END IF 
   
   END DO

   !!! Minor phase
   IF(Xol(rocktype(m)) < 1d0) THEN

   !!! converting (6,6) to (3,3,3,3) matrix
   CALL S0_2_C0(2,m,mtk,mpgpa,C0,CMP0)
   
      DO nu = 1 , size

         IF(fractvoigt > 0d0) THEN 
            !!! VOIGT
            !!! C0: single crystal stiffness tensor
            !!! Cav: aggregate stiffness tensor
            CALL Caverage(2,m,nu,C0,Cav) 
         END IF
         
         IF(fractvoigt < 1d0) THEN 
            !!! REUSS
            !!! CMP0: single crystal inverse stiffness tensor
            !!! CMP: aggregate inverse stiffness tensor
            CALL Caverage(2,m,nu,CMP0,CMP) 
         END IF
   
      END DO

   END IF

   !Reset fabric if bulk moduli are NaN
   IF(isnan(Cav(1,1,1,1)) .OR. (Cav(1,1,1,1)>400d0 .AND. rocktype(m)==1)) THEN
      print *,"Nan at",m,rocktype(m),mx1(m),mx2(m),mx3(m),mtk,mpgpa,odf(m,1),acs(1,1,1,m),acs_ens(1,1,1,m),Cav(1,1,1,1)
      acs(:,:,:,m) = acs0
      acs_ens(:,:,:,m) = acs0
      odf(m,:) = 1d0/REAL(size3**3)
      odf_ens(m,:) = odf(m,:)
      goto 10
   END IF

   !!! Average stiffness matrix
   DO i = 1 , 6 ; DO j = 1 , 6
      Voigt(i,j) = Cav(l1(i),l2(i),l1(j),l2(j))
   END DO ; END DO

   !!! Inverse of CMP
   IF(fractvoigt < 1d0) THEN
   
      DO i = 1 , 6 ; DO j = 1 , 6
         CMP6(i,j) = CMP(l1(i),l2(i),l1(j),l2(j))*mandel_scale(i,j)
      END DO ; END DO

      CALL inverse(CMP6,Reuss,6)

      DO i = 1 , 6 ; DO j = 1 , 6
         Sav1(i,j) = fractvoigt*Voigt(i,j) + (1d0-fractvoigt)*Reuss(i,j)
      END DO ; END DO

   ELSE

      Sav1 = Voigt  

   END IF

!!! Correct isotropic component of aggregate elastic tensor with that from thermodynamic database
   IF(ptmod > 0) THEN

      !Find current isotropic elastic moduli
      XE = 0d0 
      CALL V21D(Sav1,XE)
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

   END IF 
!!! End correct isotropic component of single crystal elastic tensor with that from thermodynamic database

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

      DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
         C0(i,j,k,ll) = S0(mtype,ijkl(i,j),ijkl(k,ll))
      END DO ; END DO ; END DO ; END DO
  
      IF(fractvoigt < 1d0) THEN
         S0dummy=S0(mtype,:,:)
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

      DO i = 1, 6 ; DO j = 1 , 6
         S0pt(i,j) = S0(mtype,i,j) - dSdt(i,j)*DTk + &
               dSdp(i,j)*DP + 0.5*dS0dp2(mtype,i,j)*(DP**2) + dS0dpt(mtype,i,j)*DP*DTk
      END DO ; END DO

      DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
         C0(i,j,k,ll) = S0pt(ijkl(i,j),ijkl(k,ll))
      END DO ; END DO ; END DO ; END DO
   
      IF(fractvoigt < 1d0) THEN
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

   SUBROUTINE V21D(C,X)

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

   END SUBROUTINE V21D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
