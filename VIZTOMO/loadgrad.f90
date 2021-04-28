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
!!! Loading hdf5 files                                                     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE load(dt_str4,cijkl_dir)

   USE comvar
   USE hdf5   

   IMPLICIT NONE
 
   INTEGER :: i1,i2,i3,i4,t,n1,n2,n3,nx(3),nodenum0,yy
   CHARACTER (4) :: dt_str4
   CHARACTER (500) :: filename,str
   CHARACTER (len=*) :: cijkl_dir
   CHARACTER (len(trim(cijkl_dir))) :: str1  
   DOUBLE PRECISION, DIMENSION(2,2) :: Adummy2,dV2
   DOUBLE PRECISION, DIMENSION(3,3) :: Adummy,dV
   DOUBLE PRECISION :: dum,datadb(2),lr,rr,cr,xx,yyy,zz,x2z2
   DOUBLE PRECISION ,DIMENSION(:), ALLOCATABLE :: V1,V2,V3!,Fd0,Tk0,Pa0
   INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id,dcpl,memtype !Handles
   INTEGER     ::   error  ! Error flag

   !!!WARNING
   basicstag = 1
   
   str1=trim(cijkl_dir)
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

   !Load velocity field
   nodenum0 = nx1*nx2*nx3

   yy = yinyang
   ALLOCATE(Ui(yy,3,nx1,nx2,nx3),Ux(yy,3,nx1,nx2,nx3),Dij(yy,3,3,nx1,nx2,nx3))
   ALLOCATE(V1(nodenum0*yy),V2(nodenum0*yy),V3(nodenum0*yy))
   ALLOCATE(X1(nx1),X2(nx2),X3(nx3))
   CALL loadsave_double(0,1,file_id,nodenum0*yy,H5T_NATIVE_DOUBLE,V1,'Nodes/V1',0)!Vx/VPhi
   CALL loadsave_double(0,1,file_id,nodenum0*yy,H5T_NATIVE_DOUBLE,V2,'Nodes/V2',0)!Vy/Vradial
   IF(dimensions == 3) CALL loadsave_double(0,1,file_id,nodenum0*yy,H5T_NATIVE_DOUBLE,V3,'Nodes/V3',0)!Vz/Vcolat 

   X1 = X1n
   X2 = X2n
   X3 = X3n

   !Terminate access to the file.
   CALL H5Fclose_f(file_id, error)

   !Close FORTRAN interface.
   CALL H5close_f(error)

   Ui = 0d0 ; Ux = 0d0

   DO yy = 1 , yinyang

   DO i1 = 1, nx1
      DO i2 = 1, nx2
         DO i3 = 1, nx3
         
         !Global index
         i4 = nodenum0*(yy-1) + i2 + (i1-1)*nx2 + (i3-1)*nx1*nx2

         Ui(yy,1,i1,i2,i3)=V1(i4)
         Ui(yy,2,i1,i2,i3)=V2(i4)
         IF(dimensions == 3) Ui(yy,3,i1,i2,i3)=V3(i4)

         END DO
      END DO
   END DO

   END DO

   DEALLOCATE(V1,V2)
   IF(dimensions == 3) DEALLOCATE(V3)

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
                  IF (i2 < nx2) Dij(1,1,1,i1+1,i2+1,1) = (Ux(1,1,i1+1,i2,1)-Ux(1,1,i1,i2,1))/(X1(i1+1)-X1(i1))
                  !At (0,0)
                  IF (i1 > 1) Dij(1,2,1,i1,i2,1) = (Ux(1,2,i1,i2,1)-Ux(1,2,i1-1,i2,1))/(X1(i1+1)-X1(i1-1))*2d0
                  IF (x1periodic > 0 .AND. i1 == 1) Dij(1,2,1,i1,i2,1) = (Ux(1,2,i1,i2,1)-Ux(1,2,nx1-1,i2,1))/(X1(i1+1)-X1(i1)+X1(nx1)-X1(nx1-1))*2d0
               END IF

               !dVi/dy
               IF (i2 .LT. nx2) THEN
                  !At (0,0)
                  IF (i2 > 1) Dij(1,1,2,i1,i2,1) = (Ux(1,1,i1,i2,1)-Ux(1,1,i1,i2-1,1))/(X2(i2+1)-X2(i2-1))*2d0
                  !At (+xstp/2,+ystp/2)
                  IF (i1 < nx1) Dij(1,2,2,i1+1,i2+1,1) = (Ux(1,2,i1,i2+1,1)-Ux(1,2,i1,i2,1))/(X2(i2+1)-X2(i2))
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
                   !|  dPhi/dx   dPhi/dy   dPhi/dz| = |-sin(Phi)/sin(Theta/R  cos(Phi)/sin(Theta)/R             0|
                   !|    dR/dx     dR/dy     dR/dz|   |   sin(Theta)cos(Phi)     sin(Theta)sin(Phi)    cos(Theta)|
                   !|dTheta/dx dTheta/dy dTheta/dz|   | cos(Theta)cos(Phi)/R   cos(Theta)sin(Phi)/R -sin(Theta)/R|
                   IF(yy == 1) THEN
                   
                   Adummy(1,1) =-sin(X1(i1))/sin(X3(i3))/X2(i2) 
                   Adummy(1,2) = cos(X1(i1))/sin(X3(i3))/X2(i2) 
                   Adummy(1,3) = 0d0                 

                   Adummy(2,1) = sin(X3(i3))*cos(X1(i1))
                   Adummy(2,2) = sin(X3(i3))*sin(X1(i1))
                   Adummy(2,3) = cos(X3(i3))

                   Adummy(3,1) = cos(X3(i3))*cos(X1(i1))/X2(i2) 
                   Adummy(3,2) = cos(X3(i3))*sin(X1(i1))/X2(i2) 
                   Adummy(3,3) = -sin(X3(i3))/X2(i2) 

                   !Yang grid
                   ELSE
                
                   CALL yin2yang(X1(i1),X3(i3),lr,cr)
                   rr = X2(i2)
                   !x2z2 = (sin(cr)*cos(lr))**2.0 + cos(cr)**2.0 ! (x^2 + z^2)/r^2
                   xx = rr*sin(cr)*cos(lr)
                   yyy= rr*sin(cr)*sin(lr)
                   zz = rr*cos(cr)
                   x2z2 = (xx**2.0 + zz**2.0)**0.5
                   Adummy(1,1) =-xx*zz/(ABS(xx)*x2z2**2.0)
                   Adummy(1,2) = 0d0
                   Adummy(1,3) = ABS(xx)/x2z2**2.0

                   Adummy(2,1) = sin(cr)*cos(lr)
                   Adummy(2,2) = sin(cr)*sin(lr)
                   Adummy(2,3) = cos(cr)

                   Adummy(3,1) = xx*yyy/((rr)**2.0*x2z2)
                   Adummy(3,2) =-x2z2/(rr)**2.0                    
                   Adummy(3,3) = yyy*zz/((rr)**2.0*x2z2)

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

   END IF
!!! END 3D model

   END SUBROUTINE load

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

   xx = (X1s-X1n(i1s))/(X1n(i1s+1)-X1n(i1s))
   yy = (X2s-X2n(i2s))/(X2n(i2s+1)-X2n(i2s))

   res = (1d0-xx)*(1d0-yy)*y1 + &
               (1d0-xx)*yy*y2 + &
               xx*(1d0-yy)*y3 + &
                     xx*yy*y4

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

   DOUBLE PRECISION :: X1s,X2s,x1c,x2c,x1cc,x2ccc,shiftX,shiftY
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

   if(shiftX .ne. 0) then
      if(x1periodic .ne. 0 .and. i1s== nx1-1) then
         xx=(X1s-shiftX)/((X1n(i1s+1)-X1n(i1s)+X1n(2)-X1n(1))/2)
      else
         xx=(X1s-shiftX)/((X1n(i1s+2)-X1n(i1s))/2)
      end if
   end if

   if(shiftY .ne. 0) then
      if(x2periodic .ne. 0 .and. i2s== nx2-1) then
         yy=(X2s-shiftY)/((X2n(i2s+1)-X2n(i2s)+X2n(2)-X2n(1))/2)
      else
         yy=(X2s-shiftY)/((X2n(i2s+2)-X2n(i2s))/2)
      end if
   end if
   
   if(shiftX .eq. 0) xx=(X1s-X1n(i1s))/(X1n(i1s+1)-X1n(i1s))
   if(shiftY .eq. 0) yy=(X2s-X2n(i2s))/(X2n(i2s+1)-X2n(i2s))

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

   INTEGER :: i1s,i2s,i1,i2,m
   ! indices of the UPPER-LEFT grid point closest to the considered point

   DOUBLE PRECISION :: y1,y2,y3,y4
   ! dummies for interpolation

   DOUBLE PRECISION, DIMENSION(3,3) :: ee
   ! dummy for reference strain rate calculation

   INTEGER :: nrot,tid
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
!!! subroutine UPPERLEFT2D, calculation of upper left node                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE upperleft2Dn(x10,x20,i1,i2)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: x10,x20

   INTEGER :: i1,i2

   i1 = 1; i2 = 1

   DO WHILE (X1n(i1+1) .LT. x10 .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

   DO WHILE (X2n(i2+1) .LT. x20 .AND. i2 .LT. nx2-1) ; i2 = i2+1 ; ENDDO
      
   RETURN

   END SUBROUTINE upperleft2Dn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine UPPERLEFT, calculation of upper left node                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE upperleft3Dn(x10,x20,x30,i1,i2,i3)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   DOUBLE PRECISION :: x10,x20,x30

   INTEGER :: i1,i2,i3 

   i1 = 1; i2 = 1; i3 = 1

   DO WHILE (X1n(i1+1) .LT. x10 .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

   DO WHILE (X2n(i2+1) .LT. x20 .AND. i2 .LT. nx2-1) ; i2 = i2+1 ; ENDDO
      
   DO WHILE (X3n(i3+1) .LT. x30 .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO
                 
   RETURN

   END SUBROUTINE upperleft3Dn


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

   xx = (X1s-X1n(i1s))/(X1n(i1s+1)-X1n(i1s))
   yy = (X2s-X2n(i2s))/(X2n(i2s+1)-X2n(i2s))
   zz = (X3s-X3n(i3s))/(X3n(i3s+1)-X3n(i3s))

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

   DOUBLE PRECISION :: X1s,X2s,X3s,x1c,x2c,x3c,x1cc,x2cc,x3cc,shiftX,shiftY,shiftZ
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

   if(shiftX .ne. 0) then
      if(x1periodic .ne. 0 .and. i1s== nx1-1) then
         xx=(X1s-shiftX)/((X1n(i1s+1)-X1n(i1s)+X1n(2)-X1n(1))/2)
      else
         xx=(X1s-shiftX)/((X1n(i1s+2)-X1n(i1s))/2)
      end if
   end if
   
   if(shiftY .ne. 0) then
      if(x2periodic .ne. 0 .and. i2s== nx2-1) then
         yy=(X2s-shiftY)/((X2n(i2s+1)-X2n(i2s)+X2n(2)-X2n(1))/2)
      else
         yy=(X2s-shiftY)/((X2n(i2s+2)-X2n(i2s))/2)
      end if
   end if
   
   if(shiftZ .ne. 0) then
      if(x3periodic .ne. 0 .and. i3s== nx3-1) then
         zz=(X3s-shiftZ)/((X3n(i3s+1)-X3n(i3s)+X3n(2)-X3n(1))/2)
      else
         zz=(X3s-shiftZ)/((X3n(i3s+2)-X3n(i3s))/2)
      end if
   end if
   
   if(shiftX .eq. 0) xx=(X1s-X1n(i1s))/(X1n(i1s+1)-X1n(i1s))
   if(shiftY .eq. 0) yy=(X2s-X2n(i2s))/(X2n(i2s+1)-X2n(i2s))
   if(shiftZ .eq. 0) zz=(X3s-X3n(i3s))/(X3n(i3s+1)-X3n(i3s))

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

   INTEGER :: nrot,tid
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

   CALL interpshift(X1dum,X2dum,X3dum,i1s-1,i2s-1,i3s-1,(X1(i1s)+X1(i1s-1))/2d0,(X2(i2s)+X2(i2s-1))/2d0,(X3(i3s)+X3(i3s-1))/2d0,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,1,1))

   !dVy/dy at -stepx1/2,-stepx2/2,-stepx3/2

   y1 = Dij(yy,2,2,ix(1),iy(1),iz(1)) ; y2 = Dij(yy,2,2,ix(2),iy(1),iz(1))
   y3 = Dij(yy,2,2,ix(2),iy(1),iz(2)) ; y4 = Dij(yy,2,2,ix(1),iy(1),iz(2))
   y5 = Dij(yy,2,2,ix(1),iy(2),iz(1)) ; y6 = Dij(yy,2,2,ix(2),iy(2),iz(1))
   y7 = Dij(yy,2,2,ix(2),iy(2),iz(2)) ; y8 = Dij(yy,2,2,ix(1),iy(2),iz(2))

   CALL interpshift(X1dum,X2dum,X3dum,i1s-1,i2s-1,i3s-1,(X1(i1s)+X1(i1s-1))/2d0,(X2(i2s)+X2(i2s-1))/2d0,(X3(i3s)+X3(i3s-1))/2d0,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,2,2))

   !dVz/dz
   y1 = Dij(yy,3,3,ix(1),iy(1),iz(1)) ; y2 = Dij(yy,3,3,ix(2),iy(1),iz(1))
   y3 = Dij(yy,3,3,ix(2),iy(1),iz(2)) ; y4 = Dij(yy,3,3,ix(1),iy(1),iz(2))
   y5 = Dij(yy,3,3,ix(1),iy(2),iz(1)) ; y6 = Dij(yy,3,3,ix(2),iy(2),iz(1))
   y7 = Dij(yy,3,3,ix(2),iy(2),iz(2)) ; y8 = Dij(yy,3,3,ix(1),iy(2),iz(2))

   CALL interpshift(X1dum,X2dum,X3dum,i1s-1,i2s-1,i3s-1,(X1(i1s)+X1(i1s-1))/2d0,(X2(i2s)+X2(i2s-1))/2d0,(X3(i3s)+X3(i3s-1))/2d0,y1,y2,y3,y4,y5,y6,y7,y8,l(tid,3,3))

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
