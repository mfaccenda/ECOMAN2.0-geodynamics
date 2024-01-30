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
      IF(ln_fse(m)>=ln_fse_min) THEN 
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
      IF(ln_fse(m)>=ln_fse_min) THEN
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

