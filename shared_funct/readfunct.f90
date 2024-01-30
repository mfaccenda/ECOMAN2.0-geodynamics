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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE read_par_char(unitt,str)

   IMPLICIT NONE

   INTEGER :: i,ier,unitt
   CHARACTER(500) :: str

10 READ(unitt,"(a)",advance='yes',IOSTAT=ier) str
   i=index(str,'#')
   if(i==1 .OR. LEN_TRIM(str)==0) then
      goto 10
   end if

   END SUBROUTINE read_par_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE read_par_int(unitt,intg)

   IMPLICIT NONE

   INTEGER ::i,ier,intg,unitt
   CHARACTER(500) :: str

10 READ(unitt,"(a)",advance='yes',IOSTAT=ier) str  
   i=index(str,'#')
   if(i==1 .OR. LEN_TRIM(str)==0) then
      goto 10
   else
      read(str,*) intg
   end if

   END SUBROUTINE read_par_int 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE read_par_double(unitt,dbl)

   IMPLICIT NONE

   INTEGER ::i,ier,unitt
   DOUBLE PRECISION :: dbl
   CHARACTER(500) :: str

10 READ(unitt,"(a)",advance='yes',IOSTAT=ier) str  
   i=index(str,'#')
   if(i==1 .OR. LEN_TRIM(str)==0) then
      goto 10
   else
      read(str,*) dbl
   end if

   END SUBROUTINE read_par_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE read_par_double6(unitt,dbl)

   IMPLICIT NONE

   INTEGER ::i,ier,unitt
   DOUBLE PRECISION :: dbl(6)
   CHARACTER(500) :: str

10 READ(unitt,"(a)",advance='yes',IOSTAT=ier) str  
   i=index(str,'#')
   if(i==1 .OR. LEN_TRIM(str)==0) then
      goto 10
   else
      read(str,*) dbl(:)
   end if

   END SUBROUTINE read_par_double6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

