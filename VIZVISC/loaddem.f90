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

SUBROUTINE loadtensordem

USE comvarvis 
USE hdf5

IMPLICIT NONE

INTEGER :: nx(4),error
INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id,dcpl,memtype !Handles

!Initialize FORTRAN interface.

CALL H5open_f (error)

! Open a new file using default properties.

CALL H5Fopen_f('../EXEV/DEMviscous/viscoustensordem.h5', H5F_ACC_RDONLY_F, file_id, error)

CALL loadsave_integer(0,1,file_id,4,H5T_NATIVE_INTEGER,nx(1:4),'gridnum',0)

etacontrnum = nx(1)
volnum     = nx(2)
esanum(1)  = nx(3)
esanum(2)  = nx(4)

ALLOCATE(C11(etacontrnum,volnum,esanum(1),esanum(2)))
ALLOCATE(C22(etacontrnum,volnum,esanum(1),esanum(2)))
ALLOCATE(C33(etacontrnum,volnum,esanum(1),esanum(2)))
ALLOCATE(C44(etacontrnum,volnum,esanum(1),esanum(2)))
ALLOCATE(C55(etacontrnum,volnum,esanum(1),esanum(2)))
ALLOCATE(C66(etacontrnum,volnum,esanum(1),esanum(2)))
ALLOCATE(C12(etacontrnum,volnum,esanum(1),esanum(2)))
ALLOCATE(C13(etacontrnum,volnum,esanum(1),esanum(2)))
ALLOCATE(C23(etacontrnum,volnum,esanum(1),esanum(2)))
ALLOCATE(etacontr(etacontrnum))
ALLOCATE(Vol(volnum))
ALLOCATE(esa12(esanum(1),esanum(2)))
ALLOCATE(esa23(esanum(1),esanum(2)))
ALLOCATE(Cij(etacontrnum,volnum,esanum(1),esanum(2),9))

CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C11,'C11',0)
CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C22,'C22',0)
CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C33,'C33',0)
CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C44,'C44',0)
CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C55,'C55',0)
CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C66,'C66',0)
CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C12,'C12',0)
CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C13,'C13',0)
CALL loadsave_double(0,4,file_id,nx(1:4),H5T_NATIVE_DOUBLE,C23,'C23',0)
CALL loadsave_double(0,1,file_id,nx(1),H5T_NATIVE_DOUBLE,etacontr,'nucontr',0)
!CALL loadsave_double(0,1,file_id,nx(1),H5T_NATIVE_DOUBLE,etacontr,'etacontr',0)
CALL loadsave_double(0,1,file_id,nx(2),H5T_NATIVE_DOUBLE,Vol,'Vol',0)
CALL loadsave_double(0,2,file_id,nx(3:4),H5T_NATIVE_DOUBLE,esa12,'esa12',0)
CALL loadsave_double(0,2,file_id,nx(3:4),H5T_NATIVE_DOUBLE,esa23,'esa23',0)

Cij(:,:,:,:,1) = C11
Cij(:,:,:,:,2) = C22
Cij(:,:,:,:,3) = C33
Cij(:,:,:,:,4) = C44
Cij(:,:,:,:,5) = C55
Cij(:,:,:,:,6) = C66
Cij(:,:,:,:,7) = C12
Cij(:,:,:,:,8) = C13
Cij(:,:,:,:,9) = C23

!Terminate access to the file.
CALL H5Fclose_f(file_id, error)

!Close FORTRAN interface.
CALL H5close_f(error)

RETURN

END SUBROUTINE loadtensordem
