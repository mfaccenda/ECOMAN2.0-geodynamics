
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

!-----------------------------------------------------------------------------------------
!  MODULE: DistComm
!
!> @brief Distribution Communications management
!
!> @author Giovanni Isotton
!
!> @version 1.0
!
!> @date june 2021
!
!> @par License:
!! This program is intended for internal research only and can not be distributed
!! elsewhere without authors' consent.
!-----------------------------------------------------------------------------------------

module class_DistComm

use class_precision

implicit none

include 'mpif.h'

! private variables
integer, parameter, private :: tagList  = 0 
integer, parameter, private :: tagValue = 1 

! Public variables
type, public :: t_DistComm

   ! Send list
   integer, pointer :: SendList(:)   => null()
   integer, pointer :: SendOffset(:) => null()

   ! Send communications
   integer          :: SendNcomm     = 0
   integer, pointer :: SendComm(:)   => null()
   integer, pointer :: SendStatus(:) => null()

   ! Recv list
   integer, pointer :: RecvList(:)   => null()
   integer, pointer :: RecvOffset(:) => null()

   ! Recv communications
   integer          :: RecvNcomm     = 0
   integer, pointer :: RecvComm(:)   => null()
   integer, pointer :: RecvStatus(:) => null()

end type

! Public member function/subroutines
public :: StartDistList,StartDistValueR,StartDistValueI,dlt_DistComm
public :: WaitSendComm,WaitRecvComm

!-----------------------------------------------------------------------------------------

contains

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: StartDistList
   !
   !> @brief Start to distributs Lists
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------

   subroutine StartDistList(nbuff,DistComm_var)

   implicit none

   ! Input variables
   integer, intent(in) :: nbuff

   ! Input/output variables
   type(t_DistComm), intent(inout) :: DistComm_var

   ! Local variables
   integer          :: i,j,istart,ierr
   integer          :: SendNcomm,SendBuff
   integer, pointer :: SendList(:)   => null()
   integer, pointer :: SendOffset(:) => null()
   integer, pointer :: SendComm(:)   => null()
   integer, pointer :: SendStatus(:) => null()
   integer          :: RecvNcomm,RecvBuff
   integer, pointer :: RecvList(:)   => null()
   integer, pointer :: RecvOffset(:) => null()
   integer, pointer :: RecvComm(:)   => null()
   integer, pointer :: RecvStatus(:) => null()

   ! Set handles
   SendNcomm  =  DistComm_var%SendNcomm
   SendList   => DistComm_var%SendList
   SendOffset => DistComm_var%SendOffset
   SendComm   => DistComm_var%SendComm
   SendStatus => DistComm_var%SendStatus
   RecvNcomm  =  DistComm_var%RecvNcomm
   RecvList   => DistComm_var%RecvList
   RecvOffset => DistComm_var%RecvOffset
   RecvComm   => DistComm_var%RecvComm
   RecvStatus => DistComm_var%RecvStatus

   ! Start Send Communications
   istart= 1
   do i = 1,SendNcomm
      j = SendComm(i)
      SendBuff = nbuff * ( SendOffset(i+1) - SendOffset(i) )
      call MPI_Isend(SendList(istart),SendBuff,MPI_INTEGER,j-1,tagList,MPI_COMM_WORLD,    &
      &              SendStatus(i),ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in StartDistList : MPI_Isend"
         stop
      endif
      istart= istart + SendBuff
   end do

   ! Start Recv Communications
   istart= 1
   do i = 1,RecvNcomm
      j = RecvComm(i)
      RecvBuff = nbuff * ( RecvOffset(i+1) - RecvOffset(i) )
      call MPI_Irecv(RecvList(istart),RecvBuff,MPI_INTEGER,j-1,tagList,MPI_COMM_WORLD,    &
      &              RecvStatus(i),ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in StartDistList : MPI_Irecv"
         stop
      endif
      istart = istart + RecvBuff
   end do

   end subroutine StartDistList

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: StartDistValueR
   !
   !> @brief Start to distributs real Values
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------

   subroutine StartDistValueR(nbuff,SendN,SendX,RecvN,RecvX,DistComm_var)

   implicit none

   ! Input variables
   integer,           intent(in) :: nbuff
   integer,           intent(in) :: SendN
   real(kind=double), intent(in) :: SendX(nbuff*SendN)
   integer,           intent(in) :: RecvN
   real(kind=double), intent(in) :: RecvX(nbuff*RecvN)

   ! Input/output variables
   type(t_DistComm), intent(inout) :: DistComm_var

   ! Local variables
   integer          :: i,j,istart,ierr
   integer          :: SendNcomm,SendBuff
   integer, pointer :: SendOffset(:) => null()
   integer, pointer :: SendComm(:)   => null()
   integer, pointer :: SendStatus(:) => null()
   integer          :: RecvNcomm,RecvBuff
   integer, pointer :: RecvOffset(:) => null()
   integer, pointer :: RecvComm(:)   => null()
   integer, pointer :: RecvStatus(:) => null()

   ! Set handles
   SendNcomm  =  DistComm_var%SendNcomm
   SendOffset => DistComm_var%SendOffset
   SendComm   => DistComm_var%SendComm
   SendStatus => DistComm_var%SendStatus
   RecvNcomm  =  DistComm_var%RecvNcomm
   RecvOffset => DistComm_var%RecvOffset
   RecvComm   => DistComm_var%RecvComm
   RecvStatus => DistComm_var%RecvStatus

   ! Start Send Communications
   istart= 1
   do i = 1,SendNcomm
      j = SendComm(i)
      SendBuff = nbuff * ( SendOffset(i+1) - SendOffset(i) )
      call MPI_Isend(SendX(istart),SendBuff,MPI_DOUBLE,j-1,tagValue,MPI_COMM_WORLD,       &
      &              SendStatus(i),ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in StartDistValueR : MPI_Isend"
         stop
      endif
      istart= istart + SendBuff
   end do

   ! Start Recv Communications
   istart= 1
   do i = 1,RecvNcomm
      j = RecvComm(i)
      RecvBuff = nbuff * ( RecvOffset(i+1) - RecvOffset(i) )
      call MPI_Irecv(RecvX(istart),RecvBuff,MPI_DOUBLE,j-1,tagValue,MPI_COMM_WORLD,       &
      &              RecvStatus(i),ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in StartDistValueR : MPI_Irecv"
         stop
      endif
      istart = istart + RecvBuff
   end do

   end subroutine StartDistValueR

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: StartDistValueI
   !
   !> @brief Start to distributs integer Values
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------

   subroutine StartDistValueI(nbuff,SendN,SendX,RecvN,RecvX,DistComm_var)

   implicit none

   ! Input variables
   integer, intent(in) :: nbuff
   integer, intent(in) :: SendN
   integer, intent(in) :: SendX(nbuff*SendN)
   integer, intent(in) :: RecvN
   integer, intent(in) :: RecvX(nbuff*RecvN)

   ! Input/output variables
   type(t_DistComm), intent(inout) :: DistComm_var

   ! Local variables
   integer          :: i,j,istart,ierr
   integer          :: SendNcomm,SendBuff
   integer, pointer :: SendOffset(:) => null()
   integer, pointer :: SendComm(:)   => null()
   integer, pointer :: SendStatus(:) => null()
   integer          :: RecvNcomm,RecvBuff
   integer, pointer :: RecvOffset(:) => null()
   integer, pointer :: RecvComm(:)   => null()
   integer, pointer :: RecvStatus(:) => null()

   ! Set handles
   SendNcomm  =  DistComm_var%SendNcomm
   SendOffset => DistComm_var%SendOffset
   SendComm   => DistComm_var%SendComm
   SendStatus => DistComm_var%SendStatus
   RecvNcomm  =  DistComm_var%RecvNcomm
   RecvOffset => DistComm_var%RecvOffset
   RecvComm   => DistComm_var%RecvComm
   RecvStatus => DistComm_var%RecvStatus

   ! Start Send Communications
   istart= 1
   do i = 1,SendNcomm
      j = SendComm(i)
      SendBuff = nbuff * ( SendOffset(i+1) - SendOffset(i) )
      call MPI_Isend(SendX(istart),SendBuff,MPI_INTEGER,j-1,tagValue,MPI_COMM_WORLD,       &
      &              SendStatus(i),ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in StartDistValueI : MPI_Isend"
         stop
      endif
      istart= istart + SendBuff
   end do

   ! Start Recv Communications
   istart= 1
   do i = 1,RecvNcomm
      j = RecvComm(i)
      RecvBuff = nbuff * ( RecvOffset(i+1) - RecvOffset(i) )
      call MPI_Irecv(RecvX(istart),RecvBuff,MPI_INTEGER,j-1,tagValue,MPI_COMM_WORLD,       &
      &              RecvStatus(i),ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in StartDistValueI : MPI_Irecv"
         stop
      endif
      istart = istart + RecvBuff
   end do

   end subroutine StartDistValueI

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: WaitSendComm
   !
   !> @brief Wait end of sending communications
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------

   subroutine WaitSendComm(DistComm_var)

   implicit none

   ! Input/output variables
   type(t_DistComm), intent(inout) :: DistComm_var

   ! Local variables
   integer          :: i,ierr
   integer          :: Ncomm
   integer, pointer :: Comm(:) => null()
   integer          :: SendNcomm
   integer, pointer :: SendStatus(:) => null()
   logical          :: flag
   integer          :: statusMPI(MPI_STATUS_SIZE)

   ! Set handles
   SendNcomm  =  DistComm_var%SendNcomm
   SendStatus => DistComm_var%SendStatus

   if ( SendNcomm .eq. 0 ) return

   ! initialize
   Ncomm = SendNcomm
   allocate(Comm(Ncomm),stat=ierr)
   if( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** ERROR in WaitSendComm : allocating"
      stop
   endif
   do i = 1,Ncomm
      Comm(i) = i
   end do

   ! loop and test
   i = 1
   do while (Ncomm .gt. 0)
      call MPI_Test(SendStatus(Comm(i)),flag,statusMPI,ierr)
      if ( flag ) then
         ! sended, remove from list
         Comm(i) = Comm(Ncomm)
         Ncomm = Ncomm - 1
      else
         ! not sended yet, go on
         i = i + 1
      endif
      if ( i .gt. Ncomm ) i = 1
   end do

   ! deallocate
   deallocate(Comm,stat=ierr)
   if( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** ERROR in WaitSendComm : deallocating"
      stop
   endif

   end subroutine WaitSendComm

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: WaitRecvComm
   !
   !> @brief Wait end of receiving communications
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------

   subroutine WaitRecvComm(DistComm_var)

   implicit none

   ! Input/output variables
   type(t_DistComm), intent(inout) :: DistComm_var

   ! Local variables
   integer          :: i,ierr
   integer          :: Ncomm
   integer, pointer :: Comm(:) => null()
   integer          :: RecvNcomm
   integer, pointer :: RecvStatus(:) => null()
   logical          :: flag
   integer          :: statusMPI(MPI_STATUS_SIZE)

   ! Set handles
   RecvNcomm  =  DistComm_var%RecvNcomm
   RecvStatus => DistComm_var%RecvStatus

   if ( RecvNcomm .eq. 0 ) return

   ! initialize
   Ncomm = RecvNcomm
   allocate(Comm(Ncomm),stat=ierr)
   if( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** ERROR in WaitRecvComm : allocating"
      stop
   endif
   do i = 1,Ncomm
      Comm(i) = i
   end do

   ! loop and test
   i = 1
   do while (Ncomm .gt. 0)
      call MPI_Test(RecvStatus(Comm(i)),flag,statusMPI,ierr)
      if ( flag ) then
         ! sended, remove from list
         Comm(i) = Comm(Ncomm)
         Ncomm = Ncomm - 1
      else
         ! not sended yet, go on
         i = i + 1
      endif
      if ( i .gt. Ncomm ) i = 1
   end do

   ! deallocate
   deallocate(Comm,stat=ierr)
   if( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** ERROR in WaitRecvComm : deallocating"
      stop
   endif

   end subroutine WaitRecvComm

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: dlt_DistComm
   !
   !> @brief Release data structure
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------

   subroutine dlt_DistComm(DistComm_var)

   implicit none

   ! Input/output variables
   type(t_DistComm), intent(inout) :: DistComm_var

   ! Local variables
   integer :: ierr

   ! Release memory
   if ( DistComm_var%SendNcomm .gt. 0 ) then
      deallocate(DistComm_var%SendList,DistComm_var%SendOffset,DistComm_var%SendComm,     &
      &          DistComm_var%SendStatus,stat=ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in dlt_DistComm : deallocating Send"
         stop
      endif
      DistComm_var%SendList   => null()
      DistComm_var%SendOffset => null()
      DistComm_var%SendComm   => null()
      DistComm_var%SendStatus => null()
   endif
   DistComm_var%SendNcomm = 0

   if ( DistComm_var%RecvNcomm .gt. 0 ) then
      deallocate(DistComm_var%RecvList,DistComm_var%RecvOffset,DistComm_var%RecvComm,     &
      &          DistComm_var%RecvStatus,stat=ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in dlt_DistComm : deallocating Recv"
         stop
      endif
      DistComm_var%RecvList   => null()
      DistComm_var%RecvOffset => null()
      DistComm_var%RecvComm   => null()
      DistComm_var%RecvStatus => null()
   endif
   DistComm_var%RecvNcomm = 0

   end subroutine dlt_DistComm

!-----------------------------------------------------------------------------------------

end module class_DistComm
