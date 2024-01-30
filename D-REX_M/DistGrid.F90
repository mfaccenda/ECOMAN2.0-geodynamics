
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
!  MODULE: DistGrid
!
!> @brief Structured grid distibution
!
!> @details
!! This module is used to manage structured grid among MPI ranks
!
!> @author Giovanni Isotton
!
!> @version 1.0
!
!> @date June 2021
!
!> @par License:
!! This program is intended for internal research only and can not be distributed
!! elsewhere without authors' consent.
!-----------------------------------------------------------------------------------------
module class_DistGrid

!    Distribution of MPI ranks
!
!    Cartesian grid ----------------------------------------------------------------------
!
!    Y
!    |                       Z = 1
!    ---------------------
!    | 3 | 6 | 9 | 12| 15|
!    ---------------------
!    | 2 | 5 | 8 | 11| 14|
!    ---------------------
!    | 1 | 4 | 7 | 10| 13|
!    ------------------------X
!
!    Y
!    |                       Z = 2
!    ---------------------
!    | 18| 21| 24| 27| 30|
!    ---------------------
!    | 17| 20| 23| 26| 29|
!    ---------------------
!    | 16| 19| 22| 25| 28|
!    ------------------------X
!
!    iproc1 --> X
!    iproc2 --> Y
!    iproc3 --> Z
!
!    rank 28 : iproc1 = 5; iproc2 = 1; iproc3 = 2
!
!
!    Spheric   grid ----------------------------------------------------------------------
!
!
!                      COLAT = 1
!    ---------\
!    |        /\
!    |  6    /  \
!    -------\    \  LONG
!    |     / \    \
!    |  5 /   \ 3  \
!    -----\  2 \    \
!    | 4/ 1\    \    \
!     ------------------RADIUS
!
!
!                       COLAT = 2
!    ---------\
!    |        /\
!    | 12    /  \
!    -------\    \  LONG
!    |     / \    \
!    | 11 /   \ 9  \
!    -----\  8 \    \
!    |10/ 7\    \    \
!     ------------------RADIUS
!
!    iproc1 --> LONG
!    iproc2 --> RADIUS
!    iproc3 --> COLAT
!
!    rank  6 : iproc1 = 2; irpoc2 = 3; iproc3 = 1
!

! Modules
use class_precision
use class_DistComm

implicit none

!!!!!!!!!!!!!!
!#define DEBUG
!!!!!!!!!!!!!!

! Local private parameters
real(kind=double), parameter, private :: pi = 3.141592653589793238462643383279d0

! Public function/subroutines
public :: DistGrid,RmDistGrid,DistMrkValR,DistMrkValI
public :: GetCoordSphe,GetCoordCart

! Private function/subroutines
private :: DistInfo,SetUnk
private :: DistGridCart,DistGridSphe
private :: DistMrkInd
private :: MarkLoc2GloCart,MarkLoc2GloSphe
private :: iglo2ix
private :: resizei,binsearch

! Private types

! - type for storing MPI info
type, private :: t_infoMPI
   integer :: rankMPI
   integer :: sizeMPI
   integer :: nproc1
   integer :: nproc2
   integer :: nproc3
end type

! - type for storing Computational domain
type, private :: t_CompDom
   ! - Eulerian Grid / Axis 1 (X-cart or Long)
   real(kind=double) :: x1min
   real(kind=double) :: x1max
   integer           :: nx1
   integer           :: x1periodic
   ! - Eulerian Grid / Axis 2 (Y-cart or Radial)
   real(kind=double) :: x2min
   real(kind=double) :: x2max
   integer           :: nx2
   integer           :: x2periodic
   ! - Eulerian Grid / Axis 3 (Z-cart or Colat)
   real(kind=double) :: x3min
   real(kind=double) :: x3max
   integer           :: nx3
   integer           :: x3periodic
   ! - Lagrangian Grid / Axis 1 (X-cart or Long)
   real(kind=double) :: mx1min
   real(kind=double) :: mx1max
   real(kind=double) :: mx1stp
   ! - Lagrangian Grid / Axis 2 (Y-cart or Radial)
   real(kind=double) :: mx2min
   real(kind=double) :: mx2max
   real(kind=double) :: mx2stp
   ! - Lagrangian Grid / Axis 3 (Z-cart or Colat)
   real(kind=double) :: mx3min
   real(kind=double) :: mx3max
   real(kind=double) :: mx3stp
end type

! Private variables

! - MPI info
type(t_infoMPI), private :: infoMPI

! - Computational domain
type(t_CompDom), public :: CompDom

! - Rank computational domain
type(t_CompDom), private :: RankCompDom

! - Pointers to the first global nodes
integer, allocatable, private :: ptx1(:)
integer, allocatable, private :: ptx2(:)
integer, allocatable, private :: ptx3(:)

! - Communication management
type(t_DistComm), private :: DistMrk

! Lagrangian grid
! - number of marks to print (owned + received from communications)
integer, private :: nnmrk
! - Pointers to the first global marknum (to print output)
integer, allocatable, public :: ptnmrk(:)
! - Marks to print in global reference (unsorted: as received from lists communication)
integer, allocatable, private :: indglomrk(:)

! Local Rank Lagrangian grid
! - number of local marks
integer, public :: Ranknnmrk
! - Local marks in global reference (local rank index to global index)
integer, allocatable, private :: mrkglo(:)

!-----------------------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: DistGrid
   !
   !> @brief Distributes the structured grid
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------

   subroutine DistGrid(nproc1,nproc2,nproc3,                                              &
   &                   dimensions,cartspher,                                              &
   &                   x1min,x1max,nx1,x1periodic,                                        &
   &                   x2min,x2max,nx2,x2periodic,                                        &
   &                   x3min,x3max,nx3,x3periodic,                                        &
   &                   mx1min,mx1max,mx1stp,                                              &
   &                   mx2min,mx2max,mx2stp,                                              &
   &                   mx3min,mx3max,mx3stp)

   implicit none

   ! Input variables
   integer, intent(in) :: nproc1,nproc2,nproc3 
   integer, intent(in) :: dimensions,cartspher

   ! Input/Output variables
   real(kind=double), intent(inout) :: x1min,x1max
   integer,           intent(inout) :: nx1,x1periodic
   real(kind=double), intent(inout) :: x2min,x2max
   integer,           intent(inout) :: nx2,x2periodic
   real(kind=double), intent(inout) :: x3min,x3max
   integer,           intent(inout) :: nx3,x3periodic
   real(kind=double), intent(inout) :: mx1min,mx1max,mx1stp 
   real(kind=double), intent(inout) :: mx2min,mx2max,mx2stp 
   real(kind=double), intent(inout) :: mx3min,mx3max,mx3stp 

   ! Local variables
   integer           :: rankMPI,sizeMPI,errMPI
   integer           :: i1,i2,i3
   integer           :: allmnx1,allmnx2,allmnx3
   integer           :: allmrk
   real(kind=double) :: mlongmin,mlongmax,mlongstp
   real(kind=double) :: mradiusmin,mradiusmax,mradiusstp
   real(kind=double) :: mcolatmin,mcolatmax,mcolatstp
   real(kind=double) :: longlen
   real(kind=double) :: colatlen,mcolatradiansstp,mcolat
   real(kind=double) :: mR
   integer           :: mrnum,mlongnum,mcolatnum
 
   ! Store MPI info
   call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeMPI,errMPI)
   call MPI_COMM_RANK(MPI_COMM_WORLD,rankMPI,errMPI)
   rankMPI = rankMPI + 1
   infoMPI%sizeMPI = sizeMPI 
   infoMPI%rankMPI = rankMPI 
   infoMPI%nproc1  = nproc1
   infoMPI%nproc2  = nproc2
   infoMPI%nproc3  = nproc3

   ! Check number of MPI ranks
   if ( sizeMPI .ne. nproc1 * nproc2 * nproc3 ) then
      write(*,'(a,i5)') "*** ERROR in DistGrid : sizeMPI /= nproc1 * nproc2 * nproc3",rankMPI
      stop
   endif

   ! Check dimension
   if ( dimensions .eq. 2 .and. nx3 .gt. 1 ) then
      write(*,'(a,i5)') "*** ERROR in DistGrid : dimensions == 2 and nx3 > 1",rankMPI
      stop
   endif
   if ( dimensions .eq. 3 .and. nx3 .eq. 1 ) then
      write(*,'(a,i5)') "*** ERROR in DistGrid : dimensions == 3 and nx3 == 1",rankMPI
      stop
   endif

   ! Store computational domain
   CompDom%x1min      = x1min
   CompDom%x1max      = x1max
   CompDom%nx1        = nx1
   CompDom%x1periodic = x1periodic
   CompDom%x2min      = x2min
   CompDom%x2max      = x2max
   CompDom%nx2        = nx2
   CompDom%x2periodic = x2periodic
   CompDom%x3min      = x3min
   CompDom%x3max      = x3max
   CompDom%nx3        = nx3
   CompDom%x3periodic = x3periodic
   CompDom%mx1min     = mx1min
   CompDom%mx1max     = mx1max
   CompDom%mx1stp     = mx1stp
   CompDom%mx2min     = mx2min
   CompDom%mx2max     = mx2max
   CompDom%mx2stp     = mx2stp
   CompDom%mx3min     = mx3min
   CompDom%mx3max     = mx3max
   CompDom%mx3stp     = mx3stp

   if ( sizeMPI .gt. 1 ) then

      if ( cartspher .eq. 1 ) then
         ! Cartesian
         call DistGridCart(dimensions)
      else
         ! Spherical
         call DistGridSphe(dimensions)
      endif

   else

      if ( cartspher .eq. 1 ) then

         ! Cartesian
         allmnx1 = nint((CompDom%mx1max-CompDom%mx1min)/CompDom%mx1stp)
         allmnx2 = nint((CompDom%mx2max-CompDom%mx2min)/CompDom%mx2stp)
         allmnx3 = 1
         if ( dimensions .eq. 3 ) allmnx3 = nint((CompDom%mx3max-CompDom%mx3min)/CompDom%mx3stp)
         allmrk = allmnx1 * allmnx2 * allmnx3

      else

         ! Spherical

         mlongmin   = CompDom%mx1min ! Radians 
         mlongmax   = CompDom%mx1max ! Radians
         mradiusmin = CompDom%mx2min ! Distance
         mradiusmax = CompDom%mx2max ! Distance
         mcolatmin  = CompDom%mx3min ! Radians 
         mcolatmax  = CompDom%mx3max ! Radians
         mlongstp   = CompDom%mx1stp ! Distance
         mradiusstp = CompDom%mx2stp ! Distance 
         mcolatstp  = CompDom%mx3stp ! Distance

         mrnum = floor( (mradiusmax - mradiusmin) / mradiusstp )

         ! Get number of marks for each radius
         allmrk = 0
         ! Loop over radius
         do i2 = 1,mrnum
            ! Find radius
            mR = mradiusmin + mradiusstp * ( i2 - 0.5d0 )
            ! Colatitude step in radians
            mcolatnum = 1
            mcolatradiansstp = 0.d0 
            ! Check dimensions to compute number of aggregates along colatitude arc
            if ( dimensions .eq. 3 ) then
               ! Find lenght of colatitude arc (great circle arc)
               colatlen = ( mcolatmax - mcolatmin ) * mR
               ! Number of aggregates along colatitude arc
               mcolatnum = floor ( colatlen / mcolatstp)
               ! Colatitude step in radians
               mcolatradiansstp = ( mcolatmax - mcolatmin ) / real(mcolatnum)
            endif
            ! Loop over colatitude
            do i3 = 1,mcolatnum
               ! Colatitude in radians
               if ( dimensions .eq. 3 ) then
                  mcolat = mcolatmin + mcolatradiansstp * ( i3 - 0.5d0 )
               else
                  ! Set colatitude on the equatorial plane
                  mcolat = pi/2.d0
               endif
               ! Find lenght of longitude arc 
               longlen = ( mlongmax - mlongmin ) * mR * dsin(mcolat)
               ! Number of aggregates along longitude   
               mlongnum = floor ( longlen / mlongstp )
               ! Loop over radius
               do i1 = 1, mlongnum
                  allmrk = allmrk + 1
               ! End loop over radius
               end do
            ! End loop over colatitude
            end do
         ! End loop over radius
         end do

      endif

      Ranknnmrk = allmrk
      nnmrk     = allmrk

      RankCompDom = CompDom

   endif

   ! Set rank computational domain
   x1min      = RankCompDom%x1min
   x1max      = RankCompDom%x1max
   nx1        = RankCompDom%nx1
   x1periodic = RankCompDom%x1periodic
   x2min      = RankCompDom%x2min
   x2max      = RankCompDom%x2max
   nx2        = RankCompDom%nx2
   x2periodic = RankCompDom%x2periodic
   x3min      = RankCompDom%x3min
   x3max      = RankCompDom%x3max
   nx3        = RankCompDom%nx3
   x3periodic = RankCompDom%x3periodic
   mx1min     = RankCompDom%mx1min
   mx1max     = RankCompDom%mx1max
   mx1stp     = RankCompDom%mx1stp
   mx2min     = RankCompDom%mx2min
   mx2max     = RankCompDom%mx2max
   mx2stp     = RankCompDom%mx2stp
   mx3min     = RankCompDom%mx3min
   mx3max     = RankCompDom%mx3max
   mx3stp     = RankCompDom%mx3stp

   end subroutine DistGrid

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: RmDistGrid
   !
   !> @brief Remove all
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------

   subroutine RmDistGrid

   implicit none

   ! Local variables
   integer :: ierr

   if ( infoMPI%sizeMPI .eq. 1 ) return

   deallocate(ptx1,ptx2,ptx3,stat=ierr)
   if ( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** ERROR in RmDistGrid : deallocating ptx",infoMPI%rankMPI
      stop
   endif

   call dlt_DistComm(DistMrk)

   deallocate(mrkglo,stat=ierr)
   if ( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** ERROR in RmDistGrid : deallocating mrkglo",infoMPI%rankMPI
      stop
   endif

   deallocate(ptnmrk,stat=ierr)
   if ( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** ERROR in RmDistGrid : deallocating ptnmrk",infoMPI%rankMPI
      stop
   endif

   deallocate(indglomrk,stat=ierr)
   if ( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** ERROR in RmDistGrid : deallocating indglomrk",infoMPI%rankMPI
      stop
   endif

   end subroutine RmDistGrid

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: DistGridSphe
   !
   !> @brief Distributes the structured grid in spheric coordinates
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------

   subroutine DistGridSphe(dimensions)

   implicit none

   ! Input variables
   integer, intent(in) :: dimensions

   ! Local variables
   ! - infoMPI
   integer :: rankMPI
   integer :: nproc1,nproc2,nproc3
   ! - CompDom
   integer           :: nx1,nx2,nx3
   real(kind=double) :: x1min,x2min,x3min
   real(kind=double) :: x1max,x2max,x3max
   real(kind=double) :: mx1min,mx2min,mx3min
   real(kind=double) :: mx1max,mx2max,mx3max
   real(kind=double) :: mx1stp,mx2stp,mx3stp
   ! - local
   integer :: ierr
#if defined(DEBUG)
   integer :: i
#endif
   integer              :: iproc1,iproc2,iproc3
   real(kind=double)    :: mlongmin,mlongmax,mlongstp 
   real(kind=double)    :: mradiusmin,mradiusmax,mradiusstp 
   real(kind=double)    :: mcolatmin,mcolatmax,mcolatstp
   real(kind=double)    :: mR 
   integer              :: mrnum,mlongnum,mcolatnum
   integer              :: mrnum0,mlongnum0,mcolatnum0
   real(kind=double)    :: longlen,mlongradiansstp
   real(kind=double)    :: colatlen,mcolatradiansstp,mcolat
   integer              :: i1,i2,i3
   integer              :: allmrk,m1,m2
   integer, allocatable :: nmrkrad(:),ptrad(:)

   ! Set handles
   rankMPI = infoMPI%rankMPI
   nproc1  = infoMPI%nproc1
   nproc2  = infoMPI%nproc2
   nproc3  = infoMPI%nproc3
   x1min   = CompDom%x1min
   x1max   = CompDom%x1max
   nx1     = CompDom%nx1
   x2min   = CompDom%x2min
   x2max   = CompDom%x2max
   nx2     = CompDom%nx2
   x3min   = CompDom%x3min
   x3max   = CompDom%x3max
   nx3     = CompDom%nx3
   mx1min  = CompDom%mx1min
   mx1max  = CompDom%mx1max
   mx1stp  = CompDom%mx1stp
   mx2min  = CompDom%mx2min
   mx2max  = CompDom%mx2max
   mx2stp  = CompDom%mx2stp
   mx3min  = CompDom%mx3min
   mx3max  = CompDom%mx3max
   mx3stp  = CompDom%mx3stp

   ! allocate poniters to the first node (of each axis)
   allocate(ptx1(nproc1+1),ptx2(nproc2+1),ptx3(nproc3+1),stat=ierr)
   if ( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** error in DistGridSphe : allocating ptx",rankMPI
      stop
   endif

   ! Compute Rank computational domain

   ! Eulerian grid coordinates

   ! X
   RankCompDom%x1min      = x1min
   RankCompDom%x1max      = x1max
   RankCompDom%nx1        = nx1 
   RankCompDom%x1periodic = CompDom%x1periodic 

   ! Y
   RankCompDom%x2min      = x2min
   RankCompDom%x2max      = x2max
   RankCompDom%nx2        = nx2 
   RankCompDom%x2periodic = CompDom%x2periodic 

   ! Z
   RankCompDom%x3min      = x3min
   RankCompDom%x3max      = x3max
   RankCompDom%nx3        = nx3 
   RankCompDom%x3periodic = CompDom%x3periodic 

   ! Lagrangian grid with crystal aggregate

   call iglo2ix(rankMPI,nproc1,nproc2,iproc1,iproc2,iproc3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a)') "iproc"
   write(1000000+rankMPI,'(3i10)') iproc1,iproc2,iproc3 
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Define quantities appropriate for polar coordinates 

   ! X is Longitude
   mlongmin   = mx1min ! Radians 
   mlongmax   = mx1max ! Radians
   ! Y is Radius
   mradiusmin = mx2min ! Distance
   mradiusmax = mx2max ! Distance
   ! Z is Colatitude
   mcolatmin  = mx3min ! Radians 
   mcolatmax  = mx3max ! Radians

   mlongstp   = mx1stp ! Distance
   mradiusstp = mx2stp ! Distance 
   mcolatstp  = mx3stp ! Distance

   ! Check number of proc in Y / RADIUS
   ! Number of vertical layers
   mrnum = floor( (mradiusmax - mradiusmin) / mradiusstp )
   if ( mrnum .lt. 2*nproc2 ) then
      write(*,'(a,i5)') "*** ERROR in DistGridSphe : mrnum < 2*nproc2",rankMPI
      stop
   endif
   ! Find radius
   mR = mradiusmin + mradiusstp * 0.5d0
 
   ! Check number of proc in Z / COLAT
   mcolatnum = 1 
   mcolat = pi/2.d0
   if ( dimensions .eq. 3 ) then
      ! Find lenght of latitude arc (great circle arc)
      colatlen = ( mcolatmax - mcolatmin ) * mR
      ! Number of aggregates along latitude arc
      mcolatnum = floor( colatlen / mcolatstp ) 
      ! Check number of proc in X / LONG
      ! Latitude step in radians
      mcolatradiansstp = ( mcolatmax - mcolatmin ) / real(mcolatnum)   
      ! Colatitude in radians
      mcolat = mcolatmin + mcolatradiansstp * 0.5d0
      if ( mcolatnum .lt. 2*nproc3 ) then
         write(*,'(a,i5)') "*** ERROR in DistGridSphe : mcolatnum < 2*nproc3",rankMPI
         stop
      endif
   else
      if ( mcolatnum .lt. nproc3 ) then
         write(*,'(a,i5)') "*** ERROR in DistGridSphe : mcolatnum < nproc3",rankMPI
         stop
      endif
   endif

   ! Find lenght of longitude arc
   longlen = ( mlongmax - mlongmin ) * mR * dsin(mcolat)
   ! Number of aggregates along longitude
   mlongnum = floor(longlen / mlongstp)
   if ( mlongnum .lt. 2*nproc1 ) then
      write(*,'(a,i5)') "*** ERROR in DistGridSphe : mlongnum < 2*nproc1",rankMPI
      stop
   endif

   ! Store value
   mrnum0     = mrnum
   mcolatnum0 = mcolatnum
   mlongnum0  = mlongnum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a,i10)') "mlongnum0 ",mlongnum
   write(1000000+rankMPI,'(a,i10)') "mrnum0    ",mrnum 
   write(1000000+rankMPI,'(a,i10)') "mcolatnum0",mcolatnum
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Compute ptx1 (X / LONG)
   call SetUnk(nproc1,mlongnum0,ptx1)
   
   ! Compute ptx2 (Y / RADIUS)
   ! Not a simple distribution, but a uniform ditributoins in terms of number of marks

   ! allocate scratches
   allocate(nmrkrad(mrnum),ptrad(nproc2+1),stat=ierr)

   ! Get number of marks for each radius
   allmrk = 0
   ! Loop over radius
   do i2 = 1,mrnum
      ! Find radius
      mR = mradiusmin + mradiusstp * ( i2 - 0.5d0 )
      ! Initialize number of marks
      nmrkrad(i2) = 0 
      ! Colatitude step in radians
      mcolatnum = 1
      mcolatradiansstp = 0.d0 
      ! Check dimensions to compute number of aggregates along colatitude arc
      if ( dimensions .eq. 3 ) then
         ! Find lenght of colatitude arc (great circle arc)
         colatlen = ( mcolatmax - mcolatmin ) * mR
         ! Number of aggregates along colatitude arc
         mcolatnum = floor ( colatlen / mcolatstp)
         ! Colatitude step in radians
         mcolatradiansstp = ( mcolatmax - mcolatmin ) / real(mcolatnum)
      endif
      ! Loop over colatitude
      do i3 = 1,mcolatnum
         ! Colatitude in radians
         if ( dimensions .eq. 3 ) then
            mcolat = mcolatmin + mcolatradiansstp * ( i3 - 0.5d0 )
         else
            ! Set colatitude on the equatorial plane
            mcolat = pi/2.d0
         endif
         ! Find lenght of longitude arc 
         longlen = ( mlongmax - mlongmin ) * mR * dsin(mcolat)
         ! Number of aggregates along longitude   
         mlongnum = floor ( longlen / mlongstp )
         ! Update number of marks
         nmrkrad(i2) = nmrkrad(i2) + mlongnum 
         ! Loop over radius
         do i1 = 1, mlongnum
            allmrk = allmrk + 1
         ! End loop over radius
         end do
      ! End loop over colatitude
      end do
   ! End loop over radius
   end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a,i10)') "allmrk",allmrk
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call SetUnk(nproc2,allmrk,ptrad)

   ptx2(1) = 1
   m2 = 2
   m1 = nmrkrad(1) 
   do i2 = 2,mrnum
      m1 = m1 + nmrkrad(i2)
      if ( m1 .gt. ptrad(m2) ) then
         ptx2(m2) = i2
         m2 = m2 + 1
      endif
   end do
   ptx2(nproc2+1) = mrnum + 1

   ! deallocate scratch
   deallocate(ptrad)

   ! Compute ptx3 (Z / COLAT)
   call SetUnk(nproc3,mcolatnum0,ptx3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a)') "ptx"
   do i = 1,nproc1
      write(1000000+rankMPI,'(3i10)') i,ptx1(i),ptx1(i+1)-1
   end do
   write(1000000+rankMPI,*)
   do i = 1,nproc2
      write(1000000+rankMPI,'(3i10)') i,ptx2(i),ptx2(i+1)-1
   end do
    write(1000000+rankMPI,*)
   do i = 1,nproc3
      write(1000000+rankMPI,'(3i10)') i,ptx3(i),ptx3(i+1)-1
   end do
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Compute local min/max

   ! X / LONG
   mlongradiansstp = ( mlongmax - mlongmin ) / real(mlongnum0)
   RankCompDom%mx1min = mx1min + (ptx1(iproc1)   - 1) * mlongradiansstp
   RankCompDom%mx1max = mx1min + (ptx1(iproc1+1) - 1) * mlongradiansstp
   if ( RankCompDom%mx1max .gt. mx1max ) RankCompDom%mx1max = mx1max
   RankCompDom%mx1stp = mx1stp

   ! Y / RADIUS
   RankCompDom%mx2min = mx2min + (ptx2(iproc2)   - 1) * mradiusstp
   RankCompDom%mx2max = mx2min + (ptx2(iproc2+1) - 1) * mradiusstp
   if ( RankCompDom%mx2max .gt. mx2max ) RankCompDom%mx2max = mx2max
   RankCompDom%mx2stp = mx2stp

   ! Z / COLAT
   if ( dimensions .eq. 3 ) then
      mcolatradiansstp = ( mcolatmax - mcolatmin ) / real(mcolatnum0)
      RankCompDom%mx3min = mx3min + (ptx3(iproc3)   - 1) * mcolatradiansstp
      RankCompDom%mx3max = mx3min + (ptx3(iproc3+1) - 1) * mcolatradiansstp
      if ( RankCompDom%mx3max .gt. mx3max ) RankCompDom%mx3max = mx3max
   else
      RankCompDom%mx3min = mx3min
      RankCompDom%mx3max = mx3max
   endif
   RankCompDom%mx3stp = mx3stp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a)') "mxmin - mxmax - mxstp"
   write(1000000+rankMPI,'(a,3e15.6)') "x1 ",RankCompDom%mx1min,RankCompDom%mx1max,RankCompDom%mx1stp
   write(1000000+rankMPI,'(a,3e15.6)') "x2 ",RankCompDom%mx2min,RankCompDom%mx2max,RankCompDom%mx2stp
   write(1000000+rankMPI,'(a,3e15.6)') "x3 ",RankCompDom%mx3min,RankCompDom%mx3max,RankCompDom%mx3stp
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Retrieve global reference for the Lagrangian nodes
   call MarkLoc2GloSphe(dimensions,allmrk,mrnum,nmrkrad,nproc1,nproc2,nproc3,iproc1,iproc2,iproc3)

   ! Deallocate scratch
   deallocate(nmrkrad)

   ! Index communications of the Lagrangian nodes 
   call DistMrkInd(allmrk)

   end subroutine DistGridSphe

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: MarkLoc2GloSphe
   !
   !> @brief Retrieves global marknum index for spherical grid
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------
   subroutine MarkLoc2GloSphe(dimensions,allmrk,allmrnum,nmrkrad,nproc1,nproc2,nproc3,    &
   &                          iproc1,iproc2,iproc3)

   implicit none

   ! Input variables
   integer, intent(in) :: dimensions
   integer, intent(in) :: allmrk
   integer, intent(in) :: allmrnum 
   integer, intent(in) :: nmrkrad(allmrnum) 
   integer, intent(in) :: nproc1,nproc2,nproc3 
   integer, intent(in) :: iproc1,iproc2,iproc3 

   ! Local variables
   ! - infoMPI
   integer :: rankMPI
   ! - RankCompDom
   real(kind=double) :: mx1min,mx2min,mx3min
   real(kind=double) :: mx1max,mx2max,mx3max
   real(kind=double) :: mx1stp,mx2stp,mx3stp
   ! - local
   integer           :: ierr,errMPI
   real(kind=double) :: mlongmin,mlongmax,mlongstp 
   real(kind=double) :: mradiusmin,mradiusmax,mradiusstp 
   real(kind=double) :: mcolatmin,mcolatmax,mcolatstp
   real(kind=double) :: mR 
   integer           :: mlongnum,mcolatnum
   real(kind=double) :: longlen,mlongradiansstp,mlong
   real(kind=double) :: colatlen,mcolatradiansstp,mcolat
   integer           :: i1,i2,i3
   integer           :: ii1,ii2,ii3
   integer           :: redmrk
   integer           :: imrk
   integer           :: npre1,npre2,npre3

   ! Set handles
   rankMPI = infoMPI%rankMPI

   ! Compute local number of marks

   ! Set handles
   mx1min  = RankCompDom%mx1min
   mx1max  = RankCompDom%mx1max
   mx1stp  = RankCompDom%mx1stp
   mx2min  = RankCompDom%mx2min
   mx2max  = RankCompDom%mx2max
   mx2stp  = RankCompDom%mx2stp
   mx3min  = RankCompDom%mx3min
   mx3max  = RankCompDom%mx3max
   mx3stp  = RankCompDom%mx3stp

   ! X is Longitude
   mlongmin   = mx1min ! Radians 
   mlongmax   = mx1max ! Radians
   ! Y is Radius
   mradiusmin = mx2min ! Distance
   mradiusmax = mx2max ! Distance
   ! Z is Colatitude
   mcolatmin  = mx3min ! Radians 
   mcolatmax  = mx3max ! Radians

   mlongstp   = mx1stp ! Distance
   mradiusstp = mx2stp ! Distance 
   mcolatstp  = mx3stp ! Distance

   ! Initialiaze local nnmrk
   Ranknnmrk = 0

   ! Loop over radius
   loop2a: do i2 = 1,allmrnum

      ! Find radius
      mR = CompDom%mx2min + mradiusstp * ( i2 - 0.5d0 )
 
      if ( mR .gt. mradiusmax ) exit  loop2a
      if ( mR .eq. mradiusmax .and. iproc2 .ne. nproc2 ) exit loop2a
      if ( mR .lt. mradiusmin ) cycle loop2a
 
      ! Colatitude step in radians
      mcolatnum = 1
      mcolatradiansstp = 0.d0

      ! Check dimensions to compute number of aggregates along colatitude arc
      if ( dimensions .eq. 3 ) then
         ! Find lenght of colatitude arc (great circle arc)
         colatlen = ( CompDom%mx3max - CompDom%mx3min ) * mR
         ! Number of aggregates along colatitude arc
         mcolatnum = floor ( colatlen /CompDom%mx3stp)
         ! Colatitude step in radians
         mcolatradiansstp = ( CompDom%mx3max - CompDom%mx3min ) / real(mcolatnum)
      endif

      ! Loop over colatitude
      loop3a: do i3 = 1,mcolatnum

         ! Colatitude in radians
         if ( dimensions .eq. 3 ) then
            mcolat = CompDom%mx3min + mcolatradiansstp * ( i3 - 0.5d0 )
            if ( mcolat .gt. mcolatmax ) exit loop3a
            if ( mcolat .eq. mcolatmax .and. iproc3 .ne. nproc3 ) exit loop3a
            if ( mcolat .lt. mcolatmin ) cycle loop3a 
         else
            ! Set colatitude on the equatorial plane
            mcolat = pi/2.d0
         endif

         ! Find lenght of longitude arc 
         longlen = ( CompDom%mx1max - CompDom%mx1min ) * mR * dsin(mcolat)
         ! Number of aggregates along longitude   
         mlongnum = floor ( longlen / CompDom%mx1stp )
         ! Longitude step in radians
         mlongradiansstp = ( CompDom%mx1max - CompDom%mx1min ) / real(mlongnum)

         loop1a: do i1 = 1,mlongnum
            mlong = CompDom%mx1min + mlongradiansstp*(i1-0.5d0)
            if ( mlong .gt. mlongmax ) exit loop1a
            if ( mlong .eq. mlongmax .and. iproc1 .ne. nproc1 ) exit loop1a
            if ( mlong .ge. mlongmin ) Ranknnmrk = Ranknnmrk + 1 
         end do loop1a
 
      ! End loop over colatitude
      end do loop3a

   ! End loop over radius
   end do loop2a

   ! Check marknum
   redmrk = 0
   call MPI_Allreduce(Ranknnmrk,redmrk,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,errMPI)
   if ( allmrk .ne. redmrk ) then
      write(*,'(a,i5)') "*** ERROR in MarkLoc2GloSphe : wrong mark computation",rankMPI
      write(*,'(3i5)') Ranknnmrk,redmrk,allmrk 
      stop
   endif

   ! Allocate mrkglo
   allocate(mrkglo(Ranknnmrk),stat=ierr)
   if ( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** ERROR in MarkLoc2GloSphe : allocating mrkglo",rankMPI
      stop
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a, i10)') "Rankmark",Ranknnmrk
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set global indexes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a)') "mrkglo"
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Loop over radius
   imrk = 0
   npre2 = 0
   ii2 = 0
   loop2b: do i2 = 1,allmrnum

      ! Find radius
      mR = CompDom%mx2min + mradiusstp * ( i2 - 0.5d0 )
 
      if ( mR .gt. mradiusmax ) exit loop2b
      if ( mR .eq. mradiusmax .and. iproc2 .ne. nproc2 ) exit loop2b
      if ( mR .lt. mradiusmin ) then
         ! Compute radius offset
         npre2 = npre2 + nmrkrad(i2) 
         cycle loop2b
      endif

      ii2 = ii2 + 1

      ! Colatitude step in radians
      mcolatnum = 1
      mcolatradiansstp = 0.d0

      ! Check dimensions to compute number of aggregates along colatitude arc
      if ( dimensions .eq. 3 ) then
         ! Find lenght of colatitude arc (great circle arc)
         colatlen = ( CompDom%mx3max - CompDom%mx3min ) * mR
         ! Number of aggregates along colatitude arc
         mcolatnum = floor ( colatlen /CompDom%mx3stp)
         ! Colatitude step in radians
         mcolatradiansstp = ( CompDom%mx3max - CompDom%mx3min ) / real(mcolatnum)
      endif

      ! Loop over colatitude
      npre3 = 0
      ii3 = 0
      loop3b: do i3 = 1,mcolatnum

         ! Colatitude in radians
         if ( dimensions .eq. 3 ) then
            mcolat = CompDom%mx3min + mcolatradiansstp * ( i3 - 0.5d0 )
            if ( mcolat .gt. mcolatmax ) exit loop3b
            if ( mcolat .eq. mcolatmax .and. iproc3 .ne. nproc3 ) exit loop3b
            if ( mcolat .lt. mcolatmin ) then
               ! Colatitude in radians
               mcolat = CompDom%mx3min + mcolatradiansstp * ( i3 - 0.5d0 )
               ! Find lenght of longitude arc 
               longlen = ( CompDom%mx1max - CompDom%mx1min ) * mR * dsin(mcolat)
               ! Number of aggregates along longitude   
               mlongnum = floor( longlen / CompDom%mx1stp )
               ! Colatitude offset
               npre3 = npre3 + mlongnum
               cycle loop3b
            endif
         else
            ! Set colatitude on the equatorial plane
            mcolat = pi/2.d0
         endif

         ii3 = ii3 + 1

         ! Find lenght of longitude arc 
         longlen = ( CompDom%mx1max - CompDom%mx1min ) * mR * dsin(mcolat)
         ! Number of aggregates along longitude   
         mlongnum = floor ( longlen / CompDom%mx1stp )
         ! Longitude step in radians
         mlongradiansstp = ( CompDom%mx1max - CompDom%mx1min ) / real(mlongnum)

         ! Loop over longitude
         npre1 = 0
         ii1 = 0
         loop1b: do i1 = 1,mlongnum
            mlong = CompDom%mx1min + mlongradiansstp*(i1-0.5d0)
            if ( mlong .gt. mlongmax ) exit loop1b
            if ( mlong .eq. mlongmax .and. iproc1 .ne. nproc1 ) exit loop1b
            if ( mlong .lt. mlongmin ) then
               ! Compute longitude offset
               npre1 = npre1 + 1
               cycle loop1b
            else
               ii1 = ii1 + 1
               imrk = imrk + 1
               mrkglo(imrk) = ii1 + npre1 + npre2 + npre3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
               write(1000000+rankMPI,'(5i10)') ii1,ii2,ii3,imrk,mrkglo(imrk)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            endif 
         ! End loop over longitude
         end do loop1b 
 
         ! Update colatitude offset
         if ( dimensions .eq. 3 ) then
            npre3 = npre3 + mlongnum 
         endif
 
      ! End loop over colatitude
      end do loop3b

      ! Update radius offset
      npre2 = npre2 + nmrkrad(i2) 

   ! End loop over radius
   end do loop2b

   end subroutine MarkLoc2GloSphe

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: GetCoordSphe
   !
   !> @brief Retrieves local marknum coordinates for spherical grid
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------
   subroutine GetCoordSphe(dimensions,innmrk,outmx1,outmx2,outmx3)

   implicit none

   ! Input variables
   integer, intent(in) :: dimensions
   integer, intent(in) :: innmrk

   ! Output variables
   real(kind=double), intent(out) :: outmx1(innmrk)
   real(kind=double), intent(out) :: outmx2(innmrk)
   real(kind=double), intent(out) :: outmx3(innmrk)

   ! Local variables
   ! - infoMPI
   integer :: rankMPI,nproc1,nproc2,nproc3
   ! - RankCompDom
   real(kind=double) :: mx1min,mx2min,mx3min
   real(kind=double) :: mx1max,mx2max,mx3max
   real(kind=double) :: mx1stp,mx2stp,mx3stp
   ! - local
   integer           :: iproc1,iproc2,iproc3
   real(kind=double) :: mlongmin,mlongmax,mlongstp 
   real(kind=double) :: mradiusmin,mradiusmax,mradiusstp 
   real(kind=double) :: mcolatmin,mcolatmax,mcolatstp
   real(kind=double) :: mR 
   integer           :: mrnum,mlongnum,mcolatnum
   real(kind=double) :: longlen,mlongradiansstp,mlong
   real(kind=double) :: colatlen,mcolatradiansstp,mcolat
   integer           :: i1,i2,i3
   integer           :: imrk

   ! Set handles
   rankMPI = infoMPI%rankMPI
   nproc1  = infoMPI%nproc1
   nproc2  = infoMPI%nproc2
   nproc3  = infoMPI%nproc3
   mx1min  = RankCompDom%mx1min
   mx1max  = RankCompDom%mx1max
   mx1stp  = RankCompDom%mx1stp
   mx2min  = RankCompDom%mx2min
   mx2max  = RankCompDom%mx2max
   mx2stp  = RankCompDom%mx2stp
   mx3min  = RankCompDom%mx3min
   mx3max  = RankCompDom%mx3max
   mx3stp  = RankCompDom%mx3stp

   ! Check nnmrk
   if ( innmrk .ne. Ranknnmrk ) then
      write(*,'(a,i5)') "*** ERROR in GetCoordSphe : innmrk /= Ranknnmrk",rankMPI
      stop
   endif
   
   ! Retrieve proc info  
   call iglo2ix(rankMPI,nproc1,nproc2,iproc1,iproc2,iproc3)

   ! Compute local number of marks
   ! X is Longitude
   mlongmin   = mx1min ! Radians 
   mlongmax   = mx1max ! Radians
   ! Y is Radius
   mradiusmin = mx2min ! Distance
   mradiusmax = mx2max ! Distance
   ! Z is Colatitude
   mcolatmin  = mx3min ! Radians 
   mcolatmax  = mx3max ! Radians

   mlongstp   = mx1stp ! Distance
   mradiusstp = mx2stp ! Distance 
   mcolatstp  = mx3stp ! Distance

   ! Initialiaze main loop
   imrk = 0

   ! Loop over radius
   mrnum = floor( (CompDom%mx2max - CompDom%mx2min) / mradiusstp )
   loop2a: do i2 = 1,mrnum

      ! Find radius
      mR = CompDom%mx2min + mradiusstp * ( i2 - 0.5d0 )
 
      if ( mR .gt. mradiusmax ) exit loop2a
      if ( mR .eq. mradiusmax .and. iproc2 .ne. nproc2 ) exit loop2a
      if ( mR .lt. mradiusmin ) cycle loop2a
 
      ! Colatitude step in radians
      mcolatnum = 1
      mcolatradiansstp = 0.d0

      ! Check dimensions to compute number of aggregates along colatitude arc
      if ( dimensions .eq. 3 ) then
         ! Find lenght of colatitude arc (great circle arc)
         colatlen = ( CompDom%mx3max - CompDom%mx3min ) * mR
         ! Number of aggregates along colatitude arc
         mcolatnum = floor ( colatlen /CompDom%mx3stp)
         ! Colatitude step in radians
         mcolatradiansstp = ( CompDom%mx3max - CompDom%mx3min ) / real(mcolatnum)
      endif

      ! Loop over colatitude
      loop3a: do i3 = 1,mcolatnum

         ! Colatitude in radians
         if ( dimensions .eq. 3 ) then
            mcolat = CompDom%mx3min + mcolatradiansstp * ( i3 - 0.5d0 )
            if ( mcolat .gt. mcolatmax ) exit loop3a
            if ( mcolat .eq. mcolatmax .and. iproc3 .ne. nproc3 ) exit loop3a
            if ( mcolat .lt. mcolatmin ) cycle loop3a 
         else
            ! Set colatitude on the equatorial plane
            mcolat = pi/2.d0
         endif

         ! Find lenght of longitude arc 
         longlen = ( CompDom%mx1max - CompDom%mx1min ) * mR * dsin(mcolat)
         ! Number of aggregates along longitude   
         mlongnum = floor ( longlen / CompDom%mx1stp )
         ! Longitude step in radians
         mlongradiansstp = ( CompDom%mx1max - CompDom%mx1min ) / real(mlongnum)

         loop1a: do i1 = 1,mlongnum
            mlong = CompDom%mx1min + mlongradiansstp*(i1-0.5d0)
            if ( mlong .gt. mlongmax ) exit loop1a
            if ( mlong .eq. mlongmax .and. iproc1 .ne. nproc1 ) exit loop1a
            if ( mlong .ge. mlongmin ) then
               imrk = imrk + 1
               !Spherical coordinates
               !Longitude
               outmx1(imrk) = mlong  ! mR * sin(mcolat) * cos(mlong)
               !Radial
               outmx2(imrk) = mR     ! mR * sin(mcolat) * sin(mlong)
               !Colatitude
               outmx3(imrk) = mcolat ! mR * cos(mcolat)
            endif
         end do loop1a
 
      ! End loop over colatitude
      end do loop3a

   ! End loop over radius
   end do loop2a

   end subroutine GetCoordSphe

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: DistGridCart
   !
   !> @brief Distributes the structured grid in cartesian coordinates
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------

   subroutine DistGridCart(dimensions)

   implicit none

   ! Input variables
   integer, intent(in) :: dimensions

   ! Local variables
   ! - infoMPI
   integer :: rankMPI
   integer :: nproc1,nproc2,nproc3
   ! - CompDom
   integer           :: nx1,nx2,nx3
   real(kind=double) :: x1min,x2min,x3min
   real(kind=double) :: x1max,x2max,x3max
   real(kind=double) :: mx1min,mx2min,mx3min
   real(kind=double) :: mx1max,mx2max,mx3max
   real(kind=double) :: mx1stp,mx2stp,mx3stp
   ! - local
   integer :: ierr
#if defined(DEBUG)
   integer :: i
#endif
   integer :: iproc1,iproc2,iproc3
   integer :: mnx1,mnx2,mnx3
   integer :: allmrk

   ! set handles
   rankMPI = infoMPI%rankMPI
   nproc1  = infoMPI%nproc1
   nproc2  = infoMPI%nproc2
   nproc3  = infoMPI%nproc3
   x1min   = CompDom%x1min
   x1max   = CompDom%x1max
   nx1     = CompDom%nx1
   x2min   = CompDom%x2min
   x2max   = CompDom%x2max
   nx2     = CompDom%nx2
   x3min   = CompDom%x3min
   x3max   = CompDom%x3max
   nx3     = CompDom%nx3
   mx1min  = CompDom%mx1min
   mx1max  = CompDom%mx1max
   mx1stp  = CompDom%mx1stp
   mx2min  = CompDom%mx2min
   mx2max  = CompDom%mx2max
   mx2stp  = CompDom%mx2stp
   mx3min  = CompDom%mx3min
   mx3max  = CompDom%mx3max
   mx3stp  = CompDom%mx3stp

   ! allocate poniters to the first node (of each axis)
   allocate(ptx1(nproc1+1),ptx2(nproc2+1),ptx3(nproc3+1),stat=ierr)
   if ( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** error in DistGridCart : allocating ptx",rankmpi
      stop
   endif

   ! Compute Rank computational domain

   ! Eulerian grid coordinates

   ! X
   RankCompDom%x1min      = x1min
   RankCompDom%x1max      = x1max
   RankCompDom%nx1        = nx1 
   RankCompDom%x1periodic = CompDom%x1periodic 

   ! Y
   RankCompDom%x2min      = x2min
   RankCompDom%x2max      = x2max
   RankCompDom%nx2        = nx2 
   RankCompDom%x2periodic = CompDom%x2periodic 

   ! Z
   RankCompDom%x3min      = x3min
   RankCompDom%x3max      = x3max
   RankCompDom%nx3        = nx3 
   RankCompDom%x3periodic = CompDom%x3periodic 

   ! Lagrangian grid with crystal aggregate

   ! Retrieve process indexes
   call iglo2ix(rankMPI,nproc1,nproc2,iproc1,iproc2,iproc3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a)') "iproc"
   write(1000000+rankMPI,'(3i10)') iproc1,iproc2,iproc3 
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Check number of process

   ! X
   mnx1 = nint((mx1max-mx1min)/mx1stp)
   if ( mnx1 .lt. 2*nproc1 ) then
      write(*,'(a,i5)') "*** ERROR in DistGridCart : mnx1 < 2*nproc1",rankMPI
      stop
   endif

   ! Y
   mnx2 = nint((mx2max-mx2min)/mx2stp)
   if ( mnx2 .lt. 2*nproc2 ) then
      write(*,'(a,i5)') "*** ERROR in DistGridCart : mnx2 < nproc2",rankMPI
      stop
   endif

   ! Z
   mnx3 = 1
   if ( dimensions .eq. 3 ) then
      mnx3 = nint((mx3max-mx3min)/mx3stp)
      if ( mnx3 .lt. 2*nproc3 ) then
         write(*,'(a,i5)') "*** ERROR in DistGridCart : mnx3 < 2*nproc3",rankMPI
         stop
      endif
   else
      if ( mnx3 .lt. nproc3 ) then
         write(*,'(a,i5)') "*** ERROR in DistGridCart : mnx3 < nproc3",rankMPI
         stop
      endif
   endif

   allmrk = mnx1 * mnx2 * mnx3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,*) "allmnx1",mnx1
   write(1000000+rankMPI,*) "allmnx2",mnx2
   write(1000000+rankMPI,*) "allmnx3",mnx3
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,*) "allmrk",allmrk
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Compute ptx

   ! X
   call SetUnk(nproc1,mnx1,ptx1)

   ! Y
   call SetUnk(nproc2,mnx2,ptx2)
   
   ! Z
   call SetUnk(nproc3,mnx3,ptx3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a)') "ptx"
   do i = 1,nproc1
      write(1000000+rankMPI,'(3i10)') i,ptx1(i),ptx1(i+1)-1
   end do
   write(1000000+rankMPI,*)
   do i = 1,nproc2
      write(1000000+rankMPI,'(3i10)') i,ptx2(i),ptx2(i+1)-1
   end do
    write(1000000+rankMPI,*)
   do i = 1,nproc3
      write(1000000+rankMPI,'(3i10)') i,ptx3(i),ptx3(i+1)-1
   end do
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Compute local min/max

   ! X
   RankCompDom%mx1min = mx1min + (ptx1(iproc1)   - 1) * mx1stp
   RankCompDom%mx1max = mx1min + (ptx1(iproc1+1) - 1) * mx1stp
   if ( RankCompDom%mx1max .gt. mx1max ) RankCompDom%mx1max = mx1max
   RankCompDom%mx1stp = mx1stp

   ! Y
   RankCompDom%mx2min = mx2min + (ptx2(iproc2)   - 1) * mx2stp
   RankCompDom%mx2max = mx2min + (ptx2(iproc2+1) - 1) * mx2stp
   if ( RankCompDom%mx2max .gt. mx2max ) RankCompDom%mx2max = mx2max
   RankCompDom%mx2stp = mx2stp

   ! Z
   if ( dimensions .eq. 3 ) then
      RankCompDom%mx3min = mx3min + (ptx3(iproc3)   - 1) * mx3stp
      RankCompDom%mx3max = mx3min + (ptx3(iproc3+1) - 1) * mx3stp
      if ( RankCompDom%mx3max .gt. mx3max ) RankCompDom%mx3max = mx3max
   else
      RankCompDom%mx3min = mx3min
      RankCompDom%mx3max = mx3max
   endif
   RankCompDom%mx3stp = mx3stp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a)') "mxmin - mxmax - mxstp"
   write(1000000+rankMPI,'(a,3e15.6)') "x1 ",RankCompDom%mx1min,RankCompDom%mx1max,RankCompDom%mx1stp
   write(1000000+rankMPI,'(a,3e15.6)') "x2 ",RankCompDom%mx2min,RankCompDom%mx2max,RankCompDom%mx2stp
   write(1000000+rankMPI,'(a,3e15.6)') "x3 ",RankCompDom%mx3min,RankCompDom%mx3max,RankCompDom%mx3stp
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Retrieve global reference for the Lagrangian nodes
   call MarkLoc2GloCart(dimensions,mnx1,mnx2,mnx3,nproc1,nproc2,nproc3,iproc1,iproc2,iproc3)

   ! Index communications of the Lagrangian nodes 
   call DistMrkInd(allmrk)

   end subroutine DistGridCart

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: MarkLoc2GloCart
   !
   !> @brief Retrieves global marknum index for cartesian grid
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------
   subroutine MarkLoc2GloCart(dimensions,allmnx1,allmnx2,allmnx3,nproc1,nproc2,nproc3,    &
   &                          iproc1,iproc2,iproc3)

   implicit none

   ! Input variables
   integer, intent(in) :: dimensions
   integer, intent(in) :: allmnx1,allmnx2,allmnx3
   integer, intent(in) :: nproc1,nproc2,nproc3 
   integer, intent(in) :: iproc1,iproc2,iproc3 

   ! Local variables
   ! - infoMPI
   integer :: rankMPI
   ! - RankCompDom
   real(kind=double) :: mx1min,mx2min,mx3min
   real(kind=double) :: mx1max,mx2max,mx3max
   real(kind=double) :: mx1stp,mx2stp,mx3stp
   ! - local
   integer           :: ierr,errMPI
   integer           :: i1,i2,i3
   integer           :: ii1,ii2,ii3
   integer           :: npre1,npre2,npre3 
   real(kind=double) :: xx,yy,zz 
   integer           :: allmrk,redmrk
   integer           :: imrk 

   ! Set handles
   rankMPI = infoMPI%rankMPI
   mx1min  = RankCompDom%mx1min
   mx1max  = RankCompDom%mx1max
   mx1stp  = RankCompDom%mx1stp
   mx2min  = RankCompDom%mx2min
   mx2max  = RankCompDom%mx2max
   mx2stp  = RankCompDom%mx2stp
   mx3min  = RankCompDom%mx3min
   mx3max  = RankCompDom%mx3max
   mx3stp  = RankCompDom%mx3stp

   ! Initialize local nnmrk
   Ranknnmrk = 0

   ! loop over X
   loop1a: do i1 = 1,allmnx1

      xx = CompDom%mx1min + mx1stp * ( i1 - 0.5d0 )

      if ( xx .gt. mx1max ) exit loop1a
      if ( xx .eq. mx1max .and. iproc1 .ne. nproc1 ) exit loop1a
      if ( xx .lt. mx1min ) cycle loop1a

      ! loop over Y
      loop2a: do i2 = 1,allmnx2

         yy = CompDom%mx2min + mx2stp * ( i2 - 0.5d0 )

         if ( yy .gt. mx2max ) exit loop2a
         if ( yy .eq. mx2max .and. iproc2 .ne. nproc2 ) exit loop2a
         if ( yy .lt. mx2min ) cycle loop2a

         ! loop over Z
         loop3a: do i3 = 1,allmnx3

            if (dimensions .eq. 3 ) then
               zz = CompDom%mx3min + mx3stp * ( i3 - 0.5d0 )
            else
               zz = 0.d0
            endif

            if ( zz .gt. mx3max ) exit loop3a
            if ( zz .eq. mx3max .and. iproc3 .ne. nproc3 ) exit loop3a
            if ( zz .ge. mx3min ) Ranknnmrk = Ranknnmrk + 1 

         ! end loop over Z
         end do loop3a

      ! end loop over Y
      end do loop2a

   ! end loop over X
   end do loop1a

   ! Check marknum
   allmrk = allmnx1 * allmnx2 * allmnx3
   redmrk = 0
   call MPI_Allreduce(Ranknnmrk,redmrk,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,errMPI)
   if ( allmrk .ne. redmrk ) then
      write(*,'(a,i5)') "*** ERROR in MarkLoc2GloCart : wrong mark computation",rankMPI
      write(*,'(3i5)') Ranknnmrk,redmrk,allmrk
      stop
   endif

   ! Allocate mrkglo
   allocate(mrkglo(Ranknnmrk),stat=ierr)
   if ( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** ERROR in MarkLoc2GloCart : allocating mrkglo",rankMPI
      stop
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a, i10)') "Rankmark",Ranknnmrk
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Set global indeces
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a)') "mrkglo"
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! loop over X
   imrk = 0
   npre1 = 0
   ii1 = 0
   loop1b: do i1 = 1,allmnx1

      xx = CompDom%mx1min + mx1stp * ( i1 - 0.5d0 )

      if ( xx .gt. mx1max ) exit loop1b
      if ( xx .eq. mx1max .and. iproc1 .ne. nproc1 ) exit loop1b
      if ( xx .lt. mx1min ) then
         ! compute X offset
         npre1 = npre1 + allmnx2 * allmnx3
         cycle loop1b
      endif

      ii1 = ii1 + 1

      ! loop over Y
      npre2 = 0
      ii2 = 0
      loop2b: do i2 = 1,allmnx2

         yy = CompDom%mx2min + mx2stp * ( i2 - 0.5d0 )

         if ( yy .gt. mx2max ) exit loop2b
         if ( yy .eq. mx2max .and. iproc2 .ne. nproc2 ) exit loop2b
         if ( yy .lt. mx2min ) then
            npre2 = npre2 + allmnx3
            cycle loop2b
         endif

         ii2 = ii2 + 1

         ! loop over Z
         npre3 = 0
         ii3 = 0
         loop3b: do i3 = 1,allmnx3

            if (dimensions .eq. 3 ) then
               zz = CompDom%mx3min + mx3stp * ( i3 - 0.5d0 )
            else
               zz = 0.d0
            endif

            if ( zz .gt. mx3max ) exit loop3b
            if ( zz .eq. mx3max .and. iproc3 .ne. nproc3 ) exit loop3b
            if ( zz .lt. mx3min ) then
               npre3 = npre3 + 1
               cycle loop3b
            else
               ii3 = ii3 + 1
               imrk = imrk + 1
               mrkglo(imrk) = ii3 + npre1 + npre2 + npre3 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
               write(1000000+rankMPI,'(5i10)') ii1,ii2,ii3,imrk,mrkglo(imrk)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            endif

         ! end loop over Z
         end do loop3b

         ! update X offset
         npre2 = npre2 + allmnx3

      ! end loop over Y
      end do loop2b

      ! update X offset
      npre1 = npre1 + allmnx2 * allmnx3

   ! end loop over X
   end do loop1b

   end subroutine MarkLoc2GloCart

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: GetCoordCart
   !
   !> @brief Retrieves local marknum coordinates for cartesian grid
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------
   subroutine GetCoordCart(dimensions,innmrk,outmx1,outmx2,outmx3)

   implicit none

   ! Input variables
   integer, intent(in) :: dimensions
   integer, intent(in) :: innmrk

   ! Output variables
   real(kind=double), intent(out) :: outmx1(innmrk)
   real(kind=double), intent(out) :: outmx2(innmrk)
   real(kind=double), intent(out) :: outmx3(innmrk)

   ! Local variables
   ! - infoMPI
   integer :: rankMPI,nproc1,nproc2,nproc3
   ! - RankCompDom
   real(kind=double) :: mx1min,mx2min,mx3min
   real(kind=double) :: mx1max,mx2max,mx3max
   real(kind=double) :: mx1stp,mx2stp,mx3stp
   ! - local
   integer           :: iproc1,iproc2,iproc3
   integer           :: allmnx1,allmnx2,allmnx3
   integer           :: i1,i2,i3
   integer           :: imrk
   real(kind=double) :: xx,yy,zz 

   ! Set handles
   rankMPI = infoMPI%rankMPI
   nproc1  = infoMPI%nproc1
   nproc2  = infoMPI%nproc2
   nproc3  = infoMPI%nproc3
   mx1min  = RankCompDom%mx1min
   mx1max  = RankCompDom%mx1max
   mx1stp  = RankCompDom%mx1stp
   mx2min  = RankCompDom%mx2min
   mx2max  = RankCompDom%mx2max
   mx2stp  = RankCompDom%mx2stp
   mx3min  = RankCompDom%mx3min
   mx3max  = RankCompDom%mx3max
   mx3stp  = RankCompDom%mx3stp

   ! Check nnmrk
   if ( innmrk .ne. Ranknnmrk ) then
      write(*,'(a,i5)') "*** ERROR in GetCoordCart : innmrk /= Ranknnmrk",rankMPI
      stop
   endif

   ! Retrieve global info
   allmnx1 = nint( ( CompDom%mx1max - CompDom%mx1min ) / mx1stp )
   allmnx2 = nint( ( CompDom%mx2max - CompDom%mx2min ) / mx2stp )
   allmnx3 = 1
   if ( dimensions .eq. 3 ) then
      allmnx3 = nint( ( CompDom%mx3max - CompDom%mx3min ) / mx3stp )
   endif

   ! Retrieve proc info
   call iglo2ix(rankMPI,nproc1,nproc2,iproc1,iproc2,iproc3)

   ! Set coordinates
   imrk = 0
   loop1a: do i1 = 1,allmnx1

      xx = CompDom%mx1min + mx1stp * ( i1 - 0.5d0 )

      if ( xx .gt. mx1max ) exit loop1a
      if ( xx .eq. mx1max .and. iproc1 .ne. nproc1 ) exit loop1a
      if ( xx .lt. mx1min ) cycle loop1a

      ! loop over Y
      loop2a: do i2 = 1,allmnx2

         yy = CompDom%mx2min + mx2stp * ( i2 - 0.5d0 )

         if ( yy .gt. mx2max ) exit loop2a
         if ( yy .eq. mx2max .and. iproc2 .ne. nproc2 ) exit loop2a
         if ( yy .lt. mx2min ) cycle loop2a

         ! loop over Z
         loop3a: do i3 = 1,allmnx3

            if (dimensions .eq. 3 ) then
               zz = CompDom%mx3min + mx3stp * ( i3 - 0.5d0 )
            else
               zz = 0.d0
            endif

            if ( zz .gt. mx3max ) exit loop3a
            if ( zz .eq. mx3max .and. iproc3 .ne. nproc3 ) exit loop3a
            if ( zz .ge. mx3min ) then
               imrk = imrk + 1
               outmx1(imrk) = xx
               outmx2(imrk) = yy
               outmx3(imrk) = zz
            endif 

         ! end loop over Z
         end do loop3a

      ! end loop over Y
      end do loop2a

   ! end loop over X
   end do loop1a

   end subroutine GetCoordCart

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: DistMrkInd
   !
   !> @brief Distributes indeces for Lagrangian nodes (marknum)
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------
   subroutine DistMrkInd(allmrk)

   implicit none

   ! Input variables
   integer, intent(in) :: allmrk

   ! Local types
   type :: t_list
      integer          :: nterm = 0
      integer          :: lsize = 0
      integer, pointer :: terms(:) => null()
   end type

   ! Local parameters
   integer,           parameter :: size0 = 2
   real(kind=double), parameter :: fac_resize = 1.5d0

   ! Local variables
   ! - local
   integer :: ierr
   integer :: i,j
   ! - infoMPI
   integer :: rankMPI,sizeMPI
   ! - list communication
   integer          :: lsize
   integer          :: dest
   integer          :: nterm
   integer, pointer :: terms(:) => null()
   integer          :: istart
   type(t_list)     :: ind_list(infoMPI%sizeMPI),own_list 
   integer          :: SendNterm,SendAllNterm(infoMPI%sizeMPI),SendNcomm
   integer          :: RecvNterm,RecvAllNterm(infoMPI%sizeMPI),RecvNcomm
   integer, pointer :: SendList(:) => null(),SendOffset(:) => null()
   integer, pointer :: SendComm(:) => null(),SendStatus(:) => null()
   integer, pointer :: RecvList(:) => null(),RecvOffset(:) => null()
   integer, pointer :: RecvComm(:) => null(),RecvStatus(:) => null()

   ! Set handles
   sizeMPI = infoMPI%sizeMPI
   rankMPI = infoMPI%rankMPI

   ! Initialize error flag
   ierr = 0

   ! Allocate poniters to the first global nodes
   allocate(ptnmrk(sizeMPI+1),stat=ierr)
   if ( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** ERROR in DistMrkInd : allocating ptnnum",rankMPI
      stop
   endif

   ! Retrieve distribution info
   call SetUnk(sizeMPI,allmrk,ptnmrk)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a,i10)') "allmrk",allmrk
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a)') "ptnmrk"
   do i = 1,sizeMPI
      write(1000000+rankMPI,'(3i10)') i,ptnmrk(i),ptnmrk(i+1)-1
   end do
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Initialization
   lsize = max(size0,(ptnmrk(rankMPI+1)-ptnmrk(rankMPI))/size0)

   SendNterm = 0
   do i = 1,sizeMPI
      SendAllNterm(i) = 0
      ind_list(i)%lsize = lsize
      allocate(ind_list(i)%terms(lsize),stat=ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in DistMrkInd : allocating ind_list",rankMPI
         stop
      endif
   end do
   own_list%lsize = lsize
   allocate(own_list%terms(lsize),stat=ierr)
   if( ierr .ne. 0 ) then
      write(*,'(a,i5)') "*** ERROR in DistMrkInd : allocating own_list",rankMPI
      stop
   endif

   ! Loop over nodes
   SendNcomm = 0
   do i = 1, Ranknnmrk 

      ! Retrieve destination ranks
      dest = 0
      j = 1
      do while ( ptnmrk(j) .le. mrkglo(i) )
         dest = dest + 1
         j = j + 1
      end do

     ! Update lists
     if ( dest .ne. rankMPI ) then
        ! not owned node
        ! update number of communications
        if ( ind_list(dest)%nterm .eq. 0 ) SendNcomm = SendNcomm + 1
        ! update number of terms to communicate
        ind_list(dest)%nterm = ind_list(dest)%nterm + 1
        SendNterm = SendNterm + 1
        ! store term to communicate
        if ( ind_list(dest)%nterm .lt. ind_list(dest)%lsize ) then
           ind_list(dest)%terms(ind_list(dest)%nterm) = i
        else
           lsize = int ( real ( lsize )  * fac_resize )
           ierr = resizei(lsize,ind_list(dest)%terms)
           if( ierr .ne. 0 ) then
              write(*,'(a,i10,i5)') "*** ERROR in DistMrkInd : resizing ind_list",lsize,rankMPI
              stop
           endif
           ind_list(dest)%lsize = lsize
           ind_list(dest)%terms(ind_list(dest)%nterm) = i
        endif
     else
        ! owned node
        own_list%nterm = own_list%nterm + 1
        if ( own_list%nterm .lt. own_list%lsize ) then
           own_list%terms(own_list%nterm) = i
        else
           lsize = int ( real ( lsize )  * fac_resize )
           ierr = resizei(lsize,own_list%terms)
           if( ierr .ne. 0 ) then
              write(*,'(a,i5)') "*** ERROR in DistMrkInd : resizeing own_list",rankMPI
              stop
           endif
           own_list%lsize = lsize
           own_list%terms(own_list%nterm) = i
        endif
     endif

   ! end loop over nodes
   end do
 
   ! Set-up sending buffer and remove lists
   if ( SendNcomm .gt. 0 ) then
      allocate(SendList(SendNterm),SendOffset(SendNcomm+1),                                  &
      &        SendComm(SendNcomm),SendStatus(SendNcomm),stat=ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in DistNodInd : allocating Send scratches",rankMPI
         stop
      endif
      SendNcomm = 0
      SendOffset(1) = 1
      istart = 1
      do i = 1,sizeMPI
         nterm =  ind_list(i)%nterm
         terms => ind_list(i)%terms
         do j = 1,nterm
            SendList(istart) = mrkglo(terms(j))
            istart = istart + 1
         end do
         SendAllNterm(i) = nterm
         if( nterm .gt. 0 )then
            SendNcomm = SendNcomm + 1
            SendComm(SendNcomm) = i
            SendOffset(SendNcomm+1) = SendOffset(SendNcomm) + nterm
         endif
         ind_list(i)%nterm = 0
         deallocate(ind_list(i)%terms,stat=ierr)
         if( ierr .ne. 0 ) then
            write(*,'(a,i5)') "*** ERROR in DistNodInd : deallocating ind_list",rankMPI
            stop
         endif
      end do
   else
      do i = 1,sizeMPI
         deallocate(ind_list(i)%terms,stat=ierr)
         if( ierr .ne. 0 ) then
            write(*,'(a,i5)') "*** ERROR in DistNodInd : deallocating ind_list",rankMPI
            stop
         endif
      end do
   endif
 
   ! Scatter sizes and set-up receiving buffers
   RecvNterm = 0
   RecvNcomm = 0
   do i = 1,sizeMPI
      call MPI_scatter(SendAllNterm,1,MPI_INTEGER,RecvAllNterm(i),1,MPI_INTEGER,i-1,         &
      &                MPI_COMM_WORLD,ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in DistNodInd : scattering sizes",rankMPI
         stop
      endif
       RecvNterm = RecvNterm + RecvAllNterm(i)
       if ( RecvAllNterm(i) .gt. 0 ) RecvNcomm = RecvNcomm + 1
   end do
   if ( RecvNcomm .gt. 0 ) then
      allocate(RecvList(RecvNterm),RecvOffset(RecvNcomm+1),                                  &
      &        RecvComm(RecvNcomm),RecvStatus(RecvNcomm),stat=ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in DistNodInd : allocating Recv scratches",rankMPI
         stop
      endif
      RecvNcomm = 0
      RecvOffset(1) = 1
      do i = 1,sizeMPI
         if ( RecvAllNterm(i) .gt. 0 ) then
            RecvNcomm = RecvNcomm + 1
            RecvComm(RecvNcomm) = i
            RecvOffset(RecvNcomm+1) = RecvOffset(RecvNcomm) + RecvAllNterm(i)
         endif
      end do
   endif

   ! Store ditribution communications
   DistMrk%SendNcomm  =  SendNcomm
   DistMrk%SendList   => SendList
   DistMrk%SendOffset => SendOffset
   DistMrk%SendComm   => SendComm
   DistMrk%SendStatus => SendStatus
   DistMrk%RecvNcomm  =  RecvNcomm
   DistMrk%RecvList   => RecvList
   DistMrk%RecvOffset => RecvOffset
   DistMrk%RecvComm   => RecvComm
   DistMrk%RecvStatus => RecvStatus

   ! Start communicating list
   call StartDistList(1,DistMrk)

   ! Copy own data into the GS data structure
   nnmrk = own_list%nterm + RecvNterm
   allocate(indglomrk(nnmrk),stat=ierr)
   if (ierr .ne. 0) then
      write(*,'(a,i5)') "*** ERROR in DistMrkInd : allocating indglomrk",rankMPI
      stop
   else
      terms => own_list%terms
      do i = 1,own_list%nterm
         indglomrk(i) = mrkglo(terms(i))
      end do
      deallocate(own_list%terms,stat=ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in DistMrkInd : deallocating own_list",rankMPI
         stop
      endif
   endif
 
   ! Wait end of Recv communications
   call WaitRecvComm(DistMrk)

   ! Copy received list into the GS data structure
   istart = own_list%nterm
   do i = 1,RecvNterm
      indglomrk(istart+i) = RecvList(i)
   end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a)') "indglomrk"
   write(1000000+rankMPI,'(i10)') nnmrk
   do i = 1,nnmrk
      write(1000000+rankMPI,'(2i10)') i,indglomrk(i)
   end do
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Wait end of send communications
   call WaitSendComm(DistMrk)
 
   end subroutine DistMrkInd

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: DistMrkValR
   !
   !> @brief Distributes real values for lagrangian nodes
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !     
   !> @date June 2021
   !--------------------------------------------------------------------------------------
   subroutine DistMrkValR(n1,n2,valI,valO)

   implicit none

   ! Input variables
   integer,           intent(in) :: n1
   integer,           intent(in) :: n2
   real(kind=double), intent(in) :: valI(n1,n2,Ranknnmrk)
    
   ! Output variables
   real(kind=double), intent(out) :: valO(n1,n2,nnmrk)
    
   ! Local variables
   ! - local 
   integer                        :: ierr
   integer                        :: i,ii,j,j1,j2,k,kk,kkk
   integer                        :: istr,iend
   real(kind=double), allocatable :: valscr(:)
#if defined(DEBUG)
   integer :: errMPI
#endif
   ! - infoMPI
   integer :: rankMPI,sizeMPI
   ! - list communication
   integer                        :: nown
   integer                        :: dest
   integer                        :: SendNterm,SendNcomm
   integer,           pointer     :: SendComm(:),SendOffset(:)
   integer                        :: RecvNterm,RecvNcomm
   integer,           pointer     :: RecvOffset(:)
   real(kind=double), pointer     :: SendSCR(:) => null()
   integer,           allocatable :: nterm(:)

   ! Set handles
   rankMPI       =  infoMPI%rankMPI
   sizeMPI       =  infoMPI%sizeMPI
   SendNcomm     =  DistMrk%SendNcomm
   SendComm      => DistMrk%SendComm
   SendOffset    => DistMrk%SendOffset
   RecvNcomm     =  DistMrk%RecvNcomm
   RecvOffset    => DistMrk%RecvOffset

   ! No communications
   if ( sizeMPI .eq. 1 .and. Ranknnmrk .ne. nnmrk ) then
      write(*,'(a,i5)') "*** ERROR in DistMrkValR : sizeMPI == 1 and Ranknnmrk /= nnmrk ",rankMPI
      stop
   endif
   if ( sizeMPI .eq. 1 ) then
      do i = 1,nnmrk
         do j2 = 1,n2
            do j1 = 1,n1
               valO(j1,j2,i) = valI(j1,j2,i) 
            end do
         end do
      end do
      return
   endif

   ! Allocate receiving scratch
   allocate(valscr(n1*n2*nnmrk),stat=ierr)
   if (ierr .ne. 0) then
      write(*,'(a,i5)') "*** ERROR in DistMrkValR : allocating valscr",rankMPI
      stop
   endif

   ! Allocate sending scratch
   SendNterm = 0
   if ( SendNcomm .gt. 0 ) then
      SendNterm = SendOffset(SendNcomm+1)-SendOffset(1)
      allocate(SendSCR(n1*n2*SendNterm),nterm(SendNcomm),stat=ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in DistMrkValR : allocating sending scratches",rankMPI
         return
      endif
   endif

   ! Loop over nodes
   nown = 0
   do i = 1,SendNcomm
      nterm(i) = 0
   end do
   do i = 1,Ranknnmrk

      ! Retrieve destination ranks
      dest = 0
      j = 1
      do while ( ptnmrk(j) .le. mrkglo(i) )
         dest = dest + 1
         j = j + 1
      end do

      ! Update lists
      if ( dest .ne. rankMPI ) then
         ! not owned node
         ! store value to communicate
         k = binsearch(dest,SendNcomm,SendComm)
         nterm(k) = nterm(k) + 1
         kk = (SendOffset(k)-1+nterm(k)-1)*n1*n2
         kkk = 0
         do j2 = 1,n2
            do j1 = 1,n1
               kkk = kkk + 1
               SendSCR(kk+kkk) = valI(j1,j2,i)
            end do
         end do
      else
         ! owned node
         ! store owned value
         nown = nown + 1
         kk = (nown-1)*n1*n2
         kkk = 0
         do j2 = 1,n2
            do j1 = 1,n1
               kkk = kkk + 1
               valscr(kk+kkk) = valI(j1,j2,i)
            end do
         end do
      endif

   ! end loop over nodes
   end do

   ! Start communicationg values
   RecvNterm = 0
   if ( RecvNcomm .gt. 0 ) then
      RecvNterm = RecvOffset(RecvNcomm+1)-RecvOffset(1)
   endif
   istr = nown*n1*n2+1
   iend = n1*n2*nnmrk 
   call StartDistValueR(n1*n2,SendNterm,SendSCR,RecvNterm,valscr(istr:iend),DistMrk)

   ! Wait end of communications and remove sending buffer
   call WaitRecvComm(DistMrk)
   call WaitSendComm(DistMrk)
   if ( SendNcomm .gt. 0 ) then
      deallocate(SendSCR,nterm,stat=ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in DistMrkValR : deallocating sending scratches",rankMPI
         stop
      endif
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a)') "valscr"
   write(1000000+rankMPI,'(i10)') nnmrk
   do i = 1,nnmrk
      kk = (i-1)*n1*n2
      kkk = 0
      do j2 = 1,n1
         do j1 = 1,n1
            kkk = kkk + 1
            write(1000000+rankMPI,'(3i10,e25.11)') i,j2,j1,valscr(kk+kkk)
         end do
      end do
   end do
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Sorting
   do i = 1,nnmrk
      ii = indglomrk(i) - ptnmrk(rankMPI) + 1
      kk = (i-1)*n1*n2
      kkk = 0
      do j2 = 1,n2
         do j1 = 1,n1
            kkk = kkk + 1
            valO(j1,j2,ii) = valscr(kk+kkk) 
         end do
      end do
   end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   write(1000000+rankMPI,*)
   write(1000000+rankMPI,'(a)') "valO"
   write(1000000+rankMPI,'(i10)') nnmrk
   do i = 1,nnmrk
      do j2 = 1,n2
         do j1 = 1,n1
            write(1000000+rankMPI,'(3i10,e25.11)') i,j2,j1,valO(j1,j2,i)
         end do
      end do
   end do
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(DEBUG)
   do ii = 1,sizeMPI
      if ( ii .eq. rankMPI) then
         ! open file
         if ( rankMPI .eq. 1 ) then
            open(1000000,file="DEBUG",status='replace')
         else
            open(1000000,file="DEBUG",status='old',position='append')
         endif
         do i = 1,nnmrk
            do j2 = 1,n2
               do j1 = 1,n1
!                 write(1000000,'(5i10,e25.11)') i+ptnmrk(ii)-1,ii,i,j2,j1,valO(j1,j2,i) 
                  write(1000000,'(3i10,e25.11)') i+ptnmrk(ii)-1,j2,j1,valO(j1,j2,i) 
               end do
            end do
         end do
         close(1000000)
      endif
      call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   end do
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Deallocate receiving scratch
   deallocate(valscr,stat=ierr)
   if (ierr .ne. 0) then
      write(*,'(a,i5)') "*** ERROR in DistMrkValR : deallocating valscr",rankMPI
      stop
   endif

   end subroutine DistMrkValR

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: DistMrkValI
   !
   !> @brief Distributes integer values for lagrangian nodes
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !     
   !> @date June 2021
   !--------------------------------------------------------------------------------------
   subroutine DistMrkValI(n1,n2,valI,valO)

   implicit none

   ! Input variables
   integer, intent(in) :: n1
   integer, intent(in) :: n2
   integer, intent(in) :: valI(n1,n2,Ranknnmrk)
    
   ! Output variables
   integer, intent(out) :: valO(n1,n2,nnmrk)
    
   ! Local variables
   ! - local 
   integer              :: ierr
   integer              :: i,ii,j,j1,j2,k,kk,kkk
   integer              :: istr,iend
   integer, allocatable :: valscr(:)
   ! - infoMPI
   integer :: rankMPI,sizeMPI
   ! - list communication
   integer              :: nown
   integer              :: dest
   integer              :: SendNterm,SendNcomm
   integer, pointer     :: SendComm(:),SendOffset(:)
   integer              :: RecvNterm,RecvNcomm
   integer, pointer     :: RecvOffset(:)
   integer, pointer     :: SendSCR(:) => null()
   integer, allocatable :: nterm(:)

   ! Set handles
   rankMPI       =  infoMPI%rankMPI
   sizeMPI       =  infoMPI%sizeMPI
   SendNcomm     =  DistMrk%SendNcomm
   SendComm      => DistMrk%SendComm
   SendOffset    => DistMrk%SendOffset
   RecvNcomm     =  DistMrk%RecvNcomm
   RecvOffset    => DistMrk%RecvOffset

   ! No communications
   if ( sizeMPI .eq. 1 .and. Ranknnmrk .ne. nnmrk ) then
      write(*,'(a,i5)') "*** ERROR in DistMrkValI : sizeMPI == 1 and Ranknnmrk /= nnmrk ",rankMPI
      stop
   endif
   if ( sizeMPI .eq. 1 ) then
      do i = 1,nnmrk
         do j2 = 1,n2
            do j1 = 1,n1
               valO(j1,j2,i) = valI(j1,j2,i) 
            end do
         end do
      end do
      return
   endif

   ! Allocate receiving scratch
   allocate(valscr(n1*n2*nnmrk),stat=ierr)
   if (ierr .ne. 0) then
      write(*,'(a,i5)') "*** ERROR in DistMrkValI : allocating valscr",rankMPI
      stop
   endif

   ! Allocate sending scratch
   SendNterm = 0
   if ( SendNcomm .gt. 0 ) then
      SendNterm = SendOffset(SendNcomm+1)-SendOffset(1)
      allocate(SendSCR(n1*n2*SendNterm),nterm(SendNcomm),stat=ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in DistMrkValI : allocating sending scratches",rankMPI
         return
      endif
   endif

   ! Loop over nodes
   nown = 0
   do i = 1,SendNcomm
      nterm(i) = 0
   end do
   do i = 1,Ranknnmrk

      ! Retrieve destination ranks
      dest = 0
      j = 1
      do while ( ptnmrk(j) .le. mrkglo(i) )
         dest = dest + 1
         j = j + 1
      end do

      ! Update lists
      if ( dest .ne. rankMPI ) then
         ! not owned node
         ! store value to communicate
         k = binsearch(dest,SendNcomm,SendComm)
         nterm(k) = nterm(k) + 1
         kk = (SendOffset(k)-1+nterm(k)-1)*n1*n2
         kkk = 0
         do j2 = 1,n2
            do j1 = 1,n1
               kkk = kkk + 1
               SendSCR(kk+kkk) = valI(j1,j2,i)
            end do
         end do
      else
         ! owned node
         ! store owned value
         nown = nown + 1
         kk = (nown-1)*n1*n2
         kkk = 0
         do j2 = 1,n2
            do j1 = 1,n1
               kkk = kkk + 1
               valscr(kk+kkk) = valI(j1,j2,i)
            end do
         end do
      endif

   ! end loop over nodes
   end do

   ! Start communicationg values
   RecvNterm = 0
   if ( RecvNcomm .gt. 0 ) then
      RecvNterm = RecvOffset(RecvNcomm+1)-RecvOffset(1)
   endif
   istr = nown*n1*n2+1
   iend = n1*n2*nnmrk 
   call StartDistValueI(n1*n2,SendNterm,SendSCR,RecvNterm,valscr(istr:iend),DistMrk)

   ! Wait end of communications and remove sending buffer
   call WaitRecvComm(DistMrk)
   call WaitSendComm(DistMrk)
   if ( SendNcomm .gt. 0 ) then
      deallocate(SendSCR,nterm,stat=ierr)
      if( ierr .ne. 0 ) then
         write(*,'(a,i5)') "*** ERROR in DistMrkValI : deallocating sending scratches",rankMPI
         stop
      endif
   endif

   ! Sorting
   do i = 1,nnmrk
      ii = indglomrk(i) - ptnmrk(rankMPI) + 1
      kk = (i-1)*n1*n2
      kkk = 0
      do j2 = 1,n2
         do j1 = 1,n1
            kkk = kkk + 1
            valO(j1,j2,ii) = valscr(kk+kkk) 
         end do
      end do
   end do

   ! Deallocate receiving scratch
   deallocate(valscr,stat=ierr)
   if (ierr .ne. 0) then
      write(*,'(a,i5)') "*** ERROR in DistMrkValI : deallocating valscr",rankMPI
      stop
   endif

   end subroutine DistMrkValI

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: DistInfo
   !
   !> @brief Computes distribution info
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------

   subroutine DistInfo(sizeMPI,nUnk,res,blksize)

   implicit none

   ! Input variables
   integer, intent(in) :: sizeMPI
   integer, intent(in) :: nUnk

   ! Output variables
   integer, intent(out) :: res
   integer, intent(out) :: blksize

   res = mod(nUnk,sizeMPI)
   blksize = nUnk/sizeMPI + 1

   end subroutine DistInfo

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: SetUnk
   !
   !> @brief Sets pointer to the first unknown for each MPI process
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------

   subroutine SetUnk(sizeMPI,nUnk,ptUnk)

   implicit none

   ! Input variables
   integer, intent(in) :: sizeMPI
   integer, intent(in) :: nUnk

   ! Output variables
   integer, intent(out) :: ptUnk(sizeMPI+1)

   ! Local variables
   integer :: i
   integer :: res,blksize

   ! Retrieve distribution info
   call DistInfo(sizeMPI,nUnk,res,blksize)

   ! Set pointers to the first unknown
   ptUnk(1) = 1
   do i = 1,res
      ptUnk(i+1) = ptUnk(i) + blksize
   end do
   blksize = blksize - 1
   do i = res+1,sizeMPI
     ptUnk(i+1) = ptUnk(i) + blksize
   end do

   end subroutine SetUnk

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! SUBROUTINE: iglo2ix 
   !
   !> @brief Retrieves local axis indeces from global nodes
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------

   subroutine iglo2ix(iglo,nx1,nx2,ix1,ix2,ix3) 

   implicit none

   ! Input variables
   integer, intent(in) :: iglo
   integer, intent(in) :: nx1,nx2

   ! Output variables
   integer, intent(out) :: ix1,ix2,ix3

   ! Local variables
   integer :: itmp

   ! Retrieve local indices
   ix3  = (iglo-1)/(nx2*nx1) + 1

   itmp = iglo - (ix3-1)*nx2*nx1

   ix1  = (itmp-1)/nx2 + 1

   ix2  = itmp - (ix1-1)*nx2

   end subroutine iglo2ix 

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! FUNCION: resizei 
   !
   !> @brief Resizes an integer array
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------
   function resizei(n,vec) result(info)

   implicit none

   !Input variables
   integer,          intent(in)    :: n
   integer, pointer, intent(inout) :: vec(:)

   ! Output variables
   integer :: info

   ! Local variables
   integer          :: n_old,i
   integer, pointer :: tmp(:),tmp1(:)

   info = 0
   n_old = size(vec)
   if (n_old .ne. n) then
      ! Allocate a temporary array
      allocate(tmp(n),stat=info)
      if (info.ne.0) then
         return
      endif
      ! Swap pointers
      tmp1 => vec
      vec  => tmp
      ! Copy old entries in the new array
      do i = 1,min(n,n_old)
         vec(i) = tmp1(i)
      end do
      ! Deallocate old array
      deallocate(tmp1,stat=info)
      tmp => null()
   endif

   end function resizei

!-----------------------------------------------------------------------------------------

   !--------------------------------------------------------------------------------------
   ! FUNCION: binsearch
   !
   !> @brief Retrieves the position of a term in a sorted list
   !
   !> @author Giovanni Isotton
   !
   !> @version 1.0
   !
   !> @date June 2021
   !--------------------------------------------------------------------------------------
   function binsearch(ii,n,v) result(ipos)

   implicit none

   integer, intent(in) :: ii,n,v(n)
   integer             :: ipos
   integer             :: istart,iend,imid

   ! Initialization
   ipos = 0
   istart = 1
   iend = n

   ! Check if ii is into the bounds
   if ((ii .lt. v(istart)) .or. (ii .gt. v(iend)) ) return

   do while (iend - istart .gt. 1)
      imid = (istart + iend) / 2
      if (v(imid) .gt. ii) then
         iend = imid
      else
         istart = imid
      endif
   enddo

   ! Check last two positions
   if (v(istart) .eq. ii) ipos = istart
   if (v(iend)   .eq. ii) ipos = iend

   end function binsearch

!-----------------------------------------------------------------------------------------

end module class_DistGrid
