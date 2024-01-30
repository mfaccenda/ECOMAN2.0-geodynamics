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
!!! Module of common variables                                             !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   MODULE comvar

   INTEGER :: num_layers,seismo_station_num,maxlayers=5000      
   INTEGER :: depthaxis,cartspher,dimensions,yinyang
   DOUBLE PRECISION :: dt,timesum,maxdist,minlayer,Xscale,pi,deg2rad
   CHARACTER (3) :: dt_to_str
   CHARACTER (4) :: dt_str

!!! Seismic station grid points and mantle domain

   INTEGER :: nsx1, nsx3
   DOUBLE PRECISION :: x1size, x2size, x3size 
   DOUBLE PRECISION :: stepx1, stepx3 
   DOUBLE PRECISION :: i1first,i2first,i3first
   DOUBLE PRECISION :: i1last, i2last ,i3last
   !! domain size in X1, X2 and X3 direction

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X1,X1s,X3,X3s
   ! Seismic stations coordinates 

!!! Tracers grid   

   INTEGER :: marknum !! number of markers
   INTEGER, DIMENSION(:), ALLOCATABLE :: mYY
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mx1
   ! X grid points coordinates 
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mx2
   ! Y grid points coordinates 
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mx3
   ! Z grid points coordinates 
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rhom
   ! density                    

   INTEGER, DIMENSION(3,3) :: ijkl ! tensor of indices to form Cijkl from Sij
   INTEGER, DIMENSION(6) :: l1,l2 ! tensot of indices to form Sij from Cijkl
 
!!! Elastic tensor 

   DOUBLE PRECISION, DIMENSION (:,:,:), ALLOCATABLE :: Sav,Savtmp,Savtmp2

   ! stiffness matrix

   END MODULE comvar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program read_unformatted_file

USE comvar
USE hdf5   

IMPLICIT NONE

CHARACTER(30) :: filename,arg
CHARACTER(1000) :: dummy
INTEGER :: i,j,k,i1,i2,i3,m,mode,dum_int(3)
DOUBLE PRECISION :: buf(40)
DOUBLE PRECISION :: Ax,Ay,Az,Bx,By,Bz
DOUBLE PRECISION :: datadb(3),lr,cr
DOUBLE PRECISION ,DIMENSION(:,:), ALLOCATABLE :: Xsave
INTEGER(HID_T)  :: file_id, group_id, dataspace_id, dataset_id, attr_id,dcpl,memtype !Handles
INTEGER     ::   error  ! Error flag

DO i = 1, iargc()
   CALL getarg(i, arg)
    if(i==1) read(arg,"(i1)") mode
    if(i==2) dt_str = arg
END DO

!!! Open input file stack_input.dat

OPEN(15,file="stack_input.dat")

!!! Read parameters for splitting
!Input files directory
write(*,*)
write(*,"(a)"),'********************************************************'
write(*,"(a)"),'********************************************************'
write(*,*)
write(*,"(a,a)"),' READING INPUT FILE stack_input.dat '
write(*,*)
call read_par_int(15,nsx1)
call read_par_int(15,nsx3)
call read_par_int(15,depthaxis)
!WARNING: longitude is set to 1, colatitude is set to 3, while radius is set to 2
call read_par_double(15,i1first)
call read_par_double(15,i1last)
call read_par_double(15,i2first)
call read_par_double(15,i2last)
call read_par_double(15,i3first)
call read_par_double(15,i3last)
call read_par_double(15,maxdist)
call read_par_double(15,minlayer)
call read_par_double(15,Xscale)

CLOSE(15)

pi = 3.141592653589793238462643383279

!Reading Cijkl file
filename = 'Cijkl'//dt_str//'.h5'
print *
print *,'Input file : ',filename
print *

!Initialize FORTRAN interface.

CALL H5open_f (error)

! Create a new file using default properties.

CALL H5Fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

!Load domain infos      
CALL loadsave_integer(1,1,file_id,2,H5T_NATIVE_INTEGER,dum_int(1:2),'CoordinateSystem',0)
dimensions = dum_int(1)
cartspher  = dum_int(2)
CALL loadsave_integer(1,1,file_id,3,H5T_NATIVE_INTEGER,dum_int,'Periodic',0)

yinyang = 1
if(cartspher == 2 .AND. dimensions == 3 .AND. dum_int(1) > 0 .AND. dum_int(3) > 0) yinyang = 2

CALL H5Gopen_f(file_id, "/Particles", group_id, error)

CALL loadsave_double(1,1,group_id,1,H5T_NATIVE_DOUBLE,datadb(1),'marknum',0)
marknum = datadb(1)

print *,'Marknum = ',marknum
print *

ALLOCATE(mx1(marknum),mx2(marknum),mx3(marknum),rhom(marknum),Sav(6,6,marknum),Xsave(21,marknum))

CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mx1,'mx1',0)
CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mx2,'mx2',0)
CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,mx3,'mx3',0)
CALL loadsave_double(0,1,group_id,marknum,H5T_NATIVE_DOUBLE,rhom,'rho',0)

IF(yinyang == 2) THEN

   ALLOCATE(mYY(marknum))
   CALL loadsave_integer(0,1,group_id,marknum,H5T_NATIVE_INTEGER,mYY(1:marknum),'mYY',0)

   DO m = 1, marknum
      IF(mYY(m)==2) THEN
         CALL yin2yang(mx1(m),mx3(m),lr,cr)
         mx1(m) = lr
         mx3(m) = cr
      END IF
      IF(mx1(m)<0   ) mx1(m) = 2*pi + mx1(m)
      IF(mx1(m)>2*pi) mx1(m) = mx1(m) - 2*pi
      IF(mx3(m)<0   ) mx3(m) = - mx3(m)
      IF(mx3(m)>pi  ) mx3(m) = 2*pi - mx3(m)
   END DO

END IF

!Sav(GPa)
dum_int(1)=21
dum_int(2)=marknum
CALL loadsave_double(0,2,group_id,dum_int(1:2),H5T_NATIVE_DOUBLE,Xsave,'Sav',0)

DO m = 1, marknum

    !Save elastic tensor: only 21 elastic moduli
    Sav(1,1,m) = Xsave(1,m)
    Sav(2,2,m) = Xsave(2,m)
    Sav(3,3,m) = Xsave(3,m)
    Sav(2,3,m) = Xsave(4,m)
    Sav(3,2,m) = Xsave(4,m)
    Sav(1,3,m) = Xsave(5,m)
    Sav(3,1,m) = Xsave(5,m)
    Sav(1,2,m) = Xsave(6,m)
    Sav(2,1,m) = Xsave(6,m)
    Sav(4,4,m) = Xsave(7,m)
    Sav(5,5,m) = Xsave(8,m)
    Sav(6,6,m) = Xsave(9,m)
    Sav(1,4,m) = Xsave(10,m)
    Sav(4,1,m) = Xsave(10,m)
    Sav(2,5,m) = Xsave(11,m)
    Sav(5,2,m) = Xsave(11,m)
    Sav(3,6,m) = Xsave(12,m)
    Sav(6,3,m) = Xsave(12,m)
    Sav(3,4,m) = Xsave(13,m)
    Sav(4,3,m) = Xsave(13,m)
    Sav(1,5,m) = Xsave(14,m)
    Sav(5,1,m) = Xsave(14,m)
    Sav(2,6,m) = Xsave(15,m)
    Sav(6,2,m) = Xsave(15,m)
    Sav(2,4,m) = Xsave(16,m)
    Sav(4,2,m) = Xsave(16,m)
    Sav(3,5,m) = Xsave(17,m)
    Sav(5,3,m) = Xsave(17,m)
    Sav(1,6,m) = Xsave(18,m)
    Sav(6,1,m) = Xsave(18,m)
    Sav(5,6,m) = Xsave(19,m)
    Sav(6,5,m) = Xsave(19,m)
    Sav(4,6,m) = Xsave(20,m)
    Sav(6,4,m) = Xsave(20,m)
    Sav(4,5,m) = Xsave(21,m)
    Sav(5,4,m) = Xsave(21,m)

    !Convert density from kg/m^3 to g/cm^3
    rhom(m) = rhom(m)/1000.0 

    !Cij (GPa) / rho (g/cm^3)
    Sav(:,:,m) = Sav(:,:,m) / rhom(m)

END DO

DEALLOCATE(Xsave)

CALL H5Gclose_f(group_id, error)

!Terminate access to the file.
CALL H5Fclose_f(file_id, error)

CALL H5close_f(error)
!Close FORTRAN interface.

IF(cartspher==2) THEN
   deg2rad = pi/180.0
   i1first = i1first*deg2rad
   i1last  = i1last*deg2rad
   i3first = i3first*deg2rad
   i3last  = i3last*deg2rad
END IF

!Define grid of virtual seismic stations        
IF(dimensions == 2) nsx3 = 1
ALLOCATE(X1(nsx1),X3(nsx3))

X1 = 0d0

x1size=i1last-i1first
X1(1)=i1first;
stepx1 = x1size/nsx1

!Seismic stations placed inside the box by stepx/2
X1(1)=X1(1)+stepx1/2;
DO i1 = 2, nsx1
   X1(i1) = X1(i1-1) + stepx1
END DO

x2size = i2last - i2first

if(dimensions == 2) THEN
 
   IF(cartspher == 1) THEN
      i3first = 0d0
      i3last  = 0d0
      mx3     = 0d0   
      X3(1)   = 0d0
   ELSE
      i3first = pi/2d0
      i3last  = pi/2d0
      mx3     = pi/2d0   
      X3(1)   = pi/2d0
   END IF
END IF

IF(dimensions == 3) THEN

   X3 = 0d0   

   x3size=i3last-i3first
   X3(1)=i3first 
   stepx3 = x3size/nsx3
   X3(1)=X3(1)+stepx3/2 
   DO i3 = 2, nsx3
      X3(i3) = X3(i3-1) + stepx3
   END DO

END IF

!Number of seismic stations at the surface
seismo_station_num=nsx1*nsx3

!Print grid infos
if(cartspher==1) then
   write(*,'(a,i0,a,i0)',advance='no') ' Grid with number of seismic stations: X = ',nsx1,', Z = ',nsx3
   write(*,'(a,i0)') ', Number of seismic stations = ',seismo_station_num 
   print *
   write(*,'(a)') ' Domain size is:'
   write(*,'(a,f12.4,a,f12.4,a)') ' X min = ',i1first/1d3,'; X max = ',i1last/1d3,' (km)'
   write(*,'(a,f12.4,a,f12.4,a)') ' Y min = ',i2first/1d3,'; Y max = ',i2last/1d3,' (km)'
   write(*,'(a,f12.4,a,f12.4,a)') ' Z min = ',i3first/1d3,'; Z max = ',i3last/1d3,' (km)'
   print *
end if
if(cartspher==2) then
   write(*,'(a,i0,a,i0)',advance='no') ' Grid with number of seismic stations : Long = ',nsx1,', Colat = ',nsx3
   write(*,'(a,i0)') ', Number of seismic stations = ',seismo_station_num 
   print *
   write(*,'(a)') ' Domain size is:'
   write(*,'(a,f10.2,a,f10.2,a)') ' Long   min = ',i1first/deg2rad,'; Long   max = ',i1last/deg2rad,' (deg)'
   write(*,'(a,f10.2,a,f10.2,a)') ' Radius min = ',i2first/1d3,    '; Radius max = ',i2last/1d3,' (km)'
   write(*,'(a,f10.2,a,f10.2,a)') ' Colat  min = ',i3first/deg2rad,'; Colat  max = ',i3last/deg2rad,' (deg)'
   print *
end if

write(*,*)
write(*,"(a)"),'********************************************************'
write(*,"(a)"),'********************************************************'
write(*,*)

!Define tensors of indices

ijkl(1,1) = 1 ; ijkl(1,2) = 6 ; ijkl(1,3) = 5
ijkl(2,1) = 6 ; ijkl(2,2) = 2 ; ijkl(2,3) = 4
ijkl(3,1) = 5 ; ijkl(3,2) = 4 ; ijkl(3,3) = 3

l1(1) = 1 ; l1(2) = 2 ; l1(3) = 3
l1(4) = 2 ; l1(5) = 3 ; l1(6) = 1
l2(1) = 1 ; l2(2) = 2 ; l2(3) = 3
l2(4) = 3 ; l2(5) = 1 ; l2(6) = 2
 
CALL mark2node 

STOP

end program read_unformatted_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine interp, bilinear interpolation of a function                !!!
!!! following the method given in Numerical Recipies, Press et al., p96    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE mark2node

USE comvar

IMPLICIT NONE

CHARACTER(1) :: idx1
CHARACTER(2) :: idx2
CHARACTER(3) :: idx3
CHARACTER(4) :: idx4
CHARACTER(40) :: str
INTEGER :: i1,i2,i3,i10,i20,i30,i,j,k,l,m,gi,x,y,z,s,np
INTEGER, DIMENSION(:), ALLOCATABLE :: map_idx0,map_idx
! indices of the UP-LEFT grid point closest to the considered point

DOUBLE PRECISION :: thetadist,phidist,centralangle,greatcircledist,t1,t2,dphi,m1,m3,depth1

INTEGER, DIMENSION(:), ALLOCATABLE :: depth2
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rd,rd0,depth0,depth,my,wt,rhotmp      
  
ALLOCATE(map_idx0(maxlayers),map_idx(maxlayers),rd0(maxlayers),rd(maxlayers))
ALLOCATE(depth0(maxlayers),depth(maxlayers),depth2(maxlayers),my(maxlayers))
ALLOCATE(Savtmp(6,6,maxlayers),rhotmp(maxlayers),wt(maxlayers))
ALLOCATE(X1s(seismo_station_num),X3s(seismo_station_num))

!Save surface coordinates
X1s = 0d0 ; X3s = 0d0

do i3=1,nsx3               
do i1=1,nsx1               
   s = i1 +(i3-1)*nsx1
   X1s(s)=X1(i1)
   X3s(s)=X3(i3)
end do
end do

open(88,file='Surface_coordinates.dat',status='replace')
write(88,"(10000(1PE12.4))") X1s
write(88,"(10000(1PE12.4))") X3s
close(88)

!Vertical profiles of Sav 
do i3=1,nsx3               
do i1=1,nsx1               
   s = i1 +(i3-1)*nsx1
   i = 0 ; np = 0
   Savtmp = 0d0; rd0 = 0d0; my = 0d0; map_idx0 = 0; rhotmp = 0

   !Find near markers
   DO m=1,marknum
      !Check if marker is inside vertical profile grid
      IF(isnan(Sav(1,1,m))==0) THEN
      IF(mx1(m) >= i1first .and. mx1(m) <= i1last) THEN
      IF(mx2(m) >= i2first .and. mx2(m) <= i2last) THEN
      IF(mx3(m) >= i3first .and. mx3(m) <= i3last) THEN

      !Check if marker is close to the vertical profile of the s seismic station
      m1=mx1(m) ; m3=mx3(m) ;  depth1=mx2(m)
      IF(cartspher == 1) THEN
         greatcircledist = ((m1-X1(i1))**2 + (m3-X3(i3))**2)**0.5d0
      ELSE
         !Convert colatitude into latitude
         t1 = m1 - pi/2.0
         t2 = X1s(s) - pi/2.0
         dphi = abs(m3-X3s(s)) !! this is 0 in 2D
         !Central angle
         centralangle = acos(sin(t1)*sin(t2) + cos(t1)*cos(t2)*cos(dphi))
         !Arc distance in meters along greatcircle: centralangle*radius
         greatcircledist = depth1*abs(centralangle)
      END IF

      IF(greatcircledist <= maxdist) THEN   
    
         i=i+1

         !Scaled great circle distance 
         rd0(i) = 1d0 - greatcircledist/maxdist
         !Check !!!
         if(rd0(i) < 0d0 .or. rd0(i) > 1d0) print *,'Radial distance > maximum distance for marker',m,'at',mx1(m),mx2(m),mx3(m)
 
         !Find aggregate depth relative to max. vertical extension of the mantle domain
         IF(depthaxis == 0) my(i)=ABS(i2first - depth1)
         IF(depthaxis == 1) my(i)=ABS(i2last  - depth1)
         map_idx0(i)=m 
   
      END IF
      
      END IF;END IF;END IF;END IF

   END DO

   IF(i==0) THEN
      print *,'No markers below seismic station',s,' at ',X1(i1),X3(i3)
      goto 10
   END IF
   IF(i > maxlayers) THEN
      print *,'Number of aggregates ',i,' is bigger than ',maxlayers,' below seismic station',s,' at',X1(i1),X3(i3)
      print *,'Increase the number of variable maxlayers at line 8 and recompile savstack.f90'
      stop
   END IF

   !Number of points below seismic station
   np = i
  
   !Sorting points from lowest to highest
   depth0 = 0d0 ;  map_idx = 0d0; rd = 0d0
   do j=1,np
      depth0(j) = 1d30
      !Find lowest point
      do k = 1,np
         if (my(k) <= depth0(j)) then
            depth0(j) = my(k)
            l = k
         end if
      end do
      map_idx(j) = map_idx0(l)
      rd(j) = rd0(l)
      my(l) = 1d40
   end do   
  
   !Check for markers with vertical distance < minlayer, and assign depth,Sav
   depth = 0d0; Savtmp = 0d0; wt = 0d0; rhotmp = 0d0 
   i = np !Total number of poont below seismic station
   k = 1  !Layer counter
   j = 1  !Point counter
   depth(k) = depth0(j)*rd(j)
   Savtmp(:,:,k) = Sav(:,:,map_idx(j))*rd(j)
   rhotmp(k) = rhom(map_idx(j))*rd(j) ! in g/cm^3
   wt(k) = rd(j)
   !Loop over i=np points below seismic stations
   do j=2,i   
      if(ABS(depth0(k)-depth0(j)) <= minlayer) then
        !if(s==1) print *,'Merge',k,j,depth0(k),depth0(j)
        !Merge points
        np = np - 1
        !Depth average
        depth(k) = depth(k) + depth0(j)*rd(j)
        !Sav average
        Savtmp(:,:,k) = Savtmp(:,:,k) + Sav(:,:,map_idx(j))*rd(j)
        rhotmp(k) = rhotmp(k) + rhom(map_idx(j))*rd(j) ! in g/cm^3
        wt(k) = wt(k) + rd(j)
        !If last j point, then normalize
        if(j==i) then
           depth(k) = depth(k)/wt(k)
           Savtmp(:,:,k) = Savtmp(:,:,k)/wt(k)
           rhotmp(k) = rhotmp(k)/wt(k)
        end if 
      else    
        depth(k) = depth(k)/wt(k)
        Savtmp(:,:,k) = Savtmp(:,:,k)/wt(k)
        rhotmp(k) = rhotmp(k)/wt(k)
        k = k + 1 
        !Last point
        if(j==i) then
           depth(k) = depth0(j)
           Savtmp(:,:,k)=Sav(:,:,map_idx(j))
           rhotmp(k) = rhom(map_idx(j))
           exit  
        else
           depth0(k)=depth0(j)
           depth(k) = depth0(j)*rd(j)
           Savtmp(:,:,k)=Sav(:,:,map_idx(j))*rd(j)
           rhotmp(k)=rhom(map_idx(j))*rd(j) ! in g/cm^3
           wt(k) = rd(j)
        end if 
      end if
   end do 
   
   !Unscale depth and convert to km
   depth2=NINT(depth*Xscale)

   !Saving output
   !Write file with unscaled depth
   if(s .lt. 10) then 
      write(idx1,'(i1)') s
      str = 'seismic_station_'//idx1//'.sav'
   end if
   if(s .ge. 10 .and. s .lt. 100) then
      write(idx2,'(i2)') s
      str = 'seismic_station_'//idx2//'.sav'
   end if
   if(s .ge. 100 .and. s .lt. 1000) then
      write(idx3,'(i3)') s
      str = 'seismic_station_'//idx3//'.sav'
   end if
   if(s .ge. 1000) then 
      write(idx4,'(i4)') s
      str = 'seismic_station_'//idx4//'.sav'
   end if 
   print *,str
   open(99,file=str,status='replace')
   !Saving upper right of Sav (21 matrix elements) because of tensor symmetry
   write(99,"(1000(1pi8))") depth2(1:np)
   write(99,1) rhotmp(1:np)

   !Rotate tensor such that C22 is parallel to Y
   !This is because later anicake will rotate the tensors around 
   !axis 1 by 90° and toward axis 3, and then rotate as a
   !function of the backazimuth
   IF(cartspher == 2) CALL tensorrot(np,s)

   !Save with x1 = x; x2 = y; x3 = z
   write(99,1) Savtmp(1,1,1:np)
   write(99,1) Savtmp(1,2,1:np)
   write(99,1) Savtmp(1,3,1:np)
   write(99,1) Savtmp(1,4,1:np)
   write(99,1) Savtmp(1,5,1:np)
   write(99,1) Savtmp(1,6,1:np)
   write(99,1) Savtmp(2,2,1:np)
   write(99,1) Savtmp(2,3,1:np)
   write(99,1) Savtmp(2,4,1:np)
   write(99,1) Savtmp(2,5,1:np)
   write(99,1) Savtmp(2,6,1:np)
   write(99,1) Savtmp(3,3,1:np)
   write(99,1) Savtmp(3,4,1:np)
   write(99,1) Savtmp(3,5,1:np)
   write(99,1) Savtmp(3,6,1:np)
   write(99,1) Savtmp(4,4,1:np)
   write(99,1) Savtmp(4,5,1:np)
   write(99,1) Savtmp(4,6,1:np)
   write(99,1) Savtmp(5,5,1:np)
   write(99,1) Savtmp(5,6,1:np)
   write(99,1) Savtmp(6,6,1:np)

10 close(99)

end do
end do

DEALLOCATE(map_idx0,map_idx,rd0,rd,depth0,depth,depth2,my,Savtmp,wt)

RETURN

1 FORMAT(1000(1pe14.6))


END SUBROUTINE mark2node 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   SUBROUTINE tensorrot(np,s)

   USE comvar

   IMPLICIT NONE

   INTEGER :: np,i,j,k,ll,p,q,r,ss,m,ij,s
   DOUBLE PRECISION :: phi1,theta,phi2
   DOUBLE PRECISION, DIMENSION (3,3,3,3) :: C0,Cav,Cav2
   DOUBLE PRECISION, DIMENSION(3,3) :: acs

!!! Rotate the 'np' elastic tensors such that the radial distance below the 's' seismic
!!! station becomes vertical, i.e., parallel to the Y axis. This is because SKS splitting
!!! is calculated for vertically travelling SKS waves 

   !Rotate the elastic tensor back with respect to long (X1s(s)) and colat (X3s(s)) 
   !First rotate around Z axis, then around X axis

   !If use next line, rotation to -Y
   !then mean_fazi X and Z components are cos for X and sin for Z
   phi1= 0d0 ; theta = 0 ; phi2 = -X1s(s) + pi*3.0/2.0 
   if(dimensions == 3) theta = -X3s(s) + pi/2.0

   !Transform Euler angles into direction cosine matrix
   !Z-X-Z: Euler angles in radians
   acs(1,1)=COS(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*SIN(phi2)
   acs(2,1)=COS(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*SIN(phi2)
   acs(3,1)=SIN(phi2)*SIN(theta)

   acs(1,2)=-SIN(phi2)*COS(phi1)-COS(theta)*SIN(phi1)*COS(phi2)
   acs(2,2)=-SIN(phi2)*SIN(phi1)+COS(theta)*COS(phi1)*COS(phi2)
   acs(3,2)=COS(phi2)*SIN(theta)

   acs(1,3)=SIN(theta)*SIN(phi1)
   acs(2,3)=-SIN(theta)*COS(phi1)
   acs(3,3)=COS(theta)

   DO ij = 1 , np
   
         Cav = 0d0 ; C0 = 0d0
 
         !Convert to Cijkl 4th order tensor
         DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
            C0(i,j,k,ll) = Savtmp(ijkl(i,j),ijkl(k,ll),ij)
         END DO ; END DO ; END DO ; END DO

         DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
            DO p = 1 , 3 ; DO q = 1 , 3 ; DO r = 1 , 3 ; DO ss = 1 , 3
               Cav(i,j,k,ll) = Cav(i,j,k,ll) + acs(i,p)*acs(j,q)*acs(k,r)*acs(ll,ss)*C0(p,q,r,ss)
            END DO ; END DO ; END DO ; END DO
         END DO ; END DO ; END DO ; END DO
       
        !Convert to Voigt notation
        DO i = 1 , 6 ; DO j = 1 , 6
           Savtmp(i,j,ij) = Cav(l1(i),l2(i),l1(j),l2(j))
        END DO ; END DO

   END DO

   RETURN

   END SUBROUTINE tensorrot

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

