PROGRAM fazi2stats   
 
IMPLICIT NONE

CHARACTER(len=32) :: arg
INTEGER :: ngood,i
DOUBLE PRECISION :: mean_fazi,mean_dt,mean_dazi,std_fazi,std_dt,std_dazi,meansin,meancos,da                    
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x,fazi,dt,sin_fazi,cos_fazi,dazi
DOUBLE PRECISION, PARAMETER :: pi180=0.017453292519943295769236907684886 
DO i = 1, iargc()
  CALL getarg(i, arg)
       if(i==1) read(arg,'(i3)') ngood  
!       if(i==2) read(arg,'(e3.1)') z  
END DO

!No record if 0 good measuraments
if(ngood .eq.0) then
  !write(*,*),0d0,0d0,0d0,0d0
  stop
!No statistics if only 1 good measurament
else if(ngood .eq. 1) then
  read *,i,mean_fazi,mean_dt
  write(*,*), mean_fazi,0d0,mean_dt,0d0
  stop         

else

ALLOCATE(x(ngood),fazi(ngood),dt(ngood),sin_fazi(ngood),cos_fazi(ngood),dazi(ngood))

do i=1,ngood
   read *,x(i),fazi(i),dt(i)
end do

!Convert degrees to radians and calculate sin and cos of fast azimuths
do i=1,ngood
   sin_fazi(i)=sin(2d0*fazi(i)*pi180)
   cos_fazi(i)=cos(2d0*fazi(i)*pi180)
end do

!Mean of sin(2*azi),cos(2*azi)
meansin=sum(sin_fazi)/ngood   
meancos=sum(cos_fazi)/ngood   
!Mean 2*phi azimuth in radians
mean_fazi = atan2(meansin,meancos)
!Convert phi to degrees
!do i=1,ngood
!   if(fazi(i) .gt. 90d0) fazi(i) = fazi(i) - 180d0
!end do
!mean_fazi = sum(fazi)/ngood!mean_fazi/2d0/pi180
mean_fazi = mean_fazi/2d0/pi180
!print *,mean_fazi
!Move to 0...180 range
call fix_deg_angle(mean_fazi)
if(mean_fazi .gt. 180) mean_fazi = mean_fazi -180

!Deviation from mean azimuth ( -90 < dazi < +90 )
do i=1,ngood
   call fix_deg_angle(fazi(i))
   da = fazi(i) - mean_fazi
!print *,da,fazi(i),mean_fazi
   if(da .gt. 0d0) then
      if(da .gt. 180d0) da = da - 180d0
      if(da .gt.  90d0) da = da - 180d0
   else
      if(da .lt.- 180d0) da = da + 180d0
      if(da .lt.  -90d0) da = da + 180d0
   end if
   dazi(i) = da
!print *,dazi(i)
end do

!Mean, std dev for diff_azi
!Mean of diff_azi = 0 --> stddev(dazi) = stddev(fazi)
call stats(dazi,ngood,mean_dazi,std_dazi)
std_fazi=std_dazi
!print *,mean_dazi
!Mean, std dev for dt        
call stats(dt,ngood,mean_dt,std_dt)

!Write result to output file`
write(*,*),mean_fazi,std_fazi,mean_dt,std_dt

DEALLOCATE(x,fazi,dt,sin_fazi,cos_fazi,dazi)

stop

end if


END PROGRAM


subroutine fix_deg_angle(phi)

implicit none

DOUBLE PRECISION :: phi 

do while(phi < 0 ) 
   phi = phi + 360
end do
do while(phi >360) 
   phi = phi - 360
end do

end subroutine fix_deg_angle

subroutine stats(x,ngood,ave,sdev)

implicit none

INTEGER :: i,ngood,n
DOUBLE PRECISION :: ave,sdev,var,x(ngood),s
 
s=0d0
n=0
do i=1,ngood
   if(isnan(x(i))) then
     write(*,*) 'NaN measurament for i=', i
   else
      s = s + x(i)
      n = n + 1
   end if
end do

ave=s/n

s = 0d0; sdev=0d0; var=0d0
do i=1,ngood
   if(isnan(x(i))) then
     !write(*,*) 'NaN measurament for i=', i
   else
      s = x(i)-ave
      !dev = dev + abs(s)
      var = var + s*s
   end if
end do

var = var/(n-1)
sdev = sqrt(var)

end subroutine stats
