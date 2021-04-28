PROGRAM sav2cijkl    
 
IMPLICIT NONE

CHARACTER(len=32) :: arg,fname,str
CHARACTER(len=1) :: depth1         
CHARACTER(len=2) :: depth2         
INTEGER :: i,j,k,d,beta,aniso_scale
DOUBLE PRECISION :: rho,alpha,gamma1,Sav(6,6), cijkl(3,3,3,3)
rho = 3.353 ! density in g/cm^3
!rotation angles
alpha = 0d0
beta = 0d0
gamma1 = 0d0
Sav = 0d0
DO i = 1, iargc()
  CALL getarg(i, arg)
       if(i==1) fname=arg
       if(i==2) read(arg,'(e9.1)') rho  
       if(i==3) read(arg,"(i1)") beta  
       if(i==4) read(arg,"(i1)") aniso_scale
END DO

!WRITE (*,*) fname,depth,beta,aniso_scale

rho=3.353

OPEN(51,FILE=fname,status='old')

do i=1,6
   k = i
   do j=k,6
       
      read(51,*) Sav(i,j)
      !Sav(i,j) = Sav(i,j)/rho      
      Sav(j,i) = Sav(i,j)
   
   end do
end do

close(51)


!write(*,*) Sav
! convert to cijkl
call vera_c2c4(Sav,cijkl)
!Print to Cijkl depth 
call vera_print_cijkl_ftrn(cijkl)

END PROGRAM

!
! convert a c[6,6] matrix into a c[3,3,3,3] Cijkl
! tensor
!
subroutine vera_c2c4(c,cc)   
  implicit none
  double precision ::  c(6,6),cc(3,3,3,3)
  INTEGER :: i,j,k,l 
  
  cc = 0.d0

  cc(1,1,1,1) = c(1,1)
  cc(2,2,2,2) = c(2,2)
  cc(3,3,3,3) = c(3,3)
  cc(2,3,2,3) = c(4,4)
  cc(3,2,3,2) =cc(2,3,2,3)
  cc(2,3,3,2) =cc(2,3,2,3)
  cc(3,2,2,3) =cc(2,3,2,3)
  cc(1,3,1,3) = c(5,5)
  cc(3,1,1,3) =cc(1,3,1,3)
  cc(1,3,3,1) =cc(1,3,1,3)
  cc(3,1,3,1) =cc(1,3,1,3)
  cc(1,1,2,2) = c(1,2)
  cc(2,2,1,1) =cc(1,1,2,2)
  cc(1,1,3,3) = c(1,3)
  cc(3,3,1,1) =cc(1,1,3,3)
  cc(1,1,2,3) = c(1,4)
  cc(1,1,3,2) =cc(1,1,2,3)
  cc(2,3,1,1) =cc(1,1,2,3)
  cc(3,2,1,1) =cc(1,1,2,3)
  cc(1,1,1,3) = c(1,5)
  cc(1,1,3,1) =cc(1,1,1,3)
  cc(1,3,1,1) =cc(1,1,1,3)
  cc(3,1,1,1) =cc(1,1,1,3)
  cc(1,1,1,2) = c(1,6)
  cc(1,1,2,1) =cc(1,1,1,2)
  cc(1,2,1,1) =cc(1,1,1,2)
  cc(2,1,1,1) =cc(1,1,1,2)
  cc(2,2,3,3) = c(2,3)
  cc(3,3,2,2) =cc(2,2,3,3)
  cc(2,2,2,3) = c(2,4)
  cc(2,2,3,2) =cc(2,2,2,3)
  cc(2,3,2,2) =cc(2,2,2,3)
  cc(3,2,2,2) =cc(2,2,2,3)
  cc(2,2,1,3) = c(2,5)
  cc(2,2,3,1) =cc(2,2,1,3)
  cc(1,3,2,2) =cc(2,2,1,3)
  cc(3,1,2,2) =cc(2,2,1,3)
  cc(2,2,1,2) = c(2,6)
  cc(2,2,2,1) =cc(2,2,1,2)
  cc(1,2,2,2) =cc(2,2,1,2)
  cc(2,1,2,2) =cc(2,2,1,2)
  cc(3,3,2,3) = c(3,4)
  cc(3,3,3,2) = cc(3,3,2,3)
  cc(2,3,3,3) = cc(3,3,2,3)
  cc(3,2,3,3) = cc(3,3,2,3)
  cc(3,3,1,3) = c(3,5)
  cc(3,3,3,1) = cc(3,3,1,3)
  cc(1,3,3,3) = cc(3,3,1,3)
  cc(3,1,3,3) = cc(3,3,1,3)
  cc(3,3,1,2) = c(3,6)
  cc(3,3,2,1) = cc(3,3,1,2)
  cc(1,2,3,3) = cc(3,3,1,2)
  cc(2,1,3,3) = cc(3,3,1,2)
  cc(2,3,1,3) = c(4,5)
  cc(3,2,1,3) =cc(2,3,1,3)
  cc(1,3,3,2) =cc(2,3,1,3)
  cc(1,3,2,3) =cc(2,3,1,3)
  cc(2,3,3,1) =cc(2,3,1,3)
  cc(3,2,3,1) =cc(2,3,1,3)
  cc(3,1,2,3) =cc(2,3,1,3)
  cc(3,1,3,2) =cc(2,3,1,3)
  cc(2,3,1,2) = c(4,6)
  cc(3,2,1,2) =cc(2,3,1,2)
  cc(1,2,2,3) =cc(2,3,1,2)
  cc(1,2,3,2) =cc(2,3,1,2)
  cc(2,3,2,1) =cc(2,3,1,2)
  cc(3,2,2,1) =cc(2,3,1,2)
  cc(2,1,2,3) =cc(2,3,1,2)
  cc(2,1,3,2) =cc(2,3,1,2)
  cc(1,3,1,2) = c(5,6)
  cc(3,1,1,2) =cc(1,3,1,2)
  cc(1,2,1,3) =cc(1,3,1,2)
  cc(1,2,3,1) =cc(1,3,1,2)
  cc(1,3,2,1) =cc(1,3,1,2)
  cc(3,1,2,1) =cc(1,3,1,2)
  cc(2,1,1,3) =cc(1,3,1,2)
  cc(2,1,3,1) =cc(1,3,1,2)
  cc(1,2,1,2) = c(6,6)
  cc(2,1,1,2) =cc(1,2,1,2)
  cc(1,2,2,1) =cc(1,2,1,2)
  cc(2,1,2,1) =cc(1,2,1,2)
  return
end subroutine vera_c2c4

!
! print Cijkl tensor to filep
!
subroutine vera_print_cijkl_ftrn(cij)
  implicit none
  double precision :: cij(3,3,3,3)
  integer :: i,j,k,l
  write(*,'(1p,9(e14.6,1x))') ((((cij(i,j,k,l),&
       i=1,3),j=1,3),k=1,3),l=1,3)
  
end subroutine vera_print_cijkl_ftrn
