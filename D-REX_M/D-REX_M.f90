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

   PROGRAM DREX_M

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   CHARACTER (500) :: input_dir,output_dir
   INTEGER :: t,i,m,ncyc,ncyc1,ncyc2,tid,ndt !! loop counters
   DOUBLE PRECISION :: wtime,wtime0,wtime1,wtime2,fractdisl

   wtime0 = OMP_GET_WTIME()

   !$omp parallel    
   nt = OMP_GET_NUM_THREADS()
   !$omp end parallel

   write(*,*)
   write(*,'(a,i0)'),'           NUMBER OF THREADS = ',nt
   write(*,*)

!!! Read input file and initialization
   CALL init0(input_dir,output_dir)

!!! First advect backward
   CALL backwardadvectionmain(input_dir)

!!! Set fossil fabric after backward advection
   IF(fossilfabric == 1) CALL set_fossil_fabric

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! MAIN LOOP: now move forward and compute FSE and LPO !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)
   write(*,"(a)"),' STARTING THE MAIN LOOP WITH FORWARD ADVECTION AND LPO CALCULATION'
   write(*,*)
   write(*,'(a,i0)'),' Number of input files to be processed = ',tnum
   write(*,*)
   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)

   wtime2 = OMP_GET_WTIME()

   ncyc = 0 ; ncyc1 = 0 ; ncyc2 = 0
   DO t = Tinit , Tend , Tstep 

!!! Load file with velocity components
   CALL load(t,input_dir)

   write(*,'(a,i0)'),' CALCULATE LPO FOR FILE NUMBER ',t
   write(*,*)

   IF(Tinit == Tend) THEN
      timemax = timesum + timemax
   ELSE
      timemax = -1.0
   END IF

   wtime1 = OMP_GET_WTIME()

   ndt = 0 !Number of sub-timesteps

15 ncyc  = ncyc  + 1 !Update number of cycles since the start of the forward run                       
   ncyc1 = ncyc1 + 1 !Update number of cycles since the start of this input file                        
   ncyc2 = ncyc2 + 1 !Update number of cycles since last output file  

!!! Steady-state flow (works when Tinit = Tend and timemax > 0 in input.dat)
   IF(Tinit == Tend .AND. timesum < timemax) THEN
      IF(timesum + dt > timemax) THEN
         dt = timemax - timesum
         timesum = timemax
      ELSE
         timesum = timesum + dt
      END IF
      write(*,'(a,i0,a)'),' Start cycle ',ncyc,', steady-state flow'
      write(*,*)
      write(*,'(a,1es13.6,a,1es13.6,a,1es13.6)'),'    Current time = ',timesum,' , Final time = ',timemax,' , Time to end = ',timemax-timesum
      !write(*,'(a,1f10.8,a,1f10.8,a,1f10.8)'),'    Current time = ',timesum,' , Final time = ',timemax,' , Time to end = ',timemax-timesum
      write(*,*)
   ELSE
      IF(ncyc1 == 1) timesum = timesum - dt0
      write(*,'(a,i0)'),' Start cycle ',ncyc1
      write(*,*)
      ndt = ndt + 1
      IF(numdt > 1 .AND. ndt == numdt) dt = dt0 - (numdt-1)*dt
      timesum = timesum + dt
      write(*,'(a,1es13.6,a,1es13.6)'),'    Timestep = ',dt,' Timesum = ',timesum
      !write(*,'(a,1f10.8,a,1f10.8)'),'    Timestep = ',dt,' Timesum = ',timesum
      write(*,*)
   END IF

   wtime = OMP_GET_WTIME()

   ! 2D model forward advection + LPO
   IF(dimensions == 2) CALL forwardLPOadvection2D

   ! 3D model forward advection + LPO
   IF(dimensions == 3) CALL forwardLPOadvection3D

   ! Print cycle
   OPEN(19,file='cycle.txt',access = 'append')
   write(19,'(i5)') ncyc1
   CLOSE(19)

   write(*,'(a,i0,a,i0,a,1f10.2,a)'),' INPUT FILE ',t,', CYCLE ',ncyc1,' OK (',OMP_GET_WTIME() - wtime1,' sec)'
   write(*,*)
   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)
   
!!! Check if advection timestep smaller than timestep input file
   IF(Tinit < Tend .AND. ndt < numdt) THEN
      GOTO 15
   END IF

!!! Output: compute elastic tensor and save infos into hdf5 format
   IF(ncyc2 == OutputStep .OR. (Tinit .NE. Tend .AND. t == Tend) .OR. timesum == timemax) THEN

!!! Set fossil fabric after forward advection
      IF(((Tinit .NE. Tend .AND. t == Tend) .OR. timesum == timemax) .AND. fossilfabric == 2) CALL set_fossil_fabric

      CALL stifftenzsave(t,output_dir,ncyc2)
   
      ncyc2 = 0

   END IF

!!! Steady-state flow (works when Tinit = Tend and timemax > 0 in input.dat) 
   IF(Tinit == Tend .AND. timesum < timemax) THEN
      GOTO 15
   END IF

   ncyc1 = 0

   END DO

   write(*,*)
   write(*,'(a,1f10.2,a,i0)'),' END OF FORWARD RUN (',OMP_GET_WTIME() - wtime2,' sec), TOTAL NUMBER OF CYCLES: ',ncyc
   write(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! END OF MAIN LOOP for FSE and LPO calculation               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)
   write(*,'(a,1f10.2,a)'),' ELAPSED TIME FOR THE ENTIRE RUN = ',OMP_GET_WTIME() - wtime0,' sec'
   write(*,*)
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)

   STOP

   END PROGRAM DREX_M  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Check if rocktype of aggregates is changing and reset it   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rocktypechange(m)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: m,rocktype_old

   rocktype_old = rocktype(m)

   CALL rocktypecheck(m)

   ! Reset aggregate's LPO that experience phase transformation
   IF(fabrictransformmod == 2 .AND. rocktype(m) .NE. rocktype_old) THEN
      odf(m,:) = 1d0/REAL(size3**3)
      odf_ens(m,:) = odf(m,:)
      acs(:,:,:,m) = acs0
      acs_ens(:,:,:,m) = acs0
      !!! Uncomment next line to reset also FSE
      !Fij(:,:,m) = 0d0 ; Fij(1,1,m) = 1d0 ; Fij(2,2,m) = 1d0 ; Fij(3,3,m) = 1d0
   END IF

   RETURN
 
   END SUBROUTINE rocktypechange
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Check if rocktype of aggregates is changing and reset it   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE rocktypecheck(m)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   INTEGER :: m
   DOUBLE PRECISION :: mtk,mpgpa

   IF(ptmod == 2) THEN

      CALL rhopt(m,mtk,mpgpa)

      IF(rho(m) > 3000) rocktype(m) = 1
      IF(rho(m) > 3650) rocktype(m) = 2
      IF(rho(m) > 3870) rocktype(m) = 3
      IF(rho(m) > 4150) rocktype(m) = 4
         
      IF(mtk < 300) rocktype(m) = rocktype(m) + 100

   ELSE
         
      ! Assign rocktype
      IF(mx2(m) > minx2(1) .AND. mx2(m) <= maxx2(1)) rocktype(m) = 1 ! Upper mantle (Ol + Ens)
      IF(mx2(m) > minx2(2) .AND. mx2(m) <= maxx2(2)) rocktype(m) = 2 ! Transition Zone (Wd + Grt)
      IF(mx2(m) > minx2(3) .AND. mx2(m) <= maxx2(3)) rocktype(m) = 3 ! Transition Zone (Rw + Grt)
      IF(mx2(m) > minx2(4) .AND. mx2(m) <= maxx2(4)) rocktype(m) = 4 ! Lower mantle (Pv + MgO)

   END IF
   
   RETURN
 
   END SUBROUTINE rocktypecheck
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Advect backward in time aggregates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE backwardadvectionmain(input_dir)

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   CHARACTER (500) :: input_dir
   INTEGER :: m,ncyc,ncyc1,t,ndt !! loop counters
   DOUBLE PRECISION :: wtime0,wtime1,timesumback,timemaxback,timeback

   ncyc = 0 ; ncyc1 = 0
   wtime0 = OMP_GET_WTIME()
   timeback = 0

   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)
   write(*,'(a,i0)'),' BACKWARD ADVECTION'
   write(*,*)
   write(*,'(a,i0)'),' Number of input files to be processed = ',tnum
   write(*,*)
        
   DO t = Tend, Tinit , -Tstep 

   wtime1 = OMP_GET_WTIME()

!!! Load file with velocity components
   CALL load(t,input_dir)

   IF(Tinit == Tend) THEN
      timemaxback = timesum + timemax
   ELSE
      timemaxback = -1.0
   END IF

   IF(t==Tend) timeback=timesum
   write(*,'(a,1es13.6)'),' Timesum = ',timesum
   !write(*,'(a,1f8.3)'),' Timesum = ',timesum
   write(*,*)
   
   ndt = 0

15 ncyc  = ncyc  + 1
   ncyc1 = ncyc1 + 1

!!! Steady-state flow (works when Tinit = Tend and timemax > 0 in input.dat)
   IF(Tinit == Tend .AND. timesum < timemaxback) THEN
      IF(timesum + dt > timemaxback) THEN
         dt = timemaxback - timesum
         timemaxback = timesum
      ELSE
         timemaxback = timemaxback - dt
      END IF
      write(*,'(a,i0,a)'),' Start cycle ',ncyc,', steady-state flow'
      write(*,*)
      write(*,'(a,1es13.6,a,1es13.6,a,1es13.6)'),'    Starting time = ',timesum,' , Current time = ',timemaxback,' , Time to end = ',timemaxback-timesum
      !write(*,'(a,1f10.8,a,1f10.8,a,1f10.8)'),'    Starting time = ',timesum,' , Current time = ',timemaxback,' , Time to end = ',timemaxback-timesum
      write(*,*)
   ELSE
      write(*,'(a,i0)'),' Start cycle ',ncyc1
      write(*,*)
      ndt = ndt + 1
      IF(numdt > 1 .AND. ndt == numdt) dt = dt0 - (numdt-1)*dt
      timeback = timeback - dt
      write(*,'(a,1es13.6,a,1es13.6)'),'    Timestep = ',dt,' Timesum = ',timeback
      !write(*,'(a,1f10.8,a,1f10.8)'),'    Timestep = ',dt,' Timesum = ',timeback
      write(*,*)
   END IF


   dt = -dt !Reverse time for backward advection

   ! 2D model
   IF(dimensions == 2) THEN

!!! Advection of tracers
   !$omp parallel & 
   !$omp shared(mx1,mx2,mYY,X1,X2,Ui,Tk,Pa,Fij,odf,odf_ens,acs0,acs,acs_ens,rocktype,rho,td_rho) &
   !$omp private(m) &    
   !$omp firstprivate(marknum,dt,size,size3,fabrictransformmod,ptmod,cartspher,Xol,minx2,maxx2,invdepthaxis,nx1,nx2,x1min,x2min,x1max,x2max,x1periodic,x2periodic)
   !$omp do schedule(guided,8)
   DO m = 1 , marknum
 
      IF(rocktype(m) < 100) THEN

!!! Advection of tracers

         CALL advection2D(m)

!!! Aggregate transformation
      
         IF(fabrictransformmod > 0 .AND. rocktype(m)<10) CALL rocktypechange(m)
 
      END IF

   END DO
   !$omp end do
   !$omp end parallel

   ! 3D model
   ELSE

!!! Advection of tracers
   !$omp parallel & 
   !$omp shared(mx1,mx2,mx3,mYY,X1,X2,X3,Ui,Tk,Pa,Fij,odf,odf_ens,acs0,acs,acs_ens,rocktype,rho,td_rho) &
   !$omp private(m) &    
   !$omp firstprivate(marknum,dt,size,size3,fabrictransformmod,ptmod,cartspher,Xol,minx2,maxx2,invdepthaxis,nx1,nx2,nx3,x1min,x2min,x3min,x1max,x2max,x3max,x1periodic,x2periodic,x3periodic)
   !$omp do schedule(guided,8)
   DO m = 1 , marknum
 
      IF(rocktype(m) < 100) THEN

!!! Advection of tracers

         CALL advection(m)

!!! Aggregate transformation
      
         IF(fabrictransformmod > 0 .AND. rocktype(m)<10) CALL rocktypechange(m)
 
      END IF

   END DO
   !$omp end do
   !$omp end parallel

   END IF

   dt = -dt

   write(*,'(a,1f8.2,a)'),'    Advection ok (',OMP_GET_WTIME() - wtime1,' sec)'
   write(*,*)

!!! Advection timestep smaller than timestep input file
   IF(Tinit < Tend .AND. ndt < numdt) THEN
      GOTO 15 
   END IF

!!! Steady-state flow (works when Tinit = Tend and timemax > 0 in input.dat)
   IF(Tinit == Tend .AND. timesum < timemaxback) THEN
      GOTO 15 
   END IF

   ncyc1 = 0

   END DO

   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)
   write(*,'(a,1f10.2,a)'),' BACKWARD ADVECTION OK! (',OMP_GET_WTIME() - wtime0,' sec)'
   write(*,*)
   write(*,'(a,i0)'),' TOTAL NUMBER OF CYCLES: ',ncyc
   write(*,*)
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)

   END SUBROUTINE backwardadvectionmain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Advect forward in time and compute LPO of aggregates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE forwardLPOadvection2D

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   CHARACTER (500) :: input_dir
   INTEGER :: i1,i2,m,tid
   DOUBLE PRECISION :: wtime,fractdisl

   wtime = OMP_GET_WTIME()

   !$omp parallel & 
   !$omp shared(mx1,mx2,X1,X2,Dij,Fd,Fij,e,l,epsnot,odf,odf_ens,acs,acs_ens,rocktype,rho,fractdislrock) &
   !$omp private(tid,m,i1,i2,fractdisl) &    
   !$omp firstprivate(marknum,dt,size,size3,alt,lambda,Xol,Mob,chi,tau,stressexp,acs0,fsemod,fractdislmod,uppermantlemod,x1min,x2min,x1max,x2max,nx1,nx2)
   !$omp do schedule(guided,8)
   DO m = 1 , marknum
 
!!! Compute LPO/FSE only for aggregates in the domain .OR. for upper mantle aggregates only

      IF((uppermantlemod == 0 .AND. rocktype(m) < 10) .OR. (uppermantlemod > 0 .AND. rocktype(m) == 1)) THEN

!!! Find nearest upper left node

         CALL upperleft2D(mx1(m),mx2(m),i1,i2)   

!!! Interpolate velocity gradient tensor

         tid =  OMP_GET_THREAD_NUM() + 1

         CALL gradientcalc2D(tid,mx1(m),mx2(m),i1,i2)
         IF(epsnot(tid) == 0) THEN
             write(*,'(a,i,2f10.3)'),' No strain rate for marker ',m,mx1(m),mx2(m)
             GOTO 50
         END IF

!!! Interpolate fraction of deformation accomodated by dislocation creep/anisotropic phase

         fractdisl = 1d0
         !Scaling factor for deformation accommodated by dislocation creep
         IF(fractdislmod > 0) CALL disldiff2D(m,i1,i2,fractdisl)

!!! Calculation of the LPO/FSE
      
         CALL strain(tid,m,fractdisl)

50    END IF

   END DO
   !$omp end do
   !$omp end parallel

   write(*,'(a,1f8.2,a)'),'    LPO calculation ok (',OMP_GET_WTIME() - wtime,' sec)'
   write(*,*)

   wtime = OMP_GET_WTIME()

!!! Advection of aggregates
   !$omp parallel & 
   !$omp shared(mx1,mx2,X1,X2,Ui,Fij,odf,odf_ens,acs,acs_ens,rocktype,rho,td_rho,Tk,Pa) &
   !$omp private(m) &    
   !$omp firstprivate(marknum,dt,size,size3,fabrictransformmod,ptmod,cartspher,Xol,minx2,maxx2,invdepthaxis,nx1,nx2,x1min,x2min,x1max,x2max,x1periodic,x2periodic,acs0)
   !$omp do schedule(guided,8)
   DO m = 1 , marknum
 
      IF(rocktype(m) < 100) THEN

         CALL advection2D(m)

!!! Aggregate transformation
      
         IF(fabrictransformmod > 0 .AND. rocktype(m)<10) CALL rocktypechange(m)
 
      END IF

   END DO
   !$omp end do
   !$omp end parallel

   write(*,'(a,1f8.2,a)'),'    2D Advection ok (',OMP_GET_WTIME() - wtime,' sec)'
   write(*,*)

   RETURN

   END SUBROUTINE forwardLPOadvection2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Advect forward in time and compute LPO of aggregates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE forwardLPOadvection3D

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   CHARACTER (500) :: input_dir
   INTEGER :: i1,i2,i3,m,tid
   DOUBLE PRECISION :: wtime,fractdisl

   wtime = OMP_GET_WTIME()

   !$omp parallel & 
   !$omp shared(mx1,mx2,mx3,mYY,X1,X2,X3,Dij,Ui,Fd,Fij,e,l,epsnot,max_strain,odf,odf_ens,acs,acs_ens,rocktype,rho,fractdislrock) &
   !$omp private(tid,m,i1,i2,i3,fractdisl) &    
   !$omp firstprivate(marknum,dt,size,size3,alt,lambda,Xol,Mob,chi,tau,stressexp,acs0,fsemod,fractdislmod,uppermantlemod,x1min,x2min,x3min,x1max,x2max,x3max,nx1,nx2,nx3)
   !$omp do schedule(guided,8)
   DO m = 1 , marknum

!!! Compute LPO/FSE only for aggregates in the domain .OR. for upper mantle aggregates only

      IF((uppermantlemod == 0 .AND. rocktype(m) < 10) .OR. (uppermantlemod > 0 .AND. rocktype(m) == 1)) THEN

!!! Find nearest upper left node

         CALL upperleft(mx1(m),mx2(m),mx3(m),i1,i2,i3)

!!! Interpolate velocity gradient tensor

         tid =  OMP_GET_THREAD_NUM() + 1

         CALL gradientcalc(tid,mx1(m),mx2(m),mx3(m),i1,i2,i3,mYY(m))
         IF(epsnot(tid) == 0) THEN
             write(*,'(a,i,2f10.3)'),' No strain rate for marker ',m,mx1(m),mx2(m)
             GOTO 50
         END IF

!!! Interpolate fraction of deformation accomodated by dislocation creep/anisotropic phase

         fractdisl = 1d0
         !Scaling factor for deformation accommodated by dislocation creep
         IF(fractdislmod > 0) CALL disldiff(m,i1,i2,i3,fractdisl)

!!! Calculation of the LPO/FSE

         CALL strain(tid,m,fractdisl)

50    END IF

   END DO
   !$omp end do
   !$omp end parallel

   write(*,'(a,1f8.2,a)'),'    LPO calculation ok (',OMP_GET_WTIME() - wtime,' sec)'
   write(*,*)

   wtime = OMP_GET_WTIME()

!!! Advection of aggregates
   !$omp parallel & 
   !$omp shared(mx1,mx2,mx3,mYY,X1,X2,X3,Ui,Tk,Pa,Fij,odf,odf_ens,acs0,acs,acs_ens,rocktype,rho,td_rho) &
   !$omp private(m) &    
   !$omp firstprivate(marknum,dt,size,size3,fabrictransformmod,ptmod,cartspher,Xol,minx2,maxx2,invdepthaxis,nx1,nx2,nx3,x1min,x2min,x3min,x1max,x2max,x3max,x1periodic,x2periodic,x3periodic)
   !$omp do schedule(guided,8)
   DO m = 1 , marknum

      IF(rocktype(m) < 100) THEN

         CALL advection(m)

!!! Aggregate transformation

         IF(fabrictransformmod > 0 .AND. rocktype(m)<10) CALL rocktypechange(m)

      END IF

   END DO
   !$omp end do
   !$omp end parallel

   write(*,'(a,1f8.2,a)'),'    3D Advection ok (',OMP_GET_WTIME() - wtime,' sec)'
   write(*,*)

   RETURN

   END SUBROUTINE forwardLPOadvection3D

