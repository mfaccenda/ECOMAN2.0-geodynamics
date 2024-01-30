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

   PROGRAM DREX_M

   USE comvar
   USE omp_lib
   !M3E!!!!!!!!!!!!!!
   use class_DistGrid
   !M3E!!!!!!!!!!!!!!

   IMPLICIT NONE

   CHARACTER (500) :: input_dir,output_dir
   INTEGER :: t,ncyc,ncyc1,ncyc2,ndt,nfile !! loop counters
   DOUBLE PRECISION :: wtime,wtime0,wtime1,wtime2,wtimet

   !M3E!!!!!!!!!!!!!!!!!!!!!
   integer :: errMPI,rankMPI
   !M3E!!!!!!!!!!!!!!!!!!!!!

   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Initialize MPI
   call MPI_INIT(errMPI)

   ! Get MPI info
   call MPI_COMM_RANK(MPI_COMM_WORLD,rankMPI,errMPI)
   rankMPI = rankMPI + 1
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime0 = MPI_Wtime() 

   !$omp parallel    
   nt = OMP_GET_NUM_THREADS()
   !$omp end parallel

   if ( rankMPI .eq. 1 ) then
      write(*,*)
      write(*,'(a,i0)') '           NUMBER OF THREADS = ',nt
      write(*,*)
   endif

!!! Read input file and initialization
   CALL init0(input_dir,output_dir)

!!! First advect backward
   CALL backwardadvectionmain(input_dir)

!!! Set fossil fabric after backward advection
   IF(fossilfabric == 1) CALL set_fossil_fabric

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! MAIN LOOP: now move forward and compute FSE and LPO !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if ( rankMPI .eq. 1 ) then
      write(*,"(a)") '********************************************************'
      write(*,"(a)") '********************************************************'
      write(*,*)
      write(*,"(a)") ' STARTING THE MAIN LOOP WITH FORWARD ADVECTION AND LPO CALCULATION'
      write(*,*)
      write(*,'(a,i0)') ' Number of input files to be processed = ',tnum
      write(*,*)
      write(*,"(a)") '--------------------------------------------------------'
      write(*,*)
   endif

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime2 = MPI_Wtime()

   ncyc = 0 ; ncyc1 = 0 ; ncyc2 = 0 ; nfile = 0
   DO t = Tinit , Tend , Tstep 

!!! Load file with velocity components
   CALL load(t,input_dir)

   if ( rankMPI .eq. 1 ) then
      write(*,'(a,i0)') ' CALCULATE LPO FOR FILE NUMBER ',t
      write(*,*)
   endif

   IF(Tinit == Tend) THEN
      timemax = timesum + timemax
   ELSE
      timemax = -1.0
   END IF

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime1 = MPI_Wtime()

   ndt = 0 !Number of sub-timesteps

   nfile = nfile + 1 !Update number of input file since last output file

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
      if ( rankMPI .eq. 1 ) then
         write(*,'(a,i0,a)') ' Start cycle ',ncyc,', steady-state flow'
         write(*,*)
         write(*,'(a,1es13.6,a,1es13.6,a,1es13.6)') '    Current time = ',timesum,           &
         &                                          ' , Final time = ',timemax,              &
         &                                          ' , Time to end = ',timemax-timesum
         write(*,*)
      endif
   ELSE
      IF(ncyc1 == 1) timesum = timesum - dt0
      if ( rankMPI .eq. 1 ) then
         write(*,'(a,i0)') ' Start cycle ',ncyc1
         write(*,*)
      endif
      ndt = ndt + 1
      IF(numdt > 1 .AND. ndt == numdt) dt = dt0 - (numdt-1)*dt
      timesum = timesum + dt
      if ( rankMPI .eq. 1 ) then
         write(*,'(a,1es13.6,a,1es13.6)') '    Timestep = ',dt,' Timesum = ',timesum
         write(*,*)
      endif
   END IF

!  call MPI_Barrier(MPI_COMM_WORLD,errMPI)
!  wtime = MPI_Wtime()

   ! 2D model forward advection + LPO
   IF(dimensions == 2) CALL forwardLPOadvection2D

   ! 3D model forward advection + LPO
   IF(dimensions == 3) CALL forwardLPOadvection3D

   ! Print cycle
   if ( rankMPI .eq. 1 ) then
      OPEN(19,file='cycle.txt',access = 'append')
      write(19,'(i5)') ncyc1
      CLOSE(19)
   endif

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtimet = MPI_Wtime() - wtime1
   if ( rankMPI .eq. 1 ) then
      write(*,'(a,i0,a,i0,a,1f10.2,a)') ' INPUT FILE ',t,', CYCLE ',ncyc1,' OK (',wtimet,' sec)'
      write(*,*)
      write(*,"(a)") '--------------------------------------------------------'
      write(*,*)
   endif
   
!!! Check if advection timestep smaller than timestep input file
   IF(Tinit < Tend .AND. ndt < numdt) THEN
      GOTO 15
   END IF

!!! Output: compute elastic tensor and save infos into hdf5 format
   IF( (Tinit == Tend .AND. ncyc2 == OutputStep) .OR. (Tinit .NE. Tend .AND. (nfile == OutputStep .OR. t == Tend)) .OR. timesum == timemax) THEN

!!! Set fossil fabric after forward advection
      IF(((Tinit .NE. Tend .AND. t == Tend) .OR. timesum == timemax) .AND. fossilfabric == 2) CALL set_fossil_fabric

      CALL stifftenzsave_mpi(t,output_dir,ncyc2)

      if(Tinit .NE. Tend .AND. nfile == OutputStep) nfile = 0
      if(Tinit == Tend .AND. ncyc2 == OutputStep) ncyc2 = 0

   END IF

!!! Steady-state flow (works when Tinit = Tend and timemax > 0 in input.dat) 
   IF(Tinit == Tend .AND. timesum < timemax) THEN
      GOTO 15
   END IF

   ncyc1 = 0

   END DO

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtimet = MPI_Wtime() - wtime2
   if ( rankMPI .eq. 1 ) then
      write(*,*)
      write(*,'(a,1f10.2,a,i0)') ' END OF FORWARD RUN (',wtimet,' sec), TOTAL NUMBER OF CYCLES: ',ncyc
      write(*,*)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! END OF MAIN LOOP for FSE and LPO calculation               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime = MPI_Wtime() - wtime0
   if ( rankMPI .eq. 1 ) then
      write(*,"(a)") '********************************************************'
      write(*,"(a)") '********************************************************'
      write(*,*)
      write(*,'(a,1f10.2,a)') ' ELAPSED TIME FOR THE ENTIRE RUN = ',wtime,' sec'
      write(*,*)
      write(*,"(a)") '********************************************************'
      write(*,"(a)") '********************************************************'
      write(*,*)
   endif

   !M3E!!!!!!!!!!!!!!!!!!!!!
   ! Remove DistGrid
   call RmDistGrid
   ! End MPI
   call MPI_FINALIZE(errMPI)
   !M3E!!!!!!!!!!!!!!!!!!!!!

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
      IF(sbfmod == 0) THEN
         odf(m,:) = 1d0/REAL(size3**3)
         odf_ens(m,:) = odf(m,:)
         acs(:,:,:,m) = acs0
         acs_ens(:,:,:,m) = acs0
         !!! Uncomment next line to reset also FSE
         !Fij(:,:,m) = 0d0 ; Fij(1,1,m) = 1d0 ; Fij(2,2,m) = 1d0 ; Fij(3,3,m) = 1d0
      ELSE
         !!! Uncomment next line to reset also FSE
         Fij(:,:,m) = 0d0 ; Fij(1,1,m) = 1d0 ; Fij(2,2,m) = 1d0 ; Fij(3,3,m) = 1d0
      END IF
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

      !Brd -> pPv (Oganov and Ono, Nature, 2004)
      IF(rocktype(m) == 4 .AND. mpgpa > 98.7+9.56d-3*mtk) rocktype(m) = 5
      !Brd <- pPv
      IF(rocktype(m) == 5 .AND. mpgpa < 98.7+9.56d-3*mtk) rocktype(m) = 4

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

   !M3E!!!!!!!!!!!!
   include 'mpif.h'
   !M3E!!!!!!!!!!!!

   CHARACTER (500) :: input_dir
   INTEGER :: m,ncyc,ncyc1,t,ndt !! loop counters
   INTEGER :: tid,i1,i2,i3,numstrainmax
   DOUBLE PRECISION :: wtime0,wtime1,wtime2,timeback,fractdisl
   !M3E!!!!!!!!!!!!!!!!!!!!!
   integer :: rankMPI,errMPI
   !M3E!!!!!!!!!!!!!!!!!!!!!

   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Get MPI info
   call MPI_COMM_RANK(MPI_COMM_WORLD,rankMPI,errMPI)
   rankMPI = rankMPI + 1
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ncyc = 0 ; ncyc1 = 0
   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime0 = MPI_Wtime()
   timeback = 0

   if ( rankMPI .eq. 1 ) then
      write(*,"(a)") '********************************************************'
      write(*,"(a)") '********************************************************'
      write(*,*)
      write(*,'(a,i0)') ' BACKWARD ADVECTION'
      write(*,*)
      write(*,'(a,i0)') ' Number of input files to be processed = ',tnum
      write(*,*)
   endif

   !Reset cumulated strain and time
   max_strain = 0.0; time_max_strain = 0.0

   !Main loop
   DO t = Tend, Tinit , -Tstep 

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime1 = MPI_Wtime()

!!! Load file with velocity components
   CALL load(t,input_dir)

   IF(Tinit == Tend) THEN
      timeback = timesum + timemax
   ELSE
      !Initialize backward time
      IF(t==Tend) timeback=timesum
      !timeback = timesum
   END IF

   if ( rankMPI .eq. 1 ) then
      write(*,'(a,1es13.6)') ' Timesum = ',timesum
      write(*,*)
   endif

   ndt = 0

15 ncyc  = ncyc  + 1
   ncyc1 = ncyc1 + 1

!!! Steady-state flow (works when Tinit = Tend and timemax > 0 in input.dat)
   IF(Tinit == Tend .AND. timesum < timeback) THEN
      IF(timesum + dt > timeback) THEN
         dt = timeback - timesum
         timeback = timesum
      ELSE
         timeback = timeback - dt
      END IF
      if ( rankMPI .eq. 1 ) then
         write(*,'(a,i0,a)') ' Start cycle ',ncyc,', steady-state flow'
         write(*,*)
         write(*,'(a,1es13.6,a,1es13.6,a,1es13.6)') '    Starting time = ',timesum,          &
         &                                          ' , Current time = ',timeback,        &
         &                                          ' , Time to end = ',timeback-timesum
         write(*,*)
      endif
   ELSE
      if ( rankMPI .eq. 1 ) then
         write(*,'(a,i0)') ' Start cycle ',ncyc1
         write(*,*)
      endif
      ndt = ndt + 1
      IF(numdt > 1 .AND. ndt == numdt) dt = dt0 - (numdt-1)*dt
      timeback = timeback - dt
      if ( rankMPI .eq. 1 ) then
         write(*,'(a,1es13.6,a,1es13.6)') '    Timestep = ',dt,' Timesum = ',timeback
         write(*,*)
      endif
   END IF


   dt = -dt !Reverse time for backward advection

   ! 2D model
   IF(dimensions == 2) THEN

   !$omp parallel do & 
   !$omp schedule(guided,8) &
   !$omp shared(mx1,mx2,X1,X2,Dij,Fd,Fij,e,l,epsnot,odf,odf_ens,acs,acs_ens,rocktype,rho,fractdislrock,max_strain,time_max_strain) &
   !$omp private(tid,m,i1,i2,fractdisl) &    
   !$omp firstprivate(marknum,dt,size,size3,alt,lambda,Xol,Mob,chi,tau,stressexp,acs0,strainmax) &
   !$omp firstprivate(fsemod,fractdislmod,uppermantlemod,x1min,x2min,x1max,x2max,nx1,nx2)
   DO m = 1 , marknum
 
!!! Compute LPO/FSE only for aggregates in the domain .OR. for upper mantle aggregates only

      IF(max_strain(m) < strainmax .AND. ((uppermantlemod == 0 .AND. rocktype(m) < 10) .OR. (uppermantlemod > 0 .AND. rocktype(m) == 1))) THEN

!!! Find nearest upper left node

         CALL upperleft2D(mx1(m),mx2(m),i1,i2)   

!!! Interpolate velocity gradient tensor

         tid =  OMP_GET_THREAD_NUM() + 1

         CALL gradientcalc2D(tid,mx1(m),mx2(m),i1,i2)
         IF(epsnot(tid) == 0) THEN
             write(*,'(a,i10,2f10.3)') ' No strain rate for marker ',m,mx1(m),mx2(m)
             GOTO 50
         END IF

!!! Interpolate fraction of deformation accomodated by dislocation creep/anisotropic phase

         fractdisl = 1d0
         !Scaling factor for deformation accommodated by dislocation creep
         IF(fractdislmod > 0) CALL disldiff2D(m,i1,i2,fractdisl)

!!! Update cumulative disl. creep strain and check if higher than threshold strain max
      
         max_strain(m) = max_strain(m) - dt*epsnot(tid)*fractdisl
         IF(max_strain(m) >= strainmax) time_max_strain(m) = timeback

50    END IF

   END DO
   !$omp end parallel do

!!! Advection of tracers
   !$omp parallel do &
   !$omp schedule(guided,8) &
   !$omp shared(mx1,mx2,mYY,X1,X2,Ui,Tk,Pa,Fij,odf,odf_ens,acs0,acs,acs_ens,rocktype,rho,td_rho,max_strain) &
   !$omp private(m) &
   !$omp firstprivate(marknum,dt,size,size3,fabrictransformmod,ptmod,cartspher,Xol,minx2,maxx2,strainmax) &
   !$omp firstprivate(invdepthaxis,nx1,nx2,x1min,x2min,x1max,x2max,x1periodic,x2periodic)
   DO m = 1 , marknum
 
      IF((time_max_strain(m) == timeback .OR. max_strain(m) < strainmax) .AND. rocktype(m) < 100) THEN

!!! Advection of tracers

         CALL advection2D(m)

!!! Aggregate transformation

         IF(fabrictransformmod > 0 .AND. rocktype(m)<10) CALL rocktypechange(m)
 
      END IF

   END DO
   !$omp end parallel do

   ! 3D model
   ELSE

   !$omp parallel do & 
   !$omp schedule(guided,8) &
   !$omp shared(mx1,mx2,mx3,mYY,X1,X2,X3,Dij,Ui,Fd,epsnot,rocktype,max_strain,time_max_strain) &
   !$omp private(tid,m,i1,i2,i3,fractdisl) &    
   !$omp firstprivate(marknum,dt,strainmax) &
   !$omp firstprivate(fsemod,fractdislmod,uppermantlemod,x1min,x2min,x3min,x1max,x2max,x3max,nx1,nx2,nx3)
   DO m = 1 , marknum

!!! Compute LPO/FSE only for aggregates in the domain .OR. for upper mantle aggregates only

      IF(max_strain(m) < strainmax .AND. ((uppermantlemod == 0 .AND. rocktype(m) < 10) .OR. (uppermantlemod > 0 .AND. rocktype(m) == 1))) THEN

!!! Find nearest upper left node

         CALL upperleft(mx1(m),mx2(m),mx3(m),i1,i2,i3)

!!! Interpolate velocity gradient tensor

         tid =  OMP_GET_THREAD_NUM() + 1

         CALL gradientcalc(tid,mx1(m),mx2(m),mx3(m),i1,i2,i3,mYY(m))
         IF(epsnot(tid) == 0) THEN
             write(*,'(a,i10,3f10.3)') ' No strain rate for marker ',m,mx1(m),mx2(m),mx3(m)
             GOTO 60
         END IF

!!! Interpolate fraction of deformation accomodated by dislocation creep/anisotropic phase

         fractdisl = 1d0
         !Scaling factor for deformation accommodated by dislocation creep
         IF(fractdislmod > 0) CALL disldiff(m,i1,i2,i3,fractdisl)

!!! Update cumulative disl. creep strain and check if higher than threshold strain max
      
         max_strain(m) = max_strain(m) - dt*epsnot(tid)*fractdisl
         IF(max_strain(m) >= strainmax) time_max_strain(m) = timeback

60    END IF

   END DO
   !$omp end parallel do

!!! Advection of tracers
   !$omp parallel do & 
   !$omp schedule(guided,8) &
   !$omp shared(mx1,mx2,mx3,mYY,X1,X2,X3,Ui,Tk,Pa,Fij,odf,odf_ens,acs0,acs,acs_ens,rocktype,rho,td_rho,max_strain,time_max_strain) &
   !$omp private(m) &    
   !$omp firstprivate(marknum,dt,size,size3,fabrictransformmod,ptmod,cartspher,Xol,minx2,maxx2,strainmax) &
   !$omp firstprivate(invdepthaxis,nx1,nx2,nx3,x1min,x2min,x3min,x1max,x2max,x3max,x1periodic,x2periodic,x3periodic)
   DO m = 1 , marknum
 
      IF((time_max_strain(m) == timeback .OR. max_strain(m) < strainmax) .AND. rocktype(m) < 100) THEN

!!! Advection of tracers

         CALL advection(m)

!!! Aggregate transformation
      
         IF(fabrictransformmod > 0 .AND. rocktype(m)<10) CALL rocktypechange(m)
 
      END IF

   END DO
   !$omp end parallel do

   END IF

   dt = -dt

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime2 = MPI_Wtime() - wtime1
   if ( rankMPI .eq. 1 ) then
      write(*,'(a,1f8.2,a)') '    Advection ok (',wtime2,' sec)'
      write(*,*)
   endif

!!! Check number of markers with max_strain >= strainmax 
   numstrainmax = 0
   DO m = 1 , marknum
      IF(max_strain(m) >= strainmax) numstrainmax = numstrainmax + 1
   END DO

   if ( rankMPI .eq. 1 ) then
      write(*,'(a,1i)') '    Number of aggregates that have reached strainmax = ',numstrainmax
      write(*,*)
   endif

   IF(numstrainmax == marknum) THEN
      Tinit = t
      GOTO 90
   END IF

!!! Advection timestep smaller than timestep input file
   IF(Tinit < Tend .AND. ndt < numdt) THEN
      GOTO 15 
   END IF

!!! Steady-state flow (works when Tinit = Tend and timemax > 0 in input.dat)
   IF(Tinit == Tend .AND. timesum < timeback) THEN
      GOTO 15 
   END IF

   ncyc1 = 0

   END DO

90 wtime1 = MPI_Wtime() - wtime0
   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   if (rankMPI .eq. 1 ) then
      write(*,"(a)") '--------------------------------------------------------'
      write(*,*)
      write(*,'(a,1f10.2,a)') ' BACKWARD ADVECTION OK! (',wtime1,' sec)'
      write(*,*)
      write(*,'(a,i0)') ' TOTAL NUMBER OF CYCLES: ',ncyc
      write(*,*)
      write(*,"(a)") '********************************************************'
      write(*,"(a)") '********************************************************'
      write(*,*)
   endif

   END SUBROUTINE backwardadvectionmain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Advect forward in time and compute LPO of aggregates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE forwardLPOadvection2D

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   !M3E!!!!!!!!!!!!
   include 'mpif.h'
   !M3E!!!!!!!!!!!!

   INTEGER :: i1,i2,m,tid
   DOUBLE PRECISION :: wtime,fractdisl
   !M3E!!!!!!!!!!!!!!!!!!!!!
   integer :: errMPI,rankMPI
   !M3E!!!!!!!!!!!!!!!!!!!!!

   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Get MPI info
   call MPI_COMM_RANK(MPI_COMM_WORLD,rankMPI,errMPI)
   rankMPI = rankMPI + 1
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime = MPI_Wtime()

   !$omp parallel do & 
   !$omp schedule(guided,8) &
   !$omp shared(mx1,mx2,X1,X2,Dij,Ui,Fd,Fij,e,l,epsnot,time_max_strain) &
   !$omp shared(odf,odf_ens,acs,acs_ens,rocktype,rho,fractdislrock) &
   !$omp private(tid,m,i1,i2,fractdisl) &    
   !$omp firstprivate(marknum,dt,size,size3,alt,lambda,Xol,Mob,chi,tau,stressexp,acs0) &
   !$omp firstprivate(fsemod,fractdislmod,uppermantlemod,x1min,x2min,x1max,x2max,nx1,nx2)
   DO m = 1 , marknum
 
!!! Compute LPO/FSE only for aggregates in the domain .OR. for upper mantle aggregates only

      IF(timesum >= time_max_strain(m) .AND. ((uppermantlemod == 0 .AND. rocktype(m) < 10) .OR. (uppermantlemod > 0 .AND. rocktype(m) == 1))) THEN

!!! Find nearest upper left node

         CALL upperleft2D(mx1(m),mx2(m),i1,i2)   

!!! Interpolate velocity gradient tensor

         tid =  OMP_GET_THREAD_NUM() + 1

         CALL gradientcalc2D(tid,mx1(m),mx2(m),i1,i2)
         IF(epsnot(tid) == 0) THEN
             write(*,'(a,i10,2f10.3)') ' No strain rate for marker ',m,mx1(m),mx2(m)
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
   !$omp end parallel do

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime = MPI_Wtime() - wtime

   if ( rankMPI .eq. 1 ) then
      write(*,'(a,1f8.2,a)') '    LPO calculation ok (',wtime,' sec)'
      write(*,*)
   endif

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime = MPI_Wtime()

!!! Advection of aggregates
   !$omp parallel do & 
   !$omp schedule(guided,8) &
   !$omp shared(mx1,mx2,X1,X2,Ui,Fij,odf,odf_ens,acs,acs_ens,rocktype,rho,td_rho,Tk,Pa,time_max_strain) &
   !$omp private(m) &    
   !$omp firstprivate(marknum,dt,size,size3,fabrictransformmod,ptmod,cartspher,Xol,minx2,maxx2) &
   !$omp firstprivate(invdepthaxis,nx1,nx2,x1min,x2min,x1max,x2max,x1periodic,x2periodic,acs0)
   DO m = 1 , marknum
 
      IF(timesum >= time_max_strain(m) .AND. rocktype(m) < 100) THEN

         CALL advection2D(m)

!!! Aggregate transformation
      
         IF(fabrictransformmod > 0 .AND. rocktype(m)<10) CALL rocktypechange(m)
 
      END IF

   END DO
   !$omp end parallel do

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime = MPI_Wtime() - wtime

   if ( rankMPI .eq. 1 ) then
      write(*,'(a,1f8.2,a)') '    2D Advection ok (',wtime,' sec)'
      write(*,*)
   endif

   RETURN

   END SUBROUTINE forwardLPOadvection2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Advect forward in time and compute LPO of aggregates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE forwardLPOadvection3D

   USE comvar
   USE omp_lib

   IMPLICIT NONE

   !M3E!!!!!!!!!!!!
   include 'mpif.h'
   !M3E!!!!!!!!!!!!

   INTEGER :: i1,i2,i3,m,tid
   DOUBLE PRECISION :: wtime,fractdisl
   !M3E!!!!!!!!!!!!!!!!!!!!!
   integer :: errMPI,rankMPI
   !M3E!!!!!!!!!!!!!!!!!!!!!

   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Get MPI info
   call MPI_COMM_RANK(MPI_COMM_WORLD,rankMPI,errMPI)
   rankMPI = rankMPI + 1
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime = MPI_Wtime()

   !$omp parallel do & 
   !$omp schedule(guided,8) &
   !$omp shared(mx1,mx2,mx3,mYY,X1,X2,X3,Dij,Ui,Fd,Fij,e,l,epsnot,time_max_strain) &
   !$omp shared(odf,odf_ens,acs,acs_ens,rocktype,rho,fractdislrock) &
   !$omp private(tid,m,i1,i2,i3,fractdisl) &    
   !$omp firstprivate(marknum,timesum,dt,size,size3,alt,lambda,Xol,Mob,chi,tau,stressexp,acs0) &
   !$omp firstprivate(fsemod,fractdislmod,uppermantlemod,x1min,x2min,x3min,x1max,x2max,x3max,nx1,nx2,nx3)
   DO m = 1 , marknum

!!! Compute LPO/FSE only for aggregates in the domain .OR. for upper mantle aggregates only

      IF(timesum >= time_max_strain(m) .AND. ((uppermantlemod == 0 .AND. rocktype(m) < 10) .OR. (uppermantlemod > 0 .AND. rocktype(m) == 1))) THEN

!!! Find nearest upper left node

         CALL upperleft(mx1(m),mx2(m),mx3(m),i1,i2,i3)

!!! Interpolate velocity gradient tensor

         tid =  OMP_GET_THREAD_NUM() + 1

         CALL gradientcalc(tid,mx1(m),mx2(m),mx3(m),i1,i2,i3,mYY(m))
         IF(epsnot(tid) == 0) THEN
             write(*,'(a,i10,3f10.3)') ' No strain rate for marker ',m,mx1(m),mx2(m),mx3(m)
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
   !$omp end parallel do

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime = MPI_Wtime() - wtime
   if ( rankMPI .eq. 1 ) then
      write(*,'(a,1f8.2,a)') '    LPO calculation ok (',wtime,' sec)'
      write(*,*)
   endif

   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime = MPI_Wtime()

!!! Advection of aggregates
   !$omp parallel do & 
   !$omp schedule(guided,8) &
   !$omp shared(mx1,mx2,mx3,mYY,X1,X2,X3,Ui,Tk,Pa,Fij,odf,odf_ens,acs0,acs,acs_ens,rocktype,rho,td_rho,time_max_strain) &
   !$omp private(m) &    
   !$omp firstprivate(marknum,timesum,dt,size,size3,fabrictransformmod,ptmod,cartspher,Xol,minx2,maxx2) &
   !$omp firstprivate(invdepthaxis,nx1,nx2,nx3,x1min,x2min,x3min,x1max,x2max,x3max,x1periodic,x2periodic,x3periodic)
   DO m = 1 , marknum

      IF(timesum >= time_max_strain(m) .AND. rocktype(m) < 100) THEN

         CALL advection(m)

!!! Aggregate transformation

         IF(fabrictransformmod > 0 .AND. rocktype(m)<10) CALL rocktypechange(m)

      END IF

   END DO
   !$omp end parallel do
 
   call MPI_Barrier(MPI_COMM_WORLD,errMPI)
   wtime = MPI_Wtime() - wtime
   if ( rankMPI .eq. 1 ) then
      write(*,'(a,1f8.2,a)') '    3D Advection ok (',wtime,' sec)'
      write(*,*)
   endif

   RETURN

   END SUBROUTINE forwardLPOadvection3D

