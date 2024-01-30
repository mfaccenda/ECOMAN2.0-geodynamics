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

   SUBROUTINE read_input_file(input_dir,output_dir)

   USE comvar

   IMPLICIT NONE

   !M3E!!!!!!!!!!!!
   include 'mpif.h'
   !M3E!!!!!!!!!!!!

   INTEGER :: i
   CHARACTER (len=500) :: input_dir,output_dir
   CHARACTER(200) :: arg,inputfilename
!  CHARACTER(200) :: str ! unused

   !M3E!!!!!!!!!!!!!!!!!!!!!
   integer :: errMPI,rankMPI
   !M3E!!!!!!!!!!!!!!!!!!!!!

   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Get MPI info
   call MPI_COMM_RANK(MPI_COMM_WORLD,rankMPI,errMPI)
   rankMPI = rankMPI + 1
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   DO i = 1, iargc()
      CALL getarg(i, arg)
      if(i==1) inputfilename = arg
   END DO
 
   OPEN(unit=15,file=inputfilename)

   !Input files directory
   if (rankMPI .eq. 1 ) then
      write(*,"(a)") '********************************************************'
      write(*,"(a)") '********************************************************'
      write(*,*)
      write(*,"(a,a)") 'READING INPUT FILE: ',trim(inputfilename)
      write(*,*)
   endif
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call read_par_int(15,nproc1)
   call read_par_int(15,nproc2)
   call read_par_int(15,nproc3)
   !M3E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call read_par_char(15,input_dir)
   call read_par_char(15,output_dir)
   if (rankMPI .eq. 1 ) then
      write(*,"(a,i10)") ' nproc1 :          ',nproc1
      write(*,"(a,i10)") ' nproc2 :          ',nproc2
      write(*,"(a,i10)") ' nproc3 :          ',nproc3
      write(*,*)
      write(*,"(a,a)") ' Input  directory: ',trim(input_dir)
      write(*,*)
      write(*,"(a,a)") ' Output directory: ',trim(output_dir)
      write(*,*)
   endif

   !Input files min,step,max, output frequency, max. time (in Myr) in
   !steady-state conditions
   if (rankMPI .eq. 1 ) then
      write(*,"(a)") 'INPUT - OUTPUT FILES'
      write(*,*)
   endif
   call read_par_int(15,Tinit)
   call read_par_int(15,Tstep)
   call read_par_int(15,Tend)
   call read_par_int(15,OutputStep)
   call read_par_double(15,timemax)
   !call read_par_double(15,strainmax)
   strainmax = 5.0d1
   if ( rankMPI .eq. 1 ) then
      write(*,"(a,i5)") ' Tinit :      ',Tinit
      write(*,"(a,i5)") ' Tstep :      ',Tstep
      write(*,"(a,i5)") ' Tend  :      ',Tend
      write(*,"(a,i5)") ' OutputStep : ',OutputStep
      write(*,*)
      if(Tinit == Tend) write(*,"(a,f6.2,a)")  ' Tinit = Tend --> steady-state flow for ',Timemax,' Myr'
      write(*,*)
      write(*,"(a,i5)") ' Max. strain : ',strainmax
      write(*,*)
      write(*,"(a)") '--------------------------------------------------------'
      write(*,*)
   endif

   !Define model grid
   if (rankMPI .eq. 1 ) then
      write(*,"(a)") 'DEFINE MODEL GRIDS'
      write(*,*)
   endif
   call read_par_int(15,dimensions)
   call read_par_int(15,cartspher)
   call read_par_int(15,basicstag) ;
   if (rankMPI .eq. 1 ) then
      if(dimensions==2) write(*,"(a,i5,a)") ' dimensions: ',dimensions,' --> 2D model'
      if(dimensions==3) write(*,"(a,i5,a)") ' dimensions: ',dimensions,' --> 3D model'
      if(cartspher==1) write(*,"(a,i5,a)") ' cartspher : ',cartspher,' --> cartesian coordinates system'
      if(cartspher==2 .AND. dimensions == 2) write(*,"(a,i5,a)") ' cartspher : ',cartspher,' --> polar coordinates system'
      if(cartspher==2 .AND. dimensions == 3) write(*,"(a,i5,a)") ' cartspher : ',cartspher,' --> spherical coordinates system'
   endif

   if (rankMPI .eq. 1 ) then
      write(*,*)
      write(*,"(a)") 'EULERIAN GRID'
   endif
   call read_par_double(15,x1min)   
   call read_par_double(15,x1max)   
   call read_par_int(15,nx1)
   call read_par_int(15,x1periodic)
   call read_par_double(15,x2min)
   call read_par_double(15,x2max)
   call read_par_int(15,nx2)
   call read_par_int(15,x2periodic)
   call read_par_double(15,x3min)
   call read_par_double(15,x3max)
   call read_par_int(15,nx3)
   call read_par_int(15,x3periodic)
   if (rankMPI .eq. 1) then
      write(*,*)
      write(*,"(a)") ' AXIS 1'
      if(cartspher==1) then
         write(*,"(a,f14.6,a)") ' x1min: ',x1min,' (m)'
         write(*,"(a,f14.6,a)") ' x1max: ',x1max,' (m)'
      end if
      if(cartspher==2) then
         write(*,"(a,f14.6,a)") ' x1min: ',x1min,' (degrees)'
         write(*,"(a,f14.6,a)") ' x1max: ',x1max,' (degrees)'
      end if
      write(*,"(a,i5)")    ' nx1  :          ',nx1
      write(*,"(a,i5)")    ' x1periodic:     ',x1periodic
      write(*,*)
      write(*,"(a)") ' AXIS 2'
      write(*,"(a,f14.6,a)") ' x2min: ',x2min,' (m)'
      write(*,"(a,f14.6,a)") ' x2max: ',x2max,' (m)'
      write(*,"(a,i5)")    ' nx2  :          ',nx2
      write(*,"(a,i5)")    ' x2periodic:     ',x2periodic
      write(*,*)
      write(*,"(a)") ' AXIS 3'
      if(cartspher==1) then
         write(*,"(a,f14.6,a)") ' x3min: ',x3min,' (m)'
         write(*,"(a,f14.6,a)") ' x3max: ',x3max,' (m)'
      end if
      if(cartspher==2) then
         write(*,"(a,f14.6,a)") ' x3min: ',x3min,' (degrees)'
         write(*,"(a,f14.6,a)") ' x3max: ',x3max,' (degrees)'
      end if
      write(*,"(a,i5)")    ' nx3  :          ',nx3
      write(*,"(a,i5)")    ' x3periodic:     ',x3periodic
      write(*,*)
   endif

   yinyang = 1
   if(cartspher == 2 .AND. dimensions == 3 .AND. x1periodic > 0 .AND. x3periodic > 0 ) THEN
      yinyang = 2
      if (rankMPI .eq. 1 ) then
         write(*,"(a)") ' !!!!!!! YIN-YANG GRID ACTIVE !!!!!!!!!'
         write(*,*)
      end if
   end if

   !!! Set number of nodes along 3rd axis equal to 1 if model is 2D
   if(dimensions == 2) nx3 = 1

   !Aggregate distribution
   if (rankMPI .eq. 1 ) then
      write(*,"(a)") 'LAGRANGIAN GRID'
      write(*,*)
   endif
   call read_par_double(15,mx1min)
   call read_par_double(15,mx1max)
   call read_par_double(15,mx1stp)
   call read_par_double(15,mx2min)
   call read_par_double(15,mx2max)
   call read_par_double(15,mx2stp)
   call read_par_double(15,mx3min)
   call read_par_double(15,mx3max)
   call read_par_double(15,mx3stp)
   if (rankMPI .eq. 1) then
      write(*,"(a)") ' AXIS 1'
      if(cartspher==1) then
         write(*,"(a,f14.6,a)") ' mx1min: ',mx1min,' (m)'
         write(*,"(a,f14.6,a)") ' mx1max: ',mx1max,' (m)'
      end if
      if(cartspher==2) then
         write(*,"(a,f14.6,a)") ' mx1min: ',mx1min,' (degrees)'
         write(*,"(a,f14.6,a)") ' mx1max: ',mx1max,' (degrees)'
      end if
      write(*,"(a,f14.6,a)") ' mx1stp: ',mx1stp,' (m)'
      write(*,*)
      write(*,"(a)") ' AXIS 2'
      write(*,"(a,f14.6,a)") ' mx2min: ',mx2min,' (m)'
      write(*,"(a,f14.6,a)") ' mx2max: ',mx2max,' (m)'
      write(*,"(a,f14.6,a)") ' mx2stp: ',mx2stp,' (m)'
      write(*,*)
      write(*,"(a)") ' AXIS 3'
      if(cartspher==1) then
         write(*,"(a,f14.6,a)") ' mx3min: ',mx3min,' (m)'
         write(*,"(a,f14.6,a)") ' mx3max: ',mx3max,' (m)'
      end if
      if(cartspher==2) then
         write(*,"(a,f14.6,a)") ' mx3min: ',mx3min,' (degrees)'
         write(*,"(a,f14.6,a)") ' mx3max: ',mx3max,' (degrees)'
      end if
      write(*,"(a,f14.6,a)") ' mx3stp: ',mx3stp,' (m)'
      write(*,*)
      write(*,"(a)") '--------------------------------------------------------'
      write(*,*)
   endif

   !LPO parameters
   if (rankMPI .eq. 1 ) then
      write(*,"(a)") 'LPO PARAMETERS' 
      write(*,*)
   endif
   call read_par_int(15,size3)
   if (rankMPI .eq. 1 ) then
      write(*,"(a,i5,a,i5)") ' size3:           ',size3,'  -->   Number of grains for each mineral phase: ',size3*size3*size3
      write(*,*)
   endif

   !Ol + Ens
   if (rankMPI .eq. 1 ) then
      write(*,"(a)") ' ------------------------'
      write(*,"(a)") ' UPPER MANTLE: OL + ENS'
      write(*,"(a)") ' ------------------------'
      write(*,*)
   endif
   call read_par_double(15,Xol(1))
   call read_par_double(15,minx2(1))
   call read_par_double(15,maxx2(1))
   call read_par_double(15,stressexp(1))
   call read_par_double(15,Mob(1))
   call read_par_double(15,chi(1))
   call read_par_double(15,lambda(1))
   call read_par_double(15,fractdislrock(1))
   call read_par_double(15,tau(1,1))
   call read_par_double(15,tau(1,2))
   call read_par_double(15,tau(1,3))
   call read_par_double(15,tau(1,4))
   call read_par_double(15,tau(1,5))
   call read_par_int(15,single_crystal_elastic_db(1,1))
   call read_par_int(15,single_crystal_elastic_db(1,2))
   if (rankMPI .eq. 1 ) then
      write(*,"(a,f14.6,a)") ' Xol(1): ',Xol(1),' (%)'
      write(*,"(a,f14.6,a)") ' minx2(1): ',minx2(1),' (m)'
      write(*,"(a,f14.6,a)") ' maxx2(1): ',maxx2(1),' (m)'
      write(*,"(a,f8.2)")  ' stressexp(1): ',stressexp(1)
      write(*,"(a,f14.2)") ' Mob(1): ',Mob(1)
      write(*,"(a,f14.2)") ' chi(1): ',chi(1)
      write(*,"(a,f11.2)") ' lambda(1): ',lambda(1)
      write(*,"(a,f4.2)")  ' fractdislrock(1): ',fractdislrock(1)
      write(*,"(a,es12.2,a)") ' tau(1,1): ',tau(1,1),' [100](010) Olivine '
      write(*,"(a,es12.2,a)") ' tau(1,2): ',tau(1,2),' [100](001) Olivine '
      write(*,"(a,es12.2,a)") ' tau(1,3): ',tau(1,3),' [001](010) Olivine '
      write(*,"(a,es12.2,a)") ' tau(1,4): ',tau(1,4),' [001](100) Olivine ' 
      write(*,"(a,es12.2,a)") ' tau(1,5): ',tau(1,5),' [001](100) Enstatite'
      write(*,"(a,i5)") ' single_crystal_elastic_db(1,1)  : ',single_crystal_elastic_db(1,1)
      write(*,"(a,i5)") ' single_crystal_elastic_db(1,2)  : ',single_crystal_elastic_db(1,2)
      write(*,*)
   endif

   !Wd + Grt
   if (rankMPI .eq. 1 ) then
      write(*,"(a)") ' ---------------------------------'
      write(*,"(a)") ' UPPER TRANSITION ZONE: WD + GRT'
      write(*,"(a)") ' ---------------------------------'
      write(*,*)
   endif
   call read_par_double(15,Xol(2))
   call read_par_double(15,minx2(2))
   call read_par_double(15,maxx2(2))
   call read_par_double(15,stressexp(2))
   call read_par_double(15,Mob(2))
   call read_par_double(15,chi(2))
   call read_par_double(15,lambda(2))
   call read_par_double(15,fractdislrock(2))
   call read_par_double(15,tau(2,1))
   call read_par_double(15,tau(2,2))
   call read_par_double(15,tau(2,3))
   call read_par_double(15,tau(2,4))
   call read_par_double(15,tau(2,5))
   call read_par_double(15,tau(2,6))
   call read_par_double(15,tau(2,7))
   call read_par_double(15,tau(2,8))
   call read_par_double(15,tau(2,9))
   call read_par_int(15,single_crystal_elastic_db(2,1))
   call read_par_int(15,single_crystal_elastic_db(2,2))
   if (rankMPI .eq. 1 ) then
      write(*,"(a,f14.6,a)") ' Xol(2): ',Xol(2),' (%)'
      write(*,"(a,f14.6,a)") ' minx2(2): ',minx2(2),' (m)'
      write(*,"(a,f14.6,a)") ' maxx2(2): ',maxx2(2),' (m)'
      write(*,"(a,f8.2)")  ' stressexp(2): ',stressexp(2)
      write(*,"(a,f14.2)") ' Mob(2): ',Mob(2)
      write(*,"(a,f14.2)") ' chi(2): ',chi(2)
      write(*,"(a,f11.2)") ' lambda(2): ',lambda(2)
      write(*,"(a,f4.2)")  ' fractdislrock(2): ',fractdislrock(2)
      write(*,"(a,es12.2,a)") ' tau(2,1): ',tau(2,1),' [100](001)   Wadsleyite '
      write(*,"(a,es12.2,a)") ' tau(2,2): ',tau(2,2),' [100](010)   Wadsleyite '
      write(*,"(a,es12.2,a)") ' tau(2,3): ',tau(2,3),' [100](011)   Wadsleyite '
      write(*,"(a,es12.2,a)") ' tau(2,4): ',tau(2,4),' [100](021)   Wadsleyite '
      write(*,"(a,es12.2,a)") ' tau(2,5): ',tau(2,5),' [111](10-1)  Wadsleyite '
      write(*,"(a,es12.2,a)") ' tau(2,6): ',tau(2,6),' [11_1](101)  Wadsleyite '
      write(*,"(a,es12.2,a)") ' tau(2,7): ',tau(2,7),' [1_11](10-1) Wadsleyite '
      write(*,"(a,es12.2,a)") ' tau(2,8): ',tau(2,8),' [1_1_1](101) Wadsleyite '
      write(*,"(a,es12.2,a)") ' tau(2,9): ',tau(2,9),' [001](010)   Wadsleyite '
      write(*,"(a,i5)") ' single_crystal_elastic_db(2,1)  : ',single_crystal_elastic_db(2,1)
      write(*,"(a,i5)") ' single_crystal_elastic_db(2,2)  : ',single_crystal_elastic_db(2,2)
      write(*,*)
    endif

   !Rw + Grt (no LPO parameters for lower TZ aggregates as made by isotropic crystals)
   if (rankMPI .eq. 1 ) then
      write(*,"(a)") ' -----------------------------------------------------'
      write(*,"(a)") ' LOWER TRANSITION ZONE: RW + GRT (NO LPO PARAMETERS)'
      write(*,"(a)") ' -----------------------------------------------------'
      write(*,*)
   endif
   call read_par_double(15,Xol(3))
   call read_par_double(15,minx2(3))
   call read_par_double(15,maxx2(3))
   call read_par_int(15,single_crystal_elastic_db(3,1))
   call read_par_int(15,single_crystal_elastic_db(3,2))
   if (rankMPI .eq. 1 ) then
      write(*,"(a,f14.6,a)") ' Xol(3): ',Xol(3),' (%)'
      write(*,"(a,f14.6,a)") ' minx2(3): ',minx2(3),' (m)'
      write(*,"(a,f14.6,a)") ' maxx2(3): ',maxx2(3),' (m)'
      write(*,"(a,i5)") ' single_crystal_elastic_db(3,1)  : ',single_crystal_elastic_db(3,1)
      write(*,"(a,i5)") ' single_crystal_elastic_db(3,2)  : ',single_crystal_elastic_db(3,2)
      write(*,*)
   endif

   !Brd +  MgO
   if (rankMPI .eq. 1 ) then
      write(*,"(a)") ' -------------------------'
      write(*,"(a)") ' LOWER MANTLE: BRD + MGO'
      write(*,"(a)") ' -------------------------'
      write(*,*)
   endif
   call read_par_double(15,Xol(4))
   call read_par_double(15,minx2(4))
   call read_par_double(15,maxx2(4))
   call read_par_double(15,stressexp(4))
   call read_par_double(15,Mob(4))
   call read_par_double(15,chi(4))
   call read_par_double(15,lambda(4))
   call read_par_double(15,fractdislrock(4))
   call read_par_double(15,tau(4,1))
   call read_par_double(15,tau(4,2))
   call read_par_double(15,tau(4,3))
   call read_par_double(15,tau(4,4))
   call read_par_double(15,tau(4,5))
   call read_par_double(15,tau(4,6))
   call read_par_double(15,tau(4,7))
   call read_par_double(15,tau(4,8))
   call read_par_double(15,tau(4,9))
   call read_par_double(15,tau(4,10))
   call read_par_double(15,tau(4,11))
   call read_par_double(15,tau(4,12))
   call read_par_int(15,single_crystal_elastic_db(4,1))
   call read_par_int(15,single_crystal_elastic_db(4,2))
   if (rankMPI .eq. 1 ) then
      write(*,"(a,f14.6,a)") ' Xol(4): ',Xol(4),' (%)'
      write(*,"(a,f14.6,a)") ' minx2(4): ',minx2(4),' (m)'
      write(*,"(a,f14.6,a)") ' maxx2(4): ',maxx2(4),' (m)'
      write(*,"(a,f8.2)")  ' stressexp(4): ',stressexp(4)
      write(*,"(a,f14.2)") ' Mob(4): ',Mob(4)
      write(*,"(a,f14.2)") ' chi(4): ',chi(4)
      write(*,"(a,f11.2)") ' lambda(4): ',lambda(4)
      write(*,"(a,f4.2)")  ' fractdislrock(4): ',fractdislrock(4)
      write(*,"(a,es12.2,a)") ' tau(4,1): ',tau(4,1),' [100](010)  Bridgmanite '
      write(*,"(a,es12.2,a)") ' tau(4,2): ',tau(4,2),' [100](001)  Bridgmanite '
      write(*,"(a,es12.2,a)") ' tau(4,3): ',tau(4,3),' [010](100)  Bridgmanite '
      write(*,"(a,es12.2,a)") ' tau(4,4): ',tau(4,4),' [010](001)  Bridgmanite '
      write(*,"(a,es12.2,a)") ' tau(4,5): ',tau(4,5),' [001](100)  Bridgmanite '
      write(*,"(a,es12.2,a)") ' tau(4,6): ',tau(4,6),' [001](010)  Bridgmanite '
      write(*,"(a,es12.2,a)") ' tau(4,7): ',tau(4,7),' [001](110)  Bridgmanite '
      write(*,"(a,es12.2,a)") ' tau(4,8): ',tau(4,8),' [001](-110) Bridgmanite '
      write(*,"(a,es12.2,a)") ' tau(4,9): ',tau(4,9),' [110](001)  Bridgmanite '
      write(*,"(a,es11.2,a)") ' tau(4,10): ',tau(4,10),' [-110](001) Bridgmanite '
      write(*,"(a,es11.2,a)") ' tau(4,11): ',tau(4,11),' [110](-110) Bridgmanite '
      write(*,"(a,es11.2,a)") ' tau(4,12): ',tau(4,12),' [-110](110) Bridgmanite '
      write(*,"(a,i5)") ' single_crystal_elastic_db(4,1)  : ',single_crystal_elastic_db(4,1)
      write(*,"(a,i5)") ' single_crystal_elastic_db(4,2)  : ',single_crystal_elastic_db(4,2)
      write(*,*)
      write(*,"(a)") '--------------------------------------------------------'
      write(*,*)
   endif

   !pPv +  MgO
   if (rankMPI .eq. 1 ) then
      write(*,"(a)") ' -------------------------'
      write(*,"(a)") ' LOWER MANTLE: PPV + MGO'
      write(*,"(a)") ' -------------------------'
      write(*,*)
   endif
   call read_par_double(15,Xol(5))           ;
   call read_par_double(15,stressexp(5))     ;
   call read_par_double(15,Mob(5))           ;
   call read_par_double(15,chi(5))           ;
   call read_par_double(15,lambda(5))        ;
   call read_par_double(15,fractdislrock(5)) ;
   call read_par_double(15,tau(5,1))         ;
   call read_par_double(15,tau(5,2))         ;
   call read_par_double(15,tau(5,3))         ;
   call read_par_double(15,tau(5,4))         ;
   call read_par_double(15,tau(5,5))         ;
   call read_par_double(15,tau(5,6))         ;
   call read_par_double(15,tau(5,7))         ;
   call read_par_double(15,tau(5,8))         ;
   call read_par_double(15,tau(5,9))         ;
   call read_par_double(15,tau(5,10))         ;
   call read_par_int(15,single_crystal_elastic_db(5,1))        ;
   call read_par_int(15,single_crystal_elastic_db(5,2))        ;
   if (rankMPI .eq. 1 ) then
      write(*,"(a,f8.2,a)")  ,' Xol(5):       ',Xol(5),' (%)'
      write(*,"(a,f8.2)")    ,' stressexp(5): ',stressexp(5)
      write(*,"(a,f14.2)")   ,' Mob(5): ',Mob(5)
      write(*,"(a,f14.2)")   ,' chi(5): ',chi(5)
      write(*,"(a,f11.2)")   ,' lambda(5): ',lambda(5)
      write(*,"(a,f4.2)")    ,' fractdislrock(5): ',fractdislrock(5)
      write(*,"(a,es12.2,a)"),' tau(5,1): ',tau(5,1),' [100](010)  PPv '
      write(*,"(a,es12.2,a)"),' tau(5,2): ',tau(5,2),' [100](001)  PPv '
      write(*,"(a,es12.2,a)"),' tau(5,3): ',tau(5,3),' [010](100)  PPv '
      write(*,"(a,es12.2,a)"),' tau(5,4): ',tau(5,4),' [010](001)  PPv '
      write(*,"(a,es12.2,a)"),' tau(5,5): ',tau(5,5),' [001](100)  PPv '
      write(*,"(a,es12.2,a)"),' tau(5,6): ',tau(5,6),' [001](010)  PPv '
      write(*,"(a,es12.2,a)"),' tau(5,7): ',tau(5,7),' [001](110)  PPv '
      write(*,"(a,es12.2,a)"),' tau(5,8): ',tau(5,8),' [001](-110) PPv '
      write(*,"(a,es12.2,a)"),' tau(5,9): ',tau(5,9),' [110](1-10)  PPv '
      write(*,"(a,es11.2,a)"),' tau(5,10): ',tau(5,10),' [-110](110) PPv '
      write(*,"(a,i5)")      ,' single_crystal_elastic_db(5,1)  : ',single_crystal_elastic_db(5,1)
      write(*,"(a,i5)")      ,' single_crystal_elastic_db(5,2)  : ',single_crystal_elastic_db(5,2)
      write(*,*)
   endif

   !Pre-existing fabric
   if (rankMPI .eq. 1 ) then
      write(*,"(a)") 'SET PRE-EXISTING FABRIC'       
      write(*,*)
   endif
   call read_par_int(15,fossilfabric)
   call read_par_double(15,mx1minfab)
   call read_par_double(15,mx1maxfab)
   call read_par_double(15,mx2minfab)
   call read_par_double(15,mx2maxfab)
   call read_par_double(15,mx3minfab)
   call read_par_double(15,mx3maxfab)
   if (rankMPI .eq. 1 ) then
      write(*,"(a,i5)") ' fossilfabric:       ',fossilfabric
      if(fossilfabric /= 0) then
         if(cartspher==1) then
            write(*,"(a,f14.6,a)") ' mx1minfab: ',mx1minfab,' (m)'
            write(*,"(a,f14.6,a)") ' mx1maxfab: ',mx1maxfab,' (m)'
         end if
         if(cartspher==2) then
            write(*,"(a,f14.6,a)") ' mx1minfab: ',mx1minfab,' (degrees)'
            write(*,"(a,f14.6,a)") ' mx1maxfab: ',mx1maxfab,' (degrees)'
         end if
         write(*,"(a,f14.6,a)") ' mx2minfab: ',mx2minfab,' (m)'
         write(*,"(a,f14.6,a)") ' mx2maxfab: ',mx2maxfab,' (m)'
         if(cartspher==1) then
            write(*,"(a,f14.6,a)") ' mx3minfab: ',mx3minfab,' (m)'
            write(*,"(a,f14.6,a)") ' mx3maxfab: ',mx3maxfab,' (m)'
         end if
         if(cartspher==2) then
            write(*,"(a,f14.6,a)") ' mx3minfab: ',mx3minfab,' (degrees)'
            write(*,"(a,f14.6,a)") ' mx3maxfab: ',mx3maxfab,' (degrees)'
         end if
      end if
      write(*,*)
      write(*,"(a)") '--------------------------------------------------------'
      write(*,*)
   end if
   
   !Set operating modes        
   if (rankMPI .eq. 1 ) then
      write(*,"(a)") 'SET OPERATING MODES'       
      write(*,*)
   endif
   !call read_par_int(15,sbfmod)
   !call read_par_double(15,rmax)
   sbfmod = 0
   call read_par_int(15,fsemod)
   call read_par_int(15,uppermantlemod)
   IF(sbfmod > 0) THEN
      uppermantlemod = 1
      fsemod = 1
   END IF

   if (rankMPI .eq. 1 ) then
      write(*,"(a,i5)") ' sbfmod:             ',sbfmod
      IF(sbfmod > 0) write(*,"(a,f10.2)") ' rmax:             ',rmax
      write(*,"(a,i5)") ' fsemod:             ',fsemod
      write(*,"(a,i5)") ' uppermantlemod:     ',uppermantlemod
   endif
   IF(fsemod == 0 .OR. sbfmod > 0) THEN

      call read_par_int(15,fractdislmod)
      call read_par_int(15,fabrictransformmod)
      call read_par_int(15,ptmod)
      call read_par_int(15,eosmod)
      call read_par_double(15,fractvoigt)
      fractvoigt = fractvoigt/100d0

      IF(sbfmod > 0) fractdislmod = 0

      if (rankMPI .eq. 1 ) then
         write(*,"(a,i5)") ' fractdislmod:       ',fractdislmod
         write(*,"(a,i5)") ' fabrictransformmod: ',fabrictransformmod
         write(*,"(a,i5)") ' ptmod:              ',ptmod
         write(*,"(a,i5)") ' eosmod:             ',eosmod
         write(*,"(a,f9.2,a)") ' fractvoigt:     ',fractvoigt,' (%)'
      endif

      IF(ptmod > 0 .AND. (eosmod < 1 .OR. eosmod > 5)) THEN
         write(*,"(a,i5)") '!!!No database for the selected eosmod. Set 1 <= eosmod <= 5 !!!!' 
         stop
      END IF

   ELSE

      !Reset operating modes
      fractdislmod = 0; 
      fabrictransformmod = 0; 
      ptmod = 0; 
      eosmod = 0; 
      fractvoigt = 1; 
      if (rankMPI .eq. 1 ) then
         write(*,"(a,i5)") '!!!WARNING: fsemod > 0. NOT COMPUTING LPO !!!!' 
      endif

   END IF

   if (rankMPI .eq. 1 ) then
      write(*,*)
      write(*,"(a)") '********************************************************'
      write(*,"(a)") '********************************************************'
      write(*,*)
   endif

   CLOSE(15)

   END SUBROUTINE read_input_file
