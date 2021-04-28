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

   SUBROUTINE read_input_file(input_dir,output_dir)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i
   CHARACTER (len=500) :: input_dir,output_dir
   CHARACTER(200) :: arg,inputfilename,str

   DO i = 1, iargc()
      CALL getarg(i, arg)
      if(i==1) inputfilename = arg
   END DO
 
   OPEN(unit=15,file=inputfilename)

   !Input files directory
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)
   write(*,"(a,a)"),'READING INPUT FILE: ',trim(inputfilename)
   write(*,*)
   call read_par_char(15,input_dir) ; write(*,"(a,a)"),' Input  directory: ',trim(input_dir)
   write(*,*)
   call read_par_char(15,output_dir); write(*,"(a,a)"),' Output directory: ',trim(output_dir)
   write(*,*)

   !Input files min,step,max, output frequency, max. time (in Myr) in
   !steady-state conditions
   write(*,"(a)"),'INPUT - OUTPUT FILES'
   write(*,*)
   call read_par_int(15,Tinit)      ; write(*,"(a,i5)"),' Tinit :      ',Tinit
   call read_par_int(15,Tstep)      ; write(*,"(a,i5)"),' Tstep :      ',Tstep
   call read_par_int(15,Tend)       ; write(*,"(a,i5)"),' Tend  :      ',Tend
   call read_par_int(15,OutputStep) ; write(*,"(a,i5)"),' Outputstep : ',OutputStep
   call read_par_double(15,timemax)
   write(*,*)
   if(Tinit == Tend) write(*,"(a,f6.2,a)"), ' Tinit = Tend --> steady-state flow for ',Timemax,' Myr'
   write(*,*)
   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)

   !Define model grid
   write(*,"(a)"),'DEFINE MODEL GRIDS'
   write(*,*)
   call read_par_int(15,dimensions)
   if(dimensions==2) write(*,"(a,i5,a)"),' dimensions: ',dimensions,' --> 2D model'
   if(dimensions==3) write(*,"(a,i5,a)"),' dimensions: ',dimensions,' --> 3D model'
   call read_par_int(15,cartspher)
   if(cartspher==1) write(*,"(a,i5,a)"),' cartspher : ',cartspher,' --> cartesian coordinates system'
   if(cartspher==2 .AND. dimensions == 2) write(*,"(a,i5,a)"),' cartspher : ',cartspher,' --> polar coordinates system'
   if(cartspher==2 .AND. dimensions == 3) write(*,"(a,i5,a)"),' cartspher : ',cartspher,' --> spherical coordinates system'
   call read_par_int(15,basicstag) ; write(*,"(a,i5)"),' basicstag : ',basicstag
   write(*,*)
   write(*,"(a)"),'EULERIAN GRID'
   write(*,*)
   write(*,"(a)"),' AXIS 1'
   call read_par_double(15,x1min)   
   call read_par_double(15,x1max)   
   if(cartspher==1) then
      write(*,"(a,f14.6,a)"),' x1min: ',x1min,' (m)'
      write(*,"(a,f14.6,a)"),' x1max: ',x1max,' (m)'
   end if
   if(cartspher==2) then
      write(*,"(a,f14.6,a)"),' x1min: ',x1min,' (degrees)'
      write(*,"(a,f14.6,a)"),' x1max: ',x1max,' (degrees)'
   end if
   call read_par_int(15,nx1)        ; write(*,"(a,i5)"),   ' nx1  :          ',nx1
   call read_par_int(15,x1periodic) ; write(*,"(a,i5)"),   ' x1periodic:     ',x1periodic
   write(*,*)
   write(*,"(a)"),' AXIS 2'
   call read_par_double(15,x2min)   ; write(*,"(a,f14.6,a)"),' x2min: ',x2min,' (m)'
   call read_par_double(15,x2max)   ; write(*,"(a,f14.6,a)"),' x2max: ',x2max,' (m)'
   call read_par_int(15,nx2)        ; write(*,"(a,i5)"),   ' nx2  :          ',nx2
   call read_par_int(15,x2periodic) ; write(*,"(a,i5)"),   ' x2periodic:     ',x2periodic
   write(*,*)
   write(*,"(a)"),' AXIS 3'
   call read_par_double(15,x3min)   
   call read_par_double(15,x3max)   
   if(cartspher==1) then
      write(*,"(a,f14.6,a)"),' x3min: ',x3min,' (m)'
      write(*,"(a,f14.6,a)"),' x3max: ',x3max,' (m)'
   end if
   if(cartspher==2) then
      write(*,"(a,f14.6,a)"),' x3min: ',x3min,' (degrees)'
      write(*,"(a,f14.6,a)"),' x3max: ',x3max,' (degrees)'
   end if
   call read_par_int(15,nx3)        ; write(*,"(a,i5)"),   ' nx3  :          ',nx3
   call read_par_int(15,x3periodic) ; write(*,"(a,i5)"),   ' x3periodic:     ',x3periodic
   write(*,*)

   yinyang = 1
   if(cartspher == 2 .AND. dimensions == 3 .AND. x1periodic > 0 .AND. x3periodic > 0 ) THEN
      yinyang = 2
      write(*,"(a)"),' !!!!!!! YIN-YANG GRID ACTIVE !!!!!!!!!'
      write(*,*)
   end if

   !!! Set number of nodes along 3rd axis equal to 1 if model is 2D
   if(dimensions == 2) nx3 = 1

   !Aggregate distribution
   write(*,"(a)"),'LAGRANGIAN GRID'
   write(*,*)
   write(*,"(a)"),' AXIS 1'
   call read_par_double(15,mx1min)   
   call read_par_double(15,mx1max)   
   if(cartspher==1) then
      write(*,"(a,f14.6,a)"),' mx1min: ',mx1min,' (m)'
      write(*,"(a,f14.6,a)"),' mx1max: ',mx1max,' (m)'
   end if
   if(cartspher==2) then
      write(*,"(a,f14.6,a)"),' mx1min: ',mx1min,' (degrees)'
      write(*,"(a,f14.6,a)"),' mx1max: ',mx1max,' (degrees)'
   end if
   call read_par_double(15,mx1stp)   ; write(*,"(a,f14.6,a)"),' mx1stp: ',mx1stp,' (m)'
   write(*,*)
   write(*,"(a)"),' AXIS 2'
   call read_par_double(15,mx2min)   ; write(*,"(a,f14.6,a)"),' mx2min: ',mx2min,' (m)'
   call read_par_double(15,mx2max)   ; write(*,"(a,f14.6,a)"),' mx2max: ',mx2max,' (m)'
   call read_par_double(15,mx2stp)   ; write(*,"(a,f14.6,a)"),' mx2stp: ',mx2stp,' (m)'
   write(*,*)
   write(*,"(a)"),' AXIS 3'
   call read_par_double(15,mx3min)   
   call read_par_double(15,mx3max)   
   if(cartspher==1) then
      write(*,"(a,f14.6,a)"),' mx3min: ',mx3min,' (m)'
      write(*,"(a,f14.6,a)"),' mx3max: ',mx3max,' (m)'
   end if
   if(cartspher==2) then
      write(*,"(a,f14.6,a)"),' mx3min: ',mx3min,' (degrees)'
      write(*,"(a,f14.6,a)"),' mx3max: ',mx3max,' (degrees)'
   end if
   call read_par_double(15,mx3stp)   ; write(*,"(a,f14.6,a)"),' mx3stp: ',mx3stp,' (m)'
   write(*,*)
   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)

   !LPO parameters
   write(*,"(a)"),'LPO PARAMETERS' 
   write(*,*)
   call read_par_int(15,size3)        ; write(*,"(a,i5,a,i5)"),   ' size3:           ',size3,'  -->   Number of grains for each mineral phase: ',size3*size3*size3
   write(*,*)

   !Ol + Ens
   write(*,"(a)"),' ------------------------'
   write(*,"(a)"),' UPPER MANTLE: OL + ENS'
   write(*,"(a)"),' ------------------------'
   write(*,*)
   call read_par_double(15,Xol(1))           ; write(*,"(a,f14.6,a)"),' Xol(1): ',Xol(1),' (%)'
   call read_par_double(15,minx2(1))         ; write(*,"(a,f14.6,a)"),' minx2(1): ',minx2(1),' (m)'
   call read_par_double(15,maxx2(1))         ; write(*,"(a,f14.6,a)"),' maxx2(1): ',maxx2(1),' (m)'
   call read_par_double(15,stressexp(1))     ; write(*,"(a,f8.2)") ,' stressexp(1): ',stressexp(1)
   call read_par_double(15,Mob(1))           ; write(*,"(a,f14.2)"),' Mob(1): ',Mob(1)
   call read_par_double(15,chi(1))           ; write(*,"(a,f14.2)"),' chi(1): ',chi(1)
   call read_par_double(15,lambda(1))        ; write(*,"(a,f11.2)"),' lambda(1): ',lambda(1)
   call read_par_double(15,fractdislrock(1)) ; write(*,"(a,f4.2)") ,' fractdislrock(1): ',fractdislrock(1)
   call read_par_double(15,tau(1,1))         ; write(*,"(a,es12.2,a)"),' tau(1,1): ',tau(1,1),' [100](010) Olivine '
   call read_par_double(15,tau(1,2))         ; write(*,"(a,es12.2,a)"),' tau(1,2): ',tau(1,2),' [100](001) Olivine '
   call read_par_double(15,tau(1,3))         ; write(*,"(a,es12.2,a)"),' tau(1,3): ',tau(1,3),' [001](010) Olivine '
   call read_par_double(15,tau(1,4))         ; write(*,"(a,es12.2,a)"),' tau(1,4): ',tau(1,4),' [001](100) Olivine ' 
   call read_par_double(15,tau(1,5))         ; write(*,"(a,es12.2,a)"),' tau(1,5): ',tau(1,5),' [001](100) Enstatite'
   call read_par_int(15,single_crystal_elastic_db(1,1))        ; write(*,"(a,i5)"),   ' single_crystal_elastic_db(1,1)  : ',single_crystal_elastic_db(1,1)
   call read_par_int(15,single_crystal_elastic_db(1,2))        ; write(*,"(a,i5)"),   ' single_crystal_elastic_db(1,2)  : ',single_crystal_elastic_db(1,2)
   write(*,*)

   !Wd + Grt
   write(*,"(a)"),' ---------------------------------'
   write(*,"(a)"),' UPPER TRANSITION ZONE: WD + GRT'
   write(*,"(a)"),' ---------------------------------'
   write(*,*)
   call read_par_double(15,Xol(2))           ; write(*,"(a,f14.6,a)"),' Xol(2): ',Xol(2),' (%)'
   call read_par_double(15,minx2(2))         ; write(*,"(a,f14.6,a)"),' minx2(2): ',minx2(2),' (m)'
   call read_par_double(15,maxx2(2))         ; write(*,"(a,f14.6,a)"),' maxx2(2): ',maxx2(2),' (m)'
   call read_par_double(15,stressexp(2))     ; write(*,"(a,f8.2)") ,' stressexp(2): ',stressexp(2)
   call read_par_double(15,Mob(2))           ; write(*,"(a,f14.2)"),' Mob(2): ',Mob(2)
   call read_par_double(15,chi(2))           ; write(*,"(a,f14.2)"),' chi(2): ',chi(2)
   call read_par_double(15,lambda(2))        ; write(*,"(a,f11.2)"),' lambda(2): ',lambda(2)
   call read_par_double(15,fractdislrock(2)) ; write(*,"(a,f4.2)") ,' fractdislrock(2): ',fractdislrock(2)
   call read_par_double(15,tau(2,1))         ; write(*,"(a,es12.2,a)"),' tau(2,1): ',tau(2,1),' [100](001)   Wadsleyite '
   call read_par_double(15,tau(2,2))         ; write(*,"(a,es12.2,a)"),' tau(2,2): ',tau(2,2),' [100](010)   Wadsleyite '
   call read_par_double(15,tau(2,3))         ; write(*,"(a,es12.2,a)"),' tau(2,3): ',tau(2,3),' [100](011)   Wadsleyite '
   call read_par_double(15,tau(2,4))         ; write(*,"(a,es12.2,a)"),' tau(2,4): ',tau(2,4),' [100](021)   Wadsleyite '
   call read_par_double(15,tau(2,5))         ; write(*,"(a,es12.2,a)"),' tau(2,5): ',tau(2,5),' [111](10-1)  Wadsleyite '
   call read_par_double(15,tau(2,6))         ; write(*,"(a,es12.2,a)"),' tau(2,6): ',tau(2,6),' [11_1](101)  Wadsleyite '
   call read_par_double(15,tau(2,7))         ; write(*,"(a,es12.2,a)"),' tau(2,7): ',tau(2,7),' [1_11](10-1) Wadsleyite '
   call read_par_double(15,tau(2,8))         ; write(*,"(a,es12.2,a)"),' tau(2,8): ',tau(2,8),' [1_1_1](101) Wadsleyite '
   call read_par_double(15,tau(2,9))         ; write(*,"(a,es12.2,a)"),' tau(2,9): ',tau(2,9),' [001](010)   Wadsleyite '
   call read_par_int(15,single_crystal_elastic_db(2,1))        ; write(*,"(a,i5)"),   ' single_crystal_elastic_db(2,1)  : ',single_crystal_elastic_db(2,1)
   call read_par_int(15,single_crystal_elastic_db(2,2))        ; write(*,"(a,i5)"),   ' single_crystal_elastic_db(2,2)  : ',single_crystal_elastic_db(2,2)
   write(*,*)

   !Rw + Grt (no LPO parameters for lower TZ aggregates as made by isotropic crystals)
   write(*,"(a)"),' -----------------------------------------------------'
   write(*,"(a)"),' LOWER TRANSITION ZONE: RW + GRT (NO LPO PARAMETERS)'
   write(*,"(a)"),' -----------------------------------------------------'
   write(*,*)
   call read_par_double(15,Xol(3))           ; write(*,"(a,f14.6,a)"),' Xol(3): ',Xol(3),' (%)'
   call read_par_double(15,minx2(3))         ; write(*,"(a,f14.6,a)"),' minx2(3): ',minx2(3),' (m)'
   call read_par_double(15,maxx2(3))         ; write(*,"(a,f14.6,a)"),' maxx2(3): ',maxx2(3),' (m)'
   call read_par_int(15,single_crystal_elastic_db(3,1))        ; write(*,"(a,i5)"),   ' single_crystal_elastic_db(3,1)  : ',single_crystal_elastic_db(3,1)
   call read_par_int(15,single_crystal_elastic_db(3,2))        ; write(*,"(a,i5)"),   ' single_crystal_elastic_db(3,2)  : ',single_crystal_elastic_db(3,2)
   write(*,*)

   !Brd +  MgO
   write(*,"(a)"),' -------------------------'
   write(*,"(a)"),' LOWER MANTLE: BRD + MGO'
   write(*,"(a)"),' -------------------------'
   write(*,*)
   call read_par_double(15,Xol(4))           ; write(*,"(a,f14.6,a)"),' Xol(4): ',Xol(4),' (%)'
   call read_par_double(15,minx2(4))         ; write(*,"(a,f14.6,a)"),' minx2(4): ',minx2(4),' (m)'
   call read_par_double(15,maxx2(4))         ; write(*,"(a,f14.6,a)"),' maxx2(4): ',maxx2(4),' (m)'
   call read_par_double(15,stressexp(4))     ; write(*,"(a,f8.2)") ,' stressexp(4): ',stressexp(4)
   call read_par_double(15,Mob(4))           ; write(*,"(a,f14.2)"),' Mob(4): ',Mob(4)
   call read_par_double(15,chi(4))           ; write(*,"(a,f14.2)"),' chi(4): ',chi(4)
   call read_par_double(15,lambda(4))        ; write(*,"(a,f11.2)"),' lambda(4): ',lambda(4)
   call read_par_double(15,fractdislrock(4)) ; write(*,"(a,f4.2)") ,' fractdislrock(4): ',fractdislrock(4)
   call read_par_double(15,tau(4,1))         ; write(*,"(a,es12.2,a)"),' tau(4,1): ',tau(4,1),' [100](010)  Bridgmanite '
   call read_par_double(15,tau(4,2))         ; write(*,"(a,es12.2,a)"),' tau(4,2): ',tau(4,2),' [100](001)  Bridgmanite '
   call read_par_double(15,tau(4,3))         ; write(*,"(a,es12.2,a)"),' tau(4,3): ',tau(4,3),' [010](100)  Bridgmanite '
   call read_par_double(15,tau(4,4))         ; write(*,"(a,es12.2,a)"),' tau(4,4): ',tau(4,4),' [010](001)  Bridgmanite '
   call read_par_double(15,tau(4,5))         ; write(*,"(a,es12.2,a)"),' tau(4,5): ',tau(4,5),' [001](100)  Bridgmanite '
   call read_par_double(15,tau(4,6))         ; write(*,"(a,es12.2,a)"),' tau(4,6): ',tau(4,6),' [001](010)  Bridgmanite '
   call read_par_double(15,tau(4,7))         ; write(*,"(a,es12.2,a)"),' tau(4,7): ',tau(4,7),' [001](110)  Bridgmanite '
   call read_par_double(15,tau(4,8))         ; write(*,"(a,es12.2,a)"),' tau(4,8): ',tau(4,8),' [001](-110) Bridgmanite '
   call read_par_double(15,tau(4,9))         ; write(*,"(a,es12.2,a)"),' tau(4,9): ',tau(4,9),' [110](001)  Bridgmanite '
   call read_par_double(15,tau(4,10))         ; write(*,"(a,es11.2,a)"),' tau(4,10): ',tau(4,10),' [-110](001) Bridgmanite '
   call read_par_double(15,tau(4,11))         ; write(*,"(a,es11.2,a)"),' tau(4,11): ',tau(4,11),' [110](-110) Bridgmanite '
   call read_par_double(15,tau(4,12))         ; write(*,"(a,es11.2,a)"),' tau(4,12): ',tau(4,12),' [-110](110) Bridgmanite '
   call read_par_int(15,single_crystal_elastic_db(4,1))        ; write(*,"(a,i5)"),   ' single_crystal_elastic_db(4,1)  : ',single_crystal_elastic_db(4,1)
   call read_par_int(15,single_crystal_elastic_db(4,2))        ; write(*,"(a,i5)"),   ' single_crystal_elastic_db(4,2)  : ',single_crystal_elastic_db(4,2)
   write(*,*)
   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)

   !Pre-existing fabric
   write(*,"(a)"),'SET PRE-EXISTING FABRIC'       
   write(*,*)
   call read_par_int(15,fossilfabric)        ; write(*,"(a,i5)"),' fossilfabric:       ',fossilfabric
   call read_par_double(15,mx1minfab)   
   call read_par_double(15,mx1maxfab)   
   call read_par_double(15,mx2minfab)   
   call read_par_double(15,mx2maxfab)   
   call read_par_double(15,mx3minfab)   
   call read_par_double(15,mx3maxfab)   
   if(fossilfabric /= 0) then

      if(cartspher==1) then
         write(*,"(a,f14.6,a)"),' mx1minfab: ',mx1minfab,' (m)'
         write(*,"(a,f14.6,a)"),' mx1maxfab: ',mx1maxfab,' (m)'
      end if
      if(cartspher==2) then
         write(*,"(a,f14.6,a)"),' mx1minfab: ',mx1minfab,' (degrees)'
         write(*,"(a,f14.6,a)"),' mx1maxfab: ',mx1maxfab,' (degrees)'
      end if
      write(*,"(a,f14.6,a)"),' mx2minfab: ',mx2minfab,' (m)'
      write(*,"(a,f14.6,a)"),' mx2maxfab: ',mx2maxfab,' (m)'
      if(cartspher==1) then
         write(*,"(a,f14.6,a)"),' mx3minfab: ',mx3minfab,' (m)'
         write(*,"(a,f14.6,a)"),' mx3maxfab: ',mx3maxfab,' (m)'
      end if
      if(cartspher==2) then
         write(*,"(a,f14.6,a)"),' mx3minfab: ',mx3minfab,' (degrees)'
         write(*,"(a,f14.6,a)"),' mx3maxfab: ',mx3maxfab,' (degrees)'
      end if
   end if

   write(*,*)
   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)
   
   !Set operating modes        
   write(*,"(a)"),'SET OPERATING MODES'       
   write(*,*)
   call read_par_int(15,fsemod)              ; write(*,"(a,i5)"),' fsemod:             ',fsemod
   call read_par_int(15,uppermantlemod)      ; write(*,"(a,i5)"),' uppermantlemod:     ',uppermantlemod
   IF(fsemod == 0) THEN

      call read_par_int(15,fractdislmod)        ; write(*,"(a,i5)"),' fractdislmod:       ',fractdislmod
      call read_par_int(15,fabrictransformmod)  ; write(*,"(a,i5)"),' fabrictransformmod: ',fabrictransformmod
      call read_par_int(15,ptmod)               ; write(*,"(a,i5)"),' ptmod:              ',ptmod
      call read_par_int(15,eosmod)              ; write(*,"(a,i5)"),' eosmod:             ',eosmod
      call read_par_double(15,fractvoigt)       ; write(*,"(a,f9.2,a)"),' fractvoigt:     ',fractvoigt,' (%)'
      fractvoigt = fractvoigt/100d0

      IF(ptmod > 0 .AND. (eosmod < 1 .OR. eosmod > 5)) THEN
         write(*,"(a,i5)"),'!!!No database for the selected eosmod. Set 1 <= eosmod <= 5 !!!!' 
         stop
      END IF

   ELSE

      !Reset operating modes
      fractdislmod = 0; 
      fabrictransformmod = 0; 
      ptmod = 0; 
      eosmod = 0; 
      fractvoigt = 1; 
      write(*,"(a,i5)"),'!!!WARNING: fsemod > 0. NOT COMPUTING LPO !!!!' 

   END IF

   write(*,*)
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)

   CLOSE(15)

   END SUBROUTINE read_input_file

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
