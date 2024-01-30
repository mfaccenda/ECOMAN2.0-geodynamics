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

   SUBROUTINE read_input_file

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
   call read_par_char(15,output_name); write(*,"(a,a,a)"),' Output file name: ',trim(output_name),'*.h5'
   write(*,*)

   !Input files min,step,max, output frequency, max. time (in Myr) in
   !steady-state conditions
   call read_par_double(15,strain_max); write(*,"(a,f6.2)"),' Max. strain = ',strain_max
   write(*,*)
   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)

   !Scaling factors
   call read_par_double(15,l(1,1,1))      ;
   call read_par_double(15,l(1,2,2))      ;
   call read_par_double(15,l(1,3,3))      ;
   call read_par_double(15,l(1,1,2))      ;
   call read_par_double(15,l(1,2,1))      ;
   call read_par_double(15,l(1,1,3))      ;
   call read_par_double(15,l(1,3,1))      ;
   call read_par_double(15,l(1,2,3))      ;
   call read_par_double(15,l(1,3,2))      ; 
   
   !LPO parameters
   write(*,"(a)"),'LPO PARAMETERS' 
   write(*,*)
   call read_par_int(15,size3)
   write(*,"(a,i5,a,i5)") ' size3:           ',size3,'  -->   Number of grains for each mineral phase: ',size3*size3*size3
   write(*,*)

   !Rocktype
   call read_par_int(15,rocktype(1))        ; 

   !Ol + Ens
   call read_par_double(15,Xol(1))           ; 
   call read_par_double(15,stressexp(1))     ; 
   call read_par_double(15,Mob(1))           ; 
   call read_par_double(15,chi(1))           ; 
   call read_par_double(15,lambda(1))        ; 
   call read_par_double(15,fractdislrock(1)) ; 
   call read_par_double(15,tau(1,1))         ; 
   call read_par_double(15,tau(1,2))         ; 
   call read_par_double(15,tau(1,3))         ; 
   call read_par_double(15,tau(1,4))         ; 
   call read_par_double(15,tau(1,5))         ; 
   call read_par_int(15,single_crystal_elastic_db(1,1))        ; 
   call read_par_int(15,single_crystal_elastic_db(1,2))        ; 
   if(rocktype(1) == 1) then
      write(*,"(a)")         ,' ------------------------'
      write(*,"(a)")         ,' UPPER MANTLE: OL + ENS'
      write(*,"(a)")         ,' ------------------------'
      write(*,*)
      write(*,"(a,f8.2,a)")  ,' Xol(1):       ',Xol(1),' (%)'
      write(*,"(a,f8.2)")    ,' stressexp(1): ',stressexp(1)
      write(*,"(a,f14.2)")   ,' Mob(1): ',Mob(1)
      write(*,"(a,f14.2)")   ,' chi(1): ',chi(1)
      write(*,"(a,f11.2)")   ,' lambda(1): ',lambda(1)
      write(*,"(a,f4.2)")    ,' fractdislrock(1): ',fractdislrock(1)
      write(*,"(a,es12.2,a)"),' tau(1,1): ',tau(1,1),' [100](010) Olivine '
      write(*,"(a,es12.2,a)"),' tau(1,2): ',tau(1,2),' [100](001) Olivine '
      write(*,"(a,es12.2,a)"),' tau(1,3): ',tau(1,3),' [001](010) Olivine '
      write(*,"(a,es12.2,a)"),' tau(1,4): ',tau(1,4),' [001](100) Olivine ' 
      write(*,"(a,es12.2,a)"),' tau(1,5): ',tau(1,5),' [001](100) Enstatite'
      write(*,"(a,i5)")      ,' single_crystal_elastic_db(1,1)  : ',single_crystal_elastic_db(1,1)
      write(*,"(a,i5)")      ,' single_crystal_elastic_db(1,2)  : ',single_crystal_elastic_db(1,2)
      write(*,*)
   end if

   !Wd + Grt
   call read_par_double(15,Xol(2))           ; 
   call read_par_double(15,stressexp(2))     ; 
   call read_par_double(15,Mob(2))           ; 
   call read_par_double(15,chi(2))           ; 
   call read_par_double(15,lambda(2))        ; 
   call read_par_double(15,fractdislrock(2)) ; 
   call read_par_double(15,tau(2,1))         ; 
   call read_par_double(15,tau(2,2))         ; 
   call read_par_double(15,tau(2,3))         ; 
   call read_par_double(15,tau(2,4))         ; 
   call read_par_double(15,tau(2,5))         ; 
   call read_par_double(15,tau(2,6))         ; 
   call read_par_double(15,tau(2,7))         ; 
   call read_par_double(15,tau(2,8))         ; 
   call read_par_double(15,tau(2,9))         ; 
   call read_par_int(15,single_crystal_elastic_db(2,1))        ; 
   call read_par_int(15,single_crystal_elastic_db(2,2))        ; 
   if(rocktype(1) == 2) then
      write(*,"(a)")         ,' ---------------------------------'
      write(*,"(a)")         ,' UPPER TRANSITION ZONE: WD + GRT'
      write(*,"(a)")         ,' ---------------------------------'
      write(*,*)
      write(*,"(a,f8.2,a)")  ,' Xol(2):       ',Xol(2),' (%)'
      write(*,"(a,f8.2)")    ,' stressexp(2): ',stressexp(2)
      write(*,"(a,f14.2)")   ,' Mob(2): ',Mob(2)
      write(*,"(a,f14.2)")   ,' chi(2): ',chi(2)
      write(*,"(a,f11.2)")   ,' lambda(2): ',lambda(2)
      write(*,"(a,f4.2)")    ,' fractdislrock(2): ',fractdislrock(2)
      write(*,"(a,es12.2,a)"),' tau(2,1): ',tau(2,1),' [100](001)   Wadsleyite '
      write(*,"(a,es12.2,a)"),' tau(2,2): ',tau(2,2),' [100](010)   Wadsleyite '
      write(*,"(a,es12.2,a)"),' tau(2,3): ',tau(2,3),' [100](011)   Wadsleyite '
      write(*,"(a,es12.2,a)"),' tau(2,4): ',tau(2,4),' [100](021)   Wadsleyite '
      write(*,"(a,es12.2,a)"),' tau(2,5): ',tau(2,5),' [111](10-1)  Wadsleyite '
      write(*,"(a,es12.2,a)"),' tau(2,6): ',tau(2,6),' [11_1](101)  Wadsleyite '
      write(*,"(a,es12.2,a)"),' tau(2,7): ',tau(2,7),' [1_11](10-1) Wadsleyite '
      write(*,"(a,es12.2,a)"),' tau(2,8): ',tau(2,8),' [1_1_1](101) Wadsleyite '
      write(*,"(a,es12.2,a)"),' tau(2,9): ',tau(2,9),' [001](010)   Wadsleyite '
      write(*,"(a,i5)")      ,' single_crystal_elastic_db(2,1)  : ',single_crystal_elastic_db(2,1)
      write(*,"(a,i5)")      ,' single_crystal_elastic_db(2,2)  : ',single_crystal_elastic_db(2,2)
      write(*,*)
   end if

   !Rw + Grt
   call read_par_double(15,Xol(3))           ; 
   call read_par_int(15,single_crystal_elastic_db(3,1))        ; 
   call read_par_int(15,single_crystal_elastic_db(3,2))        ; 
   if(rocktype(1) == 3) then
      write(*,"(a)")         ,' ---------------------------------'
      write(*,"(a)")         ,' LOWER TRANSITION ZONE: RW + GRT'
      write(*,"(a)")         ,' ---------------------------------'
      write(*,*)
      write(*,"(a,f8.2,a)")  ,' Xol(3):       ',Xol(3),' (%)'
      write(*,"(a,i5)")      ,' single_crystal_elastic_db(3,1)  : ',single_crystal_elastic_db(3,1)
      write(*,"(a,i5)")      ,' single_crystal_elastic_db(3,2)  : ',single_crystal_elastic_db(3,2)
      write(*,*)
   end if

   !Brd +  MgO
   call read_par_double(15,Xol(4))           ; 
   call read_par_double(15,stressexp(4))     ; 
   call read_par_double(15,Mob(4))           ; 
   call read_par_double(15,chi(4))           ; 
   call read_par_double(15,lambda(4))        ; 
   call read_par_double(15,fractdislrock(4)) ; 
   call read_par_double(15,tau(4,1))         ; 
   call read_par_double(15,tau(4,2))         ; 
   call read_par_double(15,tau(4,3))         ; 
   call read_par_double(15,tau(4,4))         ; 
   call read_par_double(15,tau(4,5))         ; 
   call read_par_double(15,tau(4,6))         ; 
   call read_par_double(15,tau(4,7))         ; 
   call read_par_double(15,tau(4,8))         ; 
   call read_par_double(15,tau(4,9))         ; 
   call read_par_double(15,tau(4,10))         ; 
   call read_par_double(15,tau(4,11))         ; 
   call read_par_double(15,tau(4,12))         ; 
   call read_par_int(15,single_crystal_elastic_db(4,1))        ; 
   call read_par_int(15,single_crystal_elastic_db(4,2))        ; 
   if(rocktype(1) == 4) then
      write(*,"(a)")         ,' -------------------------'
      write(*,"(a)")         ,' LOWER MANTLE: BRD + MGO'
      write(*,"(a)")         ,' -------------------------'
      write(*,*)
      write(*,"(a,f8.2,a)")  ,' Xol(4):       ',Xol(4),' (%)'
      write(*,"(a,f8.2)")    ,' stressexp(4): ',stressexp(4)
      write(*,"(a,f14.2)")   ,' Mob(4): ',Mob(4)
      write(*,"(a,f14.2)")   ,' chi(4): ',chi(4)
      write(*,"(a,f11.2)")   ,' lambda(4): ',lambda(4)
      write(*,"(a,f4.2)")    ,' fractdislrock(4): ',fractdislrock(4)
      write(*,"(a,es12.2,a)"),' tau(4,1): ',tau(4,1),' [100](010)  Bridgmanite '
      write(*,"(a,es12.2,a)"),' tau(4,2): ',tau(4,2),' [100](001)  Bridgmanite '
      write(*,"(a,es12.2,a)"),' tau(4,3): ',tau(4,3),' [010](100)  Bridgmanite '
      write(*,"(a,es12.2,a)"),' tau(4,4): ',tau(4,4),' [010](001)  Bridgmanite '
      write(*,"(a,es12.2,a)"),' tau(4,5): ',tau(4,5),' [001](100)  Bridgmanite '
      write(*,"(a,es12.2,a)"),' tau(4,6): ',tau(4,6),' [001](010)  Bridgmanite '
      write(*,"(a,es12.2,a)"),' tau(4,7): ',tau(4,7),' [001](110)  Bridgmanite '
      write(*,"(a,es12.2,a)"),' tau(4,8): ',tau(4,8),' [001](-110) Bridgmanite '
      write(*,"(a,es12.2,a)"),' tau(4,9): ',tau(4,9),' [110](001)  Bridgmanite '
      write(*,"(a,es11.2,a)"),' tau(4,10): ',tau(4,10),' [-110](001) Bridgmanite '
      write(*,"(a,es11.2,a)"),' tau(4,11): ',tau(4,11),' [110](-110) Bridgmanite '
      write(*,"(a,es11.2,a)"),' tau(4,12): ',tau(4,12),' [-110](110) Bridgmanite '
      write(*,"(a,i5)")      ,' single_crystal_elastic_db(4,1)  : ',single_crystal_elastic_db(4,1)
      write(*,"(a,i5)")      ,' single_crystal_elastic_db(4,2)  : ',single_crystal_elastic_db(4,2)
      write(*,*)
   end if

   !PPv +  MgO
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
   if(rocktype(1) == 5) then
      write(*,"(a)")         ,' -------------------------'
      write(*,"(a)")         ,' LOWER MANTLE: PPv + MGO'
      write(*,"(a)")         ,' -------------------------'
      write(*,*)
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
   end if

   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)

   !Set operating modes for the output       
   write(*,"(a)"),'SET OPERATING MODES FOR THE OUTPUT'       
   write(*,*)
   !call read_par_int(15,sbfmod)         ; write(*,"(a,i9)"),    ' sbfmod:      ',sbfmod
   !call read_par_double(15,rmax)        ; write(*,"(a,f9.2)")  ,' rmax:        ',rmax
   sbfmod = 0; rmax = 0.9;
   call read_par_int(15,ptmod)          ; write(*,"(a,i9)"),    ' ptmod:       ',ptmod
   call read_par_int(15,eosmod)         ; write(*,"(a,i9)"),    ' eosmod:      ',eosmod
   call read_par_double(15,mpgpa0)      ; write(*,"(a,f9.2)")  ,' pressure:    ',mpgpa0  
   call read_par_double(15,mtk0)        ; write(*,"(a,f9.2)")  ,' temperature: ',mtk0  
   call read_par_double(15,fractvoigt)  ; write(*,"(a,f9.2,a)"),' fractvoigt:  ',fractvoigt,' (%)'
   fractvoigt = fractvoigt/100d0
   IF(sbfmod > 0) fsemod = 1

   call read_par_int(15,spomod)         ; write(*,"(a,i9)"),    ' spomod:      ',spomod

   write(*,*)
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)

   CLOSE(15)

   END SUBROUTINE read_input_file

