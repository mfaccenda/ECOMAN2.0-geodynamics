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

   SUBROUTINE read_input_file(cijkl_dir,output_dir)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i,dum(3)
   CHARACTER (len=500) :: cijkl_dir,output_dir,command
   CHARACTER(200) :: arg,inputfilename,str

   DO i = 1, iargc()
      CALL getarg(i, arg)
      if(i==1) inputfilename = arg
   END DO
 
   OPEN(unit=15,file=inputfilename,status='old',action='read')

   !Input files directory
   write(*,*)
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)
   write(*,"(a,a)"),' READING INPUT FILE: ',trim(inputfilename)
   write(*,*)
   call read_par_char(15,cijkl_dir) ; write(*,"(a,a)"),' Cijkl  directory: ',trim(cijkl_dir)
   write(*,*)
   call read_par_char(15,output_dir); write(*,"(a,a)"),' Output directory: ',trim(output_dir)
   write(*,*)

   command = 'mkdir '//trim(output_dir)
   CALL SYSTEM(command)
   write(*,*)
   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)

   !Input files min,step,max, output frequency, max. time (in Myr) in
   !steady-state conditions
   write(*,"(a)"),' INPUT - OUTPUT FILES'
   write(*,*)
   call read_par_int(15,Tinit) ; write(*,"(a,i5)"),' Tinit :       ',Tinit
   call read_par_int(15,Tstep) ; write(*,"(a,i5)"),' Tstep :       ',Tstep
   call read_par_int(15,Tend)  ; write(*,"(a,i5)"),' Tend  :       ',Tend
   write(*,*)
   write(*,"(a)"),'--------------------------------------------------------'
   write(*,*)

   !Lagrangian fields
   call read_par_int(15,Lagrangian)
   call read_par_double(15,ln_fse_min)  
   call read_par_int(15,uppermantlemod)
   call read_par_int(15,rocktypemod)
   call read_par_int(15,fse3Dmod)
   call read_par_int(15,fseminmod)
   call read_par_int(15,fsemaxmod)
   call read_par_int(15,TIaxismod)
   call read_par_int(15,vpmaxmod)
   call read_par_int(15,dvsmaxmod)
   if(Lagrangian /= 0) then
      write(*,*)
      write(*,"(a)"),' LAGRANGIAN AGGREGATES VISUALIZATION IS ACTIVE'
      write(*,*)
      write(*,"(a,f6.2)"),' ln_fse_min : ',ln_fse_min
      write(*,"(a,i1)")  ,' uppermantlemod :  ',uppermantlemod
      write(*,"(a,i1)")  ,' rocktypemod :     ',rocktypemod
      write(*,"(a,i1)")  ,' fse3Dmod :        ',fse3Dmod
      write(*,"(a,i1)")  ,' fseminmod :       ',fseminmod
      write(*,"(a,i1)")  ,' fsemaxmod :       ',fsemaxmod
      write(*,"(a,i1)")  ,' TIaxismod :       ',TIaxismod
      write(*,"(a,i1)")  ,' vpmaxmod :        ',vpmaxmod
      write(*,"(a,i1)")  ,' dvsmaxmod :       ',dvsmaxmod
      write(*,*)
      write(*,"(a)"),'--------------------------------------------------------'
   end if

   !SPO modelling
   call read_par_int(15,spomod)
   
   !Eulerian fields
   call read_par_int(15,Eulerian)
   call read_par_double(15,n1first)  
   call read_par_double(15,n1last)  
   call read_par_int(15,nx11)
   call read_par_double(15,n2first)  
   call read_par_double(15,n2last)  
   call read_par_int(15,nx21)
   call read_par_double(15,n3first)  
   call read_par_double(15,n3last)  
   call read_par_int(15,nx31)

   call read_par_int(15,vpvsmod)
   call read_par_int(15,dvpvsmod)
   if(vpvsmod == 0 .and. dvpvsmod > 0) write(*,"(a)"),'WARNING: dvpvsmod is > 0, but vpvsmod is = 0, so no dVp and dVs will be printed' 
   call read_par_int(15,zoeppritzmod)

   call read_par_double(15,cosx1)
   call read_par_double(15,cosx2)
   call read_par_double(15,cosx3)

   call read_par_int(15,nx1ref)
   call read_par_int(15,nx3ref)

   call read_par_int(15,radialmod)

   call read_par_int(15,azimod)
   call read_par_double(15,aziscalex1)
   call read_par_double(15,aziscalex2)
   call read_par_double(15,aziscalex3)

   call read_par_int(15,reflectmod)
   call read_par_int(15,replicateZmod)
   call read_par_int(15,specfem3Dmod)
   call read_par_int(15,psitomomod)
   if(Eulerian /= 0) then 
      write(*,*)
      write(*,"(a)"),' EULERIAN GRIDDING AND VISUALIZATION IS ACTIVE'
      write(*,*)
      write(*,"(a,i1)"),' vpvsmod :         ',vpvsmod
      if(vpvsmod == 2) write(*,"(a,f6.2,a,f6.2,a,f6.2,a)"),' cosx1 = ',cosx1,', cosx2 = ',cosx2,', cosx3 = ',cosx3,' (deg)'
      write(*,"(a,i1)"),' dvpvsmod :        ',dvpvsmod
      if(dvpvsmod == 2)      write(*,"(a,i5)"),' nx1ref :      ',nx1ref
      if(dvpvsmod == 2)      write(*,"(a,i5)"),' nx3ref :      ',nx3ref
      write(*,"(a,i1)"),' zoepptrizmod :    ',zoeppritzmod
      write(*,"(a,i1)"),' radialmod :       ',radialmod
      write(*,"(a,i1)"),' azimod :          ',azimod
      write(*,"(a,i2)"),' reflectmod :     ',reflectmod
      write(*,"(a,i1)"),' replicateZmod :   ',replicateZmod
      write(*,"(a,i1)"),' specfem3Dmod :    ',specfem3Dmod
      write(*,"(a,i1)"),' psitomomod :      ',psitomomod
      write(*,*)
      write(*,"(a)"),'--------------------------------------------------------'
   end if

   write(*,*)
   write(*,"(a)"),'********************************************************'
   write(*,"(a)"),'********************************************************'
   write(*,*)

   CLOSE(15)

   END SUBROUTINE read_input_file
