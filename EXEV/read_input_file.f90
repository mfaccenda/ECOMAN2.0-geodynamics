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

   !Viscosity contrast and inclusions volume fraction
   write(*,"(a)"),' VISCOSITY CONTRAST AND INCLUSIONS VOLUME FRACTIONS'
   write(*,*)
   call read_par_double(15,etac1)   ; write(*,"(a,f8.4)"),' etac :   ', etac1
   call read_par_double(15,volinc1) ; write(*,"(a,f6.2)"),' volinc : ', volinc1
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
      write(*,*)
      write(*,"(a)"),'--------------------------------------------------------'
   end if

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

   call read_par_int(15,radialmod)

   call read_par_int(15,azimod)
   call read_par_double(15,aziscalex1)
   call read_par_double(15,aziscalex2)
   call read_par_double(15,aziscalex3)

   if(Eulerian /= 0) then 
      write(*,*)
      write(*,"(a)"),' EULERIAN GRIDDING AND VISUALIZATION IS ACTIVE'
      write(*,*)
      write(*,"(a,i1)"),' radialmod :       ',radialmod
      write(*,"(a,i1)"),' azimod :          ',azimod
      write(*,*)
      write(*,"(a)"),'--------------------------------------------------------'
   end if

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE read_par_double6(unitt,dbl)

   IMPLICIT NONE

   INTEGER ::i,ier,unitt
   DOUBLE PRECISION :: dbl(6)
   CHARACTER(500) :: str

10 READ(unitt,"(a)",advance='yes',IOSTAT=ier) str  
   i=index(str,'#')
   if(i==1 .OR. LEN_TRIM(str)==0) then
      goto 10
   else
      read(str,*) dbl(:)
   end if

   END SUBROUTINE read_par_double6
