!-----------------------------------------------------------------------------
!
!  -THIS CODE GENERATES A PLASMA PROFILE FILE FOR PENTA3.-
!
!
!-----------------------------------------------------------------------------
!   plasma_profiles_XXX.dat - A file containing the plasma profile
!       data, where XXX is arg7 above.  This file has the format:
!       row 1: number of radial points for which data is specified.
!       all other rows: r/a, ne, Te, Ti1, ni1, Ti2, ni2 ...
!       where the densities are in units of 10**12/cc and the 
!       temperatures in eV and i1, i2... are the ion species.
!   
!-----------------------------------------------------------------------------
!   NOTES
!    - All units are SI, except T is [eV] (unless otherwise noted)
!
!  11/29/2010 JL
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!
!+ The main program
!
Program gen_pprof_file
 
! Description: 
!  QQ
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 1.0     11/29/2010   Original code. JL 
! 
! Author(s): J. Lore
!
! Code Description:
!   Language:   Fortran 90
! 

Use penta_kind_mod
Use penta_math_routines_mod, Only : rlinspace

Implicit None

! Local variables (scalar)
Integer(iknd) :: num_rad, num_ion, ir, is
Integer(iknd) :: iu_out = 10

! Local arrays (allocatable)
Real(rknd), Allocatable :: ne(:), Te(:), ni(:,:), Ti(:,:), roa(:)
Character(Len = 20) :: fstatus,fpos,str_num

!- End of header -------------------------------------------------------------

Write(*,*) "Starting gen_pprof_file"



! Set number of radial locations at which data are specified
num_rad = 10

! Number of ion species  -- for now this only works for one
num_ion = 1

! Allocate the data arrays
Allocate(roa(num_rad),ne(num_rad),Te(num_rad))
Allocate(ni(num_rad,num_ion),Ti(num_rad,num_ion))

! Specify the radial coordinate data
roa = rlinspace(0._rknd,1._rknd,num_rad)

! Specify the data: ne (cm**-3), Te (eV), ni (cm**-3), Ti (eV)

ne = 1.d12
ni(:,1) = ne
Ti(:,1) = 11.d0
Te = Ti(:,1)


! Set specifiers for opening output files
fstatus = "unknown"
fpos = "asis"

! Open output files
Open(unit=iu_out, file="plasma_profiles_new.dat", &
  position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))



Write(str_num,*) num_rad
Write(iu_out,*) Trim(adjustl(str_num))

is = 1
Write(str_num,*) num_ion*2 + 2
Do ir = 1,num_rad

  Write(iu_out,'(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.7))') &
       roa(ir),ne(ir)/1.d12,Te(ir),ni(ir,is)/1.d12,Ti(ir,is)

EndDo
! Close output file
Close(iu_out)

! Deallocate variables
Deallocate(roa,ne,Te,ni,Ti)
Write(*,*) 'Finished'

End Program gen_pprof_file


