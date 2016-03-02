!-----------------------------------------------------------------------------
!+ Sets file I/O unit numbers
!-----------------------------------------------------------------------------
Module io_unit_spec

!
! Description:
!   This module sets the file i/o numbers to avoid conflicts.
!
! History:
! Version   Date      Comment
! -------   ----      -------
!  1.0     01/7/2009   Original Code.  JL
!  1.1     5/24/2010   Updated for PENTA3. JL 
! 
! Author(s): J. Lore 7/2009 - 9/29/2010 
!

  Implicit none
  
  Integer, parameter :: &
    iu_nl=21,           &   ! Ion parameter namelist file (input)
    iu_vmec=22,         &   ! VMEC data file (input)
    iu_pprof=23,        &   ! Plasma profile file (input)
    iu_coeff=24,        &   ! DKES coefficient file (input)
    iu_Ufile=25,        &   ! <U**2> file (input)
    iu_flux_out=10,     &   ! Fluxes vs r/a (output)
    iu_pprof_out=11,    &   ! Plasma profile check (output)
    iu_fvEr_out=12,     &   ! Fluxes vs Er (output)
    iu_flows_out=13         ! Parallel flows vs r/a (output)
End module io_unit_spec
!- End of header -------------------------------------------------------------
