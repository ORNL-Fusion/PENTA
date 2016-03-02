
!-----------------------------------------------------------------------------
!+ Module for variables read from the DKES coeff. files
!-----------------------------------------------------------------------------
Module coeff_var_pass

!
! Description:
!   This module contains the DKES coeff. data variables.
!   For definitions see subroutine read_dkes_star_files 
!    in read_input_file_mod.f90.
!   For directions on how to generate the coeff. file see penta.f90.
!
! History:
! Version     Date      Comment
! -------   --------    -------
!  1.0     01/07/2009   Original Code.  JL
!  1.1     05/24/2010   Updated for PENTA3. JL 
! 
! Author(s): J. Lore 7/2009 - 8/31/2010 
!     
  
  ! Modules used:
  Use penta_kind_mod         ! Import rknd, iknd specifications

  Implicit none

  ! Array variables
  Real(rknd), allocatable :: cmul_D11(:),     & ! cmul arrays for coeffs 
    cmul_D13(:), cmul_D31(:),                 &
    cmul_D33(:)
  Real(rknd), allocatable :: efield_D11(:),   & ! efield arrays for coeffs
    efield_D13(:), efield_D31(:),             &
    efield_D33(:)
  Real(rknd), allocatable :: D11_mat(:,:),    & ! 2D coefficient arrays
    D13_mat(:,:), D31_mat(:,:), D33_mat(:,:) 

  ! Scalar variables
  Integer(iknd) :: num_c_D11, num_e_D11,    & ! Number of efield and cmul vals
    num_c_D13, num_e_D13,                   &
    num_c_D31, num_e_D31,                   &
    num_c_D33, num_e_D33

End module coeff_var_pass
!- End of header -------------------------------------------------------------
