!
!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a BETA version of the program LIJS, which is currently
!       under development by D. A. Spong at the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a BETA version, this program is subject to change
!       and improvement without notice.
!
c	Modified 7/20/2008 to remove unused variables
c---------------------------------------------------------------------------------
c	This module contains the kind specifications for PENTA and all subroutines 
c	 and modules.
c	
	module penta_kind_mod
	integer, parameter :: rknd=SELECTED_REAL_KIND(15,300)  !added by JL (5/14/08)
	integer, parameter :: iknd=SELECTED_INT_KIND(8)        !added by JL (5/14/08)
	end module penta_kind_mod
c
c-----------------------------------------------------------------
c	Physical constants used in several routines
c
	module phys_const
	use penta_kind_mod
	!particle constants
      real(rknd),parameter :: p_mass = 1.672621637e-27_rknd !proton mass
	real(rknd),parameter :: e_mass = 9.10938215e-31_rknd  !electron mass
      real(rknd),parameter :: qq = 1.602176487e-19_rknd	  !elementary charge
	real(rknd),parameter :: mpovme=1836.15267247_rknd	  !proton to electron mass ratio
	!other constants
	real(rknd), parameter :: ep0 = 8.854187817e-12_rknd

	end module phys_const
c
c-----------------------------------------------------------------
c	Variables used only in dkes_coef and the integrand function fun123
c
      module Vblck
	use penta_kind_mod
	real(rknd), dimension(:),allocatable ::  cmullg, efldlg,
     1			cmul, efield
	real(rknd) efldlg_min, efldlg_max, efld0
	integer(iknd) :: jpass
	logical :: log_coeff,use_log_interp
      end module Vblck
c
c-----------------------------------------------------------------
c	Variable passed between routines xlim and fxlim
c
      module fxvars
	use penta_kind_mod
      real(rknd) cmlx
      end module fxvars
c
c-----------------------------------------------------------------
c	These variables are used in penta, sugama_app_c, dkes_coef and colint
c
      module lijs_transfer
	use penta_kind_mod
	real(rknd) :: ma,qa,na,nb,qb,vta,vtb,lnlambda
      end module lijs_transfer
c
c-----------------------------------------------------------------
c	These variables are only used in sugama_app_c and dkes_coef
c
	module dkes_coef_transfer
	use penta_kind_mod
	real(rknd) :: abs_Er_ovth
	integer(iknd) :: lcoef_opt
	end module dkes_coef_transfer
c
c-----------------------------------------------------------------
c	These variables are only used in penta.f and the sugama_app_c subroutine
c
      module l_components
	use penta_kind_mod
	logical :: beam,lijs_write,dkes_limit
      real(rknd) :: ion_mass, temp_i_keV, temp_e_keV,
     >   factor_m_e, factor_m_i, factor_n_e, zee,
     >  factor_n_i, inv_factor_m_e, inv_factor_m_i,
     >  factor_Ee1, factor_Ee2, den_e_mks, den_i_mks
      real(rknd), dimension(2,2) :: lam_e, lam_i
      end module L_components
c	
c-----------------------------------------------------------------
c	Variables used in calling 2-D bspline routines
c
      module V2dspline_bspline 
	use penta_kind_mod
      integer(iknd) :: nc, ne, kxord, kyord
      real(rknd), dimension(:), allocatable :: xt_e, xt_c
      real(rknd), dimension(:,:), allocatable :: c_spl
      end module V2dspline_bspline
c
c-----------------------------------------------------------------
c	These variables are used in penta.f and the rad_flux subroutine
c
      module er_solve
	use penta_kind_mod
      integer(iknd) m_er, ne_pts
	real(rknd), dimension(:), allocatable :: x_er, c_er
      end module er_solve
c
c-----------------------------------------------------------------
c	Viscosity variables used in penta.f and the sugama_app_c subroutine
c
      module viscosity_pol_tor
	use penta_kind_mod
       real(rknd), dimension(2,2) :: mnl1i, mnl2i, mnl3i,
     >     mnl1e, mnl2e, mnl3e
      logical :: pol_tor_cmpt
      end module viscosity_pol_tor