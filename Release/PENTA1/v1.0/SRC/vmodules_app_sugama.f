c
c-----------------------------------------------------------------
c
      module Vblck
      integer   itkf, ice, icmul, jefld, ncpts, kxord, kyord,
     1   ncpts4, klam1, icon1, m1_spl, istat1, isw1, i_loc
      real*8, dimension(60) :: cmul, efield
      real*8, dimension(3600) :: gamcd
      real*8, dimension(60) :: cmullg, efldlg, efldlg_nrm
      real*8, dimension(60,60) :: gamln2d
      real*8, dimension(3600) :: gamln
      real*8, dimension(60,20) :: rlij
      real*8, dimension(60) :: rnuowt
      real*8, dimension(20) :: rlamvl
      real*8, dimension(2000) :: wrk
      real*8, dimension(50) :: xcol, colspl
      real*8, dimension(54) :: rknot, cofspl
      real*8, dimension(:), allocatable :: xt1, vw1, c1_spl
      real*8 wefac,cmul0,cmin,cmax,efld0,emin,emax,one,
     1   rootm, zeff, efldlg_min, efldlg_max, pfac, col0
      end module Vblck
c
c-----------------------------------------------------------------
c
      module Vcfxs
      real*8 cmlxmn, cmlxmx
      end module Vcfxs
c
c-----------------------------------------------------------------
c
      module V2dspline
c
c-----------------------------------------------------------------
c
      integer :: px, py, ic4, ie4, icie, iwork_dim,
     1   ierr, istat
      integer, dimension(:), allocatable :: iwrksp
      real*8, dimension(:), allocatable :: lambdasp,
     1  mu, cc, wrksp
      end module V2dspline
c
c-----------------------------------------------------------------
c
      
      module lijs_transfer
      real*8 temp, mass, zee, trat, masrat, arad
      real*8 nulo,nuhi,c11,c12,c22,c31,c32,vs1,vs2,vs3
      integer lcoef_opt, klm, npts
      logical beam
      logical lijs_write
      end module lijs_transfer
c
c-----------------------------------------------------------------
c
      module l_components
      real*8 :: mass_i, mass_e, trat_i, trat_e,
     >  masrat_i, masrat_e, nu_i, nu_e, temp_i, temp_e,
     >   factor_m_e, factor_m_i, factor_n_e, nu_i0, nu_e0,
     >  factor_n_i, inv_factor_m_e, inv_factor_m_i,
     >  factor_Ee1, factor_Ee2, den_e, den_i, qq, den_e_mks, den_i_mks
      real*8, dimension(2,2) :: ll_e, ll_i, lam_e, lam_i
      end module L_components
c
c
c-----------------------------------------------------------------
c
c      module V1dspline_ssl2
c      integer :: mspl_1d, icon, nc, iswc,
c     >     ic, idim_xt, idim_vw, nmin, nmax, nn
c      real*8, dimension(:), allocatable :: xt_1d, vw_1d
c      real*8, dimension(:), allocatable :: c_spl_1d
c      real*8 :: s, err
c      end module V1dspline_ssl2
c
c-----------------------------------------------------------------
c
c      module V2dspline_ssl2
c      integer :: mspl, icon, nc, ne, iswc,
c     >     iswe, ic, ie, idim_xt, idim_vw, nmin, nmax, nn
c      real*8, dimension(:), allocatable :: xt, vw
c      real*8, dimension(:,:), allocatable :: c_spl, zint
c      real*8 :: s, err
c      end module V2dspline_ssl2
c
c
c
c-----------------------------------------------------------------
c
      module V2dspline_bspline
      integer :: nc, ne, ic, ie, nn, icon, nmin, nmax
      real*8, dimension(:), allocatable :: xt_e, xt_c
      real*8, dimension(:,:), allocatable :: c_spl, zint
      real*8 :: s, err
      end module V2dspline_bspline
c
c-----------------------------------------------------------------
c
      module er_solve
        integer m_er, icon_er, ne_pts, isw_er, i_loc_er
	real*8, dimension(:), allocatable :: c_er,
     >    er, net_flux
	real*8, dimension(:), allocatable :: x_er
c	real*8, dimension(610) :: v_er            !ssl2 option
	real*8 :: ermin, ermax, tol_er
c
	real*8, dimension(:), allocatable :: c_qi, x_qi,
     >     c_qe, x_qe
c	real*8, dimension(610) :: v_qi, v_qe      !ssl2 option
      end module er_solve
c
c
      module viscosity_pol_tor
       real*8, dimension(2,2) :: mnl1i, mnl2i, mnl3i,
     >     mnl1e, mnl2e, mnl3e
       logical pol_tor_cmpt
      end module viscosity_pol_tor
c   The following modules are only needed if the NAG options are used:
c
c      module VblckEr0
c      integer   itkf, ice, icmul, npts, ncpts,
c     1   ncpts4, lcoef_opt
c      real*8, dimension(50) :: cmul
c      real*8, dimension(2500) :: gamcd
c      real*8, dimension(50) :: cmullg
c      real*8, dimension(2500) :: gamln
c      real*8, dimension(50) :: rlij
c      real*8, dimension(50) :: rnuowt
c      real*8, dimension(2000) :: wrk
c      real*8, dimension(50) :: xcol, colspl
c      real*8, dimension(54) :: rknot, cofspl
c      real*8 cmul0,cmin,cmax,
c     1   rootm,zeff,c11,c12,c22,c31,c32
c      end module VblckEr0
c
c-----------------------------------------------------------------
c
c      module VcfxsEr0
c      real*8 cmlxmn, cmlxmx
c      end module VcfxsEr0
c
c-----------------------------------------------------------------
c
c      module VsplineEr0
c      integer icmul4,lw_cmul
c      real*8, dimension(54) :: cknots, cmspline
c      real*8, dimension(2000) :: wrk_cmul
c
c-----------------------------------------------------------------
c
c      end module VsplineEr0
