c	-PENTA-
c
c  To compile, use: 
c     build_app_c_sugama_linux script for Linux systems with the Lahey compiler
c     build_penta_macintel For Mac Intel systems with the ifort compiler
c     build_penta_ibm for IBM's with the xlf compiler
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  To run,type: xpnt [command line arguments]
c
c  where the command line arguments:
c
c  1  unique identifier text on DKES database files
c      i.e., for (m,l,n)star_lijs_***_s##, this would be
c      the text that replaces ***
c  2  Er,min
c  3  Er,max  - args 2 and 3 specify the interval in Er over which
c      the search is done for the ambipolar electric field root.
c      Er,min and Er,max are dimensionless normalized electric
c      field variables, i.e.: e<a>Er,min/max/kT_e
c	 - Note if the logical variable input_is_Er below is true then 
c	   the search interval is set by args 2 and 3 as Er (V/cm)
c  4  flux surface number
c  5  parameter that determines whether output data should be
c      written into new files (=0) or appended to
c      existing files (=1). Typically,when a script is run for a
c      sequence of flux surfaces, this is set to 0 for the first
c      surface and to 1 for subsequent surfaces
c  6  extension of the profile_data_*** file (i.e., the *** text)
c      that comes fromextracting the appropriate profiles from
c      the VMEC run using the profile_extractor.f code
c  7  central ion temperature
c  8  central density
c  9  edge density pedestal
c 10  edge temperature pedestal
c 11  central electron temperature
c 12  unique identifier on plasma profiles file, i.e., the ***
c      text in the filename plasma_profiles***.dat. This allows
c      multiple profiles to be run without having to rename files
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	NOTES
c	Equation references are denoted by parentheses and are refer to
c		Sugama, Nishimura PoP 9 4637 (2002) except where noted.
c	
c	Added Don's latest changes as of 7_7_08
c
c
	use penta_kind_mod
	use phys_const
      use lijs_transfer
      use l_components
      use er_solve
      use bspline
      use viscosity_pol_tor
      implicit none

      logical, parameter :: plots=.false.             !true requires plplot libraries
      logical, parameter :: er_fixed=.false.          !option for fixed electric field profile
      logical, parameter :: beam_on=.false.           !option for external parallel momentum source
      logical, parameter :: exact_sym=.false.         !if true, turns off Er terms in Xi1, Xe1
      logical, parameter :: search_down=.true.       !direction of Er root search (if false, then upward search is done)
      logical, parameter :: profile_scaling = .true. !Use pedestal and axis profile scaling (JL 5/22/08)
      logical, parameter :: dkes_limit_on = .false.	!Should revert to dkes limit, but seperate L* file required - where L* is D11*
      logical, parameter :: input_is_Er = .false.	    !If true, command line inputs are Er (V/cm) otherwise they are e<a>Er/kT_e
      logical, parameter :: idebug = .false.	    !More text output

      real(rknd), parameter :: k_boltz = 1.602176487e-16_rknd	!Boltzmann factor in Joules/keV
	
      real(rknd) :: ermin, ermax, tol_er
      real(rknd), dimension(:), allocatable :: c_qi, c_qe
      real(rknd), dimension(:), allocatable :: er, net_flux  
            
	real(rknd) :: mu, rmajor, abs_Er, Er_mks, temp_e_eV
      real(rknd) :: Ti_0, Te_0, r, ra, den0,temp_i_eV
      real(rknd) :: vth_e, vth_i, er_root
      real(rknd) :: tau_ee, tau_ii, arad
      real(rknd) :: dTe_dr_oTe, dTi_dr_oTi, dn_dr_on
      real(rknd) :: ion_flux,elec_flux, den_e, den_i
      real(rknd) :: xe1, xe2, xi1, xi2, xE, ep02, pi, forpi2,xi1e

      real(rknd) :: e_phi_kt_root, uprl_B, bzeta, btheta, uprl0
      real(rknd) :: chip, psip, vp, bsq, upol, utor, upol_Er_only,
     >          uprl_B_Er_only,nu_i0, nu_e0
      real(rknd) :: jbs_e, jbs_i, jbs_e_prl, jbs_tot, fit_scl,
     >    qflux_e, qflux_i, gam_e, gam_i, gi_m_ge, ped_d, ped_t,
     >    jbs_e0, jbs_i0, jbs_e_prl0, jbs_tot0, beam_source,
     >    fit0,fit1,fit2,fit3,fit4,fit5,fit6, temp_axis_elec,
     >    temp_axis_ion, den_axis, iota, chk
      real(rknd), dimension(5,5) :: lmat
      real(rknd), dimension(2,2) :: fl_mat, proj_mat, ll_e, ll_i
      real(rknd), dimension(2,2) :: mpt1i, mpt2i, mpt3i,
     >  mpt1e, mpt2e, mpt3e, inter_mat, pre_mat, post_mat
      real(rknd) :: ymin, ymax, xmin, xmax
      real(rknd), dimension(:), allocatable :: e_phi_kT,
     1   gamma_e, gamma_i, q_e, q_i, jbs,gmi_e1,
     2   gmi_e2, gmi_i1, gmi_i2, gmi_i3, gmi_Eprl
	real(rknd), dimension(:), allocatable :: gme_e1,
     2   gme_e2, gme_e3, gme_i1, gme_i2, gme_Eprl
      real(rknd), dimension(:), allocatable :: r_m, r_norm,
     >  chip_profile, psip_profile, btheta_profile,
     >  bzeta_profile, vp_profile, bsq_profile,iota_profile,
     >  r_prof, temp_prof_elec, den_prof, dn_knot_array,
     >  te_knot_array, spl_dn, spl_te,temp_prof_ion,
     >  ti_knot_array, spl_ti
      integer(iknd) :: i, ii, iroot, i_append, istat, js,
     1   js_min, js_max, j, index, np_prof, kord_prof, idrv,
	2	numargs
      character*10 :: ch_dum
      character*50 :: lijs_inp_ext, plasma_prof_file
      CHARACTER*60 :: arg1, arg2, arg3, arg4, arg5, arg6, arg7,
     1   arg8, arg9, arg10,arg11, arg12, arg13, arg14,
     2   profile_data
      CHARACTER*1 :: tb
      real(rknd), external :: rad_flux, zeroin
	integer(iknd), external :: iargc
c
c	Define constants
c
      tb = char(9)
      ep02 = ep0**2_iknd
      pi = 4._rknd*atan(1._rknd)
	forpi2 = 4._rknd*pi*pi

	!Set logical variables
      beam = beam_on
	dkes_limit = dkes_limit_on

	!Ion parameters:
      ion_mass = p_mass	   !ion mass
	mu = 1._rknd           !ion to proton mass ratio
	zee = 1._rknd          !ion charge number
c      
c     Optional fixed electric field fit:
c
      if(er_fixed) then
        fit_scl = 1.e-4_rknd
        fit0 = 0._rknd*fit_scl
        fit1 = -1.92_rknd*fit_scl
        fit2 = -115.4_rknd*fit_scl
        fit3 = 117.3_rknd*fit_scl
        fit4 = 0._rknd
        fit5 = 0._rknd
        fit6 = 0._rknd
      end if
c      er_root = fit0 + fit1*r + fit2*(r**2) + fit3*(r**3) + fit4*(r**4) + fit5*(r**5) + fit6*(r**6)
c
      Te_0 = 0._rknd            !intial value, indicating one is not provided
c
c	Get command line arguments
c
      CALL getarg(1, arg1)
      CALL getarg(2, arg2)
      CALL getarg(3, arg3)
      CALL getarg(4, arg4)
      CALL getarg(5, arg5)
      CALL getarg(6, arg6)
      CALL getarg(7, arg7)
      CALL getarg(8, arg8)
      CALL getarg(9, arg9)
      CALL getarg(10, arg10)
      CALL getarg(11, arg11)
      CALL getarg(12, arg12)
	
	numargs=iargc()
	
	if ( numargs .gt. 12 ) then
	  CALL getarg(13, arg13)
	end if
c
c	!store command line args
c
      lijs_inp_ext = trim(adjustl(arg1))
      read(arg2,'(f10.5)') ermin
      read(arg3,'(f10.5)') ermax
      read(arg4,'(i8)') js
      read(arg5,'(i1)') i_append
      read(arg7,'(f10.5)') Ti_0
      read(arg8,'(e15.7)') den0
      read(arg9,'(f10.5)') ped_d
      read(arg10,'(f10.5)') ped_t
      read(arg11,'(f10.5)') Te_0
c
c	Read beam file (if beam_on)
c
	if(beam_on) then
        open(unit=33,file="beam_params",status="old")
        read(33,*) beam_source
	  close(33)
      end if
c
c	Display logical variable data to screen (if i_append==0)
c
	if(i_append .eq. 0) then
	  write(*,*) 
	  write(*,*) "Welcome to PENTA, please note the following settings:"
	  write(*,*)
	  if ( input_is_Er == .true. ) then
	    write(*,*) 'Interpreting input range as Er (V/cm)'
	  else
	    write(*,*) 'Interpreting input range as e<a>Er/kTe'
	  end if
        if ( er_fixed == .true. ) write(*,*) 'Using fixed Er profile.'
        if ( beam_on == .true. ) write(*,*) 'Using beam data.'
        if ( exact_sym == .true. ) write(*,*) 'Assuming exact symmetry.'
        if ( search_down == .true. ) then
	    write(*,*) 'Search direction for Er root finding is down.'
	  else 
	    write(*,*) 'Search direction for Er root finding is up.'
	  end if
	  if ( profile_scaling == .true. ) then
	    write(*,*) 'Using plasma profile pedestal parameters.'
	  else
	    write(*,*) 'Not using plasma profile pedestal parameters.'
	  end if
	  if ( dkes_limit_on == .true. ) write(*,*) 'Calculating
     > quantities in dkes limit, L* file should contain D11*.'
	  if ( plots == .true. ) then
	    write(*,*) 'Writing extra information to screen and file'
	  end if
	  write(*,*)
	  write(*,*) " <r>/<a>","   Er(V/cm)","       <a>Er/kTe"
     >   ,"       upol_sug","         utor_sug"
	end if
c
c	open output files and (if i_append==0) write legends 
c
      if(i_append .eq. 0) then
        open(unit=20,file="run_sug",status="unknown")
        open(unit=23,file="fluxes_vs_r",status="unknown")
        open(unit=11,file="flux_avg_data0",status="unknown")
        open(unit=19,file="flux_avg_data1",status="unknown")
        open(unit=18,file="jbs_vs_r",status="unknown")
        open(unit=15,file="flow_viz1.dat",status="unknown")	
        write(11,'("*",/,"<r>/<a>",a1,"Er(V/cm)",a1,"<a>Er/kTe",  
     >   a1,"Upol(contra)",a1,"Utor(contra)",a1,"Uprl",a1,
     >   "Upol(contra)_Er_only")')
     >    tb, tb, tb, tb, tb, tb
	!added label to flux_avg_data1 (5/08 JL)
	  write(19,'("*",/,"<r>/<a>",a1,"xi1",a1,"uprl_B",
     >   a1,"bsq",a1,"arad",a1,"vp",a1,"btheta",a1
     >   ,"bzeta",a1,"psip",a1,"chip",a1,"xi1e")')
     >    tb, tb, tb, tb, tb, tb, tb, tb, tb
        write(18,'("*",/,"<r>/<a>",a1,"ne",a1,"Te",
     >   a1,"T_i",a1,"Jbs_e",a1,"Jbs_i",a1,"Jbs_E||",a1,"Jbs_total",
     >   a1,"Jbs_e(Er=0)",a1,"Jbs_i(Er=0)",a1,"Jbs_E||(Er=0)",a1,
     >   "Jbs_total(Er=0)",a1,"net_flux",a1,"qflux_e"
     >   ,a1,"qflux_i",a1,"L_EE")')
     >    tb, tb, tb, tb, tb, tb, tb, tb, tb, tb, tb, tb, tb, tb, tb
         !       Fixed label error (8/08) JL
        write(23,'("*",/,"<r>/<a>",a1,"Gamma_i",a1,"Gamma_e",a1	
     >     "Q_i",a1,"Q_e",a1,"Gmi_e1",a1,"Gmi_e2",a1,"Gmi_i1",
     >     a1,"Gmi_i2",a1,"Gmi_Eprl")') tb,tb,tb,tb,tb,tb,tb,tb,tb
      else if(i_append .eq. 1) then
         open(unit=20,file="run_sug",position="append",status="old")
         open(unit=23,file="fluxes_vs_r",position="append",
     >        status="old")
         open(unit=11,file="flux_avg_data0",position="append",
     >        status="old")
         open(unit=19,file="flux_avg_data1",position="append",
     >        status="old")
         open(unit=18,file="jbs_vs_r",position="append",
     >        status="old")
         open(unit=15,file="flow_viz1.dat",position="append",
     >        status="old")
      endif
      open(unit=12,file="cmul_E_requests",status="unknown")
c
c	run_sug file has a legend for each surface
c
      write(20,'("*",/,"<r>/<a>",a1,"E_r (V/cm)",a1,"<a>ePhi/kTe",a1,
     >     "gamma_i",a1,"gamma_e",a1,"Q_i",a1,"Q_e")')
     >     tb, tb, tb, tb, tb, tb
c
c    Read VMEC profile data file and select data for current flux position
c
      profile_data = "profile_data_" // trim(adjustl(arg6))
      open(unit=25,file=profile_data,status="old")
	read(25,*) js_min, js_max
	allocate(r_m(js_max), stat=istat)
      allocate(r_norm(js_max), stat=istat)
      allocate(chip_profile(js_max), stat=istat)
      allocate(psip_profile(js_max), stat=istat)
      allocate(btheta_profile(js_max), stat=istat)
      allocate(bzeta_profile(js_max), stat=istat)
      allocate(vp_profile(js_max), stat=istat)
      allocate(bsq_profile(js_max), stat=istat)
      allocate(iota_profile(js_max), stat=istat)
      read(25,*) arad !, rmajor
	rmajor = 1.0_rknd
      read(25,'(a10)') ch_dum
	!write(*,*) 'change back to i4!!!' idebug
	do j = js_min,js_max
c       read(25,'(i4,9(a1,e15.7))') index,tb,r_m(j),tb,r_norm(j),tb,
c     >    chip_profile(j),tb,psip_profile(j),tb,btheta_profile(j),
c     >    tb,bzeta_profile(j),tb,vp_profile(j),
c     >    tb,bsq_profile(j),tb,iota_profile(j)
       read(25,*) index,r_m(j),r_norm(j),
     >    chip_profile(j),psip_profile(j),btheta_profile(j),
     >    bzeta_profile(j),vp_profile(j),
     >    bsq_profile(j),iota_profile(j)
      end do
      chip = chip_profile(js); psip = psip_profile(js)
      bsq = bsq_profile(js); vp = vp_profile(js)
      btheta = btheta_profile(js); bzeta = bzeta_profile(js)
      iota = iota_profile(js)
	!NOTE that the variable 'r'=r/a=rho=sqrt(normalized toroidal flux) and 'ra'=r=rho*a
	ra = r_m(js); r = r_norm(js);	
c
c	Read plasma profiles from file   
c	
	!open plasma_profiles_check file and write legend
      if(i_append .eq. 0) then  
        open(unit=34,file="plasma_profiles_check",status="unknown")
	  write(34,116) tb, tb, tb, tb, tb, tb
 116    format("<r>/<a>",a1,"den",a1,"temp_e_keV",a1,"temp_i_keV",a1,		!fixed legend 4/2009
     >        "dndr_ovr_n",a1,"dTe_dr_oTe",a1,"dTi_dr_oTi")
      else if(i_append .eq. 1) then
        open(unit=34,file="plasma_profiles_check",
     >          position="append",status="old")
      end if
	!read r/a, n,Te,Ti from "plasma_profiles.dat" file
	plasma_prof_file="plasma_profiles"//trim(adjustl(arg12))//".dat"
      open(unit=35,file=plasma_prof_file,status="old")
      read(35,*) np_prof
      allocate(r_prof(np_prof), stat=istat)
      allocate(den_prof(np_prof), stat=istat)
      allocate(temp_prof_elec(np_prof), stat=istat)
      allocate(temp_prof_ion(np_prof), stat=istat)
	!again, here r_prof is r/a
      do j=1,np_prof
	  read(35,*) r_prof(j),den_prof(j),temp_prof_elec(j),
     >             temp_prof_ion(j)
	end do
	close(unit=35)
c
c	If logical variable "profile_scaling" is 'true' then scale the plasma
c		profile data using this information, otherwise these inputs are ignored.
c
	if ( profile_scaling )  then
	  !assume that first point is r=0
	  den_axis = den_prof(1); temp_axis_elec = temp_prof_elec(1)
	  temp_axis_ion = temp_prof_ion(1)
	  !scale profiles by axis and pedestal values
        do j=1,np_prof
          den_prof(j) = (den_prof(j)/den_axis + ped_d)/(1. + ped_d)
	    temp_prof_elec(j) =
     >    (temp_prof_elec(j)/temp_axis_elec + ped_t)/(1. + ped_t)
	    temp_prof_ion(j) =
     >    (temp_prof_ion(j)/temp_axis_ion + ped_t)/(1. + ped_t)
        end do
	else
	  !Normalize profiles using command line inputs
        den_prof = den_prof/(den0/1.e12_rknd)
	  temp_prof_elec = temp_prof_elec/(Te_0*1000._rknd)
	  temp_prof_ion = temp_prof_ion/(Ti_0*1000._rknd)
	end if
c
c	spline fit n,Te,Ti profiles and evaluate at r/a of current surface
c
	kord_prof = 3 !spline order for profile fitting
	allocate(dn_knot_array(np_prof+kord_prof), stat=istat)
      allocate(te_knot_array(np_prof+kord_prof), stat=istat)
      allocate(ti_knot_array(np_prof+kord_prof), stat=istat)
      allocate(spl_dn(np_prof), stat=istat)
      allocate(spl_te(np_prof), stat=istat)
      allocate(spl_ti(np_prof), stat=istat)
      call dbsnak(np_prof,r_prof,kord_prof,dn_knot_array)
      call dbsnak(np_prof,r_prof,kord_prof,te_knot_array)
      call dbsnak(np_prof,r_prof,kord_prof,ti_knot_array)
      call dbsint(np_prof,r_prof,den_prof,kord_prof,
     >            dn_knot_array,spl_dn)
      call dbsint(np_prof,r_prof,temp_prof_elec,kord_prof,
     >            te_knot_array,spl_te)
      call dbsint(np_prof,r_prof,temp_prof_ion,kord_prof,
     >              ti_knot_array,spl_ti)

	!evaluate spline fit at r/a of the test surface
      den_i = den0*dbsval(r,kord_prof,dn_knot_array,np_prof,spl_dn)
      den_e = den_i
      temp_i_keV = Ti_0*dbsval(r,kord_prof,ti_knot_array,np_prof,spl_ti)
      temp_e_keV = Te_0*dbsval(r,kord_prof,te_knot_array,np_prof,spl_te)
	!evaluate derivatives (dn/dr/n,dTe/dr/Te,dTi/dr/Ti) at r/a
      idrv = 1
      dn_dr_on = den0*dbsder(idrv,r,kord_prof,dn_knot_array,
     >                   np_prof,spl_dn)/(arad*den_i)
      dTe_dr_oTe = Te_0*dbsder(idrv,r,kord_prof,te_knot_array,
     >                   np_prof,spl_te)/(arad*temp_e_keV)
      dTi_dr_oTi = Ti_0*dbsder(idrv,r,kord_prof,ti_knot_array,
     >                   np_prof,spl_ti)/(arad*temp_i_keV)
	!write profile data to plasma_profile_check
      write(34,117) r,tb,den_i,tb,temp_e_keV,tb,temp_i_keV,tb,dn_dr_on,
     >              tb,dTe_dr_oTe,tb,dTi_dr_oTi
 117  format(f7.4,6(a1,e15.7))
      close(unit=34)

	!define mks densities
	den_e_mks = 1.e+06_rknd*den_e
      den_i_mks = 1.e+06_rknd*den_i
	!define T in eV
	temp_e_eV=temp_e_keV*1000._rknd
	temp_i_eV=temp_i_keV*1000._rknd
c
c	Calculate thermal velocities and loglambda 
c
      vth_i = sqrt(2._rknd*temp_i_keV*k_boltz/p_mass)  
      vth_e = sqrt(2._rknd*temp_e_keV*k_boltz/e_mass)
      if(temp_e_keV .gt. .05_rknd) then
       lnlambda = 25.3_rknd - 1.15_rknd*log10(den_e) + 
     >   2.3_rknd*log10(temp_e_eV)
      else
       lnlambda = 23.4_rknd - 1.15_rknd*log10(den_e) + 
     >   3.45_rknd*log10(temp_e_eV)
      endif
c
c   Calculate collisional times based on formule in S&N following
c    eq. (25). Also calculate collision frequencies.
c
      tau_ee = 3._rknd*sqrt(pi)*pi*ep02*e_mass**2*vth_e**3  
     >      /(den_e_mks*(qq**4)*lnlambda)
      tau_ii = 3._rknd*sqrt(pi)*pi*ep02*ion_mass**2*vth_i**3
     >      /(den_i_mks*(qq**4)*lnlambda)
	nu_e0=1._rknd/tau_ee * (3._rknd*sqrt(pi)/4._rknd)  !I'm pretty sure this (3*sqrt(pi)/4) factor is required
	nu_i0=1._rknd/tau_ii * (3._rknd*sqrt(pi)/4._rknd)  ! to get nu.  See Callen Eq. 2.86 and Don's new_coll. JL
													   ! In any case this is not used anywhere except as output
c
c   Form parallel friction coefficient matrices based on S&N's eqn. (C8):
c	The definitions of l_ij components are just before (C8)
c
      ll_e(1,1) = zee
      ll_e(1,2) = 1.5_rknd*zee
      ll_e(2,1) = ll_e(1,2)
      ll_e(2,2) = sqrt(2._rknd) + 13._rknd*zee/4._rknd
      ll_i(1,1) = 0._rknd; ll_i(1,2) = 0._rknd; ll_i(2,1) = 0._rknd
	ll_i(2,2) = sqrt(2._rknd)
      lam_e = ll_e
      lam_e(1,2) = -ll_e(1,2)
      lam_e(2,1) = -ll_e(1,2)
      lam_i = ll_i
c
c	Form scalar coefficients that are used in eqs. (C8)-(C13):
c
      factor_m_e = tau_ee/(den_e_mks*e_mass*bsq)       !in Eq C8
      factor_m_i = tau_ii/(den_i_mks*ion_mass*bsq)     !in Eq C8
      factor_n_e = tau_ee/(den_e_mks*e_mass)           !in Eq C8
      factor_n_i = tau_ii/(den_i_mks*ion_mass)         !in Eq C8
      inv_factor_m_e = den_e_mks*e_mass/(bsq*tau_ee)   !in Eq C9,C10
      inv_factor_m_i = den_i_mks*ion_mass/(bsq*tau_ii) !in Eq C9
      factor_Ee1 = den_e_mks*qq/(sqrt(bsq))            !in Eq C12
      factor_Ee2 = den_e_mks*qq*qq*tau_ee/(e_mass)     !in Eq C13
c
c	Optional screen output
c
      if(plots) then
        write(*,'("tau_ee, tau_ii = ",e15.7,3x,e15.7)') tau_ee,tau_ii
        write(*,36) factor_m_e, factor_m_i, factor_n_e, factor_n_i,
     >     inv_factor_m_e, inv_factor_m_i, factor_Ee1, factor_Ee2
  36    format("tau_ee/(den_e*m_e*<B**2>) = ", e15.7,/,
     >  "tau_ii/(den_i*m_i*<B**2>) = ", e15.7,/,
     >  "tau_ee/(den_e*m_e) = ", e15.7,/,
     >  "tau_ii/(den_i*m_i) = ", e15.7,/,
     >  "den_e*m_e*<B**2>/tau_ee = ", e15.7,/,
     >  "den_i*m_i*<B**2>/tau_ii = ", e15.7,/,
     >  "n_e*q/<B> = ", e15.7,/,
     >  "den_e*q**2*tau_ee/m_e = ", e15.7)
      end if
c
c	Allocate variables
c
      m_er = 3	  !Order of spline fit for Er, net_flux, qe, qi
      ne_pts = 200  !number of electric fields to use over search range
      iroot = 1	  !index used for root finding
      allocate(er(ne_pts), stat=istat)          !Er field array
      allocate(e_phi_kT(ne_pts), stat=istat)    !Er field array
      allocate(net_flux(ne_pts), stat=istat)    !gamma_i - gamma_e 
      allocate(gamma_e(ne_pts), stat=istat)     !elec particle flux
      allocate(gamma_i(ne_pts), stat=istat)     !ion particle flux
      
      allocate(gmi_e1(ne_pts), stat=istat)      !ion particle flux component
      allocate(gmi_e2(ne_pts), stat=istat)      !ion particle flux component
      allocate(gmi_i1(ne_pts), stat=istat)      !ion particle flux component
      allocate(gmi_i2(ne_pts), stat=istat)      !ion particle flux component
      allocate(gmi_i3(ne_pts), stat=istat)      !ion particle flux component
      allocate(gmi_Eprl(ne_pts), stat=istat)    !ion particle flux component

      allocate(gme_e1(ne_pts), stat=istat)      !electron particle flux component
      allocate(gme_e2(ne_pts), stat=istat)      !electron particle flux component
      allocate(gme_e3(ne_pts), stat=istat)      !electron particle flux component
      allocate(gme_i1(ne_pts), stat=istat)      !electron particle flux component
      allocate(gme_i2(ne_pts), stat=istat)      !electron particle flux component
      allocate(gme_Eprl(ne_pts), stat=istat)    !electron particle flux component
      
      allocate(jbs(ne_pts), stat=istat)         !bootstrap current
      allocate(q_e(ne_pts), stat=istat)         !elec heat flux
      allocate(q_i(ne_pts), stat=istat)         !ion heat flux 
      allocate(c_er(ne_pts), stat=istat)        !net_flux spline coeffs
      allocate(c_qe(ne_pts), stat=istat)        !qe spline coeffs 
      allocate(c_qi(ne_pts), stat=istat)        !qi spline coeffs 
      allocate(x_er(ne_pts+m_er), stat=istat)   !er spline knots 
c
c     Loop over electic fields for ambipolar solver. There are 4 different
c      electric field representations used here: er = E_r in Volts/cm,
c      Er_mks in Volts/m e_phi_kT = e<a>E_r/kT_e,  and abs_Er = abs(er*100) = 
c      input variable for velocity integration in Volts/m (dkes_coef subroutine).
c
      do i=1,ne_pts
c
c	Determine current electric field value, using either e<a>Er/kTe or Er (in V/cm)
c
	  if ( input_is_Er == .true. ) then
	    er(i) = ermin + 
     >            (ermax-ermin)*real(i-1,rknd)/real(ne_pts-1,rknd)
		!convert er to e_phi_kT
		e_phi_kT(i)=arad*er(i)*100._rknd/temp_e_eV
        else
	    e_phi_kT(i) = ermin + 
     >               (ermax-ermin)*real(i-1,rknd)/real(ne_pts-1,rknd)    !e*<a>*phi/k*T_e
          !convert e_phi_kT to er(volt/cm)
          er(i) = (temp_e_eV/(100._rknd*arad))*e_phi_kT(i)
	  end if
	!define Er in V/m and abs(Er)	
	  Er_mks=er(i)*100._rknd
        abs_Er = abs(Er_mks)
c 
c     Form thermodynamic forces from S&N eq. (10) and (C6):
c
        xe1 = -temp_e_keV*k_boltz*(dn_dr_on + dTe_dr_oTe
     >       + er_mks/temp_e_eV)
        xe2 = -temp_e_keV*k_boltz*dTe_dr_oTe
       
        xi1 = -temp_i_keV*k_boltz*(dn_dr_on + dTi_dr_oTi
     >       - zee*er_mks/temp_i_eV)
        xi2 = -temp_i_keV*k_boltz*dTi_dr_oTi
	
	  if(exact_sym) then
	    xe1 = -temp_e_keV*k_boltz*(dn_dr_on + dTe_dr_oTe)
          xi1 = -temp_i_keV*k_boltz*(dn_dr_on + dTi_dr_oTi)
        end if
c
c	Use beam info to define X_E else set to zero
c	
        if(beam_on) then 
	    xE = beam_source*4._rknd*(r**1.2_rknd)*(1._rknd-r*r)**2
	  else
	    xE=0._rknd
	  end if
c
c     The subroutine sugama_app_c_matrix provides transport coefficient matrices.
c         returns lmat = 5x5 transport matrix given in eqn. (C5) of S&N,   
c
	  lijs_write = .false.	  ! Writes additional screen output if true
	  pol_tor_cmpt = .false.  ! if(pol_tor_cmpt), save components needed for calculating poloidal/toroidal viscosities
        call sugama_app_c_matrix(lmat,fl_mat,lijs_inp_ext,abs_Er)
c
c     Calculate electron, ion particle and heat fluxes (actually
c      Q_e and Q_i are the heat fluxes divided by the temperature)
c      and bootstrap current based on S&N eq. (C5):
c
        gamma_e(i) = lmat(1,1)*xe1 + lmat(1,2)*xe2 + lmat(1,3)*xi1   
     >    + lmat(1,4)*xi2 + lmat(1,5)*xE
        q_e(i) = lmat(2,1)*xe1 + lmat(2,2)*xe2 + lmat(2,3)*xi1 
     >	+ lmat(2,4)*xi2 + lmat(2,5)*xE
        gamma_i(i) = lmat(3,1)*xe1 + lmat(3,2)*xe2 + lmat(3,3)*xi1 
     >    + lmat(3,4)*xi2 + lmat(3,5)*xE
        q_i(i) = lmat(4,1)*xe1 + lmat(4,2)*xe2 + lmat(4,3)*xi1 
     >    + lmat(4,4)*xi2 + lmat(4,5)*xE
        jbs(i) = lmat(5,1)*xe1 + lmat(5,2)*xe2 + lmat(5,3)*xi1
     >    + lmat(5,4)*xi2 + lmat(5,5)*xE
c
c	Ion flux components	     
c     
        gmi_e1(i) = lmat(3,1)*xe1
        gmi_e2(i) = lmat(3,2)*xe2
        gmi_i1(i) = -lmat(3,3)*temp_i_keV*k_boltz*dn_dr_on
        gmi_i2(i) = -(lmat(3,3)+lmat(3,4))*temp_i_keV*k_boltz*dTi_dr_oTi
        gmi_i3(i) = lmat(3,3)*qq*Er_mks*zee
        gmi_Eprl(i) = lmat(3,5)*xE
c     
c	Electron flux components (note not the same as the ion components)
c     
        gme_e1(i) = -lmat(1,1)*temp_e_keV*k_boltz*dn_dr_on
        gme_e2(i) = -(lmat(1,1)+lmat(1,2))*temp_e_keV*k_boltz*dTe_dr_oTe
	  gme_e3(i) = -lmat(1,1)*qq*Er_mks
        gme_i1(i) = lmat(1,3)*xi1
        gme_i2(i) = lmat(1,4)*xi2
        gme_Eprl(i) = lmat(1,5)*xE
c
c	Calculate net flux
c
        net_flux(i) = gamma_i(i) - gamma_e(i)
c
c     Write fluxes to run_sug file
c
        write(20,'(f10.6,6(a1,e15.7))') r,tb,er(i),tb,e_phi_kT(i),tb,
     >    gamma_i(i),tb, gamma_e(i), tb, q_i(i), tb, q_e(i)

      end do !i=1,ne_pts
      close(unit=12) !cmul_E_requests
c----------------------------------------------------------------------------------
c     Solve for self-consistent Er root by first spline fitting gamma_e,i vs.
c      Er and then calling zeroin.c
c
c   The following loop makes the first approximate pass at finding the
c    electric field root. It will identify only one root for further refinement.
c    Set i = ii to pick up the last root that is encountered (as er is varied from
c    ermin to ermax) or i = ne_pts - ii to pick up the first root encountered.
c

c
c	Identify location of sign change for root guess range
c
      do ii=2,ne_pts-1
       if(.not.search_down) i = ne_pts - ii
       if(search_down) i = ii
       if(net_flux(i)*net_flux(i+1) .le. 0.) iroot=i
      end do
c
c	Spline fit the net_flux and heat fluxes
c
      call dbsnak(ne_pts,er,m_er,x_er)
      call dbsint(ne_pts,er,net_flux,m_er,x_er,c_er)
      call dbsint(ne_pts,er,q_e,m_er,x_er,c_qe)
      call dbsint(ne_pts,er,q_i,m_er,x_er,c_qi)
c
c	Define search range and use functions zeroin and rad_flux to find root
c
	if ( iroot .gt. size(er) .or. iroot == 1 ) then
	  write(*,*) 'Root not found in Er search range!'
	  stop 22
	end if
	ermin = er(iroot-1)
	ermax = er(iroot+1)
	  
	tol_er = 1.e-20_rknd	!tolerance on zeroin function
      if(.not.er_fixed)
     >  er_root = zeroin(ermin,ermax,rad_flux,tol_er)
      if(er_fixed) then  !If desired, use a fixed Er profile
        er_root = fit0 + fit1*r + fit2*(r**2) + fit3*(r**3)
     >     + fit4*(r**4) + fit5*(r**5) + fit6*(r**6)
      endif
	if(plots) write(*,'("Er root = ",e15.7)') er_root
c
c	Use spline fits to calculate net_flux, qe, qi at root
c
      gi_m_ge = dbsval(er_root, m_er, x_er, ne_pts, c_er)
      qflux_i = dbsval(er_root, m_er, x_er, ne_pts, c_qi) !Not used
      qflux_e = dbsval(er_root, m_er, x_er, ne_pts, c_qe) !Not used
c
c     Calculate self-consistent flow velocities by calling sugama_app_c_matrix for the Er root
c
      pol_tor_cmpt = .true. ! save components needed for calculating poloidal/toroidal viscosities
      lijs_write = .false.  ! Writes additional screen output if true
	abs_Er = abs(er_root*100._rknd)
c
c	Calculate lmat (C5) using Er root
c
      call sugama_app_c_matrix(lmat,fl_mat,lijs_inp_ext,abs_Er)
c
c	Define thermodynamic forces using Er root (10),(C6)
c
      xi1 = -temp_i_keV*k_boltz*(dn_dr_on + dTi_dr_oTi
     >       - zee*er_root*100._rknd/temp_i_eV)
      xi1e =  qq*zee*er_root*100._rknd
      xi2 = -temp_i_keV*k_boltz*dTi_dr_oTi

      xe1 = -temp_e_keV*k_boltz*(dn_dr_on + dTe_dr_oTe
     >       + er_root*100._rknd/temp_e_eV)
      xe2 = -temp_e_keV*k_boltz*dTe_dr_oTe
c
c	Handle beams and/or exact symmetry case
c
	if(beam_on) then 
	  xE = beam_source*(1._rknd-r*r)**2
	else
	  xE = 0._rknd
	end if
      if(exact_sym) then
       xi1 = -temp_i_keV*k_boltz*(dn_dr_on + dTi_dr_oTi)
       xi1e = 0._rknd
       xe1 = -temp_e_keV*k_boltz*(dn_dr_on + dTe_dr_oTe)
      endif
	!optional screen output
      if(plots) write(*,38) xi1,xe1,xi2,xe2
 38   format("xi1 = ",e15.7,"  xe1 = ",e15.7,/,
     >   "xi2 = ",e15.7,"  xe2 = ",e15.7)
c
c	Calculate parallel flow and components using (C3) approximation
c	 Note: fl_mat= [(M_i + Lambda_I)]**-1[N_i]
c	 uprl_B = <u_||i*B>
c
      uprl_B = (fl_mat(1,1)*xi1 + fl_mat(1,2)*xi2)           !10/3/2007 corrections -QQ (minus sign?) JL
      uprl_B_Er_only = fl_mat(1,1)*xi1e                      !10/3/2007 corrections
      uprl0 = (fl_mat(1,1)*xi1 + fl_mat(1,2)*xi2)/sqrt(bsq)  !10/3/2007 corrections
c
c	Calculate components of bootstrap current from (C5)
c
      jbs_e = (lmat(5,1)*xe1 + lmat(5,2)*xe2)
      jbs_i = (lmat(5,3)*xi1 + lmat(5,4)*xi2)
      jbs_e_prl = lmat(5,5)*xE
      jbs_tot = jbs_e + jbs_i + jbs_e_prl
c
c	Calculate heat and particle fluxes from (C5)
c
      gam_e = (lmat(1,1)*xe1 + lmat(1,2)*xe2 + lmat(1,3)*xi1
     >   + lmat(1,4)*xi2 + lmat(1,5)*xE)
      qflux_e = (lmat(2,1)*xe1 + lmat(2,2)*xe2 + lmat(2,3)*xi1
     >   + lmat(2,4)*xi2 + lmat(2,5)*xE)
      gam_i = (lmat(3,1)*xe1 + lmat(3,2)*xe2 + lmat(3,3)*xi1
     >   + lmat(3,4)*xi2 + lmat(3,5)*xE)
	qflux_i = (lmat(4,1)*xe1 + lmat(4,2)*xe2 + lmat(4,3)*xi1
     >   + lmat(4,4)*xi2 + lmat(4,5)*xE)
c
c	Write fluxes and components of ion particle flux
c
      write(23,87) r, tb,gam_i,tb,gam_e,tb,qflux_i,tb,qflux_e,tb,
     >  lmat(3,1)*xe1,tb, lmat(3,2)*xe2,tb,
     >  lmat(3,3)*xi1,tb, lmat(3,4)*xi2,tb,
     >  lmat(3,5)*xE
 87   format(e15.7,9(a1,e15.7))
c
c	Write various quantities to files
c
      write(15,51) r, ra, xi1, uprl_B, bsq, arad, psip, chip
 51   format(e15.7,7(2x,e15.7))
      write(19,91) r, xi1, uprl_B, bsq, arad, vp, btheta, bzeta,
     >  psip, chip, xi1e
 91   format(e15.7,10(2x,e15.7))
c
c	Form matrix from (B1) for calculating flow components
c
      proj_mat(1,1) = 1._rknd; proj_mat(2,1) = 1._rknd
      proj_mat(1,2) = -bzeta/(bsq*chip*qq*Zee)				!Added Zee
      proj_mat(2,2) = btheta/(bsq*psip*qq*Zee)				!Added Zee
c
c	Calculate poloidal and toroidal flow components from (B1)
c
      upol = forpi2*chip*((proj_mat(1,1)*uprl_B/bsq)           !10/3/2007 corrections
     >  + proj_mat(1,2)*xi1)/vp                                !10/3/2007 corrections
      upol_Er_only = forpi2*chip*((proj_mat(1,1)               !10/3/2007 corrections
     > *uprl_B_Er_only)/bsq + proj_mat(1,2)*xi1e)/vp           !10/3/2007 corrections
      utor = forpi2*psip*((proj_mat(2,1)*uprl_B/bsq)           !10/3/2007 corrections
     >  + proj_mat(2,2)*xi1)/vp                                !10/3/2007 corrections
c
c	Optionally write to screen
c
      if(plots) then
        write(*,'("upol = ",e15.7,
     >   3x,"utor = ",e15.7,"   uprl = ", e15.7)') upol, utor, uprl0
        write(*,'("gam_e = ",e15.7,2x,"gam_i = ",e15.7)')
     >    gam_e, gam_i
        write(*,'("jbs_e = ",e15.7,2x,"jbs_i = ",e15.7,/,
     >   "gam_e/n = ",e15.7,2x,"gam_i/n = ",e15.7,/,
     >   "q_e/nT = ",e15.7,2x,"q_i/nT = ",e15.7)')
     >   jbs_e, jbs_i, gam_e/(1.e6_rknd*den_e),
     >   gam_i/(1.e+6_rknd*den_i),qflux_e/(1.e+6_rknd*den_e*temp_e_keV),
     >   qflux_i/(1.e+6_rknd*den_i*temp_i_keV)
      endif
c
c	Write to screen and file
c
	e_phi_kt_root = arad*er_root/(10._rknd*temp_e_keV) !1/10 is from 100/1000
      if(.not.plots) write(*,'(1x,f5.3,4(a1,e12.5))') r,tb,		
     >    er_root,tb,e_phi_kt_root,tb,upol,tb,utor
      write(11,'(1x,f5.3,6(a1,e20.14))') r,tb,
     >    er_root,tb,e_phi_kt_root,tb,upol,tb,
     >    utor,tb,uprl0,tb,upol_Er_only

		!
	if (idebug) then
		write(*,*) "  gam_e  ","       u_prl ","              qe "
     > ,"             qi", "          jbs"
	write(*,'(1x,e12.5,4(a1,e12.5))') ,		
     >    gam_e,tb,uprl_B,tb,qflux_e,tb,qflux_i,tb,jbs_tot
	endif
c------------------------------------------------------------------------
c   Calculate the ion and electron poloidal/toroidal viscosity
c    components, as defined in equation (B4) of Sugama, Nishimura.
c    These are written out to files as a function of radial
c    position (sqrt(flux)). The three velocity moments for
c    each species are given. These viscosities can then be
c    used in Sugama's eqn. (38) to obtain the stress tensor
c    components.
c
c
c	Open output files and write labels if (i_append==0)
c
      if(i_append .eq. 0) then
       open(unit=38,file="ion_viscosities",status="unknown")
       open(unit=39,file="electron_viscosities",status="unknown")
       open(unit=28,file="ion_viscosities_Er0",status="unknown")
       open(unit=29,file="electron_viscosities_Er0",status="unknown")
		!changed definition of nu_i (7/22/08) JL
      write(38,'("r",a1,"nu_i",a1,"er_root",a1,"mpp1i",	
     >  a1,"mpp2i",a1,"mpp3i",a1,"mpt1i",
     >  a1,"mpt2i",a1,"mpt3i",a1,"mtt1i",a1,"mtt2i",a1,"mtt3i")')
     >   tb,tb,tb,tb,tb,tb,tb,tb,tb, tb, tb
	!fixed label error here, changed definition of nu_e (7/22/08) JL
      write(39,'("r",a1,"nu_e",a1,"er_root",a1,"mpp1e",
     >  a1,"mpp2e",a1,"mpp3e",a1,"mpt1e",
     >  a1,"mpt2e",a1,"mpt3e",a1,"mtt1e",a1,"mtt2e",a1,"mtt3e")')
     >   tb,tb,tb,tb,tb,tb,tb,tb,tb, tb, tb
      else if(i_append .ne. 0) then
       open(unit=38,file="ion_viscosities",position="append",
     >      status="old")
       open(unit=39,file="electron_viscosities",position="append",
     >      status="old")
       open(unit=28,file="ion_viscosities_Er0",position="append",
     >      status="old")
       open(unit=29,file="electron_viscosities_Er0",
     >      position="append",status="old")
      end if
c
c	Calculate normalized coll frequencies
c
      nu_i0 = nu_i0*rmajor/(vth_i*abs(iota))    !nu*i
      nu_e0 = nu_e0*rmajor/(vth_e*abs(iota))    !nu*e

c
c	set up matrices from equation (B4)
c
      pre_mat(1,1) = chip*btheta/bsq
      pre_mat(1,2) = -qq*zee*psip*chip		!added zee
      pre_mat(2,1) = psip*bzeta/bsq
      pre_mat(2,2) = qq*zee*psip*chip			!added zee
      post_mat = transpose(pre_mat)
c
c	Calculate toroidal and poloidal viscosity coefficients from (B4)
c
      inter_mat = (forpi2/vp)*matmul(mnl1i,post_mat)
      mpt1i = matmul(pre_mat,inter_mat)
      inter_mat = (forpi2/vp)*matmul(mnl2i,post_mat)
      mpt2i = matmul(pre_mat,inter_mat)
      inter_mat = (forpi2/vp)*matmul(mnl3i,post_mat)
      mpt3i = matmul(pre_mat,inter_mat)
      inter_mat = (forpi2/vp)*matmul(mnl1e,post_mat)
      mpt1e = matmul(pre_mat,inter_mat)
      inter_mat = (forpi2/vp)*matmul(mnl2e,post_mat)
      mpt2e = matmul(pre_mat,inter_mat)
      inter_mat = (forpi2/vp)*matmul(mnl3e,post_mat)
      mpt3e = matmul(pre_mat,inter_mat)
      write(38,77) r,tb,nu_i0,tb,er_root,tb,mpt1i(1,1),
     >  tb,mpt2i(1,1),tb,mpt3i(1,1),
     >  tb,mpt1i(1,2),tb,mpt2i(1,2),tb,mpt3i(1,2),
     >  tb,mpt1i(2,2),tb,mpt2i(2,2),tb,mpt3i(2,2)
      write(39,77) r,tb,nu_e0,tb,er_root,tb,mpt1e(1,1),
     >  tb,mpt2e(1,1),tb,mpt3e(1,1),
     >  tb,mpt1e(1,2),tb,mpt2e(1,2),tb,mpt3e(1,2),
     >  tb,mpt1e(2,2),tb,mpt2e(2,2),tb,mpt3e(2,2)
      close(unit=38)
      close(unit=39)
  77  format(e15.7,11(a1,e15.7))
c------------------------------------------------------------------
c	Look at Er=0 case
c
	abs_Er = 1.e-4_rknd !This should be approximately 0
c
c	Calculate (C5) 5x5 matrix for Er~0
c
      pol_tor_cmpt = .true.
	lijs_write = .false.
      call sugama_app_c_matrix(lmat,fl_mat,lijs_inp_ext,abs_Er) 
c
c	Calculate thermodynamic forces for Er~0, in this case the symmetric case is identical
c
      xi1 = -temp_i_keV*k_boltz*(dn_dr_on + dTi_dr_oTi)
      xi1e =  0._rknd
      xi2 = -temp_i_keV*k_boltz*dTi_dr_oTi
	xe1 = -temp_e_keV*k_boltz*(dn_dr_on + dTe_dr_oTe)
      xe2 = -temp_e_keV*k_boltz*dTe_dr_oTe
      if(beam_on) then
	  xE = beam_source*(1._rknd-r*r)**2
	else
        xE = 0._rknd
	end if
c
c	Calculate bootstrap current from (C5)
c      
      jbs_e0 = (lmat(5,1)*xe1 + lmat(5,2)*xe2)
      jbs_i0 = (lmat(5,3)*xi1 + lmat(5,4)*xi2)
      jbs_e_prl0 = lmat(5,5)*xE
      jbs_tot0 = jbs_e + jbs_i + jbs_e_prl
      write(18,'(f10.6,15(a1,e15.7))') r,tb,den_e_mks,tb,temp_e_keV,tb,
     >  temp_i_keV,tb,jbs_e,tb,jbs_i,tb,jbs_e_prl,tb, jbs_tot, tb,
     >  jbs_e0,tb,jbs_i0,tb,jbs_e_prl0,tb, jbs_tot0, tb,
     >  gi_m_ge, tb, qflux_e, tb, qflux_i, tb, lmat(5,5)
c
c   As done above, calculate the ion and eletron poloidal/toroidal
c    viscosity components, as defined in equation (B4) of Sugama, et al.
c    for E_r ~ 0.
c
      if(i_append .eq. 0) then
        write(28,'("r",a1,"nuovth_i",a1,"er_root",a1,"mpp1i",a1,
     >    "mpp2i",a1,"mpp3i",a1,"mpt1i",
     >    a1,"mpt2i",a1,"mpt3i",a1,"mtt1i",a1,"mtt2i",a1,"mtt3i")')
     >    tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb
        write(29,'("r",a1,"nuovth_e",a1,"er_root",a1,"mpp1e",a1,
     >    "mpp2e",a1,"mpp3e",a1,"mpt1e",
     >    a1,"mpt2e",a1,"mpt3e",a1,"mtt1e",a1,"mtt2e",a1,"mtt3e")')
     >    tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb
      end if
      inter_mat = (forpi2/vp)*matmul(mnl1i,post_mat)
      mpt1i = matmul(pre_mat,inter_mat)
      inter_mat = (forpi2/vp)*matmul(mnl2i,post_mat)
      mpt2i = matmul(pre_mat,inter_mat)
      inter_mat = (forpi2/vp)*matmul(mnl3i,post_mat)
      mpt3i = matmul(pre_mat,inter_mat)
      inter_mat = (forpi2/vp)*matmul(mnl1e,post_mat)
      mpt1e = matmul(pre_mat,inter_mat)
      inter_mat = (forpi2/vp)*matmul(mnl2e,post_mat)
      mpt2e = matmul(pre_mat,inter_mat)
      inter_mat = (forpi2/vp)*matmul(mnl3e,post_mat)
      mpt3e = matmul(pre_mat,inter_mat)
      write(28,77) r,tb,nu_i0,tb,er_root,tb,mpt1i(1,1),
     >  tb,mpt2i(1,1),tb,mpt3i(1,1),
     >  tb,mpt1i(1,2),tb,mpt2i(1,2),tb,mpt3i(1,2),
     >  tb,mpt1i(2,2),tb,mpt2i(2,2),tb,mpt3i(2,2)
      write(29,77) r,tb,nu_e0,tb,er_root,tb,mpt1e(1,1),
     >  tb,mpt2e(1,1),tb,mpt3e(1,1),
     >  tb,mpt1e(1,2),tb,mpt2e(1,2),tb,mpt3e(1,2),
     >  tb,mpt1e(2,2),tb,mpt2e(2,2),tb,mpt3e(2,2)
      close(unit=28)
      close(unit=29)
      if(plots) then
c
c***** if plplot libraries are available, remove comments from following section
c
c      call plscolbg(255,255,255);
c      call plcol(15)
c      xmin = 1.2*minval(er)
c      xmax = 1.2*maxval(er)
c      ymin = min(minval(gamma_e),minval(gamma_i))
c      ymax = max(maxval(gamma_e),maxval(gamma_i))
c      call plenv(xmin,xmax,ymin,ymax,0,0)
c      call pllab('E#dr#u (Volt/cm)','#gG#di#u, #gG#de',' ')
c      call plcol(9)
c      call plmtex('t',-1._rknd,0._rknd,0,'electron flux')
c      call plline(ne_pts,er,gamma_e)
c      call plpoin(ne_pts,er,gamma_e,17)
c      call plcol(1)
c      call plline(ne_pts,er,gamma_i)
c      call plpoin(ne_pts,er,gamma_i,4)
c      call plmtex('t',-3._rknd,0._rknd,0,'ion flux')
c      call plend
c
c*****  end of plplot section
c
c	Write flux components vs. Er to file
c      
        open(unit=44,file="flux_vs_Efield",status="unknown")
         write(44,'("*",/,"Er(V/cm)",a1,"Gamma_i",a1,"Gamma_e",a1,
     >      "Q_i",a1,"Q_e",a1,"Gmi_e1",a1,"Gmi_e2",a1,"Gmi_i1",
     >      a1,"Gmi_i2",a1,"Gmi_i3",a1,"Gmi_Eprl")') tb,tb,tb,tb,
     >      tb,tb,tb,tb,tb,tb
        do i=1,ne_pts
          write(44,88) er(i), tb,gamma_i(i),tb,gamma_e(i),tb,q_i(i),
     >     tb,q_e(i),tb,gmi_e1(i), tb, gmi_e2(i), tb, gmi_i1(i),
     >     tb, gmi_i2(i), tb, gmi_i3(i), tb, gmi_Eprl(i)
        end do
 88     format(e15.7,10(a1,e15.7))
        close(unit=44)
      endif !if(plot)
c
c	Close files
c
      close(unit=11)
      close(unit=15)
      close(unit=18)
      close(unit=19)
      close(unit=20)
c      stop
      end
c
c------------------------------------------------------------------------------------
c
c	SUBROUTINES
c
c------------------------------------------------------------------------------------
      subroutine sugama_app_c_matrix(lm, fm, file_ext, abs_Er)
c
c   This subroutine calculates the transport coefficient matrix
c    given in S&N eq. (C5). It also provides matrices for poloidal
c    and toroidal viscosity coefficients and for obtaining parallel
c    flow and heat flow.
c
c    lmat = 5x5 transport matrix given in eqn. (C5) of S&N,
c    fl_mat = ion matrix needed for getting u_parallel and
c    q_parallel based on neglecting E|| and ion-electron
c    coupling effects.
c
	use penta_kind_mod
	use phys_const
      use lijs_transfer
	use dkes_coef_transfer
      use l_components
      use viscosity_pol_tor

      implicit none
      logical :: test_nw_deriv
	real(rknd) :: Ta, Tb, mb, abs_Er
      real(rknd), dimension(2,2) :: l2i, m2i, n2i, fm,
     >  l2e, m2e, n2e, temp_1_e, temp_1_i, m2e_inv, lam_e_inv,
     >  temp_2_e, temp_2_i, temp_2_e_inv,temp_e_inv, temp_i_inv,
     >  tmp_e, tmp_i, temp_ie, temp_Ee, temp_Ei, temp_EEE, ep11
      real(rknd), dimension(2,2) :: tmp1, tmp2, tmp3, tmp4,
     >  tmp5, tmp6
      real(rknd), dimension(5,5) :: lm
      real(rknd) :: cl1e, cl2e, cl3e, cm1e, cm2e,
     >  cm3e, cn1e, cn2e, cn3e,
     >  cl1i, cl2i, cl3i, cm1i, cm2i, cm3i, cn1i, cn2i, cn3i
      integer(iknd) :: i, j
      character*50 :: file_ext
      character*60 :: file_nm

      test_nw_deriv = .false.	!use Don's derivation of 3rd term in (C9) if true

c
c   S&N eq. (C8) rotation matrix:
c
      ep11(1,1)=1._rknd
	ep11(1,2)=0._rknd
	ep11(2,1)=0._rknd
	ep11(2,2)=0._rknd
c
c     Integrate electron La,Ma,Na coefficients using (Eq. 36)
c	 to get Laj,Maj,Naj thermal coefficients:
c
c	Here species "a" are electrons and "b" are ions, species' mass,
c	 charge, temperature and thermal velocity are defined below and
c	 used in subroutine dkes_coef.  Er/vth is defined for species a.
c
c	The variable lcoef_opt is 4 if the log of the coefficient file is to be taken,
c	 5 if not.
c
	ma=e_mass; qa=-qq; Ta=temp_e_keV*1000._rknd
	na=den_e_mks; vta=sqrt(2._rknd*Ta*qq/ma)
	mb=zee*ion_mass; qb=zee*qq; Tb=temp_i_keV*1000._rknd
	nb=den_i_mks; vtb=sqrt(2._rknd*Tb*qq/mb)
	abs_Er_ovth=abs_Er/vta  
	!integrate L
	file_nm = "lstar_lijs_" // file_ext	    
	lcoef_opt = 5
      call dkes_coef(cl1e,file_nm,1)	    !Le1
      call dkes_coef(cl2e,file_nm,2)		!Le2
      call dkes_coef(cl3e,file_nm,3)		!Le3
	!integrate M
        lcoef_opt = 4
        file_nm = "mstar_lijs_" // file_ext
      call dkes_coef(cm1e,file_nm,1)		!Me1
      call dkes_coef(cm2e,file_nm,2)		!Me2
      call dkes_coef(cm3e,file_nm,3)		!Me3
	!integrate N
        lcoef_opt = 5
        file_nm = "nstar_lijs_" // file_ext
      call dkes_coef(cn1e,file_nm,1)		!Ne1
      call dkes_coef(cn2e,file_nm,2)		!Ne2
      call dkes_coef(cn3e,file_nm,3)		!Ne3
c
c     Integrate ion La,Ma,Na coefficients using (Eq. 36)
c	 to get Laj,Maj,Naj thermal coefficients:
c
c	Here species "a" are ions and "b" are electrons, species' mass,
c	 charge, temperature and thermal velocity are defined below and
c	 used in subroutine dkes_coef.  Er/vth is defined for species a.
c
c	The variable lcoef_opt is 4 if the log of the coefficient file is to be taken,
c	 5 if not.
c 
	ma=zee*ion_mass;  qa=zee*qq; Ta=temp_i_keV*1000._rknd 
	na=den_i_mks
	vta=sqrt(2._rknd*Ta*qq/ma)
	mb=e_mass; Tb=temp_e_keV*1000._rknd; qb=-qq
	nb=den_i_mks
	vtb=sqrt(2._rknd*Tb*qq/mb)

	abs_Er_ovth=abs_Er/vta   

	!integrate L
        lcoef_opt = 5
        file_nm = "lstar_lijs_" // file_ext
      call dkes_coef(cl1i,file_nm,1)			!Li1
      call dkes_coef(cl2i,file_nm,2)			!Li2
      call dkes_coef(cl3i,file_nm,3)			!Li3
	!integrate M
        lcoef_opt = 4
        file_nm = "mstar_lijs_" // file_ext
      call dkes_coef(cm1i,file_nm,1)			!Mi1
      call dkes_coef(cm2i,file_nm,2)			!Mi2
      call dkes_coef(cm3i,file_nm,3)			!Mi3
	!integrate N
        lcoef_opt = 5
        file_nm = "nstar_lijs_" // file_ext
      call dkes_coef(cn1i,file_nm,1)			!Ni1
      call dkes_coef(cn2i,file_nm,2)			!Ni2
      call dkes_coef(cn3i,file_nm,3)			!Ni3
c
c     Form 2x2 matrices from (C8):
c
      n2e(1,1) = cn1e; n2e(1,2) = cn2e	!QQ (minus sign?)  JL
      n2e(2,1) = cn2e; n2e(2,2) = cn3e	!Removed minus sign 7/1/2009
      
      n2i(1,1) = cn1i; n2i(1,2) = cn2i
      n2i(2,1) = cn2i; n2i(2,2) = cn3i
      
      m2e(1,1) = cm1e; m2e(1,2) = cm2e
      m2e(2,1) = cm2e; m2e(2,2) = cm3e
      
      m2i(1,1) = cm1i; m2i(1,2) = cm2i
      m2i(2,1) = cm2i; m2i(2,2) = cm3i
      
      l2e(1,1) = cl1e; l2e(1,2) = cl2e
      l2e(2,1) = cl2e; l2e(2,2) = cl3e
      
      l2i(1,1) = cl1i; l2i(1,2) = cl2i
      l2i(2,1) = cl2i; l2i(2,2) = cl3i
c
c   For reference:
c    factor_m_e,i = tau/(den*mass*<B**2>)
c    factor_n_e,i = tau/(den*mass)
c    inv_factor_m_e,i = den*mass/(tau*<B**2>)
c    factor_Ee1 = den*e/sqrt(<B**2>)
c    factor_Ee2 = den*(e**2)*tau_ee/mass_e
c
	m2e = factor_m_e*m2e
      m2i = factor_m_i*m2i
      n2e = -factor_n_e*n2e	!QQ (minus sign?)  JL
      n2i = -factor_n_i*n2i	!QQ (minus sign?)  JL
c
c	Optionally write to screen
c
      if(lijs_write) then
        write(*,37) m2e(1,1),m2i(1,1),n2e(1,1),
     >   n2i(1,1),l2e(1,1),l2i(1,1),lam_e(1,1),lam_i(2,2)
      
  37    format("M_e_11 = ",e15.7,"  M_i_11 = ",e15.7,/,
     >  "N_e_11 = ",e15.7,"  N_i_11 = ",e15.7,/,
     >  "L_e_11 = ",e15.7,"  L_i_11 = ",e15.7,/
     >  "lam_e_11 = ",e15.7,"  lam_i_22 = ",e15.7)
	endif
c
c   if(pol_tor_cmpt) save components needed for calculating poloidal/toroidal viscosities
c	using (B4)
c
      if(pol_tor_cmpt) then
      
       mnl1i(1,1)=cm1i; mnl1i(1,2)=cn1i
       mnl1i(2,1)=cn1i; mnl1i(2,2)=cl1i
       
       mnl1e(1,1)=cm1e; mnl1e(1,2)=cn1e
       mnl1e(2,1)=cn1e; mnl1e(2,2)=cl1e
       
       mnl2i(1,1)=cm2i; mnl2i(1,2)=cn2i
       mnl2i(2,1)=cn2i; mnl2i(2,2)=cl2i
       
       mnl2e(1,1)=cm2e; mnl2e(1,2)=cn2e
       mnl2e(2,1)=cn2e; mnl2e(2,2)=cl2e
       
       mnl3i(1,1)=cm3i; mnl3i(1,2)=cn3i
       mnl3i(2,1)=cn3i; mnl3i(2,2)=cl3i
       
       mnl3e(1,1)=cm3e; mnl3e(1,2)=cn3e
       mnl3e(2,1)=cn3e; mnl3e(2,2)=cl3e
       
      end if     !if(pol_tor_cmpt)
c
c   Calculate factors used in Eq. (C9-C13)
c
	temp_1_e = m2e + lam_e
	temp_1_i = m2i + lam_i
	
      !invert some matrices
	call mat_2by2_inverse(m2e,m2e_inv)
      call mat_2by2_inverse(lam_e,lam_e_inv)
      call mat_2by2_inverse(temp_1_e, temp_e_inv)
      call mat_2by2_inverse(temp_1_i, temp_i_inv)

      temp_2_e = m2e_inv + lam_e_inv	
      call mat_2by2_inverse(temp_2_e,temp_2_e_inv)
c
c   Form fm matrix. This is (M_i + Lambda_I)**-1 N_i, which is used
c    to get U|| and Q|| from eq. (C3) taking E|| = 0 and neglecting
c    electron-ion coupling.
c	
	fm = matmul(temp_i_inv,n2i)
c
c   Calculate terms in S&N eq. (C9) for electrons:
c	Note that for electrons delta_ai=0
c
      tmp1 = matmul(temp_e_inv,n2e)
      tmp2 = matmul(n2e,tmp1)
      tmp_e = l2e - inv_factor_m_e*tmp2
c
c   Calculate terms in S&N eq. (C9) for ions:
c
      tmp1 = matmul(temp_i_inv,n2i)
      tmp2 = matmul(n2i,tmp1)
      tmp_i = l2i - inv_factor_m_i*tmp2
c
c	Optionally write to screen
c
      if(lijs_write) then
        write(*,111) l2i(1,1), l2i(1,2), l2i(2,1), l2i(2,2),
     >   tmp2(1,1), tmp2(1,2), tmp2(2,1), tmp2(2,2)
 111  format("L_ion:   ",e15.7,3x,e15.7,/,9x,e15.7,3x,e15.7,//,
     >       "ion_part2: ",e15.7,3x,e15.7,/,11x,e15.7,3x,e15.7)
      endif

      if(test_nw_deriv .eq. .false.) then
c
c   Form the third term in eq. (C9) for ions based on S&N:
c
        tmp1 = matmul(temp_i_inv,n2i)
        tmp2 = matmul(ep11,tmp1)
        tmp3 = matmul(temp_2_e_inv,tmp2)
        tmp4 = matmul(ep11,tmp3)
        tmp5 = matmul(temp_i_inv,tmp4)
        tmp6 = matmul(n2i,tmp5)
       
        if(lijs_write) then
          write(*,114) "tmp2",tmp2(1,1), tmp2(1,2), tmp2(2,1), tmp2(2,2)
          write(*,114) "tmp3",tmp3(1,1), tmp3(1,2), tmp3(2,1), tmp3(2,2)
          write(*,114) "tmp4",tmp4(1,1), tmp4(1,2), tmp4(2,1), tmp4(2,2)
          write(*,114) "tmp5",tmp5(1,1), tmp5(1,2), tmp5(2,1), tmp5(2,2)
          write(*,114) "tmp6",tmp6(1,1), tmp6(1,2), tmp6(2,1), tmp6(2,2)
 114      format(//,a4": ",e15.7,3x,e15.7,/,6x,e15.7,3x,e15.7)
		write(*,112) tmp6(1,1), tmp6(1,2), tmp6(2,1), tmp6(2,2)
 112      format(/,"ion_part3: ",e15.7,3x,e15.7,/,11x,e15.7,3x,e15.7,/)
        end if
	
	  tmp_i = tmp_i + inv_factor_m_e*tmp6
      else if(test_nw_deriv .eq. .true.) then
c
c   Form the third term in eq. (C9) for ions based on my (Don's) derivation:
c
	  tmp1 = matmul(temp_i_inv,n2i)
        tmp2 = matmul(ep11,tmp1)
        tmp3 = matmul(lam_e,tmp2)
        tmp4 = matmul(temp_e_inv,tmp3)
        tmp5 = matmul(lam_e,tmp4)
        tmp6 = matmul(ep11,tmp5)
        if(lijs_write) then
          write(*,114) "tmp2",tmp2(1,1), tmp2(1,2), tmp2(2,1), tmp2(2,2)
          write(*,114) "tmp3",tmp3(1,1), tmp3(1,2), tmp3(2,1), tmp3(2,2)
          write(*,114) "tmp4",tmp4(1,1), tmp4(1,2), tmp4(2,1), tmp4(2,2)
          write(*,114) "tmp5",tmp5(1,1), tmp5(1,2), tmp5(2,1), tmp5(2,2)
          write(*,114) "tmp6",tmp6(1,1), tmp6(1,2), tmp6(2,1), tmp6(2,2)
        end if
        tmp1(1:2,1:2) = 0._rknd; tmp2(1:2,1:2) = 0._rknd
        tmp2 = matmul(temp_i_inv,tmp6)
        tmp1 = matmul(n2i,tmp2)
        if(lijs_write) then
          write(*,112) tmp1(1,1), tmp1(1,2), tmp1(2,1), tmp1(2,2)
        end if

	  tmp_i = tmp_i - inv_factor_m_e*tmp1
      endif	!if test_nw_deriv
c
c   Begin forming (C5) 5x5 matrix from (C9) ii, ee transport coefficients
c
      if(.not.dkes_limit) then
        lm(1,1) = tmp_e(1,1); lm(1,2) = tmp_e(1,2)
        lm(2,1) = tmp_e(2,1); lm(2,2) = tmp_e(2,2)
        lm(3,3) = tmp_i(1,1); lm(3,4) = tmp_i(1,2)
        lm(4,3) = tmp_i(2,1); lm(4,4) = tmp_i(2,2)
      endif
c	
c  To revert to 'dkes limit' just use Le, Li terms
c 
      if(dkes_limit) then
        lm(1,1) = l2e(1,1)/qq**2; lm(1,2) = l2e(1,2)/qq**2
        lm(2,1) = l2e(2,1)/qq**2; lm(2,2) = l2e(2,2)/qq**2
        lm(3,3) = l2i(1,1)/qq**2; lm(3,4) = l2i(1,2)/qq**2
        lm(4,3) = l2i(2,1)/qq**2; lm(4,4) = l2i(2,2)/qq**2
      end if
c
c     Calculate off-diagonal e-i blocks using S&N eq. (C10):
c
c      temp_ie = - inv_factor_m_e*n2e*temp_e_inv*lam_e*ep11
c     >  *temp_i_inv*n2i
      tmp1 = matmul(temp_i_inv,n2i)
      tmp2 = matmul(ep11,tmp1)
      tmp3 = matmul(lam_e,tmp2)
      tmp4 = matmul(temp_e_inv,tmp3)
      tmp5 = matmul(n2e,tmp4)
      temp_ie = - inv_factor_m_e*tmp5
c	
c	Use these terms in (C5) matrix
c
      if(.not.dkes_limit) then
	!   eq. (C10) ei components = temp_ie:
       lm(1,3) = temp_ie(1,1); lm(1,4) = temp_ie(1,2)
       lm(2,3) = temp_ie(2,1); lm(2,4) = temp_ie(2,2)
	!   eq. (C10) ie components = transpose(temp_ie):
       lm(3,1) = temp_ie(1,1); lm(3,2) = temp_ie(2,1)
       lm(4,1) = temp_ie(1,2); lm(4,2) = temp_ie(2,2)
c
c     for conventional limit - no electron/ion coupling:
c
      else if(dkes_limit) then
       lm(1,3) = 0._rknd; lm(1,4) = 0._rknd
       lm(2,3) = 0._rknd; lm(2,4) = 0._rknd
       lm(3,1) = 0._rknd; lm(3,2) = 0._rknd
       lm(4,1) = 0._rknd; lm(4,2) = 0._rknd
      endif
c
c     Calculate bootstrap and electric field terms based on eqns. (C11)-(C13):
c      temp_Ee = factor_Ee1*temp_e_inv*n2e      
c      
c   First eq. (C11):
      tmp1 = matmul(temp_e_inv,n2e)
      temp_Ee = factor_Ee1*tmp1
c	Use these terms in (C5) matrix
      lm(5,1) = temp_Ee(1,1); lm(5,2) = temp_Ee(1,2)
      lm(1,5) = -temp_Ee(1,1); lm(2,5) = -temp_Ee(1,2)
c
c	Handle beam case
c
      if(beam) then
       lm(1,5) = temp_Ee(1,1); lm(2,5) = temp_Ee(1,2)
      endif
c
c   Next eq. (C12):        
c     temp_Ei = -factor_Ee1*temp_e_inv*m2e*ep11*temp_i_inv*n2i      
c
      tmp1 = matmul(temp_i_inv,n2i)
      tmp2 = matmul(ep11,tmp1)
      tmp3 = matmul(m2e,tmp2)
      tmp4 = matmul(temp_e_inv,tmp3)
      temp_Ei = -factor_Ee1*tmp4
c	Use these terms in (C5) matrix
      lm(5,3) = temp_Ei(1,1); lm(5,4) = temp_Ei(1,2)
      lm(3,5) = -temp_Ei(1,1); lm(4,5) = -temp_Ei(1,2)
c
c   Finally, eq. (C13):
c      
      temp_EEE = -factor_Ee2*(lam_e_inv - temp_e_inv)
c	Use this term in (C5) matrix
      lm(5,5) = temp_EEE(1,1)
      return
      end
c
c-----------------------------------------------------------------------
c	This subroutine inverts a 2x2 matrix, used in sugama_app_c_matrix
c
      subroutine mat_2by2_inverse(a,a_inverse)
	use penta_kind_mod
      real(rknd), DIMENSION(2,2):: a, a_inverse
      real(rknd) :: denom
      denom = a(1,1)*a(2,2) - a(1,2)*a(2,1)
c      if(abs(denom) .lt. 1.d-4) write(*,*) denom
      if(denom .eq. 0.) then
       write(*,'("error: singular matrix")')
       write(*,*) a(1,1),a(1,2)
       write(*,*) a(2,1),a(2,2)
       stop 25
      endif
      a_inverse(1,1) = a(2,2)/denom
      a_inverse(2,2) = a(1,1)/denom
      a_inverse(1,2) = -a(1,2)/denom
      a_inverse(2,1) = -a(2,1)/denom
      return
      end
c
c----------------------------------------------------------------------
c	This function is used in conjunction with 'zeroin' to find
c	 the solution to gamma_e=gamma_i
c
      function rad_flux(xx_er) Result(rad_flux_out)
      use penta_kind_mod
      use er_solve
      use bspline
      implicit none
      real(rknd) :: rad_flux_out
      real(rknd) :: xx_er
      rad_flux_out = dbsval(xx_er, m_er, x_er, ne_pts, c_er)
      return
      end function rad_flux
