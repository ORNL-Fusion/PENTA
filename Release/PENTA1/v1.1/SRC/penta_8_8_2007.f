c  To compile, use: 
c     build_app_c_sugama_linux script for Linux systems with the Lahey compiler
c     build_penta_macintel For Mac Intel systems with the ifort compiler
c     build_penta_ibm for IBM's with the xlf compiler
c
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
c  4  flux surface number
c  5  parameter that determines whether output data should be
c      written into new files (=0) or appended to
c      existing files (=1). Typically,when a script is run for a
c      sequence of flux surfaces, this is set to 0 for the first
c      surface and to 1 for subsequent surfaces
c  6  extension of the profile_data_*** file (i.e., the *** text)
c      that comes fromextracting the appropriate profiles from
c      the VMEC run using the profile_extractor.f code
c  7  ech/ich options (obsolete - only sets central electron
c      temperature if it is not set elsewhere
c  8  central ion temperature
c  9  central density
c 10  edge density pedestal
c 11  edge temperature pedestal
c 12  central electron temperature
c 13  unique identifier on plasma profiles file, i.e., the ***
c      text in the filename plasma_profiles***.dat. This allows
c      multiple profiles to be run without having to rename files
c
c
      use Vblck
      use lijs_transfer
      use l_components
      use er_solve
      use bspline
      use viscosity_pol_tor
      implicit none
      logical, parameter :: plots=.false.        !true requires plplot libraries
      logical, parameter :: er_fixed=.false.     !option for fixed electric field profile
      logical, parameter :: beam_on=.false.      !option for external parallel momentum source
      logical, parameter :: old_profiles=.FALSE. !options for den/temp profiles
      logical, parameter :: plasma_profile_data=.TRUE.  !option to read in plasma profiles
      logical, parameter :: exact_sym=.false.    !if true, turns off Er terms in Xi1, Xe1
      logical, parameter :: search_down=.true.   !direction of Er root search (if false, then upward search is done)
      logical, parameter :: new_coll = .true.   !corrected collision frequencies (9/25/2007)
	logical, parameter :: idebug = .false.	    !More text output
      real*8 :: mu, wtovth, lnlambda, z, dsdr
      real*8 :: result, Ti_0, Te_0, r, ra, den0
      real*8 :: vth_e, vth_i, dn_dr, er_root
      real*8 :: p_mass, e_mass, tau_ee, tau_ii
      real*8 :: dTe_dr_oTe, dTi_dr_oTi, dn_dr_on
      real*8 :: ion_flux,elec_flux,lambda
      real*8 :: xe1, xe2, xi1, xi2, xe, ep02, pi, forpi2, forpi, xi1e
      real*8 :: e_phi_kt_root, uprl_B, bzeta, btheta, uprl0
      real*8 :: chip, psip, vp, bsq, upol, utor, upol_Er_only,
     >          uprl_B_Er_only
      real*8 :: jbs_e, jbs_i, jbs_e_prl, jbs_tot, fit_scl,
     >    qflux_e, qflux_i, gam_e, gam_i, gi_m_ge, ped_d, ped_t,
     >    jbs_e0, jbs_i0, jbs_e_prl0, jbs_tot0, beam_source,
     >    fit0,fit1,fit2,fit3,fit4,fit5,fit6, temp_axis_elec,
     >    temp_axis_ion, den_axis, iota, chk
      real*8, dimension(5,5) :: lmat
      real*8, dimension(2,2) :: fl_mat, proj_mat
      real*8, dimension(2,2) :: mpt1i, mpt2i, mpt3i,
     >  mpt1e, mpt2e, mpt3e, inter_mat, pre_mat, post_mat
      real*8 :: ymin, ymax, xmin, xmax
      real*8, parameter :: k_boltz = 1.6e-16    !Boltzmann factor in Joules/keV
      real*8, parameter :: ep0 = 8.854e-12
      real*8, parameter :: chrg_i = 1.6e-19
      real*8, dimension(:), allocatable :: e_phi_kT,
     1   gamma_e, gamma_i, q_e, q_i, er0, jbs,gmi_e1,
     2   gmi_e2, gmi_i1, gmi_i2, gmi_i3, gmi_Eprl
      real*8, dimension(:), allocatable :: r_m, r_norm,
     >  chip_profile, psip_profile, btheta_profile,
     >  bzeta_profile, vp_profile, bsq_profile,iota_profile,
     >  r_prof, temp_prof_elec, den_prof, dn_knot_array,
     >  te_knot_array, spl_dn, spl_te,temp_prof_ion,
     >  ti_knot_array, spl_ti
      integer :: i, ii, k, iroot, i_append, istat, js,
     1   js_min, js_max, j, index, np_prof, kord_prof, idrv
      character*10 :: ch_dum
      character*50 :: lijs_inp_ext, plasma_prof_file
      CHARACTER*60 :: arg1, arg2, arg3, arg4, arg5, arg6, arg7,
     1   arg8, arg9, arg10,arg11, arg12, arg13, profile_data
      CHARACTER*1 :: tb
      real*8, EXTERNAL :: rad_flux, zeroin
      tb = char(9)
      ep02 = ep0**2
      one = 1.d0
      pi = 4.d0*atan(one)
      pfac = 0.5d0*sqrt(pi)
      beam = beam_on
      lijs_write = .false.
c   Optional fixed electric field fit:
      if(er_fixed) then
c      fit0 = 110.02d0
c      fit1 = -1663.1d0
c      fit2 = 9143.4d0
c      fit3 = -25838.d0
c      fit4 = 38809.d0
c      fit5 = -29674.d0
c      fit6 = 9109.d0
      fit_scl = 1.e-4
      fit0 = 0.d0*fit_scl
      fit1 = -1.92d0*fit_scl
      fit2 = -115.4d0*fit_scl
      fit3 = 117.3d0*fit_scl
      fit4 = 0.d0
      fit5 = 0.d0
      fit6 = 0.d0
      end if
c      er_root = fit0 + fit1*r + fit2*(r**2) + fit3*(r**3) + fit4*(r**4) + fit5*(r**5) + fit6*(r**6)
c
      Te_0 = 0.            !intial value, indicating one is not provided
      forpi2 = 4.*pi*pi
      forpi = 4.*pi
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
      CALL getarg(13, arg13)
      lijs_inp_ext = trim(adjustl(arg1))
      read(arg2,'(f10.5)') ermin
      read(arg3,'(f10.5)') ermax
      read(arg4,'(i3)') js
      read(arg5,'(i1)') i_append
      read(arg8,'(f10.5)') Ti_0
      read(arg9,'(e15.7)') den0
      read(arg10,'(f10.5)') ped_d
      read(arg11,'(f10.5)') ped_t
      read(arg12,'(f10.5)') Te_0
      if(beam_on) then
       open(unit=33,file="beam_params",status="old")
       read(33,*) beam_source
      end if
      if(i_append .eq. 0) then
         open(unit=20,file="run_sug",status="unknown")
         open(unit=23,file="fluxes_vs_r",status="unknown")
         open(unit=11,file="flux_avg_data0",status="unknown")
         open(unit=19,file="flux_avg_data1",status="unknown")
         open(unit=18,file="jbs_vs_r",status="unknown")
         open(unit=15,file="flow_viz1.dat",status="unknown")
         write(11,'("*",/,"<r>/<a>",a1,"Er(V/cm)",a1,"<a>Er/kTe",
     >    a1,"Upol(contra)",a1,"Utor(contra)",a1,"Uprl",a1,
     >    "Upol(contra)_Er_only")')
     >     tb, tb, tb, tb, tb, tb
         write(18,'("*",/,"r",a1,"ne",a1,"Te",
     >    a1,"T_i",a1,"Jbs_e",a1,"Jbs_i",a1,"Jbs_E||",a1,"Jbs_total",
     >    a1,"Jbs_e(Er=0)",a1,"Jbs_i(Er=0)",a1,"Jbs_E||(Er=0)",a1,
     >    "Jbs_total(Er=0)",
     >    a1,"net_flux",a1,"qflux_e",a1,"qflux_i")')
     >     tb, tb, tb, tb, tb, tb, tb, tb, tb, tb, tb, tb, tb, tb
         write(23,'("*",/,"<r>/<a>",a1,"Gamma_i",a1,"Gamma_e",a1,
     >      "Q_i",a1,"Q_e",a1,"Gmi_e1",a1,"Gmi_e2",a1,"Gmi_i1",
     >      a1,"Gmi_i2",a1,"Gmi_Eprl",a1,
     >      "Gme_e1",a1,"Gme_e2",a1,"Gme_i1",
     >      a1,"Gme_i2",a1,"Gme_Eprl")') tb,tb,tb,tb,tb,
     >      tb,tb,tb,tb,tb,tb,tb,tb,tb
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
      open(unit=9,file="lstar_lijs_"//lijs_inp_ext, status="old")
      read (9, *) dsdr, wtovth
      if(plots) write(*,*) dsdr, wtovth
      close(unit=9)
      if(plots) then
      call plscol0(0, 255, 255, 255)
      call plscol0(15, 0, 0, 0)
      call plssub(1,1)
      call plinit()
      endif
      write(20,'("*",/,"E_r (V/cm)",a1,"gamma_i",a1,"gamma_e",
     >     a1,"Q_i",a1,"Q_e")')
     >     tb, tb, tb, tb
      mu = 1.d0         !ion to proton mass ratio
      lnlambda = 16.    !typical Coulomb logarithm
      z = 1.            !ion charge
c   ECH parameters:
      if(arg7 .eq. "ech") then
        if(Te_0 .eq. 0.) Te_0 = 1.5
      endif
c   ICH parameters:
      if(arg7 .eq. "ich") then
        if(Te_0 .eq. 0.) Te_0 = 0.5
      endif
c
c     Ion parameters:
       mass_i = 1.d0
       masrat_i = 1.d0/1837.d0
       zeff = 1.d0
c      Electron parameters:
       mass_e = 1.d0/1837.d0
       masrat_e = 1837.d0
c
c    Read VMEC profile data file and select data for current flux position:
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
      read(25,*) arad
      read(25,'(a10)') ch_dum
      do j = js_min,js_max
c       read(25,33) index,tb,r_m(j),tb,r_norm(j),tb,
c     >    chip_profile(j),tb,psip_profile(j),tb,btheta_profile(j),
c     >    tb,bzeta_profile(j),tb,vp_profile(j),
c     >    tb,bsq_profile(j),tb,iota_profile(j)
       read(25,*) index,r_m(j),r_norm(j),
     >    chip_profile(j),psip_profile(j),btheta_profile(j),
     >    bzeta_profile(j),vp_profile(j),
     >    bsq_profile(j),iota_profile(j)
      end do
  33  format(i2,9(a1,e15.7))
      chip = chip_profile(js); psip = psip_profile(js)
      bsq = bsq_profile(js); vp = vp_profile(js)
      btheta = btheta_profile(js); bzeta = bzeta_profile(js)
      ra = r_m(js); r = r_norm(js); iota = iota_profile(js)
c
      pre_mat(1,1) = chip*btheta/bsq
      pre_mat(1,2) = -chrg_i*psip*chip
      pre_mat(2,1) = psip*bzeta/bsq
      pre_mat(2,2) = chrg_i*psip*chip
      post_mat = transpose(pre_mat)
      
      if(.not.plasma_profile_data) then
      
      den_i = den0*(ped_d + (1. - r**6)**2)/(1. + ped_d)
      den_e = den_i
      temp_i = Ti_0*(ped_t + (1. - r*r)**2)/(1. + ped_t)
      temp_e = Te_0*(ped_t + (1. - r*r)**2)/(1. + ped_t)
      dn_dr = -16.*den0*(r**5)*(1. - r**6)/(arad*(1. + ped_d))
      dTe_dr_oTe = -4.*r*(1.-r*r)/(arad*(ped_t + (1. - r*r)**2))
      dTi_dr_oTi = -4.*r*(1.-r*r)/(arad*(ped_t + (1. - r*r)**2))
c
      if(old_profiles) then
       den_i = den0*(1.-.2*r)
       den_e = den_i
       temp_i = Ti_0*(((1. - r*r)**2 + 0.6)**2)/(1.6**2)
       temp_e = Te_0*(((1. - r*r)**2 + 0.6)**2)/(1.6**2)
       dn_dr = -.2*den0/arad
       dTe_dr_oTe = -8.*r*(1. - r*r)/(((1. - r*r)**2 + 0.6)*arad)
       dTi_dr_oTi = -8.*r*(1. - r*r)/(((1. - r*r)**2 + 0.6)*arad)
      endif
c
      dn_dr_on = dn_dr/den_i
      
      end if    !.not.plasma_profile_data

      if(plasma_profile_data) then
       plasma_prof_file =
     >     "plasma_profiles"//trim(adjustl(arg13))//".dat"
       open(unit=35,file=plasma_prof_file,status="old")
       if(i_append .eq. 0) then
         open(unit=34,file="plasma_profiles_check",status="unknown")
	 write(34,116) tb, tb, tb, tb, tb, tb
 116     format("r",a1,"den",a1,"temp_e",a1,"temp_i",a1,"dndr_ovr_n"
     >          ,a1,"dTe_dr_oTe",a1,"dTi_dr_oTi")
       else if(i_append .eq. 1) then
         open(unit=34,file="plasma_profiles_check",
     >            position="append",status="old")
       end if
       kord_prof = 3
       read(35,*) np_prof
       allocate(r_prof(np_prof), stat=istat)
       allocate(den_prof(np_prof), stat=istat)
       allocate(temp_prof_elec(np_prof), stat=istat)
       allocate(temp_prof_ion(np_prof), stat=istat)
       allocate(dn_knot_array(np_prof+kord_prof), stat=istat)
       allocate(te_knot_array(np_prof+kord_prof), stat=istat)
       allocate(ti_knot_array(np_prof+kord_prof), stat=istat)
       allocate(spl_dn(np_prof), stat=istat)
       allocate(spl_te(np_prof), stat=istat)
       allocate(spl_ti(np_prof), stat=istat)
       do j=1,np_prof
        read(35,*) r_prof(j),den_prof(j),temp_prof_elec(j),
     >       temp_prof_ion(j)
       end do
        den_axis = den_prof(1); temp_axis_elec = temp_prof_elec(1)
	temp_axis_ion = temp_prof_ion(1)
       do j=1,np_prof
        den_prof(j) = (den_prof(j)/den_axis + ped_d)/(1. + ped_d)
	temp_prof_elec(j) =
     >    (temp_prof_elec(j)/temp_axis_elec + ped_t)/(1. + ped_t)
	temp_prof_ion(j) =
     >    (temp_prof_ion(j)/temp_axis_ion + ped_t)/(1. + ped_t)
       end do
       call dbsnak(np_prof,r_prof,kord_prof,dn_knot_array)
       call dbsnak(np_prof,r_prof,kord_prof,te_knot_array)
       call dbsnak(np_prof,r_prof,kord_prof,ti_knot_array)
       call dbsint(np_prof,r_prof,den_prof,kord_prof,
     >     dn_knot_array,spl_dn)
       call dbsint(np_prof,r_prof,temp_prof_elec,kord_prof,
     >     te_knot_array,spl_te)
       call dbsint(np_prof,r_prof,temp_prof_ion,kord_prof,
     >     ti_knot_array,spl_ti)
       den_i = den0*dbsval(r,kord_prof,dn_knot_array,np_prof,spl_dn)
       den_e = den_i
       temp_i = Ti_0*dbsval(r,kord_prof,ti_knot_array,np_prof,spl_ti)
       temp_e = Te_0*dbsval(r,kord_prof,te_knot_array,np_prof,spl_te)
       idrv = 1
       dn_dr = den0*dbsder(idrv,r,kord_prof,dn_knot_array,
     >                       np_prof,spl_dn)/arad
       dn_dr_on = dn_dr/den_i
       dTe_dr_oTe = Te_0*dbsder(idrv,r,kord_prof,te_knot_array,
     >                       np_prof,spl_te)/(arad*temp_e)
       dTi_dr_oTi = Ti_0*dbsder(idrv,r,kord_prof,ti_knot_array,
     >                       np_prof,spl_ti)/(arad*temp_i)
       write(34,117) r,tb,den_i,tb,temp_e,tb,temp_i,tb,dn_dr_on,tb,
     >    dTe_dr_oTe,tb,dTi_dr_oTi
 117   format(f7.4,6(a1,e15.7))
       close(unit=35); close(unit=34)
      end if           !if(plasma_profile_data)
c
      trat_i = temp_e/temp_i

      trat_e = temp_i/temp_e
      p_mass = 1.67d-27; e_mass = 9.11d-31
      qq = 1.6d-19
      vth_i = sqrt(2.*temp_i*1.6e-16/1.673e-27)
      vth_e = sqrt(2.*temp_e*1.6e-16/9.11e-31)
      if(temp_e .gt. .05) then
       lnlambda = 25.3 - 1.15*log10(den_e) + 2.3*log10(1000.*temp_e)
      else
       lnlambda = 23.4 - 1.15*log10(den_e) + 3.45*log10(1000.*temp_e)
      endif
      if(new_coll) then
       nu_i0 = ((qq*z)**4)*lnlambda*den_i*1.d+6/
     >  (forpi*(p_mass**2)*ep02*(vth_i**3))                !ion-ion
       nu_e0 = (qq**4)*lnlambda*den_e*1.d+6/
     >  (forpi*(e_mass**2)*ep02*(vth_e**3))       !elect-elect
      else if(.not.new_coll) then
       nu_i0 = 4.8e-8*(z**4)*lnlambda*den_i/
     1  (((1.e+3*temp_i)**1.5)*sqrt(mu))
       nu_e0 = 2.91e-6*den_e*lnlambda/((1.e+3*temp_e)**1.5)
      endif
      nu_i = (nu_i0/vth_i)/wtovth
      nu_e = (nu_e0/vth_e)/wtovth
c
c   Form parallel friction coefficient matrices based on S&N's eqn. (C8):
c
      zeff = 1.d0
      zee = 1.d0
      ll_e(1,1) = zee
      ll_e(1,2) = 1.5*zee
      ll_e(2,1) = ll_e(1,2)
      ll_e(2,2) = sqrt(2.d0) + 13.d0*zee/4.d0
      ll_i(2,2) = sqrt(2.d0)
      ll_i(1,1) = 0.d0; ll_i(1,2) = 0.d0; ll_i(2,1) = 0.d0
      lam_i(1:2,1:2) = ll_i(1:2,1:2)
      lam_e = ll_e
      lam_e(1,2) = -ll_e(1,2)
      lam_e(2,1) = -ll_e(1,2)
      lam_i = ll_i
c
c
      den_e_mks = 1.e+06*den_e
      den_i_mks = 1.e+06*den_i
c
c   Calculate collisional times based on formule in S&N following
c    eq. (25). Then form scalar coefficients that are used in eqs. (C8)-(C13):
c
      tau_ee = 3.*sqrt(pi)*pi*ep02*((mass_e*p_mass)**2)*(vth_e**3)
     >      /(den_e_mks*(qq**4)*lnlambda)
      tau_ii = 3.*sqrt(pi)*pi*ep02*((mass_i*p_mass)**2)*(vth_i**3)
     >      /(den_i_mks*(qq**4)*lnlambda)
      factor_m_e = tau_ee/(den_e_mks*mass_e*p_mass*bsq)
      factor_m_i = tau_ii/(den_i_mks*mass_i*p_mass*bsq)
      factor_n_e = tau_ee/(den_e_mks*mass_e*p_mass)
      factor_n_i = tau_ii/(den_i_mks*mass_i*p_mass)
      inv_factor_m_e = den_e_mks*mass_e*p_mass/(bsq*tau_ee)
      inv_factor_m_i = den_i_mks*mass_i*p_mass/(bsq*tau_ii)      
      factor_Ee1 = den_e_mks*qq/(sqrt(bsq))
      factor_Ee2 = den_e_mks*qq*qq*tau_ee/(mass_e*p_mass)
      if(plots) then
       write(*,'("tau_ee, tau_ii = ",e15.7,3x,e15.7)') tau_ee,tau_ii
       write(*,36) factor_m_e, factor_m_i, factor_n_e, factor_n_i,
     >  inv_factor_m_e, inv_factor_m_i, factor_Ee1, factor_Ee2
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
      m_er = 3
      ne_pts = 200
      iroot = 1
      allocate(er(ne_pts), stat=istat)         !Er field array
      allocate(er0(ne_pts), stat=istat)         !Er field array
      allocate(e_phi_kT(ne_pts), stat=istat)         !Er field array
      allocate(net_flux(ne_pts), stat=istat)   !gamma_i - gamma_e 
      allocate(gamma_e(ne_pts), stat=istat)        !elec particle flux
      allocate(gamma_i(ne_pts), stat=istat)        !ion particle flux
      
      allocate(gmi_e1(ne_pts), stat=istat)        !ion particle flux component
      allocate(gmi_e2(ne_pts), stat=istat)        !ion particle flux component
      allocate(gmi_i1(ne_pts), stat=istat)        !ion particle flux component
      allocate(gmi_i2(ne_pts), stat=istat)        !ion particle flux component
      allocate(gmi_i3(ne_pts), stat=istat)        !ion particle flux component
      allocate(gmi_Eprl(ne_pts), stat=istat)       !ion particle flux component
      
      allocate(jbs(ne_pts), stat=istat)        !elec heat flux
      allocate(q_e(ne_pts), stat=istat)        !elec heat flux
      allocate(q_i(ne_pts), stat=istat)        !ion heat flux 
      allocate(c_er(ne_pts), stat=istat)       !net_flux spline coeffs
      allocate(c_qe(ne_pts), stat=istat)       !qe spline coeffs 
      allocate(c_qi(ne_pts), stat=istat)       !qi spline coeffs 
      allocate(x_er(ne_pts+m_er), stat=istat)  !er spline knots 
      allocate(x_qe(ne_pts+m_er), stat=istat)  !qe spline knots 
      allocate(x_qi(ne_pts+m_er), stat=istat)  !qi spline knots 
      if(i_append .eq. 0) then
        write(*,10)
 10     format(" <r>/<a>","   Er(V/cm)",
     >   "        <a>Er/kTe","        upol_sug",
     >   "        utor_sug")
      endif
c
c   Loop over electic fields for ambipolar solver. There are 3 different
c    electric field representations used here: er = E_r in Volts/cm,
c    e_phi_kT = e<a>E_r/kT_e, and lambda = abs(e_phi_kT) = input variable for
c    lijs calculation.
c
      do i=1,ne_pts
       e_phi_kT(i) = ermin + (ermax-ermin)*real(i-1)/dble(ne_pts-1)    !e*phi/kT_elec.
c     Ion parameters:
       mass_i = 1.d0
       trat_i = temp_e/temp_i
       masrat_i = 1.d0/1837.d0
       zeff = 1.d0
c      Electron parameters:
       mass_e = 1.d0/1837.d0
       trat_e = temp_i/temp_e
       masrat_e = 1837.d0
       lambda = abs(e_phi_kT(i))
       pol_tor_cmpt = .false.
c
c   Sugama_app_c provides transport coefficient matrices.
c    lmat = transport matrix given in eqn. (C5) of S&N,
c    fl_mat = ion matrix needed for getting u_parallel and
c    q_parallel based on neglecting E|| and ion-electron
c    coupling effects.
c
       call sugama_app_c_matrix(lmat,fl_mat,lijs_inp_ext,lambda)
c  Form thermodynamic forces from S&N eq. (10):
       xe1 = -temp_e*k_boltz*(dn_dr_on + dTe_dr_oTe
     >       + e_phi_kT(i)/arad)
      if(exact_sym) then
       xe1 = -temp_e*k_boltz*(dn_dr_on + dTe_dr_oTe)
      end if
       xe2 = -temp_e*k_boltz*dTe_dr_oTe
       xe = 0.
c       if(beam_on) xe = beam_source*(1.-r*r)**2
       if(beam_on) xe = beam_source*4.*(r**1.2)*(1.-r*r)**2
       xi1 = -temp_i*k_boltz*(dn_dr_on + dTi_dr_oTi
     >       - e_phi_kT(i)*trat_i/arad)
      if(exact_sym) then
       xi1 = -temp_i*k_boltz*(dn_dr_on + dTi_dr_oTi)
      end if
       xi2 = -temp_i*k_boltz*dTi_dr_oTi
c
c     Calculate electron, ion particle and heat fluxes (actually
c      Q_e and Q_i are the heat fluxes divided by the temperature)
c      based on S&N eq. (C5):
c
       gamma_e(i) = (lmat(1,1)*xe1 + lmat(1,2)*xe2 + lmat(1,3)*xi1
     >   + lmat(1,4)*xi2 + lmat(1,5)*xe)
       q_e(i) = (lmat(2,1)*xe1 + lmat(2,2)*xe2 + lmat(2,3)*xi1
     >   + lmat(2,4)*xi2 + lmat(2,5)*xe)
       gamma_i(i) = (lmat(3,1)*xe1 + lmat(3,2)*xe2 + lmat(3,3)*xi1
     >   + lmat(3,4)*xi2 + lmat(3,5)*xe)
     
       gmi_e1(i) = lmat(3,1)*xe1
       gmi_e2(i) = lmat(3,2)*xe2
       gmi_i1(i) = -lmat(3,3)*temp_i*k_boltz*dn_dr_on
       gmi_i2(i) = -(lmat(3,3)+lmat(3,4))
     >             *temp_i*k_boltz*dTi_dr_oTi
       gmi_i3(i) = k_boltz*lmat(3,3)*e_phi_kT(i)*temp_e/arad
       gmi_Eprl(i) = lmat(3,5)*xe
c       chk = gamma_i(i) - gmi_i1(i) - gmi_i2(i) - gmi_i3(i)
c     >       - (lmat(3,1)*xe1 + lmat(3,2)*xe2)
c       write(*,*) i, chk
       q_i(i) = (lmat(4,1)*xe1 + lmat(4,2)*xe2 + lmat(4,3)*xi1
     >   + lmat(4,4)*xi2 + lmat(4,5)*xe)
       jbs(i) = (lmat(5,1)*xe1 + lmat(5,2)*xe2 + lmat(5,3)*xi1
     >   + lmat(5,4)*xi2 + lmat(5,5)*xe)
       net_flux(i) = gamma_i(i) - gamma_e(i)
c
c     convert e_phi_kT to Er(volt/cm):
c
      er(i) = (temp_e*1000./(100.*arad))*e_phi_kT(i)
      er0(i) = er(i)
      write(20,'(f10.6,6(a1,e15.7))') r,tb,er(i),tb,e_phi_kT(i),tb,
     >   gamma_i(i),tb, gamma_e(i), tb, q_i(i), tb, q_e(i)
      end do
      close(unit=12)
c
c     Solve for self-consistent Er root by first spline fitting gamma_e,i vs.
c      Er and then calling zeroin.c
c
c   The following loop makes the first approximate pass at finding the
c    electric field root. It will identify only one root for further refinement.
c    Set i = ii to pick up the last root that is encountered (as er is varied from
c    ermin to ermax) or i = ne_pts - ii to pick up the first root encountered.
c
      do ii=2,ne_pts-1
       if(.not.search_down) i = ne_pts - ii
       if(search_down) i = ii
       if(net_flux(i)*net_flux(i+1) .le. 0.) iroot=i
      end do
      call dbsnak(ne_pts,er,m_er,x_er)
      call dbsnak(ne_pts,er,m_er,x_qe)
      call dbsnak(ne_pts,er,m_er,x_qi)
      call dbsint(ne_pts,er,net_flux,m_er,x_er,c_er)
      call dbsint(ne_pts,er,q_e,m_er,x_qe,c_qe)
      call dbsint(ne_pts,er,q_i,m_er,x_qi,c_qi)
      ermin = er(iroot-1); ermax = er(iroot+1); tol_er = 1.d-20
c      if(plots) write(*,*) iroot, ermin, ermax
      if(.not.er_fixed)
     >  er_root = zeroin(ermin,ermax,rad_flux,tol_er)
      if(plots) write(*,'("Er root = ",e15.7)') er_root
      if(er_fixed) then       
        er_root = fit0 + fit1*r + fit2*(r**2) + fit3*(r**3)
     >     + fit4*(r**4) + fit5*(r**5) + fit6*(r**6)
      endif
      gi_m_ge = dbsval(er_root, m_er, x_er, ne_pts, c_er)
      qflux_i = dbsval(er_root, m_er, x_qi, ne_pts, c_qi)
      qflux_e = dbsval(er_root, m_er, x_qe, ne_pts, c_qe)
c
c     Calculate self-consistent flow velocities
c
      lambda = abs(arad*er_root/(10.*temp_e))
      pol_tor_cmpt = .true.
      lijs_write = .false.
      call sugama_app_c_matrix(lmat,fl_mat,lijs_inp_ext,lambda)
      lijs_write = .false.
      e_phi_kt_root = arad*er_root/(10.*temp_e)
      xi1 = -temp_i*k_boltz*(dn_dr_on + dTi_dr_oTi
     >       - e_phi_kt_root*trat_i/arad)
      xi1e =  temp_i*k_boltz*e_phi_kt_root*trat_i/arad
      xi2 = -temp_i*k_boltz*dTi_dr_oTi
      xe = 0.
      if(beam_on) xe = beam_source*(1.-r*r)**2
      xe1 = -temp_e*k_boltz*(dn_dr_on + dTe_dr_oTe
     >       + e_phi_kt_root/arad)
      xe2 = -temp_e*k_boltz*dTe_dr_oTe
      if(exact_sym) then
       xi1 = -temp_i*k_boltz*(dn_dr_on + dTi_dr_oTi)
       xi1e = 0.d0
       xe1 = -temp_e*k_boltz*(dn_dr_on + dTe_dr_oTe)
      endif
      if(plots) write(*,38) xi1,xe1,xi2,xe2
 38   format("xi1 = ",e15.7,"  xe1 = ",e15.7,/,
     >   "xi2 = ",e15.7,"  xe2 = ",e15.7)

c      uprl_B = (fl_mat(1,1)*xi1 + fl_mat(1,2)*xi2)/bsq
c      uprl_B_Er_only = fl_mat(1,1)*xi1e/bsq
c      uprl0 = (fl_mat(1,1)*xi1 + fl_mat(1,2)*xi2)/sqrt(bsq)
      uprl_B = (fl_mat(1,1)*xi1 + fl_mat(1,2)*xi2)           !10/3/2007 corrections
      uprl_B_Er_only = fl_mat(1,1)*xi1e                      !10/3/2007 corrections
      uprl0 = (fl_mat(1,1)*xi1 + fl_mat(1,2)*xi2)/sqrt(bsq)  !10/3/2007 corrections

      jbs_e = (lmat(5,1)*xe1 + lmat(5,2)*xe2)
      jbs_i = (lmat(5,3)*xi1 + lmat(5,4)*xi2)
      jbs_e_prl = lmat(5,5)*xe
      jbs_tot = jbs_e + jbs_i + jbs_e_prl
      qflux_e =  (lmat(2,1)*xe1 + lmat(2,2)*xe2 + lmat(2,3)*xi1
     >   + lmat(2,4)*xi2 + lmat(2,5)*xe)
      qflux_i = (lmat(4,1)*xe1 + lmat(4,2)*xe2 + lmat(4,3)*xi1
     >   + lmat(4,4)*xi2 + lmat(4,5)*xe)
      gam_e = (lmat(1,1)*xe1 + lmat(1,2)*xe2 + lmat(1,3)*xi1
     >   + lmat(1,4)*xi2 + lmat(1,5)*xe)
      gam_i = (lmat(3,1)*xe1 + lmat(3,2)*xe2 + lmat(3,3)*xi1
     >   + lmat(3,4)*xi2 + lmat(3,5)*xe)
c      gi_m_ge = gi_m_ge/sqrt(gam_e**2 + gam_i**2)
      write(23,87) r, tb,gam_i,tb,gam_e,tb,qflux_i,tb,qflux_e,tb,
     >  lmat(3,1)*xe1,tb, lmat(3,2)*xe2,tb,
     >  lmat(3,3)*xi1,tb, lmat(3,4)*xi2,tb,
     >  lmat(3,5)*xe, tb,
     >  lmat(1,1)*xe1,tb,lmat(1,2)*xe2,tb,
     >  lmat(1,3)*xi1,tb,lmat(1,4)*xi2,tb,
     >  lmat(1,5)*xe
 87   format(e15.7,14(a1,e15.7))
      write(15,51) r, ra, xi1, uprl_B, bsq, arad, psip, chip
 51   format(e15.7,7(2x,e15.7))
      write(19,91) r, xi1, uprl_B, bsq, arad, vp, btheta, bzeta,
     >  psip, chip, xi1e
 91   format(e15.7,10(2x,e15.7))
      proj_mat(1,1) = 1.; proj_mat(2,1) = 1.
      proj_mat(1,2) = -bzeta/(bsq*chip*chrg_i)
      proj_mat(2,2) = btheta/(bsq*psip*chrg_i)

c      upol = forpi2*chip*((proj_mat(1,1)*uprl_B/bsq)
c     >  + proj_mat(1,2)*xi1)/vp
c      upol_Er_only = forpi2*chip*((proj_mat(1,1)
c     > *uprl_B_Er_only/bsq) + proj_mat(1,2)*xi1e)/vp
c      utor = forpi2*psip*((proj_mat(2,1)*uprl_B/bsq)
c     >  + proj_mat(2,2)*xi1)/vp
      upol = forpi2*chip*((proj_mat(1,1)*uprl_B/bsq)           !10/3/2007 corrections
     >  + proj_mat(1,2)*xi1)/vp                            !10/3/2007 corrections
      upol_Er_only = forpi2*chip*((proj_mat(1,1)           !10/3/2007 corrections
     > *uprl_B_Er_only)/bsq + proj_mat(1,2)*xi1e)/vp           !10/3/2007 corrections
      utor = forpi2*psip*((proj_mat(2,1)*uprl_B/bsq)           !10/3/2007 corrections
     >  + proj_mat(2,2)*xi1)/vp                            !10/3/2007 corrections

      if(plots) then
        write(*,'("upol = ",e15.7,
     >   3x,"utor = ",e15.7,"   uprl = ", e15.7)') upol, utor, uprl0
        write(*,'("gam_e = ",e15.7,2x,"gam_i = ",e15.7)')
     >    gam_e, gam_i
        write(*,'("jbs_e = ",e15.7,2x,"jbs_i = ",e15.7,/,
     >   "gam_e/n = ",e15.7,2x,"gam_i/n = ",e15.7,/,
     >   "q_e/nT = ",e15.7,2x,"q_i/nT = ",e15.7)')
     >   jbs_e, jbs_i, gam_e/(1.e+6*den_e),
     >   gam_i/(1.e+6*den_i),qflux_e/(1.e+6*den_e*temp_e),
     >   qflux_i/(1.e+6*den_i*temp_i)
      endif
      if(.not.plots) write(*,'(1x,f5.3,4(a1,e12.5))') r,tb,
     >    er_root,tb,e_phi_kt_root,tb,upol,tb,utor
      write(11,'(1x,f5.3,6(a1,e20.14))') r,tb,
     >    er_root,tb,e_phi_kt_root,tb,upol,tb,
     >    utor,tb,uprl0,tb,upol_Er_only
c
      if(i_append .eq. 0) then
       open(unit=38,file="ion_viscosities",status="unknown")
       open(unit=39,file="electron_viscosities",status="unknown")
       open(unit=28,file="ion_viscosities_Er0",status="unknown")
       open(unit=29,file="electron_viscosities_Er0",status="unknown")
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
c   Calculate the ion and eletron poloidal/toroidal viscosity
c    components, as defined in equation (B4) of Sugama, et al.
c    These are written out to files as a function of radial
c    position (sqrt(flux)). The three velocity moments for
c    each species are given. These viscosities can then be
c    used in Sugama's eqn. (38) to obtain the stress tensor
c    components.
c
      nu_i0 = nu_i0*3.733/(vth_i*abs(iota))   !nu*i specific R0 and iota
      nu_e0 = nu_e0*3.733/(vth_e**abs(iota))   !nu*i specific R0 and iota
      if(i_append .eq. 0) then
      write(38,'("r",a1,"nu_i",a1,"er_root",a1,"mpp1i",
     >  a1,"mpp2i",a1,"mpp3i",a1,"mpt1i",
     >  a1,"mpt2i",a1,"mpt3i",a1,"mtt1i",a1,"mtt2i",a1,"mtt3i")')
     >   tb,tb,tb,tb,tb,tb,tb,tb,tb, tb, tb
      write(39,'("r",a1,"nu_i",a1,"er_root",a1,"mpp1e",
     >  a1,"mpp2e",a1,"mpp3e",a1,"mpt1e",
     >  a1,"mpt2e",a1,"mpt3e",a1,"mtt1e",a1,"mtt2e",a1,"mtt3e")')
     >   tb,tb,tb,tb,tb,tb,tb,tb,tb, tb, tb
      end if
      inter_mat(:,:) = 0.d0
      inter_mat = (forpi2/vp)*matmul(mnl1i,post_mat)
      mpt1i = matmul(pre_mat,inter_mat)
      inter_mat(:,:) = 0.d0
      inter_mat = (forpi2/vp)*matmul(mnl2i,post_mat)
      mpt2i = matmul(pre_mat,inter_mat)
      inter_mat(:,:) = 0.d0
      inter_mat = (forpi2/vp)*matmul(mnl3i,post_mat)
      mpt3i = matmul(pre_mat,inter_mat)
      inter_mat(:,:) = 0.d0
      inter_mat = (forpi2/vp)*matmul(mnl1e,post_mat)
      mpt1e = matmul(pre_mat,inter_mat)
      inter_mat(:,:) = 0.d0
      inter_mat = (forpi2/vp)*matmul(mnl2e,post_mat)
      mpt2e = matmul(pre_mat,inter_mat)
      inter_mat(:,:) = 0.d0
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
c
c
      lambda = 1.d-4
      pol_tor_cmpt = .true.
      call sugama_app_c_matrix(lmat,fl_mat,lijs_inp_ext,lambda)
      e_phi_kt_root = 0.d0
      xi1 = -temp_i*k_boltz*(dn_dr_on + dTi_dr_oTi
     >       - e_phi_kt_root*trat_i/arad)
      xi1e =  temp_i*k_boltz*e_phi_kt_root*trat_i/arad
      if(exact_sym) then
       xi1 = -temp_i*k_boltz*(dn_dr_on + dTi_dr_oTi)
       xi1e =  0.d0
      endif
      xi2 = -temp_i*k_boltz*dTi_dr_oTi
      xe = 0.
      if(beam_on) xe = beam_source*(1.-r*r)**2
      xe1 = -temp_e*k_boltz*(dn_dr_on + dTe_dr_oTe
     >       + e_phi_kt_root/arad)
      xe2 = -temp_e*k_boltz*dTe_dr_oTe
      jbs_e0 = (lmat(5,1)*xe1 + lmat(5,2)*xe2)
      jbs_i0 = (lmat(5,3)*xi1 + lmat(5,4)*xi2)
      jbs_e_prl0 = lmat(5,5)*xe
      jbs_tot0 = jbs_e + jbs_i + jbs_e_prl
      write(18,'(f10.6,14(a1,e15.7))') r,tb,den_e_mks,tb,temp_e,tb,
     >  temp_i,tb,jbs_e,tb,jbs_i,tb,jbs_e_prl,tb, jbs_tot, tb,
     >  jbs_e0,tb,jbs_i0,tb,jbs_e_prl0,tb, jbs_tot0, tb,
     >  gi_m_ge, tb, qflux_e, tb, qflux_i
c
c   As done above, calculate the ion and eletron poloidal/toroidal
c    viscosity components, as defined in equation (B4) of Sugama, et al.
c    for E_r ~ 0.
c
      if(i_append .eq. 0) then
      write(28,'("r",a1,"nu_i",a1,"er_root",a1,"mpp1i",a1,
     >  "mpp2i",a1,"mpp3i",a1,"mpt1i",
     >  a1,"mpt2i",a1,"mpt3i",a1,"mtt1i",a1,"mtt2i",a1,"mtt3i")')
     >   tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb
      write(29,'("r",a1,"nu_e",a1,"er_root",a1,"mpp1e",a1,
     >  "mpp2e",a1,"mpp3e",a1,"mpt1e",
     >  a1,"mpt2e",a1,"mpt3e",a1,"mtt1e",a1,"mtt2e",a1,"mtt3e")')
     >   tb,tb,tb,tb,tb,tb,tb,tb,tb,tb,tb
      end if
      inter_mat(:,:) = 0.d0
      inter_mat = (forpi2/vp)*matmul(mnl1i,post_mat)
      mpt1i = matmul(pre_mat,inter_mat)
      inter_mat(:,:) = 0.d0
      inter_mat = (forpi2/vp)*matmul(mnl2i,post_mat)
      mpt2i = matmul(pre_mat,inter_mat)
      inter_mat(:,:) = 0.d0
      inter_mat = (forpi2/vp)*matmul(mnl3i,post_mat)
      mpt3i = matmul(pre_mat,inter_mat)
      inter_mat(:,:) = 0.d0
      inter_mat = (forpi2/vp)*matmul(mnl1e,post_mat)
      mpt1e = matmul(pre_mat,inter_mat)
      inter_mat(:,:) = 0.d0
      inter_mat = (forpi2/vp)*matmul(mnl2e,post_mat)
      mpt2e = matmul(pre_mat,inter_mat)
      inter_mat(:,:) = 0.d0
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
c
c$$$      if(plots) then
c$$$      call plscolbg(255,255,255);
c$$$      call plcol(15)
c$$$      xmin = 1.2*minval(er)
c$$$      xmax = 1.2*maxval(er)
c$$$      ymin = min(minval(gamma_e),minval(gamma_i))
c$$$      ymax = max(maxval(gamma_e),maxval(gamma_i))
c$$$c      write(*,*) xmin, xmax, ymin, ymax
c$$$      call plenv(xmin,xmax,ymin,ymax,0,0)
c$$$      call pllab('E#dr#u (Volt/cm)','#gG#di#u, #gG#de',' ')
c$$$      call plcol(9)
c$$$      call plmtex('t',-1.d0,0.d0,0,'electron flux')
c$$$      call plline(ne_pts,er0,gamma_e)
c$$$      call plpoin(ne_pts,er0,gamma_e,17)
c$$$      call plcol(1)
c$$$      call plline(ne_pts,er0,gamma_i)
c$$$      call plpoin(ne_pts,er0,gamma_i,4)
c$$$      call plmtex('t',-3.d0,0.d0,0,'ion flux')
c$$$      call plend
c$$$      open(unit=44,file="flux_vs_Efield",status="unknown")
c$$$         write(44,'("*",/,"Er(V/cm)",a1,"Gamma_i",a1,"Gamma_e",a1,
c$$$     >      "Q_i",a1,"Q_e",a1,"Gmi_e1",a1,"Gmi_e2",a1,"Gmi_i1",
c$$$     >      a1,"Gmi_i2",a1,"Gmi_i3",a1,"Gmi_Eprl")') tb,tb,tb,tb,
c$$$     >      tb,tb,tb,tb,tb,tb
c$$$      do i=1,ne_pts
c$$$       write(44,88) er0(i), tb,gamma_i(i),tb,gamma_e(i),tb,q_i(i),
c$$$     >  tb,q_e(i),tb,gmi_e1(i), tb, gmi_e2(i), tb, gmi_i1(i),
c$$$     >  tb, gmi_i2(i), tb, gmi_i3(i), tb, gmi_Eprl(i)
c$$$      end do
c$$$ 88   format(e15.7,10(a1,e15.7))
c$$$      close(unit=44)
c$$$      endif
	if (idebug) then
		write(*,*) "  gam_e  ","       u_prl ","              qe "
     > ,"             qi", "          jbs"
	write(*,'(1x,e12.5,4(a1,e12.5))') ,		
     >    gam_e,tb,uprl_B,tb,qflux_e,tb,qflux_i,tb,jbs_tot
	endif
      close(unit=11)
      close(unit=15)
      close(unit=18)
      close(unit=19)
      close(unit=20)
c      stop
      end
c
c
c
      subroutine sugama_app_c_matrix(lm, fm, file_ext,
     >    lambda)
c
c   This subroutine calculates the transport coefficient matrix
c    given in S&N eq. (C5). It also provides matrices for poloidal
c    and toroidal viscosity coefficients and for obtaining parallel
c    flow and heat flow.
c
      use lijs_transfer
      use l_components
      use viscosity_pol_tor
      implicit none
      logical :: test_nw_deriv, dkes_limit
      real*8, dimension(2,2) :: l2i, m2i, n2i, fm,
     >  l2e, m2e, n2e, temp_1_e, temp_1_i, m2e_inv, lam_e_inv,
     >  temp_2_e, temp_2_i, temp_2_e_inv,temp_e_inv, temp_i_inv,
     >  tmp_e, tmp_i, temp_ie, temp_Ee, temp_Ei, temp_EEE, ep11
      real*8, dimension(2,2) :: tmp1, tmp2, tmp3, tmp4,
     >  tmp5, tmp6
      real*8, dimension(5,5) :: lm
      real*8 :: lambda, cl1e, cl2e, cl3e, cm1e, cm2e,
     >  cm3e, cn1e, cn2e, cn3e, lambda_i,
     >  cl1i, cl2i, cl3i, cm1i, cm2i, cm3i, cn1i, cn2i, cn3i
      integer :: iterate, i, j
      character*50 :: file_ext
      character*60 :: file_nm
      npts = 2
      klm = 1
      iterate = 1
c
c     Integrate electron coefficients:
c
c   S&N eq. (C8) rotation matrix:
      ep11(1,1)=1.d0;ep11(1,2)=0.d0;ep11(2,1)=0.d0;ep11(2,2)=0.d0
c
      mass = mass_e
      temp = temp_e
      trat = trat_e
      masrat = masrat_e
      nulo = nu_e
      nuhi = 1.1*nu_e
c      write(*,99) nulo,nuhi, nu_e
c  99  format(1x,"electrons: nulo = ",e15.7,4x,"nuhi = ",e15.7,
c     >   4x, "nu_e = ",e15.7)
        lcoef_opt = 5
        file_nm = "lstar_lijs_" // file_ext
        vs1 = 1.d0; vs2 = 0.d0; vs3 = 0.d0
      call dkes_coef(cl1e,file_nm,iterate,lambda)
        vs1 = 0.d0; vs2 = 1.d0; vs3 = 0.d0
      call dkes_coef(cl2e,file_nm,iterate,lambda)
        vs1 = 0.d0; vs2 = 0.d0; vs3 = 1.d0
      call dkes_coef(cl3e,file_nm,iterate,lambda)
        lcoef_opt = 4
        file_nm = "mstar_lijs_" // file_ext
        vs1 = 1.d0; vs2 = 0.d0; vs3 = 0.d0
      call dkes_coef(cm1e,file_nm,iterate,lambda)
        vs1 = 0.d0; vs2 = 1.d0; vs3 = 0.d0
      call dkes_coef(cm2e,file_nm,iterate,lambda)
        vs1 = 0.d0; vs2 = 0.d0; vs3 = 1.d0
      call dkes_coef(cm3e,file_nm,iterate,lambda)
        lcoef_opt = 5
        file_nm = "nstar_lijs_" // file_ext
        vs1 = 1.d0; vs2 = 0.d0; vs3 = 0.d0
      call dkes_coef(cn1e,file_nm,iterate,lambda)
        vs1 = 0.d0; vs2 = 1.d0; vs3 = 0.d0
      call dkes_coef(cn2e,file_nm,iterate,lambda)
        vs1 = 0.d0; vs2 = 0.d0; vs3 = 1.d0
      call dkes_coef(cn3e,file_nm,iterate,lambda)
c
c     Integrate Ion coefficients:
c
      mass = mass_i
      temp = temp_i
      trat = trat_i
      masrat = masrat_i
      nulo = nu_i
      nuhi = 1.1*nu_i
c  98  format(1x,"ions: nulo = ",e15.7,4x,"nuhi = ",e15.7,
c     >   4x, "nu_i = ",e15.7)
        lambda_i = lambda*temp_e/temp_i
        lcoef_opt = 5
        file_nm = "lstar_lijs_" // file_ext
        vs1 = 1.d0; vs2 = 0.d0; vs3 = 0.d0
      call dkes_coef(cl1i,file_nm,iterate,lambda_i)
        vs1 = 0.d0; vs2 = 1.d0; vs3 = 0.d0
      call dkes_coef(cl2i,file_nm,iterate,lambda_i)
        vs1 = 0.d0; vs2 = 0.d0; vs3 = 1.d0
      call dkes_coef(cl3i,file_nm,iterate,lambda_i)
        lcoef_opt = 4
        file_nm = "mstar_lijs_" // file_ext
        vs1 = 1.d0; vs2 = 0.d0; vs3 = 0.d0
      call dkes_coef(cm1i,file_nm,iterate,lambda_i)
        vs1 = 0.d0; vs2 = 1.d0; vs3 = 0.d0
      call dkes_coef(cm2i,file_nm,iterate,lambda_i)
        vs1 = 0.d0; vs2 = 0.d0; vs3 = 1.d0
      call dkes_coef(cm3i,file_nm,iterate,lambda_i)
        lcoef_opt = 5
        file_nm = "nstar_lijs_" // file_ext
        vs1 = 1.d0; vs2 = 0.d0; vs3 = 0.d0
      call dkes_coef(cn1i,file_nm,iterate,lambda_i)
        vs1 = 0.d0; vs2 = 1.d0; vs3 = 0.d0
      call dkes_coef(cn2i,file_nm,iterate,lambda_i)
        vs1 = 0.d0; vs2 = 0.d0; vs3 = 1.d0
      call dkes_coef(cn3i,file_nm,iterate,lambda_i)
c
c     Form 2x2 matrices:
c
      n2e(1,1) = -den_e_mks*cn1e; n2e(1,2) = -den_e_mks*cn2e
      n2e(2,1) = -den_e_mks*cn2e; n2e(2,2) = -den_e_mks*cn3e
      
      n2i(1,1) = den_i_mks*cn1i; n2i(1,2) = den_i_mks*cn2i
      n2i(2,1) = den_i_mks*cn2i; n2i(2,2) = den_i_mks*cn3i
      
      m2e(1,1) = den_e_mks*cm1e; m2e(1,2) = den_e_mks*cm2e
      m2e(2,1) = den_e_mks*cm2e; m2e(2,2) = den_e_mks*cm3e
      
      m2i(1,1) = den_i_mks*cm1i; m2i(1,2) = den_i_mks*cm2i
      m2i(2,1) = den_i_mks*cm2i; m2i(2,2) = den_i_mks*cm3i
      
      l2e(1,1) = den_e_mks*cl1e; l2e(1,2) = den_e_mks*cl2e
      l2e(2,1) = den_e_mks*cl2e; l2e(2,2) = den_e_mks*cl3e
      
      l2i(1,1) = den_i_mks*cl1i; l2i(1,2) = den_i_mks*cl2i
      l2i(2,1) = den_i_mks*cl2i; l2i(2,2) = den_i_mks*cl3i
c
c   if(pol_tor_cmpt) save components needed for calculating poloidal/toroidal
c     viscosities
      if(pol_tor_cmpt) then
      
       mnl1i(1,1)=den_i_mks*cm1i;mnl1i(1,2)=den_i_mks*cn1i
       mnl1i(2,1)=den_i_mks*cn1i;mnl1i(2,2)=den_i_mks*cl1i
       
       mnl1e(1,1)=den_e_mks*cm1e;mnl1e(1,2)=den_e_mks*cn1e
       mnl1e(2,1)=den_e_mks*cn1e;mnl1e(2,2)=den_e_mks*cl1e
       
       mnl2i(1,1)=den_i_mks*cm2i;mnl2i(1,2)=den_i_mks*cn2i
       mnl2i(2,1)=den_i_mks*cn2i;mnl2i(2,2)=den_i_mks*cl2i
       
       mnl2e(1,1)=den_e_mks*cm2e;mnl2e(1,2)=den_e_mks*cn2e
       mnl2e(2,1)=den_e_mks*cn2e;mnl2e(2,2)=den_e_mks*cl2e
       
       mnl3i(1,1)=den_i_mks*cm3i;mnl3i(1,2)=den_i_mks*cn3i
       mnl3i(2,1)=den_i_mks*cn3i;mnl3i(2,2)=den_i_mks*cl3i
       
       mnl3e(1,1)=den_e_mks*cm3e;mnl3e(1,2)=den_e_mks*cn3e
       mnl3e(2,1)=den_e_mks*cn3e;mnl3e(2,2)=den_e_mks*cl3e
       
      end if     !if(pol_tor_cmpt)
c
c   Calculate diagonal 2x2 blocks:
c
c   factor_m_e,i = tau/(den*mass*<B**2>)
c   factor_n_e,i = tau/(den*mass)
c   inv_factor_m_e,i = den*mass/(tau*<B**2>)
c   factor_Ee1 = den*e/sqrt(<B**2>)
c   factor_Ee2 = den*(e**2)*tau_ee/mass_e
c
      do i=1,2
       do j=1,2
        m2e(i,j) = factor_m_e*m2e(i,j)
        m2i(i,j) = factor_m_i*m2i(i,j)
        n2e(i,j) = -factor_n_e*n2e(i,j)
        n2i(i,j) = -factor_n_i*n2i(i,j)
c	l2e(i,j) = l2e(i,j)
c	l2i(i,j) = l2i(i,j)
	temp_1_e(i,j) = m2e(i,j) + lam_e(i,j)
	temp_1_i(i,j) = m2i(i,j) + lam_i(i,j)
	tmp1(i,j) = 0.
	tmp2(i,j) = 0.
	tmp3(i,j) = 0.
	tmp4(i,j) = 0.
	tmp5(i,j) = 0.
	tmp6(i,j) = 0.
       end do
      end do
      if(lijs_write) then
        write(*,37) m2e(1,1),m2i(1,1),n2e(1,1),
     >   n2i(1,1),l2e(1,1),l2i(1,1),lam_e(1,1),lam_i(2,2)
      endif
  37  format("M_e_11 = ",e15.7,"  M_i_11 = ",e15.7,/,
     >  "N_e_11 = ",e15.7,"  N_i_11 = ",e15.7,/,
     >  "L_e_11 = ",e15.7,"  L_i_11 = ",e15.7,/
     >  "lam_e_11 = ",e15.7,"  lam_i_22 = ",e15.7)
      call mat_2by2_inverse(m2e,m2e_inv)
      call mat_2by2_inverse(lam_e,lam_e_inv)
      do i=1,2
       do j=1,2
        temp_2_e(i,j) = m2e_inv(i,j) + lam_e_inv(i,j)
       end do
      end do
      call mat_2by2_inverse(temp_2_e,temp_2_e_inv)
      call mat_2by2_inverse(temp_1_e, temp_e_inv)
      call mat_2by2_inverse(temp_1_i, temp_i_inv)
c
c      tmp_e = l2e - inv_factor_m_e*n2e*temp_e_inv*n2e
c      tmp_i = l2i - inv_factor_m_i*n2i*temp_i_inv*n2i
c     >  + inv_factor_m_e*n2i*temp_i_inv*ep11*temp_2_e_inv
c     >  *ep11*temp_i_inv*n2i
c
c   Calculate terms in S&N eq. (C9) for electrons:
c
      tmp1 = matmul(temp_e_inv,n2e)
      tmp2 = matmul(n2e,tmp1)
      do i=1,2
       do j=1,2
        tmp_e(i,j) = l2e(i,j) - inv_factor_m_e*tmp2(i,j)
       end do
      end do
c
c   Calculate terms in S&N eq. (C9) for ions:
c
      tmp1(1:2,1:2) = 0.; tmp2(1:2,1:2) = 0.
      tmp1 = matmul(temp_i_inv,n2i)
      tmp2 = matmul(n2i,tmp1)
      do i=1,2
       do j=1,2
        tmp_i(i,j) = l2i(i,j) - inv_factor_m_i*tmp2(i,j)
       end do
      end do
      if(lijs_write) then
        write(*,111) l2i(1,1), l2i(1,2), l2i(2,1), l2i(2,2),
     >   tmp2(1,1), tmp2(1,2), tmp2(2,1), tmp2(2,2)
 111  format("L_ion:   ",e15.7,3x,e15.7,/,9x,e15.7,3x,e15.7,//,
     >       "ion_part2: ",e15.7,3x,e15.7,/,11x,e15.7,3x,e15.7)
      endif
      
      test_nw_deriv = .false.
      dkes_limit = .false.
c
c   Form fm matrix. This is (M_i + Lambda_I)**-1 N_i, which is used
c    to get U|| and Q|| from eq. (C3) taking E|| = 0 and neglecting
c    electron-ion coupling.
c
      tmp1(1:2,1:2) = 0.; tmp2(1:2,1:2) = 0.
      tmp1 = matmul(temp_i_inv,n2i)
      fm(1:2,1:2) = tmp1(1:2,1:2)
      if(test_nw_deriv .eq. .false.) then
c
c   Form the third term in eq. (C9) for ions based on S&N:
c
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
 114  format(//,a4": ",e15.7,3x,e15.7,/,6x,e15.7,3x,e15.7)
      end if
 
      if(lijs_write) then
        write(*,112) tmp6(1,1), tmp6(1,2), tmp6(2,1), tmp6(2,2)
 112  format(/,"ion_part3: ",e15.7,3x,e15.7,/,11x,e15.7,3x,e15.7,/)
      end if
      else if(test_nw_deriv .eq. .true.) then
c
c   Form the third term in eq. (C9) for ions based on my derivation:
c
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
       tmp1(1:2,1:2) = 0.; tmp2(1:2,1:2) = 0.
       tmp2 = matmul(temp_i_inv,tmp6)
       tmp1 = matmul(n2i,tmp2)
       if(lijs_write) then
        write(*,112) tmp1(1,1), tmp1(1,2), tmp1(2,1), tmp1(2,2)
       end if
      endif
c
      do i=1,2
       do j=1,2
      if(test_nw_deriv .eq. .false.) then
        tmp_i(i,j) = tmp_i(i,j) + inv_factor_m_e*tmp6(i,j)
      else if(test_nw_deriv .eq. .true.) then
        tmp_i(i,j) = tmp_i(i,j) - inv_factor_m_e*tmp1(i,j)
      end if
       end do
      end do
c   Save eq. (C9) ii, ee transport coefficients:
      if(.not.dkes_limit) then
       lm(1,1) = tmp_e(1,1); lm(1,2) = tmp_e(1,2)
       lm(2,1) = tmp_e(2,1); lm(2,2) = tmp_e(2,2)
       lm(3,3) = tmp_i(1,1); lm(3,4) = tmp_i(1,2)
       lm(4,3) = tmp_i(2,1); lm(4,4) = tmp_i(2,2)
      endif
c  Test: use just Li and Le terms:
      if(dkes_limit) then
       lm(1,1) = l2e(1,1); lm(1,2) = l2e(1,2)
       lm(2,1) = l2e(2,1); lm(2,2) = l2e(2,2)
       lm(3,3) = l2i(1,1); lm(3,4) = l2i(1,2)
       lm(4,3) = l2i(2,1); lm(4,4) = l2i(2,2)
      end if
c
c     Calculate off-diagonal e-i blocks using S&N eq. (C10):
c
c      temp_ie = - inv_factor_m_e*n2e*temp_e_inv*lam_e*ep11
c     >  *temp_i_inv*n2i
      tmp1(1:2,1:2) = 0.; tmp2(1:2,1:2) = 0.; tmp3(1:2,1:2) = 0.
      tmp4(1:2,1:2) = 0.; tmp5(1:2,1:2) = 0.; tmp6(1:2,1:2) = 0.
      tmp1 = matmul(temp_i_inv,n2i)
      tmp2 = matmul(ep11,tmp1)
      tmp3 = matmul(lam_e,tmp2)
      tmp4 = matmul(temp_e_inv,tmp3)
      tmp5 = matmul(n2e,tmp4)
      do i=1,2
       do j=1,2
        temp_ie(i,j) = - inv_factor_m_e*tmp5(i,j)
       end do
      end do
      if(.not.dkes_limit) then
c   eq. (C10) ei components = temp_ie:
       lm(1,3) = temp_ie(1,1); lm(1,4) = temp_ie(1,2)
       lm(2,3) = temp_ie(2,1); lm(2,4) = temp_ie(2,2)
c   eq. (C10) ie components = transpose(temp_ie):
       lm(3,1) = temp_ie(1,1); lm(3,2) = temp_ie(2,1)
       lm(4,1) = temp_ie(1,2); lm(4,2) = temp_ie(2,2)
      else if(dkes_limit) then
c     for conventional limit - no electron/ion coupling:
       lm(1,3) = 0.d0; lm(1,4) = 0.d0
       lm(2,3) = 0.d0; lm(2,4) = 0.d0
       lm(3,1) = 0.d0; lm(3,2) = 0.d0
       lm(4,1) = 0.d0; lm(4,2) = 0.d0
      endif
c
c     Calculate bootstrap and electric field terms based on eqns. (C11)-(C13):
c      temp_Ee = factor_Ee1*temp_e_inv*n2e      
c      
c   First eq. (C11):
      tmp1(1:2,1:2) = 0.; tmp2(1:2,1:2) = 0.; tmp3(1:2,1:2) = 0.
      tmp1 = matmul(temp_e_inv,n2e)
      do i=1,2
       do j=1,2
        temp_Ee(i,j) = factor_Ee1*tmp1(i,j)
       end do
      end do
      lm(5,1) = temp_Ee(1,1); lm(5,2) = temp_Ee(1,2)
      lm(1,5) = -temp_Ee(1,1); lm(2,5) = -temp_Ee(1,2)
      if(beam) then
       lm(1,5) = temp_Ee(1,1); lm(2,5) = temp_Ee(1,2)
      endif
c
c   Next eq. (C12):        
c     temp_Ei = -factor_Ee1*temp_e_inv*m2e*ep11*temp_i_inv*n2i      
c
      tmp1(1:2,1:2) = 0.; tmp2(1:2,1:2) = 0.; tmp3(1:2,1:2) = 0.
      tmp4(1:2,1:2) = 0.
      tmp1 = matmul(temp_i_inv,n2i)
      tmp2 = matmul(ep11,tmp1)
      tmp3 = matmul(m2e,tmp2)
      tmp4 = matmul(temp_e_inv,tmp3)
      do i=1,2
       do j=1,2
        temp_Ei(i,j) = -factor_Ee1*tmp4(i,j)
       end do
      end do
      lm(5,3) = temp_Ei(1,1); lm(5,4) = temp_Ei(1,2)
      lm(3,5) = -temp_Ei(1,1); lm(4,5) = -temp_Ei(1,2)
c
c   Finally, eq. (C13):
c      
      do i=1,2
       do j=1,2
        temp_EEE(i,j) = -factor_Ee2*(lam_e_inv(i,j) - temp_e_inv(i,j))
       end do
      end do
      lm(5,5) = temp_EEE(1,1)
      
      return
      end
c
c
c
      subroutine mat_2by2_inverse(a,a_inverse)
      real*8, DIMENSION(2,2):: a, a_inverse
      real*8 :: denom
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
c
      real*8 function rad_flux(xx_er)
      use er_solve
      use bspline
      implicit none
      integer i
      real*8 :: xx_er, yval
      rad_flux = dbsval(xx_er, m_er, x_er, ne_pts, c_er)
      return
      end function rad_flux
c
