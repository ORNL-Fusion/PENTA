!
!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a BETA version of the program dkes_data_prep1.f, which is currently
!       under development by D. A. Spong at the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a BETA version, this program is subject to change
!       and improvement without notice.
!
!
!      This code reads the temp.ext files produced by running DKES
!       through the xrun shell and processes their data to produce
!       tables of transport and viscosity coefficients
!       as a function of collisionality (visc.ext files). The
!       viscosities are calculated using the analysis of
!       Sugama, et al., Phys. Plasmas vol. 9 (2002), pg. 4637]
!
!       compile with:
!        lf95 -o xprep1 dkes_data_prep1.f
!
!    Command line arguments:
!    xprep1 temp_file.ext profile_data_file surf_inner surf_outer delta_surf ftrapped_file(optional)
!         where temp_file.ext is just the filename extension of the temp.*
!         file produced by xrun that you want to process.
!        profile_data_file = name of profile file for this case (this has been made
!         prior to running xprep1 by running xpro (from profile_extractor.f) on the
!         wout file that was used in the DKES runs
!
!        surf_inner = innermost surface used
!        delta_surf = increment in surface number
!        surf_outer = outermost surface used
!
        implicit none
        Integer, parameter :: rknd = selected_real_kind(12,300) 
        INTEGER, PARAMETER :: mxpts=1000
	LOGICAL, PARAMETER :: itok = .TRUE.
	LOGICAL, PARAMETER :: Bscale_correction = .TRUE.
	REAL*8, PARAMETER :: e_ovr_c = 1.6e-19      !MKS units
	REAL*8, PARAMETER :: four_pi_sq = 39.47841762
	REAL*8, DIMENSION(mxpts):: cmul, efield, L11, L31, L33, chip,
     >    psip, btheta, bzeta, mstar, lstar, nstar, mpp, mpt, mtt,
     >    mpp_0, mpt_0, mtt_0, vp, lstar_0, nstar_0, lstar_0_E0, alpha
        INTEGER, DIMENSION(mxpts):: itype
        REAL*8, DIMENSION(2,2):: amat, bmat, cmat, tmp, visc
	REAL*8 :: ftrapped(1000)
	REAL*8 :: efield_lore(100)
	REAL(rknd) :: L11_lore(100,20), L31_lore(100,20),
     >                L33_lore(100,20)
	REAL*8, DIMENSION(:), ALLOCATABLE :: chip_s,
     >   psip_s, btheta_s, bzeta_s, vp_s, bsq_avg_s, iota_s
        REAL*8 :: weov,wtov,L11m,L11p,L31m,L31p,L33m,L33p,scal11,
     >    scal13,scal33,rsds_max,denom,coef,c1,c2,c3,c4,c5,c6,
     >    utilde_sq, bsq_flx_avg, merge_factor, reduct, extrap_factor,
     >    dum1, dum2, dum3, dum4, dum5, dum6, dum7, dum8, dum9, dum10,
     >    dum11, merge_factor1, l31_lownu_E0, aminor, rr, ra,s0,skip1,
     >    skip2,skip3,skip4,skip5,skip6,skip7,denom_asympt,arg_tanh,
     >    l33_sptzr, Rmajor, b0
        INTEGER :: i,ipts,nunin=15,nunout=16, num_extrap_pts,
     >     nun_merge_file=17, num_E0_merge_pts, iloc1, iloc2, ii,
     >     surf0, outer_surf, delta_surf, isurf, ie, ierr,
     >     start_surf, end_surf, numargs, iargc, nc, ne, iout,
     >     ioffset
        CHARACTER*1 :: chardum
        CHARACTER*1 :: tb, t0
        CHARACTER*60 :: file_in, file_out, file_nm, utilde_file
        CHARACTER*60 :: file_out_lore1, file_out_lore2, file_out_lore3
	CHARACTER*60 :: arg1, arg2, arg3, arg4, arg5, arg6, skp_line
	CHARACTER*20 :: dr, srf
	CHARACTER*10 :: dr_num
	LOGICAL :: merge_hi_nu, extrap_lo_nu, l31_extrap_E0
	LOGICAL :: dkes_limit
c
c	dkes_limit is an option that allows the more conventional
c       (non-momentum conserving) form of the particle fluxes to
c       be tested. Setting dkes_limit = .true. reverts to this
c       limit. This variable must also separately be set equal to .true.
c       within the PENTA code in order to go to this limit in the 
c       fluxes. Both codes must be recompiled and re-run after changing
c       dkes_limit. Setting dkes_limit = .false. retains all of the
c       terms from the Sugama analysis.
c
	dkes_limit = .false.
	numargs = iargc()
        CALL getarg(1, arg1)
        CALL getarg(2, arg2)
        CALL getarg(3, arg3)
        CALL getarg(4, arg4)
        CALL getarg(5, arg5)
        if(numargs .eq. 6) CALL getarg(6, arg6)
	read(arg3,'(i4)') surf0
	read(arg4,'(i4)') outer_surf
	read(arg5,'(i4)') delta_surf
	l31_extrap_E0 = .FALSE.
c
c
      utilde_file = "utilde_vs_s_" // trim(adjustl(arg1))
      open(unit=6,file=trim(adjustl(arg2)),status="unknown")
      open(unit=22,file=utilde_file,status="unknown")
      read(6,*) start_surf, end_surf
      read(6,*) aminor, Rmajor
      read(6,*) chardum
c
       allocate (chip_s(end_surf), stat=ierr)
       allocate (psip_s(end_surf), stat=ierr)
       allocate (btheta_s(end_surf), stat=ierr)
       allocate (bzeta_s(end_surf), stat=ierr)
       allocate (vp_s(end_surf), stat=ierr)
       allocate (bsq_avg_s(end_surf), stat=ierr)
       allocate (iota_s(end_surf), stat=ierr)
c
      do i=start_surf,end_surf
       read(6,33) ii,tb,rr,tb,ra,tb,chip_s(i),tb,
     >   psip_s(i),tb,btheta_s(i),tb,bzeta_s(i),tb,vp_s(i),
     >   tb,bsq_avg_s(i),tb,iota_s(i), b0
      end do
  33  format(i4,10(a1,e15.7))
      close(6)
c
      if (numargs .eq. 6) then
       open(unit=33,file=arg6,status="old")
        read(33,*) skp_line
        do i=1,1000
        read(33,*,END=51) s0,skip1,skip2,ftrapped(i),skip3,skip4,
     >     skip5,skip6,skip7
c        write(*,*) i,s0,ftrapped(i)
        end do
  51    continue
       close(unit=33)
       end if  !(numargs .eq. 6)
c	
c      Surface loop
c
       do isurf=surf0,outer_surf,delta_surf
       write(dr_num, '(i4)') isurf
       dr = "s" // trim(adjustl(dr_num)) // "/"
       srf = "s" // trim(adjustl(dr_num))
       file_nm = "temp." // trim(adjustl(arg1)) // "_" //
     >   trim(adjustl(srf)) // "_E"     
       file_out_lore1 = trim(adjustl(dr)) // "D11_star_" // 
     >   trim(adjustl(arg1)) // "_" // trim(adjustl(srf))
       file_out_lore2 = trim(adjustl(dr)) // "D13_star_" //
     >    trim(adjustl(arg1)) // "_" // trim(adjustl(srf))
       file_out_lore3 = trim(adjustl(dr)) // "D33_star_" //
     >    trim(adjustl(arg1)) // "_" // trim(adjustl(srf))
       open(unit=37,file=file_out_lore1,status="unknown")
       open(unit=38,file=file_out_lore2,status="unknown")
       open(unit=39,file=file_out_lore3,status="unknown")
       
       bsq_flx_avg = bsq_avg_s(isurf)
c
c      Electric field loop
c
        do ie = 1, 15
c
c       Initialize variables to zero:
c
c        do i=1,mxpts
c	 cmul(i)=0.;efield(i)=0.; L11(i)=0.; L31(i)=0.; L33(i)=0.
c         chip(i)=0.; psip(i)=0.; btheta(i)=0.; bzeta(i)=0.; mstar(i)=0.
c         lstar(i)=0.; nstar(i)=0.; mpp(i)=0.; mpt(i)=0.; mtt(i)=0.
c         mpp_0(i)=0.; mpt_0(i)=0.; mtt_0(i)=0.; vp(i)=0.; lstar_0(i)=0.
c         nstar_0(i)=0.; lstar_0_E0(i)=0.; alpha(i)=0.
c         itype(i) = 0
c	END DO
c         amat(1:2,1:2) = 0.; bmat(1:2,1:2) = 0.
c	 cmat(1:2,1:2) = 0.; tmp(1:2,1:2) = 0.; visc(1:2,1:2) = 0.
c
c
        if (ie .eq. 1) then
	  merge_hi_nu = .FALSE.
	  extrap_lo_nu = .TRUE.
	else if (ie .gt. 1) then
	  merge_hi_nu = .TRUE.
	  extrap_lo_nu = .TRUE.
	endif
	if(ie .eq. 1) then
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm))//"0.00000"
	else if(ie .eq. 2) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "1m10"
	else if(ie .eq. 3) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "1m9"
	else if(ie .eq. 4) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "1m8"
	else if(ie .eq. 5) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "1m7"
	else if(ie .eq. 6) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "1m6"
	else if(ie .eq. 7) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "1m5"
	else if(ie .eq. 8) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "0.00010"
	else if(ie .eq. 9) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "0.00030"
	else if(ie .eq. 10) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "0.00100"
	else if(ie .eq. 11) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "0.00300"
	else if(ie .eq. 12) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "0.01000"
	else if(ie .eq. 13) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "0.03000"
	else if(ie .eq. 14) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "0.10000"
	else if(ie .eq. 15) then 
	 file_in = trim(adjustl(dr)) //
     >    trim(adjustl(file_nm)) // "0.30000"
	endif
	iloc1 = index(file_in, "t")
	iloc2 = index(file_in, ".")
	file_out = trim(adjustl(file_in(1:iloc1-1))) //
     >   "visc." // trim(adjustl(file_in(iloc2+1:60)))
        open(unit=nunin,file=file_in,status="old")
        open(unit=nunout,file=file_out,status="unknown")
        tb = char(9)
        write(nunout,'("*",/,"cmul",a1,"efield",a1,"mpp",
     >   a1,"-mpt",a1,"mtt",a1,"mstar",a1,"lstar",
     >   a1,"nstar",a1,"L11",a1,"L31",a1,"L33")')
     >   tb,tb,tb,tb,tb,tb,tb,tb,tb,tb
	read(nunin,98) chardum
	read(nunin,98) chardum
       DO i=1,mxpts           ! collisionality loop
       READ(NUNIN,99,END=501) cmul(i),tb,efield(i),tb,weov,tb,wtov,
     >  tb,L11m,tb,L11p,tb,L31m,tb,L31p,tb,L33m,tb,L33p,tb,scal11,tb,
     >  scal13,tb,scal33,tb,rsds_max,tb,chip(i),tb,psip(i),tb,
     >  btheta(i),tb,bzeta(i),tb,vp(i)
       efield_lore(ie) = efield(i)
       L11(i) = 0.5*(L11m + L11p)
       L31(i) = 0.5*(L31m + L31p)
       L33(i) = 0.5*(L33m + L33p)
c
c   The following choices are useful for testing DKES convergence:
c       L11(i) = L11m
c       L31(i) = L31m
c       L33(i) = L33m
c       L11(i) = L11p
c       L31(i) = L31p
c       L33(i) = L33p
c       L33(i) = min(L33m,L33p)
c
       if(Bscale_correction) then
        L31(i) = 0.5*(L31m + L31p)*sqrt(bsq_flx_avg)
        L33(i) = 0.5*(L33m + L33p)*bsq_flx_avg
c        L31(i) = 0.5*(L31m + L31p)*sqrt(bsq_flx_avg)
c	 l33_sptzr = 2.d0/(3.d0*cmul(i))
c        L33(i) = (l33_sptzr - 0.5*(L33m + L33p))*bsq_flx_avg
       endif
       ipts=i
       END DO
       write(*,*) ' ipts=mxpts; may have missed some lines'
 501   continue
c
c
       if(ie .eq. 1) then
         utilde_sq = 1.5*L11(1)/cmul(1)
	 l31_lownu_E0 = L31(ipts)
c         write(*,'(i3,2x,"Utilde_sq = ",e12.5)') isurf, utilde_sq
	 write(22,'(i3,2x,e15.7)') isurf, utilde_sq
       end if
c
c      If extrap_lo_nu = .true., perform extrapolations of L11, L31, L33
c       down to extrap_factor*cmul(ipts) using the rules given below:
c
       if(extrap_lo_nu) then
        num_extrap_pts = 10                                                 !Case specific
        extrap_factor = 1.e-3                                               !Case specific
        reduct = exp((log(extrap_factor))/real(num_extrap_pts))
c      Define scaling type and nu-scaling of last calculated point:
       if(efield(1) .le. 0.1*cmul(ipts)) then
         itype(ipts) = 1
         alpha(ipts) = -1.
       else if(efield(1) .gt. 0.1*cmul(ipts) .and.
     >           efield(1) .le. 10.*cmul(ipts)) then
         itype(ipts) = 2
         alpha(ipts) = 1. - log10(10.*cmul(ipts)/(efield(1)))
       else if(efield(1) .gt. 10.*cmul(ipts)) then
         itype(ipts) = 3
         alpha(ipts) = 1.
       endif
c      Use asymptotic 1/nu, nu, or smooth connection formula for L11
       do i=ipts + 1,num_extrap_pts + ipts
        cmul(i) = cmul(i-1)*reduct
	 if(efield(1) .le. 0.1*cmul(i)) then
	  if(.not.itok) L11(i) = L11(i-1)*cmul(i-1)/cmul(i)
	  if(itok) L11(i) = L11(i-1)*cmul(i)/cmul(i-1)
	  itype(i) = 1
	 else if(efield(1) .gt. 0.1*cmul(i) .and.
     >           efield(1) .le. 10.*cmul(i)) then
          alpha(i) = 1. - log10(10.*cmul(i)/(efield(1)))
           if(.not.itok) L11(i)=L11(i-1)*((cmul(i)/cmul(i-1))**alpha(i))
           if(itok) L11(i)=L11(i-1)*(cmul(i)/cmul(i-1))
	  itype(i) = 2
	 else if(efield(1) .gt. 10.*cmul(i)) then
	  L11(i) = L11(i-1)*cmul(i)/cmul(i-1)
	  itype(i) = 3
	 endif
	if(.not.l31_extrap_E0) L31(i) = L31(ipts)
	if(l31_extrap_E0) L31(i) = l31_lownu_E0
	L33(i) = L33(ipts)*cmul(ipts)/cmul(i)
        psip(i) = psip(ipts)
        chip(i) = chip(ipts)
        btheta(i) = btheta(ipts)
        bzeta(i) = bzeta(ipts)
        vp(i) = vp(ipts)
        efield(i) = efield(ipts)
       end do
       ipts = ipts + num_extrap_pts
       endif
c
c
       call smooth_5pt(L11,ipts)
       call smooth_5pt(L31,ipts)
       call smooth_5pt(L33,ipts)
c
       merge_factor = 1.    !default
       DO i=1,ipts
        if(merge_hi_nu) then
         merge_factor = 1. - (tanh(4.0*(log(cmul(i))+2.)) + 1.)/2.
	endif
        denom = 1. - 1.5*(cmul(i)*L33(i)/bsq_flx_avg)
c
c     The following issue has been corrected by the Bscale_correction section:
c     ****For devices at low magnetic fields (B .le. 0.5T) the factor
c      denom can pass through zero and go negative. If a trapped
c      fraction input file has been provided, replace the DKES-based
c      denom factor with one based on the asymptotic version of the
c      D33 coefficient [derived from Sugama's eqn. (45)].****
c
c	if(numargs .eq. 6) then
c	  denom_asympt = 1.d0 - ftrapped(isurf)
c	  arg_tanh = (denom - denom_asympt)/denom_asympt
c	  if(denom .lt. denom_asympt) then
c	  denom = denom_asympt*(1. - (tanh(arg_tanh) + 1.)/2.)
c     >          + denom*(tanh(arg_tanh) + 1.)/2.
c          end if
c	  if(denom .lt. denom_asympt)
c     >       denom = 0.707*sqrt(denom**2 + denom_asympt**2)
c	end if
	mstar(i) = (cmul(i)**2)*L33(i)/denom
	lstar(i) = ((L11(i) - (2.*cmul(i)/3.)*utilde_sq)*merge_factor
     >    + 1.5*cmul(i)*(L31(i)**2)/(bsq_flx_avg*denom))
     >    /(e_ovr_c**2)
	lstar_0(i) = (L11(i) - (2.*cmul(i)/3.)*utilde_sq)*merge_factor
     >    + 1.5*cmul(i)*(L31(i)**2)/(bsq_flx_avg*denom)
        if(dkes_limit) then
	  lstar(i) = L11(i)/(e_ovr_c**2)
	  lstar_0(i) = L11(i)
	end if
        if(merge_hi_nu) then
         merge_factor1 = 1. - (tanh(6.0*(log(cmul(i))+2.3)) + 1.)/2.
	 lstar(i) = lstar(i)*merge_factor1
     >            + lstar_0_E0(i)*(1. - merge_factor1)
         lstar_0(i) = lstar_0(i)*merge_factor1
     >            + lstar_0_E0(i)*(1. - merge_factor1)*(e_ovr_c**2)
c
c        Invoke lower limits on L*
c
c         if(lstar_0(i) .le. 0.) lstar_0(i) = 1.e-12
c         if(lstar(i) .le. 0.) lstar(i) = 1.e-12/(e_ovr_c**2)
        endif
        nstar(i) = cmul(i)*L31(i)/(denom*e_ovr_c)   !minus sign to agree with Sugama example
        nstar_0(i) = cmul(i)*L31(i)/denom           !minus sign to agree with Sugama example
c
c      Symmetric limits for lstar, nstar:
c
        if(itok) then                                                       !Case specific
c	lstar(i) = (nstar(i)**2)/mstar(i)
	lstar(i) = mstar(i)*((bzeta(i)/bsq_flx_avg)**2)/
     >    ((e_ovr_c*chip(i))**2)
	nstar(i) = -mstar(i)*(bzeta(i)/bsq_flx_avg)/
     >    (e_ovr_c*chip(i))
	lstar_0(i) = mstar(i)*((bzeta(i)/bsq_flx_avg)**2)/
     >    ((chip(i))**2)
	nstar_0(i) = -mstar(i)*(bzeta(i)/bsq_flx_avg)/
     >    (chip(i))
        endif
c
c        Invoke lower limits on M*, N*
c
	if(mstar(i) .le. 0.) mstar(i) = 1.e-12
c
	amat(1,1) = chip(i)*btheta(i)/bsq_flx_avg
	amat(2,1) = psip(i)*bzeta(i)/bsq_flx_avg
	amat(1,2) = -e_ovr_c*psip(i)*chip(i)
	amat(2,2) = e_ovr_c*psip(i)*chip(i)
c
	bmat(1,1) = mstar(i)
	bmat(2,1) = nstar(i)
	bmat(1,2) = nstar(i)
	bmat(2,2) = lstar(i)
c
	cmat(1,1) = amat(1,1)
	cmat(2,1) = amat(1,2)
	cmat(1,2) = amat(2,1)
	cmat(2,2) = amat(2,2)
c
	tmp = matmul(bmat,cmat)
	visc = matmul(amat,tmp)
c
	mpp(i) = four_pi_sq*visc(1,1)/vp(i)
	mpt(i) = -four_pi_sq*visc(1,2)/vp(i)
	mtt(i) = four_pi_sq*visc(2,2)/vp(i)
c
        coef = (four_pi_sq/vp(i))*((psip(i)*chip(i))**2)
	c1 = ((btheta(i)/bsq_flx_avg)**2)/(psip(i)**2)
	c2 = 2.*btheta(i)/(bsq_flx_avg*psip(i))
	c3 = (btheta(i)*bzeta(i)/bsq_flx_avg)/(psip(i)*chip(i))
	c4 = btheta(i)/(bsq_flx_avg*psip(i))
     >     - bzeta(i)/(bsq_flx_avg*chip(i))
	c5 = ((bzeta(i)/bsq_flx_avg)**2)/(chip(i)**2)
	c6 = 2.*bzeta(i)/(bsq_flx_avg*chip(i))
        mpp_0(i) = coef*(c1*mstar(i) - c2*nstar_0(i) + lstar_0(i))
        mpt_0(i) = -coef*(c3*mstar(i) + c4*nstar_0(i) - lstar_0(i))
        mtt_0(i) = coef*(c5*mstar(i) + c6*nstar_0(i) + lstar_0(i))
        tb = char(9)
        WRITE(NUNOUT,97) cmul(i),tb,efield(i),tb,mpp_0(i),tb,
     >  mpt_0(i),tb,mtt_0(i),tb,mstar(i),tb,
     >  lstar(i),tb,nstar(i),tb,L11(i),tb,L31(i),tb,L33(i),tb,wtov
        if(ie .eq. 1) then
	  lstar_0_E0(i) = lstar(i)
        endif
       L11_lore(i,ie) = L11(i)
       L31_lore(i,ie) = L31(i)
       L33_lore(i,ie) = L33(i)
       END DO    ! i=1,ipts - Collisionality loop
       END DO    ! ie = 1,15 - Electric field loop

       ne = 14; nc = ipts
       do iout=37,39
       write(iout,*) nc, ne
       do i=1,ipts
        ioffset = ipts - i + 1
        write(iout,*) cmul(ioffset)
       end do
       do ie=1,ne
        write(iout,*) efield_lore(ie)
       end do
       do ie=1,ne
        do i=1,ipts
         ioffset = ipts - i + 1
	 if(iout .eq. 37) then
	  write(iout,*) L11_lore(ioffset,ie)
	 else if(iout .eq. 38) then
	  write(iout,*) L31_lore(ioffset,ie)
	 else if(iout .eq. 39) then
	  write(iout,*) L33_lore(ioffset,ie)
	 end if
        end do
       end do
       end do  !iout=37,39
       
       END DO    ! isurf = surf0,outer_surf,delta_surf - Flux surface loop
   98 FORMAT(a1)
   99 FORMAT(18(e12.5,a1),e12.5)
   97 FORMAT(12(e12.5,a1),e12.5)
      end
c
c
c
      subroutine smooth_5pt(y,npts)
      real*8, dimension(npts) :: y
      integer npts,i,istat
      real*8, dimension(5) :: s5
      real*8, dimension(:), allocatable :: temp
      allocate(temp(npts), stat=istat)
      s5(1) = 0.2d0; s5(2) = 0.5d0; s5(3) = 1.d0
      s5(4) = s5(2); s5(5) = s5(1)
      do i=1,npts
       temp(i) = y(i)
      end do
      do i=3,npts-2
       y(i) = (s5(1)*temp(i-2) + s5(2)*temp(i-1)
     >       + s5(3)*temp(i) + s5(4)*temp(i+1)
     >       + s5(5)*temp(i+2))/2.4d0
      end do
      deallocate(temp)
      return
      end
