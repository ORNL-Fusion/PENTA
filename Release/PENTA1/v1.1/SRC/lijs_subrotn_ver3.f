!
!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a BETA version of the program LIJS, which is currently
!       under development by D. A. Spong at the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a BETA version, this program is subject to change
!       and improvement without notice.
!
!
      subroutine dkes_coef(result,arg1,iterate,lambda)
c
c    To compile this code on an IBM RS/6000 workstation, type:
c
c      xlf vmodules.f
c      xlf -o xlijs -lnag lijs_f90.f vmodules.o
c
c    To compile on a Linux workstation, using the Lahey compiler:
c
c      lf95 -c vmodules.f
c      lf95 -c lijs_f90.f
c      lf95 -o xlijs lijs_f90.o vmodules.o /usr/local/lff95/lib/libnag.a
c
c
c   This is the subroutine version of the LIJS.FOR code used to do the
c   energy integrations of DKES transport coefficients.
c
c   This version is a modification of lijs_f90.f code
c   with the interactive request/response functions commented out.  The data
c   that was entered in interactively in lijs_f90.f is now either hard-coded in,
c   passed in through the argument list, or passed in through the lijs_transfer
c   f90 module.  The purpose of this code is to have a subroutine that can
c   be repetitively called by an external driver routine.
c
c   This version has been updated by Don Spong
c   from the original version of van Rij in the following areas:
c
c   - converted to Fortran-90 style
c   - all real variables are made double precision to avoid round-off and
c       be consistent with our NAG library
c   - accept/type statements used for interactive input changed to read/write
c   - uppercase changed to lower case
c   - amax1, alog, amin1, etc. intrinsics changed to new versions
c   - NAG function/subroutine names updated
c   - The defunct 2D spline NAG function e01ace (no longer available) was converted
c      over to a call to e01daf (produces 2D spline fit) followed by a call to
c        e02def (evaluates 2D spline fit)
c   -  zeroin function replaced by newer version
c   - Checks and write statements put inside the fun123 function called
c      by the integrator to test when data is requested outside the range
c      of the 2D spline fits.
c   - Some of the interactive input options have been
c      commented out/wired around for testing purposes.
c   - Replaced call to NAG s15aef function for error function with call to
c      fortran erf(derf) function.  Commented out use of power series approx.
c      to error function for small arguments.
c   - Removed subroutines to make DISSPLA plots and commented out
c       calls to them
c   - Merged together L11,L12,L22 and L31,L32 calculations into the same
c       code (originally separated into the lijs.f and l3j.f codes)
c   - Made the input filename prefix a command line argument
c   - Added options for the integration of the M1,M2,M3,L1,L2,L3,N1,N2,and N3
c     viscosity coeficients as defined in eqn (36) of the paper: H. Sugama,
c     S. Nishimura, Physics of Plasmas, vol. 9, pg. 4637, 2002.
c
c     Modified on 9/26/2007 to remove comment lines that are no longer needed
c       (mostly involving previously used spline and integration routines)
c
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      USE Vblck
      USE Vcfxs
      use V2dspline
      use V2dspline_bspline
      use lijs_transfer
      USE bspline
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer, dimension(260) :: iwrk
      integer :: ice0, lw, liw, npts0, npts1, klam, ncpts0, irplt, iraw
     1   , i, j, k, itabu, iprt, ijkbot, ijktop, n, ifail, ifpltg
     2   , ifpltd, l, unit, l312, j312, iterate, j4
      real*8::lamax,lamin,numin,numax,lambda
      real*8, dimension(50) :: xclo, xchi, cmulx, weovth
      real*8, dimension(3) :: xlit1, xbig1
      real*8, dimension(2) :: xlit3, xbig3
      real*8, dimension(3) :: xlit4, xbig4
      real*8 xtol, epsabs, epsrel,
     1   dsdr, wtovth, vth, efact, bvtom1, coefl, xlit, xbig, 
     2   rlbig, rlsml, cmullo, cmulhi, xminc, xmaxc, xmine, xmaxe, xmin
     3   , xmax, facmin, facmax, alo, da, cbot, ctop, ebot, etop, xbot, 
     4   xtop, abser, blnk, cminlg, cmaxlg, eminlg, emaxlg, result
      real*8 ctest,etest,gamtest,aqc8_flag
      character*50  status
      integer nargs, numargs, numchars
      character arg1*50,lijs_in*60,lijs_out*60,
     1   lijs_old*60,lijs_data*60,lijs_error*60
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real*8, EXTERNAL :: fun123, fxclo, fxchi, colint
C-----------------------------------------------
c
c  Using monoenergetic data generated by the "dkes" code
c    "lijs.f" computes linear combinations of:
c
c  lcoef_opt = 1
c  the energy-integrated diffusion coefficients l(i,j) [1 <= i,j <= 2]
c
c  lcoef_opt = 3
c  the bootstrap/flow coefficients l(i,j) [i=3, j=1,2]
c
c  lcoef_opt = 4
c  the momentum conserving transport/viscous coefficients
c  M1, M2, M3 or N1, N2, N3, or L1, L2, L3
c
c  the user must provide the input file "filename_prefix.inp", which contains
c  lists of the collisionalities and electric fields, followed by the
c  monoenergetic diffusion coefficients.  For lcoef_opt = 1, these
c  should be the DKES L11 coefficients.  For lcoef_opt = 3, these
c  should be the DKES L31 coefficients. For lcoef_opt = 4, these
c  should be either the M*, N*, or L* coefficients obtained from running
c  DKES followed by the post-processing code viscosity.f
c
c  "lijs.f" creates the output files "filename_prefix.out", the tabulated-data
c  file "filename_prefix.data", and the file "filename_prefix.err" list of
c  data requested by the integrator which is outside the database
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  there can be no more than "ice" (ice<=ice0) values of cmul or efield
c  if "ice0" is increased, the dimensions of cmul, efield,
c  gamcd, cmullg, efldlg, gamln, xx, work, am, and dee must be increased
c  accordingly
 
      data ice0, lw, liw/60, 2000, 260/
 
c  l(i,j) is evaluated at "npts" (npts<=npts0) equally spaced values of
c  log(nu/wtr) (if the user enters 0, the value "npts1" is used)
c  if "npts0" is increased, the dimensions of rlij, rnuowt, xclo, xchi,
c  and cmulx must be increased accordingly
 
c  up to "klam" values of e*phi/ti can be selected
c  if "klam" is increased, the dimensions of rlij and rlamvl must be
c  increased accordingly
 
c  the collision integral used in the calculation of the l(i,j) integral
c  is computed by means of a cubic-spline interpolation of data defined
c  at "ncpts0" points: if "ncpts0" is increased, the dimensions of xcol,
c  colspl, rknot, and cofspl must be increased accordingly
 
      data npts0, npts1, klam, ncpts0/50, 30, 20, 50/
 
c  integration limits constraints
 
      data xlit1/ 0.3466d0, 0.4691d0, 0.5550d0/
      data xbig1/ 3.373d0, 3.835d0, 4.290d0/
      data xlit3/ 0.265d0, 0.4d0/
      data xbig3/ 3.14d0, 3.7d0/
      data xlit4/ 0.35d0, 0.28d0, 0.22d0/
      data xbig4/ 3.4d0, 3.9d0, 4.3d0/
      data xtol/1.d-5/
 
c  epsabs (epsrel) is the absolute (relative) integration parameter for
c  the "nag" adaptive quadrature routine "d01ajf"
 
      data epsabs, epsrel/2.d-6, 2.d-4/
c      data epsabs, epsrel/1.d-5, 1.d-3/
 
      itkf = 0
      klam1 = klam
c
c   arg1 is the input file prefix
c
       lijs_in = trim(arg1)
       lijs_out = trim(trim(arg1) // '.out')
       lijs_old = trim(trim(arg1) // '.old')
       lijs_error = trim(trim(arg1) // '.err')
       lijs_data = trim(trim(arg1) // '.data')
c      write(*,'(8(e12.5))') c11, c12, c22, c31, c32,
c     1   nulo, temp, lambda
c  option to only plot previously calculated l(i,j)
        irplt = 0
	iraw = 0
c  read input data from user-provided file "filename_prefix.inp", and enter
c  additional physical parameters
 
      if(iterate .eq. 1) open(unit=7,file=lijs_in,status='old')
 
c  temp   = ti (kev)
c  mass   = mi (u)
c  zee    = zi
c  trat   = ti'/ti [for i - i' scattering]
c  masrat = mi'/mi [for i - i' scattering]
c  zeff   = zi'*zi'*ni'*lii'/zi*zi*ni*lii [for i - i' scattering]
c  arad   = a (m) = - phi/phi'
c 
c  Ion parameters:
c         mass = 1.d0
c         zee = 1.d0
c         trat = 1000.d0/200.d0
c         masrat = 1.d0
c         zeff = 1.d0
c         arad = 0.33d0
c  Electron parameters:
c         mass = 1.d0/1837.d0
c         zee = 1.d0
c         trat = 200.d0/1000.d0
c         masrat = 1837.d0
c         zeff = 1.d0
c         arad = 0.33d0
c
      if(iterate .eq. 1)read (7, *) dsdr, wtovth          !Read from filename_prefix.inp
  100 format(e10.2)
 
c  icmul = number of cmul values
c  jefld = number of efield values
 
      if(iterate .eq. 1) read (7, *) icmul, jefld       !Read from filename_prefix.inp
  200 format(2i5)
 
      ice = max0(icmul,jefld)
      if (ice > ice0) then
         write (*, 7) ice, ice0
    7    format(' dimension(s) exceeded : ice, ice0 = ',2i5)
         stop 11
 
c  must have cmul(i+1) > cmul(i) > 0
 
      endif
      if(iterate .eq. 1) then
      do i = 1, icmul
         read (7, *) cmul(i)         !Read from filename_prefix.inp
  300    format(e11.3,e11.3)
         cmullg(i) = log(cmul(i))
      end do
      cmin = cmul(1)
      cmax = cmul(icmul)
      endif
 
c  must have efield(j+1) > efield(j) > 0

      if(iterate .eq. 1) then 
      do j = 1, jefld
c         read (7, 300) efield(j), weovth(j)  !weovth not needed
         read (7, *) efield(j)       !Read from filename_prefix.inp
         efldlg(j) = log(efield(j))  !fit negative of log due to bsplib convention
      end do
      emin = efield(1)
      emax = efield(jefld)
      endif
c
c
c   Determine if L(1,1), L(1,2), L(2,2) coefficients,
c   L(3,1), L(3,2) coefficients, or M1,M2,M3,N1,N2,N3,L1,L2,L3
c   coefficients are desired [note: L(3,1) = -L(1,3), and
c   L(3,2) = -L(2,3)].
c
c         lcoef_opt = 1
c 
c  lijs.for will compute either L = c11*L(1,1) + c12*L(1,2) + c22*L(2,2)
c                  or c31*L(3,1) + c32*L(3,2)
c
         if(lcoef_opt .eq. 1) then
          ijkbot = 5
          if (c11 == 0.d0) ijkbot = 7
          if (c11==0.d0 .and. c12==0.d0) ijkbot = 9
          ijktop = 9
          if (c22 == 0.d0) ijktop = 7
          if (c22==0.d0 .and. c12==0.d0) ijktop = 5
          xlit = xlit1(ijkbot/2-1)
          xbig = xbig1(ijktop/2-1)
         else if(lcoef_opt .eq. 3) then 
          xlit = xlit3(1); l312 = 1
          if(c31 == 0.d0) then
            xlit = xlit3(2); l312 = 2
          endif
          xbig = xbig3(2); l312 = 2
          if(c32 == 0.d0) then
            xbig = xbig3(1); l312 = 1
	  endif
          j312 = 2 + 2*l312
         else if(lcoef_opt .ge. 4) then 
          xlit = xlit4(1);j4=0
          if(vs1 == 0.d0) then
             xlit = xlit4(2);j4=1
          endif
          if(vs1 == 0.d0 .and. vs2 == 0.d0) then
             xlit = xlit4(3);j4=2
          endif
          xbig = xbig4(3); j4 = 2
          if(vs3 == 0.d0) then
             xbig = xbig4(2);j4=1
          endif
          if(vs3 == 0.d0 .and. vs2 == 0.d0) then
             xbig = xbig4(1);j4=0
          endif
	 endif
c
      if(iterate .eq. 1)then      
      k = 0
      do j = 1, jefld
         do i = 1, icmul
            k = k + 1
            read (7, *) gamcd(k)              !Read from filename_prefix.inp
            if(lcoef_opt .eq. 1) gamln2d(i,j) = log(gamcd(k))
            if(lcoef_opt .eq. 3) gamln2d(i,j) = gamcd(k)
            if(lcoef_opt .eq. 4) gamln2d(i,j) = log(gamcd(k))
            if(lcoef_opt .eq. 5) gamln2d(i,j) = gamcd(k)
         end do
      end do
c      close(unit=7)
c
      do j=1,jefld
       do i=1,icmul
        k = jefld*(i-1) + j
        gamln(k) = gamln2d(i,j)
       end do
      end do
c      
c  close input file 
      close(unit=7)              !Finished reading from filename_prefix.inp
      endif
c
      nc = icmul
      ne = jefld
      allocate(c_spl(nc,ne), stat=istat)
      allocate(zint(nc,ne), stat=istat)
      do i=1,nc
       do j=1,ne
        zint(i,j) = gamln2d(i,j)
       end do
      end do
      efldlg_min = minval(efldlg)
      efldlg_max = maxval(efldlg)
      do j=1,ne
       efldlg_nrm(j) = (efldlg(j)-efldlg_min)/(efldlg_max-efldlg_min)
      end do
      kxord = 2; kyord = 2
      allocate(xt_c(nc+kxord), stat=istat)
      allocate(xt_e(ne+kyord), stat=istat)
      call dbsnak(ne,efldlg_nrm,kyord,xt_e)
      call dbsnak(nc,cmullg,kxord,xt_c)
      call dbs2in(nc,cmullg,ne,efldlg_nrm,zint,nc,kxord,
     >  kyord,xt_c,xt_e,c_spl)
c 
         rootm = sqrt(masrat/trat)
	 col0 = one
c	 col0 = colint(one)     !not needed for nu_0 that is used (9/27/07)
c	 write(*,*) col0
         vth = 4.392825d+05*sqrt(temp/mass)
         efact = 2.276439d-3*sqrt(temp*mass)/dsdr
         bvtom1 = 2.d0*efact/zee
         if(lcoef_opt .eq. 1) coefl=vth*bvtom1*bvtom1/sqrt(atan(1.d0))
         if(lcoef_opt .eq. 3) coefl=vth*bvtom1/sqrt(atan(1.d0))
         if(lcoef_opt .ge. 4) coefl=1.*vth*mass*1.6726d-27
     1                              /sqrt(atan(1.d0))
         efact = efact/arad
         numin = cmin/wtovth
         numax = cmax/wtovth
         lamin = emin/efact
         lamax = emax/efact
c         wefac = .5d0*wtovth*lamin/weovth(1)    !not needed
 
         if (npts == 0) npts = npts1
         npts = min0(npts0,npts) 
   18    continue
 
c  choose e*phi/ti values for l(i,j) calculation
 
         rlbig = -1.d0
         rlsml = 1.d20
   21    continue
         rlamvl(klm) = lambda
         rlbig = max(rlbig,lambda)
         rlsml = min(rlsml,lambda)
 
c  integration limits as determined by cmul data
 
   24    continue
         cmullo = wtovth*nulo
         cmulhi = wtovth*nuhi
         cmlxmn = cmulhi/cmax
         cmlxmx = cmullo/cmin
         xminc = fxclo(xtol)
         xmaxc = fxchi(xtol)
 
c  integration limits as determined by efield data
          xmine = efact*rlbig/emax
          xmaxe = efact*rlsml/emin
c  net integration limits
          xmin = max(xminc,xmine,xlit)
          xmax = min(xmaxc,xmaxe,xbig)
c 
c  integration limits diagnostics, including option to override
c  constraints

         if(lcoef_opt .eq. 1) then
           facmin = exp((-xmin**2))*xmin**ijkbot
           facmax = exp((-xmax**2))*xmax**ijktop
         else if(lcoef_opt .eq. 3) then
           facmin = exp((-xmin**2))*xmin**j312
           facmax = exp((-xmax**2))*xmax**j312
         else if(lcoef_opt .ge. 4) then
           facmin = exp((-xmin**2))*((xmin**2 - 2.5)**j4)*xmin**5
           facmax = exp((-xmax**2))*((xmax**2 - 2.5)**j4)*xmax**5
         endif
         if (xmin==xlit .and. xmax==xbig) go to 26
 
c  set up cubic spline of collision integral
 
   26    continue
c
c     Write out lower/upper bounds of database
c
      cminlg = cmullg(1)
      cmaxlg = cmullg(icmul)
      eminlg = efldlg(1)
      emaxlg = efldlg(jefld)
c
c
 
c  set up equally spaced nu/wtr mesh
 
         alo = log(nulo)
         da = (log(nuhi) - alo)/real(npts - 1)
         do n = 1, npts
            rnuowt(n) = exp(alo + (real(n - 1))*da)
            cmulx(n) = wtovth*rnuowt(n)
            cmlxmx = cmulx(n)/cmin
            cmlxmn = cmulx(n)/cmax
            xclo(n) = max(xmin,fxclo(xtol))
            xchi(n) = min(xmax,fxchi(xtol))
         end do
 
c  cmul and efield limit values used in integration
 
         cbot = colint(xmax)*cmullo
         ctop = colint(xmin)*cmulhi
         ebot = efact*rlsml/xmax
         etop = efact*rlbig/xmin
 
c  compute l(i,j) using "d01ajf", and optionally type to screen
 
         ifail = 0							!NAG option
         do k = 1, klm
            efld0 = efact*rlamvl(k)
             xmine = efld0/emax
             xmaxe = efld0/emin
            do n = 1, npts
               cmul0 = cmulx(n)
                xbot = max(xclo(n),xmine)
                xtop = min(xchi(n),xmaxe)
		
c		xbot=sqrt(1.d-6)
c	    xtop=sqrt(125.)
		call quanc8(fun123, xbot, xtop, epsabs, epsrel,
     >             s, err, nn, aqc8_flag)

               if(abs(aqc8_flag) .gt. 1.e-16) go to 30
               rlij(n,k) = s
               rlij(n,k) = coefl*rlij(n,k)
            end do
            go to 32
 
   30       continue
            write(*,
     >      '("velocity integration failed: aqc8_flag = ",e15.7)')
     >      aqc8_flag
            stop 14
   32       continue
         end do
        result = rlij(1,1)
c
c  Deallocate storage used for 2d spline fits
c
c      write(*,'("1")')
      if(allocated(iwrksp))deallocate(iwrksp, stat=istat)
      if(allocated(lambdasp))deallocate(lambdasp, stat=istat)
      if(allocated(mu))deallocate(mu, stat=istat)
      if(allocated(cc))deallocate(cc, stat=istat)
      if(allocated(wrksp))deallocate(wrksp, stat=istat)
      if(allocated(zint))deallocate(zint, stat=istat)
      if(allocated(xt1))deallocate(xt1, stat=istat)
      if(allocated(c1_spl))deallocate(c1_spl, stat=istat)
      
c      write(*,*) allocated(c_spl),allocated(xt),allocated(vw)
      
      if(allocated(c_spl))deallocate(c_spl, stat=istat)
      if(allocated(xt_e))deallocate(xt_e, stat=istat)
      if(allocated(xt_c))deallocate(xt_c, stat=istat)

      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function fxclo (xtol)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real*8 xtol
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      real*8 :: xzer, xvaln, fvaln, xvalo, fvalo
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real*8, EXTERNAL :: fxmin, zeroin
C-----------------------------------------------
 
c  calculate xminc using the routines fxmin and "zeroin"
 
 
      data xzer/1.d-8/
 
      xvaln = 1.d0
      fvaln = fxmin(xvaln)
      if (fvaln == 0.d0) go to 2
      xvalo = xvaln
      fvalo = fvaln
      xvaln = max(xzer,xvalo - .1d0)
      fvaln = fxmin(xvaln)
      if (fvaln == 0.d0) go to 2
      do while(fvalo*fvaln > 0.d0)
         xvalo = xvaln
         fvalo = fvaln
         xvaln = max(xzer,xvalo - .1d0)
         fvaln = fxmin(xvaln)
         if (fvaln == 0.d0) go to 2
      end do
      fxclo = zeroin(xvaln,xvalo,fxmin,xtol)
      go to 3
    2 continue
      fxclo = xvaln
 
    3 continue
      return 
      end function fxclo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function fxmin (xval)
      USE Vcfxs
      implicit none
      real*8 xval
      real*8, EXTERNAL :: colint
 
c  used in conjunction with the routine "zeroin" to calculate the
c  value of "xval" for which "fxmin" vanishes
 
 
      fxmin = 1.d0 - cmlxmn*colint(xval)
 
      return 
      end function fxmin
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function fxchi (xtol)
      implicit none
      real*8 xtol
      real*8 :: xvaln, fvaln, xvalo, fvalo
      real*8, EXTERNAL :: fxmax, zeroin

c  calculate xmaxc using the routines fxmax and "zeroin"
 
      xvaln = 1.d0
      fvaln = fxmax(1.d0)
      if (fvaln == 0.d0) go to 2
      xvalo = xvaln
      fvalo = fvaln
      xvaln = xvalo + .5d0
      fvaln = fxmax(xvaln)
      if (fvaln == 0.d0) go to 2
      do while(fvalo*fvaln > 0.d0)
         xvalo = xvaln
         fvalo = fvaln
         xvaln = xvalo + .5d0
         fvaln = fxmax(xvaln)
         if (fvaln == 0.d0) go to 2
      end do
      fxchi = zeroin(xvalo,xvaln,fxmax,xtol)
      go to 3
    2 continue
      fxchi = xvaln
 
    3 continue
      return 
      end function fxchi
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function fxmax (xval)
      USE Vcfxs
      implicit none
      real*8 xval
      real*8, EXTERNAL :: colint
 
c  used in conjunction with the routine "zeroin" to calculate the
c  value of "xval" for which "fxmax" vanishes
 
      fxmax = 1.d0 - cmlxmx*colint(xval)
      return 
      end function fxmax
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function colint (xval)
      USE Vblck
      implicit none
      real*8 xval
      real*8 y1, y2, a1, a2, colint1, colint2, xv,
     >  y1sq, y2sq
	real*8,external :: derf
      xv = xval
      y1 = xv; y2 = rootm*xv
      y1sq = y1*y1; y2sq = y2*y2
      a1 = 0.5d0/y1sq; a2 = 0.5d0/y2sq
      colint1 = pfac*derf(y1)*(1.d0 - a1)/y1 + a1*exp(-y1sq)
      colint2 = pfac*derf(y2)*(1.d0 - a2)/y2 + a2*exp(-y2sq)
      colint = (colint1 + zeff*rootm*colint2)/(col0*xv*xv*xv)
      return 
      end function colint
      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function fun123 (xval)
C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
      USE Vblck
      use V2dspline
      use V2dspline_bspline
      use lijs_transfer
      USE bspline
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real*8 xval
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer :: ifail, m, ii
      real*8 :: cval, c, e, xsq, c_out, e_out, e_nrm
      real*8, dimension(1) :: a1, b1, ff
      real*8, EXTERNAL :: colint
C-----------------------------------------------
 
c  used in conjunction with routine "quanc8" to calculate the
c  integrated diffusion coefficient
c 
      cval = colint(xval)
      c = log(cmul0*cval)
      e = log(efld0/xval)
      c_out =  cmul0*cval
      e_out =  efld0/xval   
c
      if(e .lt. efldlg(1)) e = efldlg(1)       
      if(c .ge. cmullg(1) .and. c .le. cmullg(icmul) .and.
     1   e .ge. efldlg(1) .and. e .le. efldlg(jefld)) then
        e_nrm = (e - efldlg_min)/(efldlg_max - efldlg_min)
        ff(1) = dbs2vl(c,e_nrm,kxord,kyord,xt_c,xt_e,nc,ne,c_spl)
c
        xsq = xval*xval
        if(lcoef_opt .eq. 1) then
          fun123=exp(ff(1) - xsq)*(c11+(c12+c22*xsq)*xsq)*xval**5
        else if(lcoef_opt .eq. 3) then
          fun123=ff(1)*exp(-xsq)*(c31+c32*xsq)*xsq*xsq
        else if(lcoef_opt .eq. 4) then
          fun123=exp(ff(1) - xsq)*2.*xval*
     1      (vs1 + (vs2 + vs3*(xsq-2.5))*(xsq-2.5))*xval**4
        else if(lcoef_opt .eq. 5) then
          fun123=ff(1)*exp(-xsq)*2.*xval*
     1      (vs1 + (vs2 + vs3*(xsq-2.5))*(xsq-2.5))*xval**4
        endif
        else
         fun123 = 0.d0   !By default, the integrand is set to zero
                        !if cmul and efield are outside of the supplied data.
			!Would be better to replace this with an extrapolation
			!in the future.
      endif
c        write(12,'(e15.7,4(2x,e15.7),2x,i3)') e_out,c_out,xval,
c     >   ff(1),fun123,lcoef_opt
      return 
      end function fun123
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c-----------------------------------------------------------------------
      real*8 function zeroin(ax,bx,f,tol)
      implicit none
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      real*8 ax, bx, f, tol
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      integer,save :: jfirst
      real*8 :: a, b, c, d, e, eps, fa, fb, fc, tol1,
     1   xm, p, q, r, s
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL f
C-----------------------------------------------
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0)
c
c
c  output..
c
c  zeroin abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
c
c  compute eps, the relative machine precision
c
      eps = 1.0
c   10 eps = eps/2.0
c      tol1 = 1.0 + eps
c      if (tol1 .gt. 1.0) go to 10
      eps = epsilon(eps)
      tol1 = 1.0 + eps
c
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (abs(fc) .ge. abs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0d0*eps*abs(b) + 0.5d0*tol
      xm = .5d0*(c - b)
      if (abs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0d0) go to 90
c
c is bisection necessary
c
      if (abs(e) .lt. tol1) go to 70
      if (abs(fa) .le. abs(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0d0*xm*s
      q = 1.0d0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
      q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
c
c adjust signs
c
   60 if (p .gt. 0.0d0) q = -q
      p = abs(p)
c
c is interpolation acceptable
c
      if ((2.0d0*p) .ge. (3.0d0*xm*q - abs(tol1*q))) go to 70
      if (p .ge. abs(0.5d0*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (abs(d) .gt. tol1) b = b + d
      if (abs(d) .le. tol1) b = b + sign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/abs(fc))) .gt. 0.0d0) go to 20
      go to 30
c
c done
c
   90 zeroin = b
      return
      end function zeroin
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c	This function calculates the error function for input xin
c
	function derf(xin) Result(derf_out)
	use penta_kind_mod
      implicit none
	integer(iknd) :: k
	real(rknd) :: eps, pi, xin2, r, er, c0
	real(rknd), intent(in) :: xin
        real(rknd) :: derf_out

      eps=1.e-15_rknd
	pi=4._rknd*atan(1._rknd)
      xin2=xin*xin
	if (dabs(xin).lt. 3.5_rknd) then
		er=1._rknd
		r=1._rknd
		do k=1,50
			r=r*xin2/(k+0.5_rknd)
			er=er+r
			if (dabs(r).le.dabs(er)*eps) exit
		enddo
		c0=2._rknd/dsqrt(pi)*xin*dexp(-xin2)
		derf_out=c0*er
	else
		er=1._rknd
		r=1._rknd
		do k=1,12
			r=-r*(k-0.5_rknd)/xin2
			er=er+r
		enddo
		c0=dexp(-xin2)/(dabs(xin)*dsqrt(pi))
		derf_out=1.0D0-c0*er
		if (xin.lt.0.0) derf_out=-derf_out
	endif
	end function derf