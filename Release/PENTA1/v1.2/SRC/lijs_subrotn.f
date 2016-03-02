      subroutine dkes_coef(result1,arg1,jval)

c   This subroutine version performs the energy integrations of the 
c	transport coefficients La, Ma, Na as defined in equation (36)
c	of [1].
c
c	INPUTS:
c		arg1 - File name in which the normalized monoenergetic coefficients are stored.
c				The file must contain lists of the collisionalities and electric fields, 
c				followed by the monoenergetic diffusion coefficients.  The file should
c				also begin with two numbers (wtov, dsdr) no longer used but maintained
c				for legacy reasons.
c		lcoef_opt - If 4 the log of the coefficients is taken, if 5 the coefficients
c				 	are read as is.
c		jval - The value of "j" in eq. (36).
c		ma,qa,na,nb,qb,vta,vtb,lnlambda - Particle species parameters as used in eqs (25,36).
c		abs_Er_ovth - |Er|/vth for species a in Volts*second/m^2
c
c     OUTPUTS
c		result1 - The integrated coefficient.
c
c	NOTES
c		-Two logical variables can be set below, see beginning of the main 
c           code for details.
c		-A function version of the error function has been appended to this file.
c		  This can be removed or commented out if a library version exists.
c		-The Lstar file is actually L*/e^2 to make the unnormalization process simpler.
c		-The spline order for interpolating the monoenergetic coefficient matricies
c		  are set below before the first dbsnak call.
c
c	[1] All references to equations are referring to [Sugama and Nishimura
c	        PoP 9 4637, (2002)]
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C-----------------------------------------------
C   M o d u l e s 
C-----------------------------------------------
	use penta_kind_mod
	use phys_const
      use Vblck
	use lijs_transfer
	use dkes_coef_transfer
      use V2dspline_bspline
      use bspline
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
	integer(iknd) :: jval
	real(rknd) :: result1
	character arg1*50

	integer(iknd) :: i, j, istat, nofun
	real(rknd) xtol, epsabs, epsrel, junk, coefl, xlit, xbig,
     1   xminc, xmaxc, xmine, xmaxe, xbot, xtop, aqc8_flag,
	2   s, err, gamtmp, cmin, cmax, emin, emax, pi, ktop, kbot
	real(rknd), dimension(3) :: xlit4, xbig4
	real(rknd), dimension(:),allocatable ::  efldlg_nrm
	real(rknd), dimension(:,:),allocatable :: gamln2d
      character :: lijs_in*60
	logical :: use_smart_integral_limits
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      real(rknd), EXTERNAL :: fun123, xlim
	real(rknd), EXTERNAL :: fun123x
C-----------------------------------------------
c
c     Set constants and logical variables
c
	pi=4._rknd*atan(1._rknd)
	use_smart_integral_limits =.true.  !If true the integration limits are checked against the cmul data and adjusted accordingly
	use_log_interp = .true.		   !If true the interpolation is performed using log(cmul) and log(efield)
c
c	Set default energy integration limits
c
      data xlit4/ 0.05_rknd, 0.28_rknd, 0.22_rknd/  
      data xbig4/ 3._rknd, 3.9_rknd, 4.3_rknd/     
c
c	Set tolerance on calculating the limits from the cmul values
c
      data xtol/1.e-5_rknd/						
c
c	Set tolerances for integration
c	  epsabs (epsrel) is the absolute (relative) integration parameter for
c	  the adaptive quadrature integration routine "quanc8"
c
      data epsabs, epsrel/2.e-5_rknd, 2.e-3_rknd/   								  
c
c     Read data from the file indicated by arg1
c
      lijs_in = trim(arg1)
      open(unit=7,file=lijs_in,status='old')
	read (7, *) junk,junk    !wtov, dsdr no longer needed    
      read (7, *) nc, ne       !number of cmul and efield vals
c
c	Allocate variables
c
	allocate(gamln2d(nc,ne), stat=istat) 	
	allocate(cmul(nc),stat=istat)
	allocate(cmullg(nc),stat=istat)
	allocate(efield(ne),stat=istat)
	allocate(efldlg(ne),stat=istat)
	allocate(efldlg_nrm(ne),stat=istat)
	allocate(c_spl(nc,ne), stat=istat)
c
c	Read cmul (nu/v) values assuming cmul(i+1) > cmul(i) > 0 and define 
c		log and limits of cmul values
c
      do i = 1, nc 
        read (7, *) cmul(i)
      end do
	cmullg = log(cmul)
      cmin = cmul(1)
      cmax = cmul(nc)
c
c	Read efield values assuming efield(i+1) > efield(i) > 0
c      
      do j = 1, ne 
        read (7, *) efield(j)       
      end do
	efldlg = log(efield)	
      emin = efield(1)
      emax = efield(ne)      
c
c	Loop over efield and cmul and read coefficients into 2-D array
c
      do j = 1, ne
        do i = 1, nc
          read (7, *) gamtmp
          if(lcoef_opt .eq. 4) gamln2d(i,j) = log(gamtmp)
          if(lcoef_opt .eq. 5) gamln2d(i,j) = gamtmp
        end do
      end do
      close(unit=7) !close input file    
c
c	Define normalized log(efield) to be used in spline fitting
c	
      efldlg_min = minval(efldlg)
      efldlg_max = maxval(efldlg)
      do j=1,ne
        efldlg_nrm(j) = (efldlg(j)-efldlg_min)/(efldlg_max-efldlg_min)
      end do
c
c	Calculate 2D spline coefficients for 'x'=log(cmul)
c	 and 'y'=normalized log(efield)
c
      kxord = 2; kyord = 2  !spline orders (1 is nearest, 2 linear, 3 spline)
      allocate(xt_c(nc+kxord), stat=istat)
      allocate(xt_e(ne+kyord), stat=istat)
	if ( use_log_interp ) then
        call dbsnak(ne,efldlg_nrm,kyord,xt_e)	!compute 'not-a-knot' sequence
        call dbsnak(nc,cmullg,kxord,xt_c)		!compute 'not-a-knot' sequence
        call dbs2in(nc,cmullg,ne,efldlg_nrm,gamln2d,nc,kxord, !2D spline fit
     >            kyord,xt_c,xt_e,c_spl)	
	else
	  call dbsnak(ne,efield,kyord,xt_e)	    !compute 'not-a-knot' sequence
        call dbsnak(nc,cmul,kxord,xt_c)		!compute 'not-a-knot' sequence
        call dbs2in(nc,cmul,ne,efield,gamln2d,nc,kxord, !2D spline fit
     >            kyord,xt_c,xt_e,c_spl)	
	end if
c
c	Determine integration limits
c
	if (jval == 1) then                 !  Default values
        xlit = xlit4(1)
	  xbig = xbig4(1)
      elseif (jval == 2) then
        xlit = xlit4(2)
	  xbig = xbig4(2)
      elseif (jval ==3) then
        xlit = xlit4(3)
	  xbig = xbig4(3)
      else 
	  write(*,*) 'bad jval = ',jval
	  stop
	endif
	xbot=xlit
	xtop=xbig
c
c	Optionally use integral limits set by the cmul and efield data.
c     These are the normalized velocities [x] at which cmul or efield 
c	  would go out of range for a thermal particle.
c
	if (use_smart_integral_limits) then         
	  xminc = xlim(xtol,cmax,1._rknd,-0.1_rknd) !  Integration limits as determined by cmul data.
        xmaxc = xlim(xtol,cmin,1._rknd,0.5_rknd)  
        xmine=abs_Er_ovth/emax					!  Integration limits as determined by efield data.
	  xmaxe=abs_Er_ovth/emin

        xbot = max(xminc,xlit,xmine) !  Use most conservative limits
        xtop = min(xmaxc,xbig,xmaxe)
	endif
c
c  compute l(i,j) using "quanc8", and optionally type to screen - NOTE Changed integration to K instead of x 8/2007
c 
      efld0 = abs_Er_ovth   !set variables for passing to fun123
	log_coeff=.false.; if ( lcoef_opt == 4 ) log_coeff=.true.  !if log_coeff then take exp(coeff) when computing integral
	jpass=jval	!j value for computing integral (see Eq. 36)

	kbot=xbot**2_iknd !actual integration limits (in K, normalized energy)
	ktop=xtop**2_iknd
	
c	kbot=1.e-6_rknd; xbot=sqrt(kbot)
c	ktop=10._rknd; xtop=sqrt(ktop)
	call quanc8(fun123, kbot, ktop, epsabs, epsrel,
     >             s, err, nofun, aqc8_flag)
c	call quanc8(fun123x, xbot, xtop, epsabs, epsrel,
c     >             s, err, nofun, aqc8_flag)

	!check for integrator error
      if(abs(aqc8_flag) .gt. 1.e-16_rknd) then
	  write(*,'("velocity integration failed: aqc8_flag = ",f15.7)')
     >        aqc8_flag
	  write(*,*) 'component was: ',arg1
	  write(*,*) 'jval was ', jval
	  write(*,*) 'Er/v was ',abs_Er_ovth
	  write(*,*) 'number of function evals: ',nofun
	  write(*,*) 'integration limits: ',kbot,ktop
        stop 14
	endif
c
c	Set output
c
	!coefl=4._rknd*vta*ma/sqrt(pi)	!The 1/e^2 (or c^2/e^2 in cgs) for the L* -> L is included in the definintion of L*
	coefl=na*2._rknd*vta*ma/sqrt(pi)	!Changed from 4.0 to 2.0 JL (7/2008) and added density
      result1 = coefl*s	
c
c  Deallocate storage
c      
      if(allocated(c_spl)) deallocate(c_spl, stat=istat)
      if(allocated(xt_e)) deallocate(xt_e, stat=istat)
      if(allocated(xt_c)) deallocate(xt_c, stat=istat)
	if(allocated(gamln2d)) deallocate(gamln2d, stat=istat)
	if(allocated(cmul)) deallocate(cmul, stat=istat)
	if(allocated(cmullg)) deallocate(cmullg, stat=istat)
	if(allocated(efield)) deallocate(efield, stat=istat)
	if(allocated(efldlg)) deallocate(efldlg, stat=istat)
	if(allocated(efldlg_nrm)) deallocate(efldlg_nrm, stat=istat)
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ADDITIONAL FUNCTIONS
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c--------------------------------------------------------------------------
c  This function calculates the (normalized) velocity x at which the collision
c	frequency defined in Eq. 25 is equal to the collision frequency 'cval'.
c
c  The search range starts at 'xstart' and proceeds in steps 'xstep' until
c	the function f(x)=1-nu(x)/cval changes sign.  The bounds determined
c	in x are then used with the function "zeroin" to find the x value
c	at which f(x)=0 within the tolerance 'xtol'
c
c  The lower limit in x is set by the parameter xzer below.
c
      function xlim (xtol,cval,xstart,xstep) Result(xlim_out)
	use penta_kind_mod
	use fxvars
      implicit none
      real(rknd) xtol,cval,xstart,xstep
      real(rknd) :: xval_old, xval_new, fval_old, fval_new
      real(rknd), EXTERNAL :: zeroin,fxlim
	real(rknd),parameter :: xzer = 1.e-8_rknd
	integer(iknd) :: icount
	integer(iknd),parameter :: imax = 1000_iknd
        real(rknd) :: xlim_out
	cmlx=cval	!cmlx is passed to the function fxlim via module fxvars
      xval_new = xstart	!starting xval
	icount=0_iknd
      fval_new = fxlim(xval_new) !evaluate f(x)
      if (fval_new /= 0._rknd) then !check if the first value was right

	  fval_old = fval_new	
        xval_new = max(xzer,xval_new + xstep)   !increase x by xstep
        fval_new = fxlim(xval_new)              !evaluate f(x) at new x
        do while(fval_old*fval_new > 0._rknd .and. icount < imax)	  !wait for a sign change in f(x) or 1000 iterations
		icount=icount+1
          xval_old = xval_new
          fval_old = fval_new
          xval_new = max(xzer,xval_old + xstep) !increase x by xstep
          fval_new = fxlim(xval_new)			  !evaluate f(x) at new x	
        end do

	  if (fval_new /= 0._rknd) then                     !after sign change use the x range to 
		xval_new = zeroin(xval_new,xval_old,fxlim,xtol)	!find the zero within tolerance xtol
	  endif
	endif
	if (icount >= imax) then
		write(*,*) 'Maximum function evaluations exceeded in xlim'
		stop
	endif
	!stop
      xlim_out = xval_new

      return 
      end function xlim
c
c--------------------------------------------------------------------------
c  This function is used by the functions "xlim" and "zeroin" to find the 
c	value 'xval' at which nu(x)=cmlx.  The value cmlx is passed via the 
c	module fxvars
c
      function fxlim (xval) Result(fxlim_out)
      USE fxvars
	use penta_kind_mod
      implicit none
      real(rknd) xval
      real(rknd), EXTERNAL :: colint
      real(rknd) :: fxlim_out
 
      fxlim_out = 1._rknd - colint(xval)/cmlx

      return 
      end function fxlim
c
c--------------------------------------------------------------------------
c  This function calculates the energy-dependent collision frequency as defined
c	in Eq. 25.  
c
c  The input xa=v/vth_a is passed by argument, all other inputs are passed
c	by the module lijs_transfer.  These inputs are na,qa,nb,qb,
c	ma and lnlambda.  The subscripts a,b refer to particle species.
c
c  The output is the collision frequency over v  (nu/v)
c	
      function colint (xa) Result(colint_out)
	USE penta_kind_mod
	use phys_const
	use lijs_transfer
      implicit none
	real(rknd) :: xa
        real(rknd) :: colint_out

      real(rknd) :: pi,part_aa,part_ab,xb,colint0,xa2,xb2
	real(rknd),external :: derf

      pi = 4._rknd*atan(1._rknd)

	xb=xa*(vta/vtb)

	xa2=xa*xa
	xb2=xb*xb

	part_aa=na*qa**2*( (2._rknd*xa2-1._rknd)*derf(xa) + 
     >        xa*2._rknd/sqrt(pi)*dexp(-xa2) )/(2._rknd*xa2)
	part_ab=nb*qb**2*( (2._rknd*xb2-1._rknd)*derf(xb) + 
     >        xb*2._rknd/sqrt(pi)*dexp(-xb2) )/(2._rknd*xb2)
	
	colint0=4._rknd*pi*qa**2*lnlambda
     >   /(ma**2*(vta*xa)**3)/(4._rknd*pi*ep0)**2
	colint_out=colint0*(part_aa+part_ab)

	colint_out=colint_out/(xa*vta)

      return 
      end function colint
c
c--------------------------------------------------------------------------
c  This function returns the integrand used in conjunction with routine 
c	"quanc8" to calculate the integrated diffusion coefficient.
c
c	Use either log(cmul) and log(efield) in the interpolation or just cmul,efield
c	  using logical variable use_log_interp (set in dkes_coef)
c
      function fun123 (Kval) Result(fun123_out)
      USE Vblck
	use penta_kind_mod
      use V2dspline_bspline
	USE bspline
      use phys_const

      implicit none
      real(rknd) :: Kval
      real(rknd) :: fun123_out

      real(rknd) :: c, e, xsq, e_nrm, ff, emin, emax, cmin, cmax
	real(rknd) :: xval
      real(rknd), EXTERNAL :: colint
c
c	K (the integration variable, normalized energy) is defined as
c	  x**2, where x=v/vt (vt is thermal velocity).
c
	xval=sqrt(Kval)
c
c	Calculate cmul and efield values for this velocity
c
	if ( use_log_interp ) then
        c = log(colint(xval))
        e = log(efld0/xval)
	  emin=efldlg(1)
	  emax=efldlg(ne)
	  cmin=cmullg(1)
	  cmax=cmullg(nc)
	else 
        c = colint(xval)
        e = efld0/xval
	  emin=efield(1)
	  emax=efield(ne)
	  cmin=cmul(1)
	  cmax=cmul(nc)
	end if
c
c	Check for cmul or efield values out of range.  If efield is below the smallest database
c	  efield then use the smallest database efield.
c	Then look up coefficient at efield (e) and cmul (c) values from spline fits if
c	  the cmul and efield values are in range
c	The integrand is defined as in Eq. (36), but with K**2 instead of sqrt(K) to
c	  account for the K**1.5 dependence of the La,Ma,Na coefficients (when converting
c	  from L*,M*,N*)
c
      if(e .lt. emin) e = emin    
      if(c .ge. cmin .and. c .le. cmax .and.
     1       e .ge. emin .and. e .le. emax) then

	  if ( use_log_interp ) then
	    e_nrm = (e - efldlg_min)/(efldlg_max - efldlg_min)
	    ff = dbs2vl(c,e_nrm,kxord,kyord,xt_c,xt_e,nc,ne,c_spl)
	  else
          ff = dbs2vl(c,e,kxord,kyord,xt_c,xt_e,nc,ne,c_spl)
	  end if

        if( log_coeff ) then
            fun123_out=exp(ff - Kval)*Kval**2*(Kval - 2.5_rknd)  !the K^2 instead of sqrt(K) accounts for K^1.5 dependance of the coefficients
     1		**(jpass-1_iknd)
        else 
            fun123_out=ff*exp(-Kval)*Kval**2*(Kval - 2.5_rknd)
     1		**(jpass-1_iknd)
        endif
      else !By default, the integrand is set to zero if the request is out of range
        fun123_out = 0._rknd   
	  !write(*,*) 'efield,cmul,Kval= ',e,c,Kval
      endif
	!uncomment this to write to file
c        write(12,'(e15.7,4(2x,e15.7),2x,i3)') e_out,c_out,xval,
c     >   ff(1),fun123_out,lcoef_opt
      return 
      end function fun123
c
c-----------------------------------------------------------------------
c
c  a zero of the function  f(x)  is computed in the interval ax,bx
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
      function zeroin(ax,bx,f,tol) Result(zeroin_out)
	use penta_kind_mod
      implicit none

      real(rknd) ax, bx, f, tol
      real(rknd) :: zeroin_out

      integer(iknd),save :: jfirst
      real(rknd) :: a, b, c, d, e, eps, fa, fb, fc, tol1,
     1   xm, p, q, r, s
      EXTERNAL f
C-----------------------------------------------

c
c  compute eps, the relative machine precision
c
      eps = epsilon(1._rknd)
      tol1 = 1._rknd + eps
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
   40 tol1 = 2._rknd*eps*abs(b) + 0.5_rknd*tol
      xm = .5_rknd*(c - b)
      if (abs(xm) .le. tol1) go to 90
      if (fb .eq. 0._rknd) go to 90
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
      p = 2._rknd*xm*s
      q = 1._rknd - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2._rknd*xm*q*(q - r) - (b - a)*(r - 1._rknd))
      q = (q - 1._rknd)*(r - 1._rknd)*(s - 1._rknd)
c
c adjust signs
c
   60 if (p .gt. 0._rknd) q = -q
      p = abs(p)
c
c is interpolation acceptable
c
      if ((2._rknd*p) .ge. (3._rknd*xm*q - abs(tol1*q))) go to 70
      if (p .ge. abs(0.5_rknd*e*q)) go to 70
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
      if ((fb*(fc/abs(fc))) .gt. 0._rknd) go to 20
      go to 30
c
c done
c
   90 zeroin_out = b
      return
      end function zeroin
c-----------------------------------------------------------------------
c	This function calculates the error function for input xin
c
	function derf(xin) Result(derf_out)
	use penta_kind_mod
      implicit none
	integer(iknd) :: k
        real(rknd) :: derf_out
	real(rknd) :: eps, pi, xin2, r, er, c0
	real(rknd), intent(in) :: xin

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
