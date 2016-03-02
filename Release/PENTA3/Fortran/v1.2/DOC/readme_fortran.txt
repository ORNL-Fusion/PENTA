Instructions for "PENTA3" created by J. Lore, 1/4/2011.

These instructions are for the Fortran version of PENTA. PENTA can be run from the 
command line or from a shell script file (ex: run_penta). See the header to 
penta.f90 for additional details on input files and command line arguments.

Essentially, PENTA is run for a given surface (for which equilibrium data and DKES
coefficient data must exist). PENTA then loops over a given range of Er values and
calculates the parallel flow moments and the radial particle fluxes. The roots of 
the ambipolarity constraint are then evaluated and the flows and fluxes are evaluated
at these roots. Clearly then, the Er range must be appropriate (such that roots are found)
and must have enough test points to accurately evaluate the roots (and not skip over multiple
roots). Note that the output files will accumulate rows based on the number of ion species, the
number of Sonine coefficients (additional flows) and the number of Er roots (additional
ambipolar fluxes, for example).


1) RUNNING PENTA
	PENTA assumes that the necessary input files have been generated ahead of 
	time and are located in a single "run folder". The program is then executed
	in the run directory.  Here the input files are described. See the header
	to penta.f90 for the command line arguments.

	a) DKES data files
		For each run surface three files are required, containing the D11
		D13, and D33 coefficients.  These coefficients should be given
		in a "normalized" form (D11*, D13*, D33*) defined such that the Dxx 
		coefficients in the following equations are in units of [m**2/s]:

		D11=D11_star * (m**2 * v**3 / 2 * q**2) 
		D13=D13_star * (m * v**2 / 2 * q) / B_0
		D33=D33_star * (v / 2) / B_0**2
	
		where m,v,q are the particle mass, velocity and charge, respectively
		and B_0 is the reference field strength (in Tesla). Note that based
		on the above definitions the "B field correction" should be carried
		out _before_ writing the DKES data files!

		The files should be named:
		
		D11_star_xxx_s#,D13_star_xxx_s#,D33_star_xxx_s#

		where xxx is a run identifier set in the script and # is the booz_xform
		(VMEC) surface index (2-ns).

		IMPORTANT:  The D33 coefficient should contain the collisional 
				"Spitzer" contribution.  If not, be sure to set
				the proper flag in the run script.
		            The cmul and efield values should be monotonically increasing.
			    The same cmul and efield values should be used for each coefficient
			        for a given surface, but may differ by surface.

		The file is formatted as follows: 
			Row 1:  NC NE
			Rows 2 thru NC+1:  CMUL values
			Rows NC+2 thru (NC+NE+1): EFIELD values
			Remaining rows: Coefficient value corresponding to nested loop 
					over efield, cmul

			Here NC is the number of CMUL values, NE the number of EFIELD
			values, EFIELD=Er/v, CMUL=nu/v.

			The file is read as:
			     ! Read cmul values
			        Do ic = 1, nc 
				  Read (iu_coeff, *) cmul_vec(ic)
				Enddo
			     ! Read efield values  
			        Do je = 1, ne 
				  Read (iu_coeff, *) efield_vec(je)     
				Enddo
						    
			     ! Loop over efield and cmul and read coefficients into 2-D array
			        Do je = 1, ne
			          Do ic = 1, nc
	  		            Read (iu_coeff, *) coeff_tmp
	  		            coef2d(ic,je) = coeff_tmp
 			          Enddo
			        Enddo

	b) ion_params file
		This file contains the information about the ion species.  The file 
		is formatted as a Fortran namelist.  A simple example is:

		&ion_params
		num_ion_species=3
		Z_ion_init=1.0,2.0,6.0
		miomp_init=1.0,2.0,12.0
		/

		Here three ion species are included, the arrry Z_ion_init specifies the 
		charge number and miomp_init the mass normalized to a proton mass. (I.e.
		this plasma contains atomic hydrogen, deuterium and fully stripped carbon).

		Note that the plasma profile file (see below) must reflect the information
		in this file.

	c) plasma profile file
		This file contains information about the density and temperature profiles
		of each plasma species.		

		Name: "plasma_profiles_zzz.dat", where zzz is a profile identifier s
			specified in the run script.

		Format:
		Row 1: npts (number of radial points where data are defined)
		Rows 2-npts+1: r/a  ne  Te  n_i1  T_i1  n_i2  T_i2 ...

		Here the densities are in units of 10^12*cm^-3 (e.g., a plasma density of 1e12 cm^-3
		would be input to the file as 1.000) and the temperatures are in eV.

		IMPORTANT:  The order of the ion density and temperature is different than
				older versions of PENTA.

		The file is read as:
		  Read(iu_pprof,*) np_prof
		  ! Loop over surfaces
		  Do jind=1,np_prof
		    Read(iu_pprof,*) roa_prof(jind),ne_prof(jind),Te_prof(jind),ion_pprof_data
		    ! Break apart ion data
		    Do ispec=1,nis
		      tmp_ind=(ispec-1)*2+1
		      ni_prof(jind,ispec)=ion_pprof_data(tmp_ind)
		      Ti_prof(jind,ispec)=ion_pprof_data(tmp_ind+1)
		    Enddo
		  Enddo

	d) VMEC data file
		This file contains geometry information (probably) calculated from VMEC.

		Name:  profile_data_xxx  (where xxx is the run ident.)
		
		Format:
		Row 1:  S_first  S_last  (the first and last VMEC surface included)
		Row 2:  a   R0		 (major and minor radii in meters)
		Row 3: legend (any text string)
		Remaining rows:  j  r  r/a  chip  psip  btheta  bzeta  dV_dr  <B**2>  iota  B_0

		where:  j is the booz_xform surface index
			r is the average minor radius in meters (not used in current version of PENTA)
			r/a is the sqare root of the normalized toroidal flux
			chip is the radial derivative of the poloidal flux
			psip is the radial derivative of the toroidal flux
			btheta is the poloidal component of B in boozer coords (not used in current version of PENTA)
			bzeta is the toroidal component of B in boozer coords (not used in current version of PENTA)
			dV_dr is the radial derivative of the surface volume
			<B**2> is the flux surface averaged B**2 in Tesla
			iota is the rotational transform (actually iota-bar)
			B_0 is the reference magnetic field (B_0,0 on the surface) in Tesla

	e) run_params
		This file containing the namelist "run_params" which contains
		run settings.  

		An example is given as:

		&run_params
		input_is_Er    = .true.
		log_interp     = .true.
		use_quanc8     = .false.
		read_U2_file   = .true.
		Add_Spitzer_to_D33 = .true.                             
		num_Er_test = 100
		numKsteps   = 1000
		kord_pprof  = 3
		keord       = 2
		kcord       = 2
		Kmin   = 1.d-5
		Kmax   = 10.d
		epsabs = 1.d-8
		epsrel = 1.d-6
		Method  = 'SN'
		/		

		Logical ::    &
		input_is_Er    = .true.,     & ! If true, Er range is (V/cm) else e<a>Er/kT_e
		log_interp     = .true.,     & ! If true, log. interp. of DKES coeffs is used
		use_quanc8     = .false.,    & ! If false, rect. approx. to convolution used
		read_U2_file   = .true.,     & ! If false <U**2> is calculated from D11*
		Add_Spitzer_to_D33 = .true.    ! If true collisional portion of D33* is added
				               !  else it is assumed to be included in-file  
		Integer(iknd) ::    &
		num_Er_test = 20_iknd,        & ! Number of Er points in search range
		numKsteps   = 1000_iknd,      & ! Number of K points (linear) for convolution 
			                        !  (used if use_quanc8=.false.)
		kord_pprof  = 3_iknd,         & ! Spline order for plasma profile fitting 
			                        !  (and <U**2> file)
		keord       = 2_iknd,         & ! Spline order for DKES coeff fitting (efield)
		kcord       = 2_iknd            ! Spline order for DKES coeff fitting (cmul)
		
		Real(rknd) ::    &
		Kmin   = 1.e-5_rknd,         & ! Minimum K in energy convolution
		Kmax   = 10._rknd,           & ! Maximum K in energy convolution
		epsabs = 1.e-8_rknd,         & ! Absolute tolerance for quanc8  
					       !  (used if use_quanc8=.true.)
		epsrel = 1.e-6_rknd            ! Relative tolerance for quanc8  
					       !  (used if use_quanc8=.true.)
		Character(Len=10) ::     &
		Method  = 'DKES'                  ! Which algorithm to use.  Options are
						  !  'T'    = Taguchi
						  !  'SN'   = Sugama-Nishimura
						  !  'MBT'  = Maassberg-Beidler-Turkin
						  !  'DKES' = Direct energy convolution

		MORE TO BE WRITTEN


	f) Utilde2_profile - Contains the quantity <U**2>, where U is the
	     PS flow function as defined by Sugama and Nishimura.  The
	     first row is the number of points, then r/a points and the
	     corresponding <U**2> value.  Note that if read_U2_file is .false.	
	     then <U**2> is calculated from D11 at high collisionality.

2) PENTA output

	a) flows_vs_roa
	b) fluxes_vs_Er
	c) fluxes_vs_roa
	d) plasma_profiles_check

	-- see penta.f90 header
	
	
