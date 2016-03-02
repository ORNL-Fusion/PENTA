Instructions for "PENTA3" created by Jeremy Lore, 5/19/2010.

These instructions are for the Matlab version of PENTA.  Typical runs will use
the script "run_PENTA.m".  Many options are set from this script and should have 
sufficient explanation.  


1) RUNNING PENTA
	PENTA assumes that the necessary input files have been generated ahead of 
	time and are located in a single "run folder".  The path to this folder is
	specified in the run script.  Here the input files are described.

	a) DKES data files
		For each run surface three files are required, containing the D11
		D13, and D33 coefficients.  These coefficients should be 
		in a "normalized" form (D11*, D13*, D33*) defined as:

		D11=D11_star * (m**2 * v**3 / 2 * q**2) 
		D13=D13_star * (m * v**2 / 2 * q) / B_0
		D33=D33_star * (v / 2) / B_0**2
	
		where m,v,q are the particle mass, velocity and charge, respectively
		and B_0 is the reference field strength (in Tesla).  

		The files should be named:
		
		D11_star_xxx_s#,D13_star_xxx_s#,D33_star_xxx_s#

		where xxx is a run identifier set in the script and # is the booz_xform
		(VMEC) surface index (2-ns).

		IMPORTANT:  The D33 coefficient should contain the collisional 
				"Spitzer" contribution.  If not, be sure to set
				the proper flag in the run script.

		The file is formatted as follows: 
			Row 1:  NC NE
			Rows 2-NC+1:  CMUL values
			Rows NC+2-(NC+NE+1): EFIELD values
			Remaining rows: Coefficient value corresponding to nested loop 
					over efield, cmul

			Here NC is the number of CMUL values, NE the number of EFIELD
			values, EFIELD=Er/v, CMUL=nu/v.

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

		Here the densities are in units of cm^-3 and the temperatures are in eV.

		IMPORTANT:  The order of the ion density and temperature is different than
				older versions of PENTA.

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

2) PENTA output
	The primary output is from the function "surf_loop_PENTA" where some postprocessing
	has already been performed.  The output quantities follow, with parentheses indicating
	the array indices.
		is:  surface index where 1 corresponds to the first surface in the run_surfs array, etc
		ie:  number of ambipolar roots.  Note that the roots are ordered by increasing Er.
		Sind: Sonine expansion index.  So the first value of a given species is the 0 index of the Sonine expansion, etc
		i_Er_search:  Index of the Er search values used
		ispec:  Index of the species.

		Notes:  Curly brackets {} indicate a cell array.
		Notes:  All units are mks unless noted

	roa_vals(is):  the r/a value of each run surface
	Er_roots{is}(ie):  The ambipolar roots in V/m.
	Gammas_ambi{is}(ispec,ie):  The ambipolar particle fluxes.  First Smax+1 entries are for electrons, then ion species 1, etc
	Flows_ambi{is}(Sind,ie):   The ambipolar particle flows.  First Smax+1 entries are for electrons, then ion species 1, etc
	gamma_e_vs_Er{is}(i_Er_search):  Electron particle flux vs Er
	gamma_i_vs_Er{is}(ispec,i_Er_search):  Ion particle flux vs Er
	Er_test_vals(i_Er_search): Er values used in the search in V/m.
	pprof_info_allsurf:  plasma profile data structure.  Should be self explanatory from variable names
	Flow_vs_Er{is}(Sind,i_Er_search):  Similar to flux vs Er above but with format of Flows_ambi
	Er_ambi(is,ie):  Ambipolar Er root matrix.  Here the width is determined by the surface with the max number of roots,
			  all other surfaces use NaN when not enough roots exist.
	gamma_e(is,ie):  Same format as Er_ambi but with electron particle flux
	gamma_i(is,ie,ispec):  Same format as Er_ambi but with ion particle flux (and species index)
	QoTs_ambi, QoT_e,QoT_i -- See Gammas_ambi, gamma_e, and gamma_i, respectively
	
	