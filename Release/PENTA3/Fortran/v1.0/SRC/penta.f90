!-----------------------------------------------------------------------------
!
!  -THIS IS A BETA VERSION OF THE PENTA CODE-
!
! -- PENTA3 v1.0
!
!-----------------------------------------------------------------------------
!
!  PENTA calculates the neoclassical parallel flows, radial particle and 
!   energy fluxes, and the radial electric field for a surface given the
!   plasma profiles (n, T), surface geometry information (from VMEC) and
!   the monoenergetic transport coefficients (from DKES)
!
! To run,type: PENTA3 [command line arguments]
!
! where the command line arguments are:
!
! 1  unique identifier text on DKES database files
!     i.e., for D11_star_***, this would be
!     the text that replaces *** for the specific surface (ex: hsx_s10)
! 2  Er,min
! 3  Er,max  - args 2 and 3 specify the interval in Er over which
!     the search is done for the ambipolar electric field root.
!    Er,min and Er,max are dimensionless normalized electric
!     field variables, i.e.: e<a>Er,min/max/kT_e
!   - Note if the logical variable input_is_Er below is true then 
!     the search interval is set by args 2 and 3 as Er in (V/cm)
! 4  flux surface number (index from VMEC, typically)
! 5  parameter that determines whether output data should be
!     written into new files (=0) or appended to
!     existing files (=1). Typically,when a script is run for a
!     sequence of flux surfaces, this is set to 0 for the first
!     surface and to 1 for subsequent surfaces
! 6  extension of the profile_data_*** file (i.e., the *** text)
! 7  unique identifier on plasma profiles file, i.e., the ***
!    text in the filename plasma_profiles_***.dat. This allows
!    multiple profiles to be run without having to rename files
! 8  <B*E||> in T*V/m  [real] (Used if a parallel electric field is applied
!     or determined from equilibrium calculations)
! 9  Smax -- The upper limit on the summation of Laguerre polynomial terms.
!     Since the summation is from 0 to Smax the number of terms used is
!     Smax+1.  This affects the number of parallel flow moments calculated and
!     output.
!
!  Input files:
!
!   [D11,D13,D33]_star_lijs_*** - where *** is arg1 above.
!       These files contain the values of efield and cmul and the
!       corresponding normalized monoenergetic coefficients from DKES.
!
!   ion_params - A file containing the namelist "ion_params" which
!       defines the variables num_ion_species, Z_ion_init and 
!       miomp_init which are the number of ion species (integer),
!       the corresponding ion charge numbers (real, array) and 
!       the ion to proton mass ratio (real, array) for each 
!       non-electron species.
!  
!  run_params - A file containing the namelist "run_params" which contains
!    run settings.  These settings are described below.
!
!   profile_data_*** - where *** is arg1 above.  Contains geometry
!       variables, most likely from a VMEC run.
!
!   plasma_profiles_XXX.dat - A file containing the plasma profile
!       data, where XXX is arg7 above.  This file has the format:
!       row 1: number of radial points for which data is specified.
!       all other rows: r/a, ne, Te, Ti1, ni1, Ti2, ni2 ...
!       where the densities are in units of 10**12/cc and the 
!       temperatures in eV and i1, i2... are the ion species.
!
!   Utilde2_profile - Contains the quantity <U**2>, where U is the
!     PS flow function as defined by Sugama and Nishimura.  The
!     first row is the number of points, then r/a points and the
!     corresponding <U**2> value.  Note that if read_U2_file is /false.
!     then <U**2> is calculated from D11 at high collisionality
!      
!   
!  Output files:
!
!   "fluxes_vs_roa", "fluxes_vs_Er", "flows_vs_roa", "plasma_profiles_check"
!   -QQ- TO BE WRITTEN
!
!   
!-----------------------------------------------------------------------------
!   NOTES
!    - based on 
!         Sugama, Nishimura PoP 9 4637 (2002), 
!         Sugama, Nishimura PoP 15, 042502 (2008),
!         Maassberg, et al, PoP 15, 072504 (2009),
!         Taguchi, Phys. Fluids B, 4 3638 (1992)
!         Spong, PoP 12, 056114 (2005)
!         Note that the code is written from PhD thesis of J. Lore,
!          from which notation typically follows.
!    - All units are SI, except T is [eV] (unless otherwise noted)
!    - Some code used from orignal PENTA code, written by Don Spong and 
!        modified by J. Lore
!
!  7/2009-9/29/2010 JL
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!
!+ The main PENTA program
!
Program penta3
 
! Description: 
!  PENTA calculates the neoclassical parallel flows, radial particle and 
!   energy fluxes, and the radial electric field for a surface given the
!   plasma profiles (n, T), surface geometry information (from VMEC) and
!   the monoenergetic transport coefficients (from DKES).  See above for
!   input/output and references.
!
! History:
! Version   Date      Comment
! -------   ----      -------
! ...     .........   Older versions did not follow this scheme
! 1.0     5/24/2010   Original code for "PENTA3". JL 
! 
! Author(s): J. Lore 7/2009 - 9/8/2010
!            D. Spong - pre 7/2009
!
! Code Description:
!   Language:   Fortran 90 with use of Fortran 2003 COMMAND_ARGUMENT_COUNT
!   Standard:   Based on "European Standards for Writing and Documenting
!               Exchangeable Fortran 90 Code"  v1.1
! 
! External subroutines used:
!   define_friction_coeffs: Contained in file penta_subroutines.f90
!   find_Er_roots: Contained in file penta_subroutines.f90
!   form_Xvec: Contained in file penta_subroutines.f90
!   fit_coeffs: Contained in file penta_subroutines.f90

! Modules used:
Use penta_kind_mod                     ! Import rknd, iknd specifications
Use io_unit_spec, Only :             &
  iu_nl,                             & ! Ion parameter namelist i/o unit #
  iu_flux_out,                       & ! flux vs r/a i/o unit #
  iu_pprof_out,                      & ! plasma profile check i/o unit #
  iu_fvEr_out,                       & ! Flux vs Er i/o unit #
  iu_flows_out                         ! flows vs r/a i/o unit #
Use read_input_file_mod, Only :      &
  ! Imported Subroutines
  read_vmec_file,                    & ! Reads VMEC data file
  read_pprof_file,                   & ! Reads plasma profile data file
  read_dkes_star_files,              & ! Reads DKES data files
  read_Utilde2_file                    ! Reads Utilde2 file (<U**2> data)
Use vmec_var_pass , Only :            & ! VMEC variables
  ! Imported scalar variables
    arad,  roa_surf, Bsq, B0
Use pprof_pass, Only :               & ! Plasma profile variables
  ! Imported scalar variables
    ne, Te, dnedr, dTedr,            &
  ! Imported array variables (1D)
    ni, Ti, dnidr, dTidr
Use coeff_var_pass, Only :           & ! DKES coeff. variables
  ! Imported scalar variables
    num_c_D11, num_e_D11,            &
    num_c_D13, num_e_D13,            &
    num_c_D31, num_e_D31,            &
    num_c_D33, num_e_D33,            &
  ! Imported array variables (1D)
    cmul_D11,cmul_D13,cmul_D31,      & 
    cmul_D33,efield_D11,             &
    efield_D13,efield_D31,           &
    efield_D33,                      &
  ! Imported array variables (2D)
    D11_mat,D13_mat,D31_mat,         &
    D33_mat
Use phys_const, Only :               & ! Physical constants
  ! Imported parameters (scalar)
    p_mass,                          & ! proton mass
    e_mass,                          & ! electron mass
    elem_charge                        ! elementary charge
Use penta_functions_mod, Only :      & ! Functions
  ! Imported functions
     calc_flows_T,                   & ! Taguchi method of parallel flows
     calc_fluxes_T,                  & ! Taguchi method of radial part. fluxes
     calc_QoTs_T,                    & ! Taguchi method of radial energy fluxes
     calc_flows_SN,                  & ! S-N method of parallel flows
     calc_fluxes_SN,                 & ! S-N method of radial part. fluxes
     calc_QoTs_SN,                   & ! S-N method of radial energy fluxes
     calc_fluxes_MBT,                & ! M-B-T method of radial part. fluxes
     calc_QoTs_MBT,                  & ! M-B-T method of radial energy fluxes
     calc_flows_DKES,                & ! DKES method of parallel flows
     calc_fluxes_DKES,               & ! DKES method of radial part. fluxes
     calc_QoTs_DKES                    ! DKES method of radial part. fluxes

Use penta_math_routines_mod, Only :  & 
  ! Imported functions
     rlinspace                         ! Creates linearly spaced arrays

Implicit None

! Default values for variables set in run_params namelist file:
Logical ::    &
  input_is_Er    = .true.,     & ! If true, Er range is (V/cm) else e<a>Er/kT_e
  log_interp     = .true.,     & ! If true, log. interp. of DKES coeffs is used
  use_quanc8     = .false.,    & ! If false, rect. approx. to convolution used
  read_U2_file   = .true.,     & ! If false <U**2> is calculated from D11*
  Add_Spitzer_to_D33 = .true.    ! If true collisional portion of D33* is added
                                 !  else it is assumed to be included in-file  
Integer(iknd) ::    &
  num_Er_test = 20_iknd,        & ! Number of Er points in search range
  numKsteps   = 1000_iknd,     & ! Number of K points (linear) for convolution 
                                 !  (used if use_quanc8=.false.)
  kord_pprof  = 3_iknd,        & ! Spline order for plasma profile fitting 
                                 !  (and <U**2> file)
  keord       = 2_iknd,        & ! Spline order for DKES coeff fitting (efield)
  kcord       = 2_iknd         ! Spline order for DKES coeff fitting (cmul)
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

! Local parameters used to set max array sizes
Integer(iknd), Parameter :: num_roots_max = 10_iknd
Integer(iknd), Parameter :: num_ion_max = 20_iknd

! Local variables (scalar)
Integer(iknd) :: numargs          ! Number of command line args.
Integer(iknd) :: js               ! Booz_xform (VMEC) index of the surface
Integer(iknd) :: i_append         ! Set to 1 to append to output files
Integer(iknd) :: num_species      ! Total number of plasma species
Integer(iknd) :: num_ion_species  ! Number of ion species
Integer(iknd) :: Smax             ! Max index of Sonine poly. expansion
Integer(iknd) :: iocheck          ! Used to check for file open errors
Integer(iknd) :: ie               ! Er loop index
Integer(iknd) :: ind_X, ind_A     ! Thermodynamic force vector indices
Integer(iknd) :: ispec1           ! Primary species loop index
Integer(iknd) :: min_ind          ! Index for finding Er = 0
Integer(iknd) :: iroot            ! Ambipolar root index
Integer(iknd) :: num_roots        ! Number of ambipolar roots
Real(rknd)    :: Er_min, Er_max   ! Er search range min and max
Real(rknd)    :: B_Eprl           ! Input parallel electric field
Real(rknd)    :: U2               ! <U**2> Related to Pfirsch-Schlüter factor
Real(rknd)    :: G2               ! <G**2> Related to Pfirsch-Schlüter factor
Real(rknd)    :: vth_e            ! Electron thermal velocity
Real(rknd)    :: loglambda        ! Coulomb logarithm
Real(rknd)    :: Er_test, abs_Er  ! Current value of Er and |Er| used in loop
Real(rknd)    :: min_Er           ! Min value of Er_test_vals
Real(rknd)    :: eaEr_o_kTe       ! Normalized Er (used as output)
! Local variables (array)
Character(Len=100) ::           & ! Command line args
  arg1, arg2, arg3, arg4,       & 
  arg5, arg6, arg7, arg8,       &
  arg9  
Character(Len=100) :: coeff_ext   ! Identifier for DKES coeff. data files
Character(Len=100) :: run_ident   ! Identifier for VMEC data file
Character(Len=20)  :: pprof_char  ! Identifier for plasma profile file
Character(Len=100) :: fpos        ! Status for writing to files (position)
Character(Len=100) :: fstatus     ! Status for writing to files
character(Len=100) :: str_num     ! Used for converting numbers to strings
Real(rknd)         ::           &
  Z_ion_init(num_ion_max),      & ! Ion charge numbers
  miomp_init(num_ion_max)         ! Ion mass ratios
Real(rknd) ::  cmins(5),        & ! Coefficient axes limits
  cmaxes(5), emins(5), emaxes(5)

! Local allocatable arrays (1D)
Real(rknd), Allocatable ::      &
  ion_mass(:),                  & ! Ion masses
  Z_ion(:),                     & ! Ion charges
  vth_i(:),                     & ! Ion thermal velocities
  charges(:),                   & ! All species charges
  dens(:),                      & ! All species densities
  masses(:),                    & ! All species masses
  temps(:),                     & ! All species temperatures (eV)
  vths(:),                      & ! All species thermal velocities  
  dTdrs(:),                     & ! All species temperature gradients (eV/m)
  dndrs(:),                     & ! All species density gradients (1/m^2)
  Xvec(:), Avec(:),             & ! Thermodynamic force vectors
  Flows(:),                     & ! Parallel flow moment array
  Gammas(:),                    & ! Radial flux array
  xt_c_D11(:),xt_c_D13(:),      & ! Spline knots for DKES coeffs
  xt_c_D31(:),xt_c_D33(:),      & 
  xt_e_D11(:),xt_e_D13(:),      &
  xt_e_D31(:),xt_e_D33(:),      &
  xt_c_Dex(:),xt_e_Dex(:),      &
  xt_c_DUa(:),xt_e_DUa(:),      &
  xt_c_Drat(:),xt_e_Drat(:),    &
  Gamma_e_vs_Er(:),             & ! Electron particle flux vs Er
  Er_test_vals(:)                 ! Er values used in ambipolar search
    

! Local allocatable arrays (2D)
Real(rknd), Allocatable ::      &
  lmat(:,:),                    & ! Classical friction coefficients
  Dspl_D11(:,:),Dspl_D13(:,:),  & ! DKES coeff spline coeffs
  Dspl_D31(:,:),Dspl_D33(:,:),  & ! DKES coeff spline coeffs
  Dspl_Dex(:,:),                & ! Extra radial flux coefficient for SN & T
  Dspl_DUa(:,:),                & ! Extra viscosity coefficient for SN method
  Dspl_Drat(:,:),               & ! Extra coefficient for SN method
  cmesh_D11(:,:),               & ! Repeated 2D array of cmul_D11
  Gamma_i_vs_Er(:,:)              ! Ion particle flux vs Er
Real(rknd), Allocatable ::      & 
  Flows_ambi(:,:),              & ! Parallel flow moments (species,roots)
  Gammas_ambi(:,:),             & ! Radial particle fluxes (species,roots)
  QoTs_ambi(:,:)                  ! Radial energy fluxes (species,roots)

! Interface blocks QQ
Real(rknd) :: Er_roots(num_roots_max)   !QQ size!
Interface
  Subroutine find_Er_roots(gamma_e,gamma_i,Er_test_vals,Z_ion, &
    num_Er_test,num_ion_species,Er_roots,num_roots)
    Use penta_kind_mod
    Real(rknd),    Intent(in)  :: gamma_e(num_Er_test)  
    Real(rknd),    Intent(in)  :: gamma_i(num_Er_test,num_ion_species)
    Real(rknd),    Intent(in)  :: Er_test_vals(num_Er_test)
    Real(rknd),    Intent(in)  :: Z_ion(num_ion_species)
    Integer(iknd), Intent(in)  :: num_Er_test
    Integer(iknd), Intent(in)  :: num_ion_species
    Real(rknd),    Intent(out) :: Er_roots(:) 
    Integer(iknd), Intent(out) :: num_roots 
  End Subroutine find_Er_roots
End Interface

! Namelist files
Namelist / ion_params / num_ion_species, Z_ion_init, miomp_init
Namelist / run_params / input_is_Er, log_interp, use_quanc8, read_U2_file, &
  Add_Spitzer_to_D33, num_Er_test, numKsteps, kord_pprof, keord, kcord,    &
  Kmin, Kmax, epsabs, epsrel, Method

!- End of header -------------------------------------------------------------

!-----------------------------------------------------------------------------
![1.0] Initialize: read command line args, allocate space, read input files
!                  display intro text, calculate thermal velocities and
!                  log(lambda)
!-----------------------------------------------------------------------------

! Read namelist file for ion parameters
Open(iu_nl,file="ion_params",status="old",form="formatted",iostat=iocheck)
If ( iocheck /= 0 ) Then
  Write(*,*) 'Error opening ion_params namelist file'
  Stop 'Exiting: I/O Error in penta.f90'
Endif
Read(iu_nl,nml=ion_params)
Close(iu_nl)

! Read namelist file for run parameters
Open(iu_nl,file="run_params",status="old",form="formatted",iostat=iocheck)
If ( iocheck /= 0 ) Then
  Write(*,*) 'Error opening run_params namelist file'
  Stop 'Exiting: I/O Error in penta.f90'
Endif
Read(iu_nl,nml=run_params)
Close(iu_nl)

! Get command line arguments
numargs=command_argument_count()
If ( numargs /= 9 ) Then
  Write(*,*) 'Incorrect number of input arguments, see penta.f90 for details'
  Stop 'Exiting: Input arguments error in penta.f90'
Endif
Call Getarg(1, arg1)
Call Getarg(2, arg2)
Call Getarg(3, arg3)
Call Getarg(4, arg4)
Call Getarg(5, arg5)
Call Getarg(6, arg6)
Call Getarg(7, arg7)
Call Getarg(8, arg8)
Call Getarg(9, arg9)

! Store command line args
coeff_ext = Trim(Adjustl(arg1))
Read(arg2,*) Er_min
Read(arg3,*) Er_max
Read(arg4,*) js
Read(arg5,*) i_append
run_ident = Trim(Adjustl(arg6))
pprof_char = Trim(Adjustl(arg7))
Read(arg8,*) B_Eprl
Read(arg9,*) Smax
  
! Allocate variables according to number of ion species defined
num_species = num_ion_species + 1_iknd
Allocate(ni(num_ion_species),Ti(num_ion_species))          ! Ion profile info
Allocate(dnidr(num_ion_species),dTidr(num_ion_species))
Allocate(Z_ion(num_ion_species),ion_mass(num_ion_species)) ! Ion parameters
Allocate(vth_i(num_ion_species))
Allocate(charges(num_species))                             ! Parameters for
Allocate(dens(num_species))                                !  all species
Allocate(vths(num_species)) 
Allocate(masses(num_species))
Allocate(Temps(num_species))
Allocate(dTdrs(num_species))
Allocate(dndrs(num_species))
Allocate(lmat((Smax+1)*num_species,(Smax+1)*num_species))  ! Clas. fric.coeffs
Allocate(Xvec(num_species*2+1),Avec(num_species*3))        ! Thermo. Force vecs.
Allocate(Flows((Smax+1)*num_species))                      ! Prl flow moments
Allocate(Gammas(num_species))                              ! Rad fluxes
Allocate(Gamma_i_vs_Er(num_Er_test,num_ion_species))       ! Ion flux vs Er
Allocate(Gamma_e_vs_Er(num_Er_test))                       ! Electron flux vs Er
Allocate(Er_test_vals(num_Er_test))                        ! Er to loop over

! Read input files
Call read_vmec_file(js,run_ident)
Call read_pprof_file(pprof_char,num_ion_species,roa_surf,arad,kord_pprof)
Call read_dkes_star_files(coeff_ext,Add_Spitzer_to_D33)

! Allocate DKES coefficients arrays
Allocate(xt_c_D11(num_c_D11 + kcord))
Allocate(xt_c_D13(num_c_D13 + kcord))
Allocate(xt_c_D31(num_c_D31 + kcord))
Allocate(xt_c_D33(num_c_D33 + kcord))
Allocate(xt_c_Dex(num_c_D11 + kcord))
Allocate(xt_c_Drat(num_c_D11 + kcord))
Allocate(xt_c_DUa(num_c_D33 + kcord))
Allocate(xt_e_D11(num_e_D11 + keord))
Allocate(xt_e_D13(num_e_D13 + keord))
Allocate(xt_e_D31(num_e_D31 + keord))
Allocate(xt_e_D33(num_e_D33 + keord))
Allocate(xt_e_Dex(num_e_D11 + keord))
Allocate(xt_e_Drat(num_e_D11 + keord))
Allocate(xt_e_DUa(num_e_D33 + keord))
Allocate(Dspl_D11(num_c_D11,num_e_D11))
Allocate(Dspl_D13(num_c_D13,num_e_D13))
Allocate(Dspl_D31(num_c_D31,num_e_D31))
Allocate(Dspl_D33(num_c_D33,num_e_D33))
Allocate(Dspl_Dex(num_c_D11,num_e_D11))
Allocate(Dspl_Drat(num_c_D11,num_e_D11)) !QQ size 
Allocate(Dspl_DUa(num_c_D11,num_e_D11)) !QQ size 
Allocate(cmesh_D11(num_c_D11,num_e_D11))

! Optionally read file containing <U**2> info.  Else this is 
! calculated from the D11* coefficient at high nu/v and Er=0.
If ( read_U2_file ) Then
  Call read_Utilde2_file(roa_surf,U2,kord_pprof)
Else
  U2=1.5*D11_mat(num_c_D11,1)/cmul_D11(num_c_D11);
Endif

! Change Er test range to V/cm if necessary
If ( input_is_Er .EQV. .false.)  Then
  Er_min = Er_min * Te / arad
  Er_max = Er_max * Te / arad
Endif

! Assign ion parameters
Z_ion    = Z_ion_init(1:num_ion_species)
ion_mass = miomp_init(1:num_ion_species) * p_mass

! Calculate fitting parameters to the D##* coefficients
Call fit_coeffs(cmul_D11,efield_D11,num_c_D11,num_e_D11,D11_mat,log_interp, &
  kcord,keord,xt_c_D11,xt_e_D11,Dspl_D11,cmins(1),cmaxes(1),emins(1),emaxes(1))
Call fit_coeffs(cmul_D13,efield_D13,num_c_D13,num_e_D13,D13_mat,log_interp, &
  kcord,keord,xt_c_D13,xt_e_D13,Dspl_D13,cmins(2),cmaxes(2),emins(2),emaxes(2))
Call fit_coeffs(cmul_D31,efield_D31,num_c_D31,num_e_D31,D31_mat,log_interp, &
  kcord,keord,xt_c_D31,xt_e_D31,Dspl_D31,cmins(3),cmaxes(3),emins(3),emaxes(3))
Call fit_coeffs(cmul_D33,efield_D33,num_c_D33,num_e_D33,D33_mat,log_interp, &
  kcord,keord,xt_c_D33,xt_e_D33,Dspl_D33,cmins(4),cmaxes(4),emins(4),emaxes(4))

! Display run data to screen (if i_append==0)
If ( i_append == 0 ) Then
  Write(*,*) 
  Write(*,*) "Welcome to PENTA3, please note the following settings:"
  Write(*,*)
  Write(*,'(a,i3)') ' Number of ion species: ',num_ion_species
  If ( input_is_Er .EQV. .true. ) Then
    Write(*,*) 'Interpreting input range as Er (V/cm)'
  Else
    Write(*,*) 'Interpreting input range as e<a>Er/kTe'
  Endif
  If ( log_interp .EQV. .true. ) Then
    Write(*,*) 'Performing logarithmic interpolation in Er, cmul'
  Else
    Write(*,*) 'Performing linear interpolation in Er,cmul'
  Endif
  If ( use_quanc8 .EQV. .true. ) Then
    Write(*,'(a,2(a,e10.4))')                           &
      ' Using quanc8 integrator with tolerances: ',     &
      'abs: ',epsabs,' rel: ', epsrel
  Else
    Write(*,'(a,i6,a)') ' Using ',numKsteps,            &
      ' point integral approximation'
  Endif
  Write(*,'(a,2(" ",e15.4))') ' K range on convolution integral: ', &
    Kmin, Kmax
  If ( Add_Spitzer_to_D33 .EQV. .true. ) Then
    Write(*,*) 'Adding collisional (Spitzer) portion to D33 coefficient'
  Endif
  Write(*,'(a,i2)') ' Number of terms in Sonine expansion: ', Smax+1
  Select Case (Method)
    Case ('T')
      Write(*,*) 'Using Taguchi Method'
    Case ('SN')
      Write(*,*) 'Using Sugama-Nishimura Method'
    Case ('MBT')
      Write(*,*) 'Using Maassberg-Beidler-Turkin Method'
    Case ('DKES')
      Write(*,*) 'Using DKES Method'
    Case Default
      Write(*,'(3a)') ' Error: ''', Trim(Adjustl(Method)), &
        ''' is not a valid Method'
      Stop 'Error: Exiting, method select error in penta.f90 (1)'
  EndSelect
      
  Write(*,*)
  Write(*,*) " <r>/<a>","   Er roots (V/cm)"
Endif ! intro text

! Set specifiers for opening output files
If (i_append == 0) Then
  fstatus = "unknown"
  fpos = "asis"
ElseIf (i_append == 1) Then
  fstatus = "old"
  fpos = "append"
Else
  Write(*,*) 'Bad value for i_append (0 or 1 expected)'
  Stop 'Error: Exiting, i_append error in penta.f90'
EndIf

! Open output files
Open(unit=iu_flux_out, file="fluxes_vs_roa", &
  position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))
Open(unit=iu_pprof_out, file="plasma_profiles_check",  &
  position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))
Open(unit=iu_fvEr_out, file="fluxes_vs_Er",  &
  position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))
Open(unit=iu_flows_out, file="flows_vs_roa",  &
  position=Trim(Adjustl(fpos)),status=Trim(Adjustl(fstatus)))


! Write legends (if i_append = 0) for most files
If ( i_append == 0 ) Then
  ! Fluxes vs r/a
  Write(iu_flux_out,'("*",/,"r/a    Er(V/cm)    e<a>Er/kTe    ",  &
    & "Gamma_e     Q_e/T_e     Gamma_i      Q_i/T_i")')
  ! Flows vs r/a
  Write(iu_flows_out,'("*",/,"r/a   Er(V/cm)    e<a>Er/kTe    ",  &
    & " <B*u_||ke>/<B**2>    <B*u_||ki>/<B**2>")')
  ! Plasma profile check
  Write(iu_pprof_out,'("*",/,"r/a    Te    ne    dnedr   dTedr ", &
    & "  Ti   ni    dnidr    dTidr")')
EndIf

! Legend for fluxes vs Er is written for each surface
Write(iu_fvEr_out,'("*",/,"r/a   Er(V/cm)   Gamma_e   Gamma_i")')


!!$Te=549.157d0  !QQ
!!$dTedr=-2.6309d4
!!$ne=3.8418d18
!!$dnedr=-2.4528d19
!!$Ti=60.741d0
!!$dTidr=-259.45d0
!!$ni=3.8418d18
!!$dnidr=-2.4528d19
!!$write(*,*) 'assuming profiles!!!!!!!!!!!'   !QQ

! Calculate thermal velocities 
vth_i = Dsqrt(2._rknd*Ti*elem_charge/ion_mass)
vth_e = Dsqrt(2._rknd*Te*elem_charge/e_mass)

! Calculate Coulomb logarithm
If ( Te > 50._rknd ) Then
  loglambda = 25.3_rknd - 1.15_rknd*Dlog10(ne/1.e6_rknd) + 2.3_rknd*Dlog10(Te)
Else
  loglambda = 23.4_rknd - 1.15_rknd*Dlog10(ne/1.e6_rknd) + 3.45_rknd*Dlog10(Te)
Endif

! Assign arrays of parameters for all species (charge, mass, n, T, v_th, dTdr)
charges=elem_charge*(/-1._rknd,Z_ion/)
dens=(/ne, ni/)
masses=(/e_mass, ion_mass/)
Temps=(/Te,Ti/)
vths=(/vth_e,vth_i/)
dTdrs=(/dTedr,dTidr/)
dndrs=(/dnedr,dnidr/)

! Define matrix of friction coefficients (lmat)
Call define_friction_coeffs(masses,charges,vths,Temps,dens,loglambda, &
                            num_species,Smax,lmat)

! Fit radial transport coefficients specific to different methods
SelectCase (Method)
  Case ('T')
    cmesh_D11 = Spread(cmul_D11,2,num_e_D11)
    ! Calculate the D11 coefficient minus the P-S contribution
    ! Also, do not allow for negative coefficients QQ
    Call fit_coeffs(cmul_D11,efield_D11,num_c_D11,num_e_D11, &
      Max(D11_mat-(2._rknd/3._rknd)*cmesh_D11*U2,0._rknd), &
      log_interp,kcord,keord,xt_c_Dex,xt_e_Dex,Dspl_Dex,     &
      cmins(5),cmaxes(5),emins(5),emaxes(5))
  Case ('SN')
    ! Calculate fits to D31*/D33*  ! QQ add check for all cmul, efield, etc equal
    Call fit_coeffs(cmul_D33,efield_D33,num_c_D33,num_e_D33, &
      D31_mat/D33_mat, &
      log_interp,kcord,keord,xt_c_Drat,xt_e_Drat,Dspl_Drat,     &
      cmins(5),cmaxes(5),emins(5),emaxes(5))

    cmesh_D11 = Spread(cmul_D11,2,num_e_D11)
    ! Calculate coefficient for Ua term
    Call fit_coeffs(cmul_D33,efield_D33,num_c_D33,num_e_D33, &
      (2._rknd*Bsq/(3._rknd*cmesh_D11) - D33_mat)*cmesh_D11/D33_mat, &
      log_interp,kcord,keord,xt_c_DUa,xt_e_DUa,Dspl_DUa,     &
      cmins(5),cmaxes(5),emins(5),emaxes(5))

    ! Calculate coefficient for radial flux 
    ! Also, do not allow for negative coefficients QQ
    Call fit_coeffs(cmul_D11,efield_D11,num_c_D11,num_e_D11, &
      Max(D11_mat-(2._rknd/3._rknd)*cmesh_D11*U2+D31_mat**2/D33_mat/B0**2, &
      0._rknd),log_interp,kcord,keord,xt_c_Dex,xt_e_Dex,Dspl_Dex,     &
      cmins(5),cmaxes(5),emins(5),emaxes(5))

  Case ('MBT')
    G2 = (Bsq/2._rknd)*U2
  Case ('DKES')
  Case Default
    Write(*,'(3a)') ' Error: ''', Trim(Adjustl(Method)), &
      ''' is not a valid Method'
    Stop 'Error: Exiting, method select error in penta.f90 (2)'
EndSelect

!      !calculate Spitzer conductivity QQ
!      sigma_S=ne**2_iknd*elem_charge**2_iknd *  
!     1   lmat(2,2)/( lmat(1,1)*lmat(2,2) + (lmat(1,2))**2_iknd)

!-----------------------------------------------------------------------------
![2.0] Efield loop: Loop over test values of Er and calculate flows and fluxes
!-----------------------------------------------------------------------------

! Define array of Er values to test [V/m]
Er_test_vals = rlinspace(Er_min,Er_max,num_Er_test)*100._rknd

! Check for Er=0, doesn't work for log interpolation
min_Er = Minval(Dabs(Er_test_vals),DIM=1) 
If ((log_interp .EQV. .true. ) .AND. ( Dabs(min_Er) <= elem_charge ))  Then
  min_ind = Minloc(Dabs(Er_test_vals),DIM=1)
  If ( min_ind == Num_Er_test ) Then 
    Er_test_vals(min_ind) = Er_test_vals(min_ind - 1)/2._rknd
  Else
    Er_test_vals(min_ind) = Er_test_vals(min_ind + 1)/2._rknd
  EndIf
  Write(*,'(a,i4,a,f10.3)') 'Cannot use Er=0 with log_interp, using Er(',  &
     min_ind, ') = ', Er_test_vals(min_ind)
EndIf

! Loop over Er to get fluxes as a function of Er
Do ie = 1,num_Er_test

  Er_test = Er_test_vals(ie)
  abs_Er = Dabs(Er_test)

  ! Form thermodynamic force vector (Xvec)
  Call form_Xvec(Er_test,Z_ion,B_Eprl,num_ion_species,Xvec)

  ! Form alternate thermodynamic force vector (Avec)
  Do ispec1 = 1,num_species
    ind_X = (ispec1-1)*2 + 1
    ind_A = (ispec1-1)*3 + 1

    Avec(ind_A)   = -Xvec(ind_X) / (Temps(ispec1)*elem_charge) &
         - 2.5_rknd*dTdrs(ispec1)/Temps(ispec1)
    Avec(ind_A+1)   = -Xvec(ind_X+1) / (Temps(ispec1)*elem_charge)
    Avec(ind_A+2)   = Xvec(num_species*2+1)*charges(ispec1) &
         * B0/(Temps(ispec1)*elem_charge*Dsqrt(Bsq))
  Enddo

  ! Select the appropriate algorithm and calculate the flows and fluxes
  SelectCase (Method)

    Case ('T')
      ! Calculate array of parallel flow moments
      Flows = calc_flows_T(num_species,Smax,abs_Er,Temps,dens,vths,charges,  &
         masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,      &
         cmins,cmaxes,emins,emaxes,xt_c_D31,xt_e_D31,Dspl_D31,xt_c_D33,      &
         xt_e_D33,Dspl_D33,num_c_D31,num_e_D31,num_c_D33,num_e_D33,kcord,    &
         keord,Avec,Bsq,lmat)

      Gammas = calc_fluxes_T(num_species,Smax,abs_Er,Temps,dens,vths,charges, &
        masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,           &
        cmins,cmaxes,emins,emaxes,xt_c_D31,xt_e_D31,Dspl_D31,xt_c_Dex,        &
        xt_e_Dex,Dspl_Dex,num_c_D31,num_e_D31,num_c_D11,num_e_D11,kcord,      &
        keord,Avec,lmat,Flows,U2)                                                

    Case ('SN')
      Flows = calc_flows_SN(num_species,Smax,abs_Er,Temps,dens,vths,charges, &
         masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,      &
         cmins,cmaxes,emins,emaxes,xt_c_Drat,xt_e_Drat,Dspl_Drat,xt_c_DUa,   &
         xt_e_DUa,Dspl_DUa,num_c_D31,num_e_D31,num_c_D33,num_e_D33,kcord,    &
         keord,Avec,lmat)                                                

      Gammas = calc_fluxes_SN(num_species,Smax,abs_Er,Temps,dens,vths,charges, &
        masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,            &
        cmins,cmaxes,emins,emaxes,xt_c_Drat,xt_e_Drat,Dspl_Drat,xt_c_Dex,      &
        xt_e_Dex,Dspl_Dex,num_c_D31,num_e_D31,num_c_D11,num_e_D11,kcord,       &
        keord,Avec,Bsq,B0,lmat,Flows,U2)      
    
    Case ('MBT')
      ! Flow methods are the same for T and MBT
      Flows = calc_flows_T(num_species,Smax,abs_Er,Temps,dens,vths,charges,  &
         masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,      &
         cmins,cmaxes,emins,emaxes,xt_c_D31,xt_e_D31,Dspl_D31,xt_c_D33,      &
         xt_e_D33,Dspl_D33,num_c_D31,num_e_D31,num_c_D33,num_e_D33,kcord,    &
         keord,Avec,Bsq,lmat)

      Gammas = calc_fluxes_MBT(num_species,Smax,abs_Er,Temps,dens,vths,charges, &
        masses,dTdrs,dndrs,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,  &
        cmins,cmaxes,emins,emaxes,xt_c_D11,xt_e_D11,Dspl_D11,xt_c_D31,        &
        xt_e_D31,Dspl_D31,num_c_D11,num_e_D11,num_c_D31,num_e_D31,kcord,      &
        keord,Avec,lmat,Flows,G2,B0)   

    Case ('DKES')
      Flows = calc_flows_DKES(num_species,Smax,abs_Er,Temps,dens,vths,charges,  &
         masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,      &
         cmins,cmaxes,emins,emaxes,xt_c_D31,xt_e_D31,Dspl_D31,xt_c_D33,      &
         xt_e_D33,Dspl_D33,num_c_D31,num_e_D31,num_c_D33,num_e_D33,kcord,    &
         keord,Avec)
      Gammas = calc_fluxes_DKES(num_species,abs_Er,Temps,dens,vths,charges, &
        masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,  &
        cmins,cmaxes,emins,emaxes,xt_c_D11,xt_e_D11,Dspl_D11,xt_c_D31,        &
        xt_e_D31,Dspl_D31,num_c_D11,num_e_D11,num_c_D31,num_e_D31,kcord,      &
        keord,Avec,B0)   

    Case Default
      Write(*,'(3a)') ' Error: ''', Trim(Adjustl(Method)), &
        ''' is not a valid Method'
      Stop 'Error: Exiting, method select error in penta.f90 (3)'
  EndSelect

  Gamma_e_vs_Er(ie)   = Gammas(1)
  Gamma_i_vs_Er(ie,:) = Gammas(2:num_species)

  ! Write fluxes vs Er
  Write(str_num,*) num_ion_species + 2  ! Convert num to string
  Write(iu_fvEr_out,'(f7.4,' // trim(adjustl(str_num)) // '(" ",e15.7))') &
    roa_surf,Er_test/100._rknd,Gamma_e_vs_Er(ie),Gamma_i_vs_Er(ie,:)

Enddo !efield loop


!-----------------------------------------------------------------------------
![3.0] Find ambipolar roots and evaluate quantities at these roots.
!-----------------------------------------------------------------------------

! Check for only one Er test value -- this is then used to evaluate the ambipolar fluxes QQ
!If ( num_Er_test  == 1 ) Then
!  Er_roots = Er_test_vals

! Find the ambipolar root(s) from gamma_e = sum(Z*gamma_i)
Call find_Er_roots(gamma_e_vs_Er,gamma_i_vs_Er,Er_test_vals,Z_ion, &
  num_Er_test,num_ion_species,Er_roots,num_roots)

! Allocate arrays according to number of ambipolar roots
Allocate(Flows_ambi((Smax+1)*num_species,num_roots)) ! Parallel flow moments
Allocate(Gammas_ambi(num_species,num_roots))         ! Rad. particle fluxes
Allocate(QoTs_ambi(num_species,num_roots))           ! Rad. energy fluxes

! Evaluate fluxes and flows at the ambipolar Er
Do iroot = 1_iknd, num_roots

  Er_test = Er_roots(iroot)
  abs_Er = Dabs(Er_test)

  ! Form thermodynamic force vector (Xvec)
  Call form_Xvec(Er_test,Z_ion,B_Eprl,num_ion_species,Xvec)

  ! Form alternate thermodynamic force vector (Avec)
  Do ispec1 = 1,num_species
    ind_X = (ispec1-1)*2 + 1
    ind_A = (ispec1-1)*3 + 1

    Avec(ind_A)   = -Xvec(ind_X) / (Temps(ispec1)*elem_charge) &
         - 2.5_rknd*dTdrs(ispec1)/Temps(ispec1)
    Avec(ind_A+1)   = -Xvec(ind_X+1) / (Temps(ispec1)*elem_charge)
    Avec(ind_A+2)   = Xvec(num_species*2+1)*charges(ispec1) &
         * B0/(Temps(ispec1)*elem_charge*Dsqrt(Bsq))
  Enddo

  ! Select the appropriate algorithm and calculate the flows and fluxes QQ formatting below!
  SelectCase (Method)

    Case ('T')
      ! Calculate array of parallel flow moments
      Flows_ambi(:,iroot) = calc_flows_T(num_species,Smax,abs_Er,Temps,dens, &
        vths,charges,masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,     &
        log_interp,cmins,cmaxes,emins,emaxes,xt_c_D31,xt_e_D31,Dspl_D31,     &
        xt_c_D33,xt_e_D33,Dspl_D33,num_c_D31,num_e_D31,num_c_D33,num_e_D33,  &
        kcord,keord,Avec,Bsq,lmat)
      ! Calculate array of radial particle fluxes
      Gammas_ambi(:,iroot) = calc_fluxes_T(num_species,Smax,abs_Er,Temps,dens, &
        vths,charges,masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,          &
        log_interp,cmins,cmaxes,emins,emaxes,xt_c_D31,xt_e_D31,Dspl_D31,       &
        xt_c_Dex,xt_e_Dex,Dspl_Dex,num_c_D31,num_e_D31,num_c_D11,num_e_D11,    &
        kcord,keord,Avec,lmat,Flows,U2)  
      ! Calculate array of radial energy fluxes
      QoTs_ambi(:,iroot) = calc_QoTs_T(num_species,Smax,abs_Er,Temps,dens, &
        vths,charges,masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,          &
        log_interp,cmins,cmaxes,emins,emaxes,xt_c_D31,xt_e_D31,Dspl_D31,       &
        xt_c_Dex,xt_e_Dex,Dspl_Dex,num_c_D31,num_e_D31,num_c_D11,num_e_D11,    &
        kcord,keord,Avec,lmat,Flows,U2)                                                

    Case ('SN')
      ! Calculate array of parallel flow moments
      Flows_ambi(:,iroot) = calc_flows_SN(num_species,Smax,abs_Er,Temps,dens, &
         vths,charges, &
         masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,      &
         cmins,cmaxes,emins,emaxes,xt_c_Drat,xt_e_Drat,Dspl_Drat,xt_c_DUa,   &
         xt_e_DUa,Dspl_DUa,num_c_D31,num_e_D31,num_c_D33,num_e_D33,kcord,    &
         keord,Avec,lmat)                                                
      ! Calculate array of radial particle fluxes
      Gammas_ambi(:,iroot) = calc_fluxes_SN(num_species,Smax,abs_Er,Temps,dens,&
        vths,charges, &
        masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,            &
        cmins,cmaxes,emins,emaxes,xt_c_Drat,xt_e_Drat,Dspl_Drat,xt_c_Dex,      &
        xt_e_Dex,Dspl_Dex,num_c_D31,num_e_D31,num_c_D11,num_e_D11,kcord,       &
        keord,Avec,Bsq,B0,lmat,Flows,U2)      
      ! Calculate array of radial energy fluxes
      QoTs_ambi(:,iroot) = calc_QoTs_SN(num_species,Smax,abs_Er,Temps,dens,&
        vths,charges, &
        masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,            &
        cmins,cmaxes,emins,emaxes,xt_c_Drat,xt_e_Drat,Dspl_Drat,xt_c_Dex,      &
        xt_e_Dex,Dspl_Dex,num_c_D31,num_e_D31,num_c_D11,num_e_D11,kcord,       &
        keord,Avec,Bsq,B0,lmat,Flows,U2)  

    Case ('MBT')
      ! Calculate array of parallel flow moments 
        ! Note: Flow methods are the same for T and MBT
      Flows_ambi(:,iroot) = calc_flows_T(num_species,Smax,abs_Er,Temps,dens, &
        vths,charges,masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,&
        cmins,cmaxes,emins,emaxes,xt_c_D31,xt_e_D31,Dspl_D31,xt_c_D33,      &
        xt_e_D33,Dspl_D33,num_c_D31,num_e_D31,num_c_D33,num_e_D33,kcord,    &
        keord,Avec,Bsq,lmat)
      ! Calculate array of radial particle fluxes
      Gammas_ambi(:,iroot) = calc_fluxes_MBT(num_species,Smax,abs_Er,Temps, &
        dens,vths,charges, &
        masses,dTdrs,dndrs,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,&
        cmins,cmaxes,emins,emaxes,xt_c_D11,xt_e_D11,Dspl_D11,xt_c_D31,         &
        xt_e_D31,Dspl_D31,num_c_D11,num_e_D11,num_c_D31,num_e_D31,kcord,       &
        keord,Avec,lmat,Flows,G2,B0)   
      ! Calculate array of radial energy fluxes
      QoTs_ambi(:,iroot) = calc_QoTs_MBT(num_species,Smax,abs_Er,Temps, &
        dens,vths,charges, &
        masses,dTdrs,dndrs,loglambda,use_quanc8,Kmin,Kmax,numKsteps,log_interp,&
        cmins,cmaxes,emins,emaxes,xt_c_D11,xt_e_D11,Dspl_D11,xt_c_D31,         &
        xt_e_D31,Dspl_D31,num_c_D11,num_e_D11,num_c_D31,num_e_D31,kcord,       &
        keord,Avec,lmat,Flows,G2,B0)   

    Case ('DKES')
      ! Calculate array of parallel flow moments 
      Flows_ambi(:,iroot) = calc_flows_DKES(num_species,Smax,abs_Er,Temps, &
       dens,vths,charges,  &
         masses,loglambda,B0,use_quanc8,Kmin,Kmax,numKsteps,log_interp,      &
         cmins,cmaxes,emins,emaxes,xt_c_D31,xt_e_D31,Dspl_D31,xt_c_D33,      &
         xt_e_D33,Dspl_D33,num_c_D31,num_e_D31,num_c_D33,num_e_D33,kcord,    &
         keord,Avec)
      ! Calculate array of radial particle fluxes
      Gammas_ambi(:,iroot) = calc_fluxes_DKES(num_species,abs_Er,Temps,dens, &
        vths,charges,masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,        &
        log_interp,cmins,cmaxes,emins,emaxes,xt_c_D11,xt_e_D11,Dspl_D11,     &
        xt_c_D31,xt_e_D31,Dspl_D31,num_c_D11,num_e_D11,num_c_D31,num_e_D31,  &
        kcord,keord,Avec,B0)  
      ! Calculate array of radial energy fluxes
      QoTs_ambi(:,iroot) = calc_QoTs_DKES(num_species,abs_Er,Temps,dens, &
        vths,charges,masses,loglambda,use_quanc8,Kmin,Kmax,numKsteps,        &
        log_interp,cmins,cmaxes,emins,emaxes,xt_c_D11,xt_e_D11,Dspl_D11,     &
        xt_c_D31,xt_e_D31,Dspl_D31,num_c_D11,num_e_D11,num_c_D31,num_e_D31,  &
        kcord,keord,Avec,B0)  

    Case Default
      Write(*,'(3a)') ' Error: ''', Trim(Adjustl(Method)), &
        ''' is not a valid Method'
      Stop 'Error: Exiting, method select error in penta.f90 (4)'
  EndSelect

EndDo ! Ambipolar root loop

!-----------------------------------------------------------------------------
![4.0] Write output files
!-----------------------------------------------------------------------------

! Loop over ambipolar Er for writing output files
Do iroot = 1_iknd, num_roots

  Er_test = Er_roots(iroot)
  eaEr_o_kTe = arad*Er_test/Te

  ! Write fluxes to file "fluxes_vs_roa"
  Write(str_num,*) 2*num_species + 2
  Write(iu_flux_out,'(f7.3,' // Trim(Adjustl(str_num)) // '(" ",e15.7))') &
    roa_surf,Er_test/100._rknd,eaEr_o_kTe,Gammas_ambi(1,iroot),  &
    QoTs_ambi(1,iroot),Gammas_ambi(2:num_species,iroot),  &
    QoTs_ambi(2:num_species,iroot)

  ! Write flows to file "flows_vs_roa"
  Write(str_num,*) (Smax+1)*num_species + 2
  Write(iu_flows_out,'(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.7))') &
    roa_surf,Er_test/100._rknd,eaEr_o_kTe,Flows_ambi(:,iroot)

EndDo ! Ambipolar root loop

! Write plasma profile information to "plasma_profiles_check"
Write(str_num,*) 4*num_species
Write(iu_pprof_out,'(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.7))') & 
  roa_surf,Te,ne,dnedr,dTedr,Ti,ni,dnidr,dTidr

!      !calculate classical parallel electric current  QQ
!      J_E_cl=sigma_S*Xvec(num_species*2+1)

!        !calculate total parallel electric current
!        J_E_tot_parts=dens*charges*B_uprl*dsqrt(Bsq)
!        J_E_tot=sum(J_E_tot_parts)

!        !calculate Bootstrap current
!        J_bs_ambi(iroot)=J_E_tot-J_E_cl

!        

! QQ write file with number of roots per surface!

! Write screen output
write(str_num,*) num_roots
write(*,'(f7.3,' // trim(adjustl(str_num)) // '(" ",e15.4))') & 
  roa_surf,er_roots(1:num_roots)/100._rknd

!-----------------------------------------------------------------------------
![4.0] Cleanup -- Deallocate variables and close files
!-----------------------------------------------------------------------------
! Deallocate variables
Deallocate(ni,Ti,dnidr,dTidr)                     ! Ion profile info
Deallocate(Z_ion,ion_mass)                        ! Ion parameters
Deallocate(cmul_D11,efield_D11)                   ! DKES coeff info
Deallocate(cmul_D13,efield_D13)
Deallocate(cmul_D31,efield_D31)
Deallocate(cmul_D33,efield_D33)
Deallocate(D11_mat,D13_mat,D31_mat,D33_mat)
Deallocate(Temps,masses,vths,charges,dens,dTdrs,dndrs)  ! All species parameters
Deallocate(lmat)                                  ! Class. fric. coeffs.
Deallocate(Xvec,Avec)                             ! Thermo. Force Vecs
Deallocate(Flows)                                 ! Prl flows
Deallocate(Gammas)                                ! Rad fluxes
Deallocate(xt_c_D11,xt_c_D13,xt_c_D31,xt_c_D33)   ! DKES spline fitting info
Deallocate(xt_c_Dex,xt_c_Drat,xt_c_DUa)
Deallocate(xt_e_D11,xt_e_D13,xt_e_D31,xt_e_D33)
Deallocate(xt_e_Dex)
Deallocate(xt_e_Drat)
Deallocate(xt_e_DUa)
Deallocate(Dspl_D11)
Deallocate(Dspl_D13)
Deallocate(Dspl_D31)
Deallocate(Dspl_D33)
Deallocate(Dspl_Dex)
Deallocate(Dspl_Drat)
Deallocate(Dspl_DUa)
Deallocate(cmesh_D11)     ! Spread array of cmul data
Deallocate(Gamma_i_vs_Er) ! Ion particle flux vs Er
Deallocate(Gamma_e_vs_Er) ! Electron particle flux vs Er
Deallocate(Gammas_ambi)   ! Ambipolar quantities
Deallocate(QoTs_ambi)
Deallocate(Flows_ambi)
Deallocate(Er_test_vals)  ! Er to loop over

! Close output files
Close(iu_flux_out)
Close(iu_pprof_out)
Close(iu_fvEr_out)
Close(iu_flows_out)

End program penta3

