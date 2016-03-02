%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This file is a run script for the Matlab version of PENTA.  The surface
%   information (from VMEC), plasma profile information and directory
%   information is loaded or defined and the main PENTA function called for
%   each surface.  This version has been written for multiple ion species.
%
%   Check all settings between the double lines of "--" below before
%   running.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   OUTPUTS:
%           - note "is" is a surface index, "ie" is an ambipolar root index
%               "ispec" is the ion species index
%            roa_vals(is)         - r/a values for the run surfaces
%            Er_ambi(is,ie)       - Ambipolar Er roots in (V/m)
%            J_bs(is,ie)          - Bootstrap current density in (A/m^2)
%            gamma_e(is,ie)       - electron particle flux in (/m^2/s^2)
%            gamma_i(is,ie,ispec) - electron particle flux in (/m^2/s^2)
%            q_e(is,ie)           - electron, ion heat flux in (J/m^2/s)
%            J_E_e(is,ie)         - electron component of total parallel
%                                    electric current (A/m^2)
%            J_E_i(is,ie)         - electron component of total parallel
%                                    electric current (A/m^2)
%            J_E_cl(is,ie)        - classical || electric current (A/m^2)
%            roa_vals(is)         - r/a values of the surfaces used
%
%            FFmat{is}{ie}(:,:) - thermal transport coefficient
%               matrix relating the thermodynamic drives to the particle
%               flux, heat flux/temperature and the bootstrap current, as
%               defined in eq. (C5) of [1].
%            X_vec{is}{ie}(:) - Thermodynamic force vector
%
%   see:
%       [1] H. Sugama and S. Nishimura, PoP 9, 4637 (2002).
%       My talks
%
%   6/2009 JL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Set up directory info
%%%%%       -data_path is the full path to the input files
%%%%%       -pprof_char is an identifier labeling the plasma profile file,
%%%%%         i.e. "A" indicates the file is named "plasma_profilesA.dat"
%%%%%       -run_ident is an identifier labeling the coefficient and VMEC
%%%%%         data files, i.e. "hsx" indicates "profile_data_hsx",
%%%%%         "nstar_lijs_hsx_s20", etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_path='D:/Transport/PENTA/tests/7_2009_PENTA_Matlab_impurity_tests/qhs1T_parabola';
run_ident='qhs1T';

pprof_char='2';     %single character labeling plasma profile file ('z')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Set up run info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_surfs=[2:10 20:40:200];     %VMEC surfaces for which to run PENTA, note "2" is the first usable surface.
run_surfs=[10];     %VMEC surfaces for which to run PENTA, note "2" is the first usable surface.
Er_min=-200;            %Minimum Er range value
Er_max=250;             %Maximum Er range value
input_is_Er = true;     %if true Er search range is V/cm, else e<a>*phi/kTe
plot_pprof=0;           %plot input profiles
postproc_plots=0;       %Make plots vs r/a after running

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Set ion parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%These are now set in the ion_params file
% p_mass=1.672621637e-27;             %proton mass
% num_ion_species=2;                  %number of ion species
% ion_mass=p_mass*[1 12];             %row vector containing ion masses
% Z_ion=[1.0 6];                      %row vector containing ion charge nums

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Set machine parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arad=.116;      %<a> of device in meters, set to [] to use value read in VMEC profile file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Define parallel electric field as a function of r/a of each surface
%%%%%    - to match eq. (C6) this must be in the form of <B*E_||>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B_Eprl=0*ones(1,length(run_surfs));

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%Call surface loop and data loading for main PENTA function
[roa_vals,Er_ambi,gamma_e,q_e,gamma_i,q_i,J_bs,J_E_e,J_E_i,J_E_cl,FFmat,X_vec]=...
    surf_loop_PENTA(data_path,run_ident,pprof_char,run_surfs,Er_min,Er_max,...
    input_is_Er,plot_pprof,postproc_plots,...
    arad,B_Eprl);
