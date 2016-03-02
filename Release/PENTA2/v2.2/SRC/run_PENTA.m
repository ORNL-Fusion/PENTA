function [roa_vals,Er_ambi,gamma_e,q_e,gamma_i,q_i,J_bs,J_E_e,J_E_i,J_E_cl,FFmat,X_vec,BUprl_e,BUprl_i,gamma_e_zero_Er_root]=run_PENTA(pprof_char_opt,Smax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function is a run script for the Matlab version of PENTA. The 
%   surface information (from VMEC), plasma profile information and 
%   directory information is loaded or defined and the main PENTA function 
%   called for each surface.  This version has been written for multiple 
%   ion species.
%
%   Check all settings between the double lines of "--" below before
%   running.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   INPUTS:
%           pprof_char_opt        -option for setting the pprof_char (see
%                                   below), if not set or left empty ([])
%                                   the value below will be used
%           Smax                  -Order of the Sonine polynomial
%                                  expansion.
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
%            J_E_i(is,ie,ispec)   - Ion components of total parallel
%                                    electric current (A/m^2)
%            J_E_cl(is,ie)        - classical || electric current (A/m^2)
%            roa_vals(is)         - r/a values of the surfaces used
%            BUprl_e(is,ie)       - Electron parallel flows <B*u_||e>
%            BUprl_i(is,ie,ispec) - Ion parallel flows <B*u_||e>
%
%            FFmat{is}{ie}(:,:) - thermal transport coefficient
%               matrix relating the thermodynamic drives to the particle
%               flux and heat flux/temperature.
%            X_vec{is}{ie}(:) - Thermodynamic force vector
%
%            gamma_e_zero_Er_root - Radial electric field roots where
%                                   electron particle flux = 0
%
%   see:
%       [1] H. Sugama and S. Nishimura, PoP 9, 4637 (2002).
%       [2] H. Sugama and S. Nishimura, PoP 15, 042502 (2008).
%       My talks
%
%   6/2009 JL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NOTES for me:
%
% -better screen out
% -save variables opt
% -other Er input

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


data_path='D:/Transport/PENTA/penta_runs/PENTA_I_runs/qhs1T_parabola';
% data_path='D:/Transport/PENTA/penta_runs/PENTA_I_runs/ideal_tokamak';

run_ident='qhs1T';
% run_ident='tok';

if nargin < 1 || isempty(pprof_char_opt)
    pprof_char='a';     %single character labeling plasma profile file
%     pprof_char='1';     %single character labeling plasma profile file
else
    pprof_char=pprof_char_opt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Set up run info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_surfs=[2:2:21 40:20:200];     %VMEC surfaces for which to run PENTA, note "2" is the first usable surface.
run_surfs=[10];     %VMEC surfaces for which to run PENTA, note "2" is the first usable surface.
Er_min=-200;            %Minimum Er range value
Er_max=500;             %Maximum Er range value
input_is_Er = true;     %if true Er search range is V/cm, else e<a>*phi/kTe
plot_pprof=0;           %plot input profiles
postproc_plots=0;       %Make plots vs r/a after running (beta)
num_Er_pts=200;         %number of Er points to use in test range
plot_flux=1;            %plot fluxes vs roa for each surface
plotit_Er_find=0;       %plot root finding results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Set ion parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%These are set in the ion_params file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Set machine parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arad=.116;      %<a> of device in meters, set to [] to use value read in VMEC profile file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Define parallel electric field as a function of r/a of each surface
%%%%%    - to match eq. (C6) this must be in the form of <B*E_||>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B_Eprl=0*ones(1,length(run_surfs));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Extra run parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Maximum order in the Sonine expansion (number of terms - 1) -- Smax >= 1
if nargin < 2 || isempty(Smax)
    Smax=2;
end

%use logarithmic interpolation in x,y (cmul,efield)
log_interp=1;

%use quad to integrate or rectangular approx
use_quad=0;

%use Ji and Held formulae for friction coefficients
use_Ji=1;

%find root where gamma_e = 0
find_gamma_e_root=0;

%run in DKES limit
DKES_limit=0;

%if 1 load LMN files directly, if 0 load D* files and convert
%note if this is zero, <U**2> is also calculated internally
load_LMN_files=0;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%Sonine expansion not necessary in DKES limit
if DKES_limit==1
    Smax=1;
end

if Smax < 1
    error('Smax must be >= 1')
end

%Call surface loop and data loading for main PENTA function
[roa_vals,Er_ambi,gamma_e,q_e,gamma_i,q_i,J_bs,J_E_e,J_E_i,J_E_cl,FFmat,X_vec,BUprl_e,BUprl_i,gamma_e_zero_Er_root]=...
    surf_loop_PENTA(data_path,run_ident,pprof_char,run_surfs,Er_min,Er_max,...
    input_is_Er,plot_pprof,postproc_plots,...
    arad,B_Eprl,num_Er_pts,plot_flux,plotit_Er_find,Smax,log_interp,use_quad,use_Ji,find_gamma_e_root,DKES_limit,load_LMN_files);
