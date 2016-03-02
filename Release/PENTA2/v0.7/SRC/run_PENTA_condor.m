function [roa_vals,Er_ambi,gamma_e,q_e,gamma_i,q_i,J_bs,J_E_e,J_E_i,J_E_cl,FFmat,X_vec]=run_PENTA_condor(pprof_char_opt)
%See file run_PENTA for details.  Modified to run on UW Condor.
%
%7/17/2009 JL

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

if nargin < 1 
    pprof_char='2';     %single character labeling plasma profile file
else
    pprof_char=pprof_char_opt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Set up run info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_surfs=[2:10 20:40:200];     %VMEC surfaces for which to run PENTA, note "2" is the first usable surface.
Er_min=-200;            %Minimum Er range value
Er_max=250;             %Maximum Er range value
input_is_Er = true;     %if true Er search range is V/cm, else e<a>*phi/kTe
plot_pprof=0;           %plot input profiles
postproc_plots=0;       %Make plots vs r/a after running

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Set ion parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%These must now be set in the ion_params file

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
