%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function is a run script for the Matlab version of PENTA. 
%
%   Note that many legacy settings from the older versions of PENTA have 
%   been removed.
%
%   Check all settings between the double lines of "--" below before
%   running.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   For instructions see the DOC directory
%
%   References:
%       [1] H. Sugama and S. Nishimura, PoP 9, 4637 (2002).
%       [2] H. Sugama and S. Nishimura, PoP 15, 042502 (2008).
%       [3] H. Maassberg, C.D. Beidler, and Y. Turkin, PoP 16, 072504
%       (2009).
%       My talks, Chapter 3 of my Thesis.
%
%   v3.1
%
%   6/2009-6/8/2010 JL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Set up directory info
%%%%%       - data_path is the full path to the input files for this run
%%%%%       - pprof_char is an identifier labeling the plasma profile file,
%%%%%          e.g. "test" indicates the file is named 
%%%%%            "plasma_profiles_test.dat"
%%%%%       - run_ident is an identifier labeling VMEC data file and the
%%%%%           files containing the DKES coefficients.  E.g. "hsx" 
%%%%%           indicates "profile_data_hsx", D11_star_hsx_s10, etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_path='D:\Transport\PENTA\Release\PENTA3\v1.1\TESTS\qhs1T';

run_ident='qhs1T';

pprof_ident='a';  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Set up run info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_surfs=[10,20];      %Booz_xform (VMEC) surfaces for which to run PENTA, 
                        %  note "2" is the first usable surface.
Er_min=-100;            %Minimum Er search range value (in V/cm)
Er_max=550;             %Maximum Er search range value (in V/cm)
num_Er_pts=50;          %Number of Er points to use in test range
Method='DKES';               %'SN': Use S-N method, 'MBT': Use M-B-T method, 
                        %'T': Use Taguchi method, 'DKES': DKES
Smax=2;                 %Maximum order in the Sonine expansion 
                        %  (number of terms - 1) -- Smax >= 1
plot_pprof=0;           %Plot input plasma profiles
plot_flux=1;            %Plot particle flux vs Er for each surface
plot_flow=0;            %Plot particle flow vs Er for each surface
plotit_Er_find=0;       %Plot root finding results  

%There is probably no reason to change these settings
Add_Spitzer_to_D33=1;   %Account for "Spitzer" portion of D33
log_interp=1;           %use logarithmic interpolation in x,y (cmul,efield)
use_Ji=1;               %use Ji and Held formulae for friction coefficients
load_U2_file=0;         %Load U2 data from Utilde2_profile.dat file, else
                        %  it is calculated from D11 (BETA)
use_quad=0;             %use quad to integrate or rectangular approx                        
calc_inductive_flow=0;  %If this flag is set, PENTA is run normally until
                        %  the fluxes and flows are evaluated at the ambi-
                        %  polar roots.  At this point the density and
                        %  temperature gradients are set to zero and the
                        %  quantity <BE_||>/<B**2>^1/2 is set to 1.  The
                        %  flows (and fluxes) returned are then the
                        %  inductive flows for a unity valued excitation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Set ion parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%These are set in the ion_params file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Set machine parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arad=.116;      %<a> of device in meters, set to [] to use 
                %  value set in VMEC profile file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Define parallel electric field as a function of r/a of each surface
%%%%%    - this must be in the form of <dot(B,E_||)>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B_Eprl=0*ones(1,length(run_surfs));

%NOT FULLY IMPLEMENTED!!!
find_gamma_e_root=0;    %find root where gamma_e = 0
input_is_Er = true;     %If true Er search range is V/cm, else e<a>*phi/kTe
use_fixed_Er_profile=0;
fixed_Er_profile_path='';

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%Set up a structure for the path info
path_info.data_path=data_path;
path_info.run_ident=run_ident;
path_info.pprof_ident=pprof_ident;

%Set up flags structure
flags.input_is_Er=input_is_Er;
flags.plot_pprof=plot_pprof;
flags.plot_flux=plot_flux;
flags.plot_flow=plot_flow;
flags.plotit_Er_find=plotit_Er_find;
flags.log_interp=log_interp;
flags.use_quad=use_quad;
flags.use_Ji=use_Ji;
flags.find_gamma_e_root=find_gamma_e_root;
flags.use_fixed_Er_profile=use_fixed_Er_profile;
flags.Method=Method;
flags.Add_Spitzer_to_D33=Add_Spitzer_to_D33;
flags.load_U2_file=load_U2_file;
flags.calc_inductive_flow=calc_inductive_flow;

if Smax < 1
    error('Smax must be >= 1')
end

%Set spline knots for fixed Er profile (BETA)
if use_fixed_Er_profile
    pp_fixed_Er_profile=pp_fixed_Er_data;
else
    pp_fixed_Er_profile='';
end

%start a counter
Tstart=tic; 

%Call surface loop and data loading for main function
[roa_vals,Er_roots,Gammas_ambi,Flows_ambi,gamma_e_vs_Er,....
    gamma_i_vs_Er,Er_test_vals,pprof_info_allsurf,Flow_vs_Er,...
    Er_ambi,gamma_e,gamma_i,QoTs_ambi,QoT_e,QoT_i]=...
    surf_loop_PENTA(path_info,run_surfs,Er_min,...
    Er_max,arad,B_Eprl,num_Er_pts,Smax,pp_fixed_Er_profile,flags);

disp(['PENTA2 run finished in ' num2str(toc(Tstart)) ' seconds.'])
