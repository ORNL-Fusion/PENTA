function [roa_vals,Er_roots,Gammas_ambi,Flows_ambi,gamma_e_vs_Er,...
    gamma_i_vs_Er,Er_test_vals,pprof_info_allsurf,Flow_vs_Er,...
    Er_ambi,gamma_e,gamma_i,QoTs_ambi,QoT_e,QoT_i]=...
    surf_loop_PENTA(path_info,run_surfs,Er_min,Er_max,...
    arad,B_Eprl,num_Er_pts,Smax,pp_fixed_Er_profile,flags)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function is the surface loop for calling the main program.  
%Data is loaded by calling "load_input_files", run text is displayed 
%to the screen and some output postprocessing is performed.
%6/2009-5/19/2010 JL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Constants
p_mass=1.672621637e-27;             %proton mass

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load input files
%   VMEC data file, ion parameter file, plasma profile file, and <U**2>
%   file (optional).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[VMEC_info,ion_params,pprof_info_allsurf,U2_data] = ...
    load_input_files(flags,path_info,arad);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Intro text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('---------------------------------------------------------------');
fprintf('\n\nWelcome to PENTA, please note the following settings:\n\n');
fprintf('Data path: %s\n',path_info.data_path);
fprintf('Device <R0>, <a>: %5.3f  %5.3f \n',VMEC_info.Rmajor_vmec,...
    VMEC_info.arad);
fprintf('Number of ion species: %3i \n',ion_params.num_ion_species);
fprintf('Sonine order: %3i \n',Smax);
for ispec=1:ion_params.num_ion_species
    fprintf('Ion species %3i: mi/mp=%5.3f  Zi=%5.3f \n',...
        ispec,ion_params.ion_mass(ispec)./p_mass, ion_params.Z_ion(ispec));
end
if flags.input_is_Er
    fprintf('Interpreting search range as Er\n');
else
    fprintf('Interpreting search range as e<a>*Phi/kTe\n');
    error('Not yet implemented (e<a>Phi/kTe)')
end
if flags.Add_Spitzer_to_D33
    fprintf('Including Spitzer portion in D_{33}*\n');
end
switch flags.Method
    case 'T',
        fprintf('Using T Method v2\n\n');        
    case 'SN',
        fprintf('Using S-N Method v2\n\n');        
    case 'DKES',
        fprintf('Using DKES Method\n\n');
    case 'MBT',
        fprintf('Using M-B-T Method v2\n\n');
    otherwise,
        error('Unknown Method')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numsurfs=length(run_surfs);
num_ion_species=ion_params.num_ion_species;
[init(1:numsurfs).surf]=deal(NaN);             %Matlab is a bit tricky with
[gamma_e_vs_Er{1:numsurfs}]=deal(init.surf);   %initializing cell arrays
[gamma_i_vs_Er{1:numsurfs}]=deal(init.surf);   %but this method works well.
[Flows_ambi{1:numsurfs}]=deal(init.surf);
[Gammas_ambi{1:numsurfs}]=deal(init.surf);
[Er_roots{1:numsurfs}]=deal(init.surf);
[Flow_vs_Er{1:numsurfs}]=deal(init.surf);
[QoTs_ambi{1:numsurfs}]=deal(init.surf);
roa_vals=NaN(1,numsurfs);
numroots_surf=NaN(1,numsurfs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop over surfaces and call main PENTA function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for isurf = 1:numsurfs
    %get index of vmec surface
    vmec_index=find(VMEC_info.js_vmec==run_surfs(isurf));
    
    %get r/a for this surface
    roa_surf=VMEC_info.roa_vmec(vmec_index);
    
    %get local plasma profile info and compile into structure for this
    %surface
    pprof_info.Te_surf=ppval(pprof_info_allsurf.pp_Te,roa_surf);
    pprof_info.ne_surf=ppval(pprof_info_allsurf.pp_ne,roa_surf);
    pprof_info.dnedr_surf=...
        ppval(pprof_info_allsurf.pp_ne_prime,roa_surf)/arad;
    pprof_info.dTedr_surf=...
        ppval(pprof_info_allsurf.pp_Te_prime,roa_surf)/arad;
    pprof_info.roa_surf=roa_surf;    
    for ispec=1:num_ion_species
        pprof_info.Ti_surf(ispec)=...
            ppval(pprof_info_allsurf.pp_Ti(ispec,:),roa_surf);
        pprof_info.ni_surf(ispec)=...
            ppval(pprof_info_allsurf.pp_ni(ispec,:),roa_surf);
        pprof_info.dTidr_surf(ispec)=...
            ppval(pprof_info_allsurf.pp_Ti_prime(ispec,:),roa_surf)/arad;
        pprof_info.dnidr_surf(ispec)=...
            ppval(pprof_info_allsurf.pp_ni_prime(ispec,:),roa_surf)/arad;    
    end
     
    %get local VMEC info and compile into structure
    surf_VMEC_info.chip=VMEC_info.chip_vmec(vmec_index);
    surf_VMEC_info.psip=VMEC_info.psip_vmec(vmec_index);
    surf_VMEC_info.Bsq=VMEC_info.Bsq_vmec(vmec_index);
    surf_VMEC_info.vp=VMEC_info.vp_vmec(vmec_index);
    surf_VMEC_info.Btheta=VMEC_info.Btheta_vmec(vmec_index);
    surf_VMEC_info.Bzeta=VMEC_info.Bzeta_vmec(vmec_index);
    surf_VMEC_info.iota=VMEC_info.iota_vmec(vmec_index);
    surf_VMEC_info.vmec_jval=VMEC_info.js_vmec(vmec_index);
    surf_VMEC_info.B0=VMEC_info.B0_vmec(vmec_index);
    
    %define local U2
    if flags.load_U2_file
        U2_surf=interp1(U2_data(:,1),U2_data(:,2),roa_surf);
    else
        U2_surf=[];  %U2 will be calculated from D11
    end
    
    %define search range
    if flags.use_fixed_Er_profile
        Er_range(1)=ppval(pp_fixed_Er_profile,roa_surf);
        Er_range(2)=ppval(pp_fixed_Er_profile,roa_surf);
        num_Er_pts=1;
    else
        Er_range=[Er_min,Er_max];
    end
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %call PENTA main function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Er_test_vals,gamma_e_vs_Er{isurf},gamma_i_vs_Er{isurf},...
        Flows_ambi{isurf},Gammas_ambi{isurf},Er_roots{isurf},...
        Flow_vs_Er{isurf},QoTs_ambi{isurf}]=...
        PENTA_main(pprof_info,surf_VMEC_info,Er_range,...
        ion_params,path_info.run_ident,...
        path_info.data_path,B_Eprl(isurf),...
        U2_surf,num_Er_pts,Smax,flags);
    
    roa_vals(isurf)=roa_surf;
    
    %print some outputs to screen
    disp(['r/a = ' num2str(roa_vals(isurf)) ...
        ' Er = ' num2str(Er_roots{isurf}/100) ' V/cm'])

    %Save plasma profile variables used (for comparison to input profiles)
    pprof_info_allsurf.Te_used(isurf)=pprof_info.Te_surf;
    pprof_info_allsurf.Ti_used(isurf,:)=pprof_info.Ti_surf;
    pprof_info_allsurf.ne_used(isurf)=pprof_info.ne_surf;
    pprof_info_allsurf.ni_used(isurf,:)=pprof_info.ni_surf;
    pprof_info_allsurf.dTedr_used(isurf)=pprof_info.dTedr_surf;
    pprof_info_allsurf.dTidr_used(isurf,:)=pprof_info.dTidr_surf;
    pprof_info_allsurf.dnedr_used(isurf)=pprof_info.dnedr_surf;
    pprof_info_allsurf.dnidr_used(isurf,:)=pprof_info.dnidr_surf;
    
    %save the number of ambipolar roots at each surface
    numroots_surf(isurf)=length(Er_roots{isurf}(:));
    
end

% Loop over ambipolar roots and define quantities
numroots=max(numroots_surf);
for is=1:numsurfs
    for ie=1:numroots
        if length(Er_roots{is}(:))>=ie
            Er_ambi(is,ie)=Er_roots{is}(ie);
            gamma_e(is,ie)=Gammas_ambi{is}(1,ie);
            QoT_e(is,ie)=QoTs_ambi{is}(1,ie);
            
            for ispec=1:num_ion_species
                gamma_i(is,ie,ispec)=Gammas_ambi{is}(ispec+1,ie);
                QoT_i(is,ie,ispec)=QoTs_ambi{is}(ispec+1,ie);
            end
        else  %Nonexistant roots will be set to NaN
            Er_ambi(is,ie)=NaN;
            gamma_e(is,ie)=NaN;
            gamma_i(is,ie)=NaN;
            for ispec=1:num_ion_species
                gamma_i(is,ie,ispec)=NaN;
                QoT_i(is,ie,ispec)=NaN;
            end
        end
    end
end