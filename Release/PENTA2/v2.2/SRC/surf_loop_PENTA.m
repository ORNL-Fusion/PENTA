function [roa_vals,Er_ambi,gamma_e,q_e,gamma_i,q_i,J_bs,J_E_e,J_E_i,J_E_cl,FFmat,X_vec,BUprl_e,BUprl_i,gamma_e_zero_Er_root]=surf_loop_PENTA(data_path,run_ident,pprof_char,run_surfs,...
    Er_min,Er_max,input_is_Er,plot_pprof,postproc_plots,...
    arad,B_Eprl,num_Er_pts,plot_flux,plotit_Er_find,Smax,log_interp,use_quad,use_Ji,find_gamma_e_root,DKES_limit,load_LMN_files)

save_pprof_data=0;
more_text_out=0;

p_mass=1.672621637e-27;             %proton mass

%load ion parameters from file
ion_file=[data_path '/ion_params'];
fid_ion=fopen(ion_file);
ion_data=textscan(fid_ion,'%s');
id1=char(ion_data{1}(2));
id2=char(ion_data{1}(3));
id3=char(ion_data{1}(4));
num_ion_species=str2num(id1(find(id1=='=')+1:end));
Z_ion=str2num(id2(find(id2=='=')+1:end));
miomp=str2num(id3(find(id3=='=')+1:end));
fclose(fid_ion);

ion_mass=miomp*p_mass;

%check if ion parameters defined
if length(ion_mass)~= num_ion_species || length(Z_ion)~= num_ion_species
    error('Variables ion_mass and Z_ion must have length==num_ion_species')
end

%Compile ion parameters
ion_params.ion_mass=ion_mass;
ion_params.Z_ion=Z_ion;

%Read plasma profile file
[pp_ne,pp_Te,pp_ne_prime,pp_Te_prime,pp_Ti,pp_Ti_prime,pp_ni,pp_ni_prime]=...
    read_pprof_file(data_path,pprof_char,num_ion_species,plot_pprof);

%read VMEC profile info
data1_vmec=dlmread([data_path '/profile_data_' run_ident],'',[0,0,1,1]);
data2_vmec=dlmread([data_path '/profile_data_' run_ident],'',3,0);
%use predefined <a> if given above
if isempty(arad)
    arad=data1_vmec(2,1);    %<a>
end
Rmajor_vmec=data1_vmec(2,2); %<R0>
js_vmec=data2_vmec(:,1);     %VMEC surface #
r_vmec=data2_vmec(:,2);      %<r>
roa_vmec=data2_vmec(:,3);    %sqrt(norm. tor. flux)
chip_vmec=data2_vmec(:,4);   %deriv of poloidal flux
psip_vmec=data2_vmec(:,5);   %deriv of toroidal flux
btheta_vmec=data2_vmec(:,6); %Boozer pol field
bzeta_vmec=data2_vmec(:,7);  %Boozer tor field
vp_vmec=data2_vmec(:,8);     %volume derivative
bsq_vmec=data2_vmec(:,9);    %<B**2>
iota_vmec=data2_vmec(:,10);  %iota

%read U^2 file
U2_data=dlmread([data_path '/Utilde2_profile'],'',1,0);

%Intro text
fprintf('---------------------------------------------------------------');
fprintf('\n\nWelcome to PENTA, please note the following settings:\n\n');

fprintf('Data path: %s\n',data_path);
fprintf('Device <R0>, <a>: %5.3f  %5.3f \n',Rmajor_vmec, arad);
fprintf('Number of ion species: %3i \n',num_ion_species);
fprintf('Sonine order: %3i \n',Smax);
for ispec=1:num_ion_species
    fprintf('Ion species %3i: mi/mp=%5.3f  Zi=%5.3f \n',...
        ispec,ion_mass(ispec)./p_mass, Z_ion(ispec));
end
if input_is_Er
    fprintf('Interpreting search range as Er\n');
else
    fprintf('Interpreting search range as e<a>*Phi/kTe\n');
    error('I haven''t implemented this yet (e<a>Phi/kTe)')
end
if DKES_limit
    fprintf('Using DKES limit\n\n');
end

%loop over surfaces and call PENTA function
numsurfs=length(run_surfs);
for isurf = 1:numsurfs
    %get index of vmec surface
    vmec_index=find(js_vmec==run_surfs(isurf)); %(should be run_surfs(isurf)-1)
    
    %get r/a for this surface
    roa_surf=roa_vmec(vmec_index);
    
    %Write screen output
%     disp(['Working on surface ' num2str(isurf) ' of ' num2str(numsurfs) ', r/a = ' num2str(roa_surf)]);
    
    %get local plasma profile info and compile into structure
    pprof_info.Te_surf=ppval(pp_Te,roa_surf);
    pprof_info.ne_surf=ppval(pp_ne,roa_surf);
    pprof_info.dnedr_surf=ppval(pp_ne_prime,roa_surf)/arad;
    pprof_info.dTedr_surf=ppval(pp_Te_prime,roa_surf)/arad;
    pprof_info.roa_surf=roa_surf;
    
    for ispec=1:num_ion_species
        pprof_info.Ti_surf(ispec)=ppval(pp_Ti(ispec,:),roa_surf);
        pprof_info.ni_surf(ispec)=ppval(pp_ni(ispec,:),roa_surf);
        pprof_info.dTidr_surf(ispec)=ppval(pp_Ti_prime(ispec,:),roa_surf)/arad;
        pprof_info.dnidr_surf(ispec)=ppval(pp_ni_prime(ispec,:),roa_surf)/arad;    
    end
     
    %get local VMEC info and compile into structure
    surf_VMEC_info.chip=chip_vmec(vmec_index);
    surf_VMEC_info.psip=psip_vmec(vmec_index);
    surf_VMEC_info.bsq=bsq_vmec(vmec_index);
    surf_VMEC_info.vp=vp_vmec(vmec_index);
    surf_VMEC_info.btheta=btheta_vmec(vmec_index);
    surf_VMEC_info.bzeta=bzeta_vmec(vmec_index);
    surf_VMEC_info.iota=iota_vmec(vmec_index);
    surf_VMEC_info.vmec_jval=js_vmec(vmec_index);
    
    %define local U2
    U2_surf=interp1(U2_data(:,1),U2_data(:,2),roa_surf);
    
    %define search range
    Er_range=[Er_min,Er_max];
    
    %call PENTA
    [gamma_e_ambi{isurf},q_e_ambi{isurf},gamma_i_ambi{isurf},q_i_ambi{isurf},...
        J_bs_ambi{isurf},Er_roots{isurf},FFmat{isurf},X_vec{isurf},...
        J_E_tote{isurf},J_E_toti{isurf},J_E_cl(isurf),sigma_NC2(isurf),...
        B_uprle{isurf},B_uprli{isurf},gamma_e_zero_Er_root(isurf)]=...
        PENTA(pprof_info,surf_VMEC_info,Er_range,ion_params,run_ident,data_path,B_Eprl(isurf),...
        num_ion_species,U2_surf,num_Er_pts,plot_flux,plotit_Er_find,Smax,log_interp,use_quad,use_Ji,find_gamma_e_root,DKES_limit,load_LMN_files);
    roa_vals(isurf)=roa_surf;
%     %print some outputs to screen
    disp(['r/a = ' num2str(roa_vals(isurf)) ' Er = ' num2str(Er_roots{isurf}/100) ' V/cm'])
    
    if more_text_out == 1
      disp(['g_e = ' num2str(gamma_e_ambi{isurf})])
    disp(['q_e = ' num2str(q_e_ambi{isurf})])
    for ispec=1:num_ion_species
        disp(['g_i for species ' num2str(ispec) ' = ' num2str(gamma_i_ambi{isurf}(ispec,:))])
        disp(['q_i for species ' num2str(ispec) ' = ' num2str(q_i_ambi{isurf}(ispec,:))])
    end
    disp(['J_bs = ' num2str(J_bs_ambi{isurf})])
    disp(['sig_NC = ' num2str(sigma_NC2(isurf))])
    disp(['BUprls = ' num2str(BUprls{isurf})])
    end

    %plasma profile variables for saving
    Te(isurf)=pprof_info.Te_surf;
    Ti(isurf,:)=pprof_info.Ti_surf;
    ne(isurf)=pprof_info.ne_surf;
    ni(isurf,:)=pprof_info.ni_surf;
    dTedr(isurf)=pprof_info.dTedr_surf;
    dTidr(isurf,:)=pprof_info.dTidr_surf;
    dnedr(isurf)=pprof_info.dnedr_surf;
    dnidr(isurf,:)=pprof_info.dnidr_surf;
    
end

% Loop over ambipolar roots and define quantities
for is=1:numsurfs
    for ie=1:3
        if length(Er_roots{is}(:))>=ie
            Er_ambi(is,ie)=Er_roots{is}(ie);
            J_bs(is,ie)=J_bs_ambi{is}(ie);
            gamma_e(is,ie)=gamma_e_ambi{is}(ie);
            q_e(is,ie)=q_e_ambi{is}(ie);
            J_E_e(is,ie)=J_E_tote{is}(ie);
            BUprl_e(is,ie)=B_uprle{is}(ie);
            
            for ispec=1:num_ion_species
                gamma_i(is,ie,ispec)=gamma_i_ambi{is}(ispec,ie);
                q_i(is,ie,ispec)=q_i_ambi{is}(ispec,ie);
                J_E_i(is,ie,ispec)=J_E_toti{is}(ispec,ie);
                BUprl_i(is,ie,ispec)=B_uprli{is}(ispec,ie);
            end
        else
            Er_ambi(is,ie)=NaN;
            J_bs(is,ie)=NaN;
            gamma_e(is,ie)=NaN;
            q_e(is,ie)=NaN;
            J_E_e(is,ie)=NaN;
            BUprl_e(is,ie)=NaN;
            for ispec=1:num_ion_species
                gamma_i(is,ie,ispec)=NaN;
                q_i(is,ie,ispec)=NaN;
                J_E_i(is,ie,ispec)=NaN;
                BUprl_i(is,ie,ispec)=NaN;
            end
        end
    end
end


if postproc_plots
%postprocessing
figure;hold on; box on;
plot(roa_vals,Er_ambi/100,'o')
xlabel('r/a');ylabel('E_r (V/cm)')

figure;hold on; box on;
plot(roa_vals,J_bs,'o-')
xlabel('r/a');ylabel('J_bs (A/m^2)')

figure;hold on; box on;
plot(roa_vals,gamma_e/1e20,'ro-')
colors=['b','m','g'];
for ispec=1:num_ion_species
    plot(roa_vals,gamma_i(:,:,ispec)/1e20,[colors(mod(ispec,length(colors))+1) 'o-'])
end
xlabel('r/a');ylabel('\Gamma')
end

if save_pprof_data
    save([data_path '/pprof_check.mat'],'roa_vals','Te','Ti','ne',...
        'ni','dTedr','dTidr','dnedr','dnidr');
end
