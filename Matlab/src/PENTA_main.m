function [Er_test_vals,gamma_e_vs_Er,gamma_i_vs_Er,...
    Flows_ambi,Gammas_ambi,Er_roots,Flow_vs_Er,QoTs_ambi]=...
    PENTA_main(pprof_info,surf_VMEC_info,Er_range,ion_params,run_ident,...
    data_path,B_Eprl,U2,num_Er_pts,Smax,flags)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function defines the main PENTA program.  
%
% 6/2009 - 6/8/2010 JL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%expand flags structure
plot_flux=flags.plot_flux;
plot_flow=flags.plot_flow;
plotit_Er_find=flags.plotit_Er_find;
use_Ji=flags.use_Ji;
Add_Spitzer_to_D33=flags.Add_Spitzer_to_D33;
log_interp=flags.log_interp;
calc_inductive_flow=flags.calc_inductive_flow;

%constants
elem_charge = 1.602176487e-19;
e_mass = 9.10938215e-31;

%expand ion params
ion_mass=ion_params.ion_mass;
Z_ion=ion_params.Z_ion;
num_ion_species=ion_params.num_ion_species;

%total number of particle species
num_species=1+num_ion_species;

%Define Er test range in V/m
num_Er_test=num_Er_pts;
Er_test_vals=linspace(Er_range(1),Er_range(2),num_Er_test).*100;

%expand plasma profile variables
Te=pprof_info.Te_surf;
Ti=pprof_info.Ti_surf;
ne=pprof_info.ne_surf;
ni=pprof_info.ni_surf;
dTedr=pprof_info.dTedr_surf;
dTidr=pprof_info.dTidr_surf;
dnedr=pprof_info.dnedr_surf;
dnidr=pprof_info.dnidr_surf;
roa_surf=pprof_info.roa_surf;

%expand VMEC profile data
Bsq=surf_VMEC_info.Bsq;
vmec_jval=surf_VMEC_info.vmec_jval;
B0_vmec=surf_VMEC_info.B0;

%Calculate thermal velocities and loglambda 
vth_i = sqrt(2*Ti*elem_charge./ion_mass);
vth_e = sqrt(2*Te*elem_charge/e_mass);
if Te > 50
    loglambda = 25.3 - 1.15*log10(ne/1e6) + 2.3*log10(Te);
else
    loglambda = 23.4 - 1.15*log10(ne/1e6) + 3.45*log10(Te);
end

%load files containing D* coefficients (11,13,31,33)
[cmul_array,efield_array,D11_star_slice,...
    D31_star_slice,D13_star_slice,D33_star_slice]=...
    load_DKES_star_files(data_path,run_ident,vmec_jval,Add_Spitzer_to_D33,Bsq);

%define <U**2>
if isempty(U2)
    U2=1.5*D11_star_slice(end,1)./cmul_array(end);
end

%compile arrays for all species
masses=[e_mass ion_mass];
charges=elem_charge*[-1 Z_ion];
v_ths=[vth_e vth_i];
Ts=[Te,Ti];
ns=[ne,ni];
dTdrs=[dTedr,dTidr];
dndrs=[dnedr,dnidr];

%Define classical friction coefficients
lmat=define_friction_coeffs(masses,charges,v_ths,ns,...
    Ts,loglambda,Smax,use_Ji);

%Predefine variables used for interpolation
[cmesh,emesh]=meshgrid(cmul_array,efield_array);

log_D11_star_slice=log(D11_star_slice);
log_D33_star_slice=log(D33_star_slice);

%Loop over Er to get fluxes as a function of Er
for ie=1:num_Er_test
    
    %test radial electric field
    Er_test=Er_test_vals(ie);   
    
    %check for Er=0 exactly when using log(Er/v) for interpolation
    if Er_test==0 && log_interp
        Er_test=Er_test_vals(ie+1)/2;
        Er_test_vals(ie)=Er_test;
        disp(['ln(Er/v) not valid for Er=0, using Er=' num2str(Er_test) ' V/m'])
    end

    %define |Er|
    abs_Er=abs(Er_test);
    
    %Define "A" version of force vector
    Avec=form_Avec(Er_test,Ts,dTdrs,ns,dndrs,charges,B_Eprl,Bsq,B0_vmec);
    
    switch flags.Method
        case {'T' , 'MBT'},
            %Calculate flows and particle fluxes using Taguchi method
            Flows=calc_flows_Taguchi(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,log_D33_star_slice,Smax,flags,Avec,lmat,Bsq,B0_vmec);
            [Gammas]=calc_Gamma_MBT(num_species,abs_Er,ns,Ts,v_ths,charges,masses,dTdrs,dndrs,loglambda,...
                 cmesh,emesh,D11_star_slice,log_D11_star_slice,D31_star_slice,Smax,B0_vmec,flags,Avec,lmat,Flows,U2);  
            [QoTs]=calc_QoT_MBT(num_species,abs_Er,ns,Ts,v_ths,charges,masses,dTdrs,dndrs,loglambda,...
                 cmesh,emesh,D11_star_slice,log_D11_star_slice,D31_star_slice,Smax,B0_vmec,flags,Avec,lmat,Flows,U2);               
        case 'SN',
            %Calculate flows and particle fluxes using S-N method
            Flows=calc_flows_SN(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,D33_star_slice,Smax,flags,Avec,lmat,Bsq,B0_vmec);      
            [Gammas]=calc_Gamma_SN(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D11_star_slice,D31_star_slice,D33_star_slice,log_D11_star_slice, ...
                Smax,flags,Avec,Bsq,U2,Flows,lmat,dTdrs,dndrs);     
            [QoTs]=calc_QoT_SN(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D11_star_slice,D31_star_slice,D33_star_slice,log_D11_star_slice, ...
                Smax,flags,Avec,Bsq,U2,Flows,lmat,dTdrs,dndrs);                    
        case 'DKES',
            %Calculate flows and particle fluxes using DKES method                     
            Flows=calc_flows_DKES(num_species,abs_Er,ns,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,D33_star_slice,Smax,flags,Avec,B0_vmec);
            [Gammas,QoTs]=calc_flux_DKES(num_species,abs_Er,ns,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D11_star_slice,D13_star_slice,flags,Avec,B0_vmec);            
        otherwise,
            error('Unknown Method')
    end
    
    %Save the flows and particle fluxes as a function of Er
    Flow_vs_Er(:,ie)=Flows;
    gamma_e_vs_Er(ie)=Gammas(1);
    gamma_i_vs_Er(:,ie)=Gammas(2:end);
    QoT_e_vs_Er(ie)=QoTs(1);
    QoT_i_vs_Er(:,ie)=QoTs(2:end);    
end

colors=['m','b','g','k'];
if plot_flux
    figure; hold on; box on;
    plot(Er_test_vals/100,gamma_e_vs_Er/10^20,'r-')
    for ispec=1:num_ion_species
        plot(Er_test_vals/100,gamma_i_vs_Er(ispec,:)/10^20,[colors(mod(ispec,length(colors))+1) '-'])
    end
    xlabel('E_r (V/cm)');ylabel('\Gamma (10^{20}m^{-2}s^{-1})')
    legend('Electron','Ions')
    title(['particle flux for r/a = ',num2str(roa_surf)])

    figure; hold on; box on;
    plot(Er_test_vals/100,QoT_e_vs_Er,'r-')
    for ispec=1:num_ion_species
        plot(Er_test_vals/100,QoT_i_vs_Er(ispec,:),[colors(mod(ispec,length(colors))+1) '-'])
    end
    xlabel('E_r (V/cm)');ylabel('Q/T')
    legend('Electron','Ions')
    title(['QoT flux for r/a = ',num2str(roa_surf)])    
end

if plot_flow    
    %plot particle flows
    figure; hold on; box on;
    plot(Er_test_vals/100,Flow_vs_Er(1,:)/1e3,'r-')
    for ispec=1:num_ion_species
        plot(Er_test_vals/100,Flow_vs_Er(Smax+2,:)/1e3,[colors(mod(ispec,length(colors))+1) '-'])
    end
    xlabel('E_r (V/cm)');ylabel('<Bu_{||}>/<B^2> (km/s)')
    legend('Electron','Ions')
    title(['particle flows for r/a = ',num2str(roa_surf)])
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the Ambipolar roots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(Er_test_vals)==1
    Er_roots=Er_test_vals;
else
    %Find the ambipolar roots from gamma_e == sum(Z_ion*gamma_i)
    [Er_roots]=find_Er_roots_PENTA(gamma_e_vs_Er,gamma_i_vs_Er,Er_test_vals,Z_ion,roa_surf,plotit_Er_find,0);
end

num_Er_roots=length(Er_roots);

%Evaluate at ambipolar root(s)
for iroot=1:num_Er_roots
    
    Er_test=Er_roots(iroot);
    abs_Er=abs(Er_test);
    Avec=form_Avec(Er_test,Ts,dTdrs,ns,dndrs,charges,B_Eprl,Bsq,B0_vmec);    

    if calc_inductive_flow
        Avec(1:end) = 0;
        Avec(3:3:end) = 1;
    end
         
    switch flags.Method
        case {'T' , 'MBT'},
            %Calculate flows and particle fluxes using Taguchi method
            Flows_ambi(:,iroot)=calc_flows_Taguchi(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,log_D33_star_slice,Smax,flags,Avec,lmat,Bsq,B0_vmec);         
            [Gammas_ambi(:,iroot)]=calc_Gamma_MBT(num_species,abs_Er,ns,Ts,v_ths,charges,masses,dTdrs,dndrs,loglambda,...
                 cmesh,emesh,D11_star_slice,log_D11_star_slice,D31_star_slice,Smax,B0_vmec,flags,Avec,lmat,Flows_ambi(:,iroot),U2);         
            [QoTs_ambi(:,iroot)]=calc_QoT_MBT(num_species,abs_Er,ns,Ts,v_ths,charges,masses,dTdrs,dndrs,loglambda,...
                 cmesh,emesh,D11_star_slice,log_D11_star_slice,D31_star_slice,Smax,B0_vmec,flags,Avec,lmat,Flows_ambi(:,iroot),U2);               
        case 'SN',
            %Calculate flows and particle fluxes using S-N method
            Flows_ambi(:,iroot)=calc_flows_SN(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,D33_star_slice,Smax,flags,Avec,lmat,Bsq,B0_vmec);          
            [Gammas_ambi(:,iroot)]=calc_Gamma_SN(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D11_star_slice,D31_star_slice,D33_star_slice,log_D11_star_slice, ...
                Smax,flags,Avec,Bsq,U2,Flows_ambi(:,iroot),lmat,dTdrs,dndrs);              
            [QoTs_ambi(:,iroot)]=calc_QoT_SN(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D11_star_slice,D31_star_slice,D33_star_slice,log_D11_star_slice, ...
                Smax,flags,Avec,Bsq,U2,Flows_ambi(:,iroot),lmat,dTdrs,dndrs);                     
        case 'DKES',
            %Calculate flows and particle fluxes using DKES method                     
            Flows_ambi(:,iroot)=calc_flows_DKES(num_species,abs_Er,ns,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,D33_star_slice,Smax,flags,Avec,B0_vmec);            
            [Gammas_ambi(:,iroot),QoTs_ambi(:,iroot)]=calc_flux_DKES(num_species,abs_Er,ns,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D11_star_slice,D13_star_slice,flags,Avec,B0_vmec);                     
        otherwise,
            error('Unknown Method')
    end    

end %iroot loop




% 
% if find_gamma_e_root
%     [gamma_e_zero_Er_root]=find_Er_roots_PENTA(gamma_e,gamma_i,Er_test_vals,Z_ion,roa_surf,plotit_Er_find,find_gamma_e_root);
% else
%     gamma_e_zero_Er_root=NaN;
% end