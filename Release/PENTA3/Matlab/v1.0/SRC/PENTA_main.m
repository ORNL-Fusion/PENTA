function [Er_test_vals,gamma_e_vs_Er,gamma_i_vs_Er,...
    Flows_ambi,Gammas_ambi,Er_roots,Flow_vs_Er,QoTs_ambi]=...
    PENTA_main(pprof_info,surf_VMEC_info,Er_range,ion_params,run_ident,...
    data_path,B_Eprl,U2,num_Er_pts,Smax,flags)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function defines the main PENTA program.  
%
% 6/2009 - 5/19/2010 JL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%expand flags structure
plot_flux=flags.plot_flux;
plot_flow=flags.plot_flow;
plotit_Er_find=flags.plotit_Er_find;
use_Ji=flags.use_Ji;
Add_Spitzer_to_D33=flags.Add_Spitzer_to_D33;
log_interp=flags.log_interp;

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
    load_DKES_star_files(data_path,run_ident,vmec_jval,Add_Spitzer_to_D33);

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
    
    %Form thermodynamic force vector [X(e)1 X(e)2 X(i1)1 X(i1)2 .... XE]
    Xvec=form_Xvec(Er_test,Te,dTedr,ne,dnedr,Ti,dTidr,ni,dnidr,Z_ion,B_Eprl,Bsq);
    
    %Also define "A" version of force vector
    for ispec1=1:num_species
        ind_X=(ispec1-1)*2+1;
        ind_A=(ispec1-1)*3+1;
        
        Avec(ind_A)   = -Xvec(ind_X  )/(Ts(ispec1)*elem_charge) - (5/2)*dTdrs(ispec1)/Ts(ispec1);
        Avec(ind_A+1) = -Xvec(ind_X+1)/(Ts(ispec1)*elem_charge);
        Avec(ind_A+2) =  Xvec(end)*charges(ispec1)*B0_vmec/(Ts(ispec1)*elem_charge*sqrt(Bsq));
    end
    
    
    switch flags.Method
        case 'T',
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate flows and particle fluxes using Taguchi method v2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

            Flows=calc_flows_Taguchi_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,D33_star_slice,Smax,flags,Avec,lmat,Bsq,B0_vmec);
            Gammas=calc_flux_Taguchi_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D11_star_slice,D31_star_slice,Smax,flags,Avec,U2,Flows,lmat,B0_vmec);              
        case 'SN',
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate flows and particle fluxes using S-N method v2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

            Flows=calc_flows_SN_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,D33_star_slice,Smax,flags,Avec,lmat,Bsq,B0_vmec);
            [Gammas,QoTs,Gamma_U,Gamma_L,Gamma_PS]=calc_flux_SN_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D11_star_slice,D31_star_slice,D33_star_slice,Smax,flags,Avec,Bsq,U2,Flows,lmat,B0_vmec);     
             gamma_e_vs_Er_L(ie)=Gamma_L(1);
             gamma_i_vs_Er_L(:,ie)=Gamma_L(2:end);
             gamma_e_vs_Er_U(ie)=Gamma_U(1);
             gamma_i_vs_Er_U(:,ie)=Gamma_U(2:end);
             gamma_e_vs_Er_PS(ie)=Gamma_PS(1);
             gamma_i_vs_Er_PS(:,ie)=Gamma_PS(2:end);    
        case 'DKES',
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate flows and particle fluxes using DKES method
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
            
            Flows=calc_flows_DKES(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,Smax,flags,Avec,B0_vmec);
            Gammas=calc_flux_DKES(num_species,abs_Er,ns,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D11_star_slice,D13_star_slice,flags,Avec,B0_vmec); 
        case 'MBT',
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate flows and particle fluxes using M-B-T method v2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Flows=calc_flows_Taguchi_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,D33_star_slice,Smax,flags,Avec,lmat,Bsq,B0_vmec);            
            G2=(Bsq/2)*U2;
            [Gammas,QoTs,DKES_flux,Uflux,Uflux_fric,PS_flux]=calc_flux_MBT_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,dTdrs,dndrs,loglambda,...
                 cmesh,emesh,D11_star_slice,D31_star_slice,D13_star_slice,Smax,B0_vmec,flags,Avec,lmat,Flows,G2);    
             gamma_e_vs_Er_DKES(ie)=DKES_flux(1);
             gamma_i_vs_Er_DKES(:,ie)=DKES_flux(2:end);
             gamma_e_vs_Er_Uflux(ie)=Uflux(1);
             gamma_i_vs_Er_Uflux(:,ie)=Uflux(2:end);
             gamma_e_vs_Er_Uflux_fric(ie)=Uflux_fric(1);
             gamma_i_vs_Er_Uflux_fric(:,ie)=Uflux_fric(2:end);
             gamma_e_vs_Er_PS_flux(ie)=PS_flux(1);
             gamma_i_vs_Er_PS_flux(:,ie)=PS_flux(2:end);             

        otherwise,
            error('Unknown Method')
    end
    
    %Save the flows and particle fluxes as a function of Er
    Flow_vs_Er(:,ie)=Flows;
    gamma_e_vs_Er(ie)=Gammas(1);
    gamma_i_vs_Er(:,ie)=Gammas(2:end);
    
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
    
    %Form thermodynamic force vector [X(e)1 X(e)2 X(i1)1 X(i1)2 .... XE]
    Xvec=form_Xvec(Er_test,Te,dTedr,ne,dnedr,Ti,dTidr,ni,dnidr,Z_ion,B_Eprl,Bsq);
    
    %Convert Xvec to Avec
    for ispec1=1:num_species
        ind_X=(ispec1-1)*2+1;
        ind_A=(ispec1-1)*3+1;
        
        Avec(ind_A)   = -Xvec(ind_X)/(Ts(ispec1)*elem_charge) - (5/2)*dTdrs(ispec1)/Ts(ispec1);
        Avec(ind_A+1) = -Xvec(ind_X+1)/(Ts(ispec1)*elem_charge);
        Avec(ind_A+2) = Xvec(end)*charges(ispec1)*B0_vmec/(Ts(ispec1)*elem_charge*sqrt(Bsq));
    end
    
    switch flags.Method
        case 'T',
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate flows and particle fluxes using Taguchi method v2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

            Flows_ambi(:,iroot)=calc_flows_Taguchi_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,D33_star_slice,Smax,flags,Avec,lmat,Bsq,B0_vmec);
            [Gammas_ambi(:,iroot),QoTs_ambi(:,iroot)]=calc_flux_Taguchi_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D11_star_slice,D31_star_slice,Smax,flags,Avec,U2,Flows_ambi(:,iroot),lmat,B0_vmec);              
        case 'SN',
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate flows and particle fluxes using S-N method v2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

            Flows_ambi(:,iroot)=calc_flows_SN_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,D33_star_slice,Smax,flags,Avec,lmat,Bsq,B0_vmec);
            [Gammas_ambi(:,iroot),QoTs_ambi(:,iroot)]=calc_flux_SN_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D11_star_slice,D31_star_slice,D33_star_slice,Smax,flags,Avec,Bsq,U2,Flows_ambi(:,iroot),lmat,B0_vmec);            
        case 'DKES',
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate flows and particle fluxes using DKES method
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

            Flows_ambi(:,iroot)=calc_flows_DKES(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,Smax,flags,Avec,B0_vmec);
            [Gammas_ambi(:,iroot),QoTs_ambi(:,iroot)]=calc_flux_DKES(num_species,abs_Er,ns,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D11_star_slice,D13_star_slice,flags,Avec,B0_vmec);               
        case 'MBT',
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate flows and particle fluxes using M-B-T method v2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Flows_ambi(:,iroot)=calc_flows_Taguchi_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
                cmesh,emesh,D31_star_slice,D33_star_slice,Smax,flags,Avec,lmat,Bsq,B0_vmec);            
            
            G2=(Bsq/2)*U2;
            [Gammas_ambi(:,iroot),QoTs_ambi(:,iroot)]=calc_flux_MBT_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,dTdrs,dndrs,loglambda,...
                 cmesh,emesh,D11_star_slice,D31_star_slice,D13_star_slice,Smax,B0_vmec,flags,Avec,lmat,Flows_ambi(:,iroot),G2);    

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