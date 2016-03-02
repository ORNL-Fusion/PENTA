function [gamma_e_ambi,q_e_ambi,gamma_i_ambi,q_i_ambi,J_bs_ambi,Er_roots,FFmat_out,X_vec_out,J_E_tote,J_E_toti_parts,J_E_cl]=...
    PENTA(pprof_info,surf_VMEC_info,Er_range,ion_params,run_ident,data_path,B_Eprl,num_ion_species,U2,num_Er_pts,plot_flux,plotit_Er_find)
%This function defines the main PENTA program.  
%
% 6/2009 JL


% plot_flux=1;      %particle fluxes vs. Er are plotted for each surface
% plotit_Er_find=1; %plot root finding results

%constants
elem_charge = 1.602176487e-19;
e_mass = 9.10938215e-31;


%Define Er test range in V/m
num_Er_test=num_Er_pts;
Er_test_vals=linspace(Er_range(1),Er_range(2),num_Er_test).*100;

%break out plasma profile variables
Te=pprof_info.Te_surf;
Ti=pprof_info.Ti_surf;
ne=pprof_info.ne_surf;
ni=pprof_info.ni_surf;
dTedr=pprof_info.dTedr_surf;
dTidr=pprof_info.dTidr_surf;
dnedr=pprof_info.dnedr_surf;
dnidr=pprof_info.dnidr_surf;
roa_surf=pprof_info.roa_surf;

%break out VMEC profile data
Bsq=surf_VMEC_info.bsq;
vmec_jval=surf_VMEC_info.vmec_jval;
chip=surf_VMEC_info.chip;
psip=surf_VMEC_info.psip;
vol_p=surf_VMEC_info.vp;
btheta=surf_VMEC_info.btheta;
bzeta=surf_VMEC_info.bzeta;
iota=surf_VMEC_info.iota;

%break out ion params
ion_mass=ion_params.ion_mass;
Z_ion=ion_params.Z_ion;

%
%	Calculate thermal velocities and loglambda 
%
vth_i = sqrt(2*Ti*elem_charge./ion_mass);
vth_e = sqrt(2*Te*elem_charge/e_mass);
if Te > 50
    loglambda = 25.3 - 1.15*log10(ne/1e6) + 2.3*log10(Te);
else
    loglambda = 23.4 - 1.15*log10(ne/1e6) + 3.45*log10(Te);
end

%%%%% Load L*,M*,N* files for this surface
[cmul_mat,efield_mat,coef2d_mat]=load_LMN_star_files(data_path,run_ident,vmec_jval);


%Define friction coefficients
lmat=define_friction_coeffs(ion_mass,Z_ion,vth_e,vth_i,ne,ni,Te,Ti,loglambda);

%calculate Spitzer conductivity
sigma_S=ne^2*elem_charge^2 *  lmat{1}(2,2)/( lmat{1}(1,1)*lmat{1}(2,2) + (lmat{1}(1,2))^2);

%Loop over Er to get fluxes as a function of Er
for ie=1:num_Er_test
    
    %test radial electric field
    Er_test=Er_test_vals(ie);   
    if Er_test==0
        Er_test=Er_test_vals(ie+1)/2;
        disp(['Can''t use Er=0, using Er=' num2str(Er_test)])
        disp('could make this switch only when log interp on')
    end

    %Form thermodynamic force vector
    Xvec=form_Xvec(Er_test,Te,dTedr,ne,dnedr,Ti,dTidr,ni,dnidr,Z_ion,B_Eprl,Bsq);
    
    %define thermal coefficient matrices (La,Ma,Na)
    [Lmat,Mmat,Nmat]=define_thermal_mats(Er_test,ne,vth_e,Z_ion,ion_mass,ni,vth_i,loglambda,cmul_mat,efield_mat,coef2d_mat);

    %define flux-force matrix
    FFmat=define_flux_force_mat(num_ion_species,lmat,Lmat,Mmat,Nmat,Bsq,ne,ni,Z_ion);
    
    %calculate PS fluxes
    [gamma_PS,q_PS]=calculate_PS_fluxes(Z_ion,num_ion_species,lmat,Xvec,U2);
    
    gamma_e(ie)=sum(FFmat(1,:).*Xvec.')+gamma_PS(2);
    q_e(ie)=(sum(FFmat(2,:).*Xvec.')+q_PS(1)).*elem_charge*Te;
    
    for ispec=1:num_ion_species
        gamma_i(ispec,ie)=sum(FFmat(3+2*(ispec-1),:).*Xvec.')+gamma_PS(ispec);
        q_i(ispec,ie)=(sum(FFmat(4+2*(ispec-1),:).*Xvec.')+q_PS(ispec)).*elem_charge*Ti(ispec);
    end
 
end


%calculate classical parallel electrical current
J_E_cl=sigma_S*Xvec(end);

colors=['b','m','g'];
if plot_flux
    figure; hold on; box on;
    plot(Er_test_vals/100,gamma_e/10^20,'r--')
    for ispec=1:num_ion_species
        plot(Er_test_vals/100,gamma_i(ispec,:)/10^20,[colors(mod(ispec,length(colors))+1) '-'])
    end
    xlabel('E_r (V/cm)');ylabel('10^{20} \Gamma')
    legend('Electron','Ions')
    title(['particle flux for r/a = ',num2str(roa_surf)])    
    
    figure; hold on; box on;
    plot(Er_test_vals/100,q_e,'r--')
    for ispec=1:num_ion_species
        plot(Er_test_vals/100,q_i(ispec,:),[colors(mod(ispec,length(colors))+1) '-'])
    end
    xlabel('E_r (V/cm)');ylabel('q')
    legend('Electron','Ions')
    title(['heat flux for r/a = ',num2str(roa_surf)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Find the ambipolar roots from gamma_e == sum(Z_ion*gamma_i)
[Er_roots]=find_Er_roots_PENTA(gamma_e,gamma_i,Er_test_vals,Z_ion,roa_surf,plotit_Er_find);

%Evaluate drives, Lmat, fluxes at ambipolar root(s)
for iroot=1:length(Er_roots)
    
    Er_test=Er_roots(iroot);
    
    
    %Form thermodynamic force vector
    Xvec=form_Xvec(Er_test,Te,dTedr,ne,dnedr,Ti,dTidr,ni,dnidr,Z_ion,B_Eprl,Bsq);
    
    %define thermal coefficient matrices (La,Ma,Na)
    [Lmat,Mmat,Nmat]=define_thermal_mats(Er_test,ne,vth_e,Z_ion,ion_mass,ni,vth_i,loglambda,cmul_mat,efield_mat,coef2d_mat);

    %define flux-force matrix
    FFmat=define_flux_force_mat(num_ion_species,lmat,Lmat,Mmat,Nmat,Bsq,ne,ni,Z_ion);
        
    %calculate PS fluxes
    [gamma_PS,q_PS]=calculate_PS_fluxes(Z_ion,num_ion_species,lmat,Xvec,U2);
    
    gamma_e_ambi(iroot)=sum(FFmat(1,:).*Xvec.')+gamma_PS(1);
    q_e_ambi(iroot)=(sum(FFmat(2,:).*Xvec.')+q_PS(1)).*elem_charge*Te;
    
    for ispec=1:num_ion_species
        gamma_i_ambi(ispec,iroot)=sum(FFmat(3+2*(ispec-1),:).*Xvec.')+gamma_PS(ispec);
        q_i_ambi(ispec,iroot)=(sum(FFmat(4+2*(ispec-1),:).*Xvec.')+gamma_PS(ispec)).*elem_charge*Ti(ispec);
    end
    
    %Now calculate the parallel flows for each species
    UeandQe=Bsq*inv(Nmat{1})*(-Lmat{1}*Xvec(1:2) + [gamma_e_ambi(iroot);q_e_ambi(iroot)/elem_charge/Te]);
    B_uprle(iroot)=UeandQe(1);
    B_qprle(iroot)=UeandQe(2)*5/2*Te*elem_charge*ne;
    
    %Ion parallel flows
    for ispec=1:num_ion_species
        UiandQi=Bsq*inv(Nmat{1+ispec})*(-Lmat{1+ispec}*Xvec(3+2*(ispec-1):4+2*(ispec-1)) + [gamma_i_ambi(ispec,iroot);q_i_ambi(ispec,iroot)/elem_charge/Ti(ispec)]);
        B_uprli(ispec,iroot)=UiandQi(1);
        B_qprli(ispec,iroot)=UiandQi(2)*5/2*elem_charge*Ti(ispec)*ni(ispec);
        
        %calculate toroidal and poloidal flow components
        u_theta_i(ispec,iroot)=chip*4*pi^2/vol_p*(B_uprli(ispec,iroot)/Bsq - Xvec(3+2*(ispec-1))*bzeta/(Z_ion(ispec)*elem_charge*chip*Bsq));
        u_zeta_i(ispec,iroot)= psip*4*pi^2/vol_p*(B_uprli(ispec,iroot)/Bsq + Xvec(3+2*(ispec-1))*btheta/(Z_ion(ispec)*elem_charge*psip*Bsq));
        
    end
    
    %calculate toroidal and poloidal flow components
    u_theta_e(iroot)=chip*4*pi^2/vol_p*(B_uprle(iroot)/Bsq + Xvec(1)*bzeta /(elem_charge*chip*Bsq));
    u_zeta_e(iroot) =psip*4*pi^2/vol_p*(B_uprle(iroot)/Bsq - Xvec(1)*btheta/(elem_charge*psip*Bsq));
    
    
    %Calculate total parallel electric current
    %ions
    for ispec=1:num_ion_species
        J_E_toti_parts(ispec,iroot)=ni(ispec)*elem_charge*Z_ion(ispec)*B_uprli(ispec,iroot)/sqrt(Bsq);
    end
    %electrons
    J_E_tote(iroot)=-ne*elem_charge*B_uprle(iroot)/sqrt(Bsq);
    %ions
    J_E_toti(iroot)=sum(J_E_toti_parts(:,iroot));
    %total
    J_E_tot(iroot)=J_E_tote(iroot)+J_E_toti(iroot);
    
    %Calculate bootstrap current
    J_bs_ambi(iroot)=J_E_tot(iroot)-J_E_cl;
    

    
    %complile force-friction matrix and force vector for output
    FFmat_out{iroot}=FFmat;
    X_vec_out{iroot}=Xvec;
    
end
