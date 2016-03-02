function [gamma_e_ambi,q_e_ambi,gamma_i_ambi,q_i_ambi,J_bs_ambi,...
    Er_roots,FFmat_out,X_vec_out,J_E_tote,J_E_toti_parts,J_E_cl,...
    sigma_NC2,B_uprle,B_uprli,gamma_e_zero_Er_root]=...
    ...
    PENTA(pprof_info,surf_VMEC_info,Er_range,ion_params,run_ident,...
    data_path,B_Eprl,num_ion_species,U2,num_Er_pts,plot_flux,...
    plotit_Er_find,Smax,log_interp,use_quad,use_Ji,find_gamma_e_root,DKES_limit,load_LMN_files)
%This function defines the main PENTA program.  
%
% 6/2009 JL

%constants
elem_charge = 1.602176487e-19;
e_mass = 9.10938215e-31;

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
Bsq=surf_VMEC_info.bsq;
vmec_jval=surf_VMEC_info.vmec_jval;
chip=surf_VMEC_info.chip;
psip=surf_VMEC_info.psip;
vol_p=surf_VMEC_info.vp;
btheta=surf_VMEC_info.btheta;
bzeta=surf_VMEC_info.bzeta;
iota=surf_VMEC_info.iota;

%expand ion params
ion_mass=ion_params.ion_mass;
Z_ion=ion_params.Z_ion;

%Calculate thermal velocities and loglambda 
vth_i = sqrt(2*Ti*elem_charge./ion_mass);
vth_e = sqrt(2*Te*elem_charge/e_mass);
if Te > 50
    loglambda = 25.3 - 1.15*log10(ne/1e6) + 2.3*log10(Te);
else
    loglambda = 23.4 - 1.15*log10(ne/1e6) + 3.45*log10(Te);
end

%Either load LMN files directly, or load DKES files and convert
if load_LMN_files == 1
    [cmul_mat,efield_mat,coef2d_mat]=load_LMN_star_files(data_path,run_ident,vmec_jval);
else
    %load files containing D* coefficients
    [cmul_mat,efield_mat,coef2d_mat]=load_DKES_star_files(data_path,run_ident,vmec_jval);
    
    %if DKES limit coef2d_mat just contains the D11_star coefficient
    if ~DKES_limit
        [cmul_mat,efield_mat,coef2d_mat,U2]=convert_D_to_LMN(cmul_mat,efield_mat,coef2d_mat,U2,Bsq);
    end
end


%compile arrays for all species
masses=[e_mass ion_mass];
charges=elem_charge*[-1 Z_ion];
v_ths=[vth_e vth_i];
Ts=[Te,Ti];
ns=[ne,ni];

lmat=define_friction_coeffs(masses,charges,v_ths,ns,Ts,loglambda,Smax,use_Ji);

%calculate Spitzer conductivity
sigma_S=-ne^2*elem_charge^2 * lmat{1}(2,2)/( lmat{1}(1,1)*lmat{1}(2,2) - (lmat{1}(1,2))^2);


%Loop over Er to get fluxes as a function of Er
for ie=1:num_Er_test
    
    %test radial electric field
    Er_test=Er_test_vals(ie);   
    
    %check for Er=0 exactly (causes DIV0)
    if Er_test==0
        Er_test=Er_test_vals(ie+1)/2;
        disp(['Can''t use Er=0, using Er=' num2str(Er_test)])
        disp('could make this switch only when log interp on')
    end

    %Form thermodynamic force vector [X(e)1 X(e)2 X(i1)1 X(i1)2 .... XE]
    Xvec=form_Xvec(Er_test,Te,dTedr,ne,dnedr,Ti,dTidr,ni,dnidr,Z_ion,B_Eprl,Bsq);
    
    %define thermal coefficient matrices (La,Ma,Na)
    [Lmat,Mmat,Nmat]=define_thermal_mats(Er_test,ns,Ts,v_ths,charges,masses,loglambda,cmul_mat,efield_mat,coef2d_mat,Smax,log_interp,use_quad,DKES_limit);
    
    if DKES_limit
        FFmat=zeros(2*num_species,length(Xvec));
        for ispec=1:num_species
            index=1+(ispec-1)*2;
            FFmat(index:index+1,index:index+1)=Lmat{ispec};
        end
        
    else
        %define flux-force matrix
        FFmat=define_flux_force_mat(num_ion_species,lmat,Lmat,Mmat,Nmat,Bsq,ne,ni,Z_ion,Smax);
    end
    %calculate PS fluxes
    [gamma_PS,qoT_PS]=calculate_PS_fluxes(Z_ion,num_ion_species,lmat,Xvec,U2);
    
    gamma_e(ie)=sum(FFmat(1,:).*Xvec.') + gamma_PS(1);
    q_e(ie)=-( sum(FFmat(2,:).*Xvec.') + qoT_PS(1) ).*elem_charge*Te;   
    
    for ispec=1:num_ion_species
        ind1=Smax+2+(ispec-1)*2+(ispec-1)*(Smax-1);
        gamma_i(ispec,ie)=sum(FFmat(ind1,:).*Xvec.') + gamma_PS(ispec+1);
        q_i(ispec,ie)=-( sum(FFmat(ind1+1,:).*Xvec.') + qoT_PS(ispec+1) ).*elem_charge*Ti(ispec);
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
    
%     figure; hold on; box on;
%     plot(Er_test_vals/100,q_e,'r--')
%     for ispec=1:num_ion_species
%         plot(Er_test_vals/100,q_i(ispec,:),[colors(mod(ispec,length(colors))+1) '-'])
%     end
%     xlabel('E_r (V/cm)');ylabel('q')
%     legend('Electron','Ions')
%     title(['heat flux for r/a = ',num2str(roa_surf)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Find the ambipolar roots from gamma_e == sum(Z_ion*gamma_i)
[Er_roots]=find_Er_roots_PENTA(gamma_e,gamma_i,Er_test_vals,Z_ion,roa_surf,plotit_Er_find,0);


%Evaluate drives, Lmat, fluxes at ambipolar root(s)
for iroot=1:length(Er_roots)
    
    Er_test=Er_roots(iroot);
    
    %Form thermodynamic force vector
    Xvec=form_Xvec(Er_test,Te,dTedr,ne,dnedr,Ti,dTidr,ni,dnidr,Z_ion,B_Eprl,Bsq);
    
    %define thermal coefficient matrices (La,Ma,Na)
    [Lmat,Mmat,Nmat]=define_thermal_mats(Er_test,ns,Ts,v_ths,charges,masses,loglambda,cmul_mat,efield_mat,coef2d_mat,Smax,log_interp,use_quad,DKES_limit);
    
    if DKES_limit
        FFmat=zeros(2*num_species,length(Xvec));
        for ispec=1:num_species
            index=1+(ispec-1)*2;
            FFmat(index:index+1,index:index+1)=Lmat{ispec};
        end
        
        %set some outputs to NaN
        J_bs_ambi(iroot)=NaN;
        J_E_tote(iroot)=NaN;
        J_E_toti_parts(ispec,iroot)=NaN;
        J_E_cl=NaN;
        sigma_NC2=NaN;
        B_uprle(iroot)=NaN;
        B_uprli(ispec,iroot)=NaN;
        
    else
        
        %calculate parallel flows
        
        Umat2=zeros((Smax+1)*num_species,2*num_species);  %extra col gets added on below
        
        for ispec1=1:num_species    %row
            for ispec2=1:num_species  %column
                if ispec1==ispec2
                    Umat1_add=Mmat{ispec1}./Bsq-lmat{ispec1,ispec2};
                    
                    Ntmp=[Nmat{ispec1}(1:Smax+1,1), -Nmat{ispec1}(1:Smax+1,2)];
                    ind1=1+(ispec1-1)*2+(ispec1-1)*(Smax-1);
                    ind2=1+(ispec2-1)*2;
                    Umat2(ind1:ind1+Smax,ind2:ind2+1)=-Ntmp;
                    
                else
                    Umat1_add=-lmat{ispec1,ispec2};
                end
                ind1=1+(ispec1-1)*2+(ispec1-1)*(Smax-1);
                ind2=1+(ispec2-1)*2+(ispec2-1)*(Smax-1);
                Umat1(ind1:ind1+Smax,ind2:ind2+Smax)=Umat1_add;
            end
        end
        
        %add last column of Umat2
        Umat2_end(1:Smax+1,1)=[-ne*elem_charge*sqrt(Bsq); zeros(Smax,1);];
        for ispec1=1:num_species-1
            ind1=Smax+2+(ispec1-1)*2+(ispec1-1)*(Smax-1);
            Umat2_end(ind1:ind1+Smax,1)=[ni(ispec1)*Z_ion(ispec1)*elem_charge*sqrt(Bsq); zeros(Smax,1);];
        end
        
        Umat2=[Umat2 Umat2_end];
        
        UFmat=(inv(Umat1)*Umat2);
        Uvec=UFmat*Xvec;
        
        %define flux-force matrix
        FFmat=define_flux_force_mat(num_ion_species,lmat,Lmat,Mmat,Nmat,Bsq,ne,ni,Z_ion,Smax);
        
        
        
        B_uprle(iroot)=Uvec(1);
        B_qprle(iroot)=-Uvec(2)*5/2*Te*elem_charge*ne;
        
        %NC parallel conductivity
        sigma_NC_e=[1,zeros(1,Smax)]*(charges(1)*ns(1)*(Bsq*inv(Nmat{1})*FFmat(1:Smax+1,end))/sqrt(Bsq));
        
        UE=inv(Umat1)*Umat2_end;
        s_NC_e=charges(1)*ns(1)*UE(1)/sqrt(Bsq);
        
        %Ion parallel flows (this uses just the bn fluxes)
        for ispec=1:num_ion_species
            ind1=Smax+2+(ispec-1)*2+(ispec-1)*(Smax-1);
            B_uprli(ispec,iroot)=Uvec(ind1);
            B_qprli(ispec,iroot)=Uvec(ind1+1)*5/2*elem_charge*Ti(ispec)*ni(ispec);
            
            %calculate toroidal and poloidal flow components
            u_theta_i(ispec,iroot)=chip*4*pi^2/vol_p*(B_uprli(ispec,iroot)/Bsq - Xvec(3+2*(ispec-1))*bzeta/(Z_ion(ispec)*elem_charge*chip*Bsq));
            u_zeta_i(ispec,iroot)= psip*4*pi^2/vol_p*(B_uprli(ispec,iroot)/Bsq + Xvec(3+2*(ispec-1))*btheta/(Z_ion(ispec)*elem_charge*psip*Bsq));
            
            %NC parallel conductivity
            ind1=Smax+2+(ispec-1)*2+(ispec-1)*(Smax-1);
            sigma_NC_i(ispec)=[1,zeros(1,Smax)]*(charges(ispec+1)*ns(1+ispec)*(Bsq*inv(Nmat{1+ispec})*FFmat(ind1:ind1+Smax,end))/sqrt(Bsq));
            s_NC_i(ispec)=charges(1+ispec)*ns(1+ispec)*UE(ind1)/sqrt(Bsq);
        end
        
        %total parallel conductivity
        sigma_NC=sigma_NC_e+sum(sigma_NC_i);
        sigma_NC2=s_NC_e+sum(s_NC_i);
        
        %calculate toroidal and poloidal flow components
        u_theta_e(iroot)=chip*4*pi^2/vol_p*(B_uprle(iroot)/Bsq + Xvec(1)*bzeta /(elem_charge*chip*Bsq));
        u_zeta_e(iroot) =psip*4*pi^2/vol_p*(B_uprle(iroot)/Bsq - Xvec(1)*btheta/(elem_charge*psip*Bsq));
        
        
        %Calculate total parallel electric current
        %ions
        for ispec=1:num_ion_species
            J_E_toti_parts(ispec,iroot)=ns(1+ispec)*charges(1+ispec)*B_uprli(ispec,iroot)/sqrt(Bsq);
        end
        %electrons
        J_E_tote(iroot)=ns(1)*charges(1)*B_uprle(iroot)/sqrt(Bsq);
        %ions
        J_E_toti(iroot)=sum(J_E_toti_parts(:,iroot));
        %total
        J_E_tot(iroot)=J_E_tote(iroot)+J_E_toti(iroot);
        
        %Calculate bootstrap current
        J_bs_ambi(iroot)=J_E_tot(iroot)-J_E_cl-sigma_NC*Xvec(end);
        
    end %DKES_limit
    %calculate PS fluxes
    [gamma_PS,qoT_PS]=calculate_PS_fluxes(Z_ion,num_ion_species,lmat,Xvec,U2);
    
    %total fluxes
    gamma_e_ambi(iroot)=sum(FFmat(1,:).*Xvec.') + +gamma_PS(1);
    q_e_ambi(iroot)=-(sum(FFmat(2,:).*Xvec.') + qoT_PS(1)).*elem_charge*Te;
    
    for ispec=1:num_ion_species
        ind1=Smax+2+(ispec-1)*2+(ispec-1)*(Smax-1);
        gamma_i_ambi(ispec,iroot)=sum(FFmat(ind1,:).*Xvec.')+ gamma_PS(ispec+1);
        q_i_ambi(ispec,iroot)=-(sum(FFmat(ind1+1,:).*Xvec.') + qoT_PS(ispec+1)).*elem_charge*Ti(ispec);
    end
    
    %complile force-friction matrix and force vector for output
    FFmat_out{iroot}=FFmat;
    X_vec_out{iroot}=Xvec;
    
end

if find_gamma_e_root
    [gamma_e_zero_Er_root]=find_Er_roots_PENTA(gamma_e,gamma_i,Er_test_vals,Z_ion,roa_surf,plotit_Er_find,find_gamma_e_root);
else
    gamma_e_zero_Er_root=NaN;
end