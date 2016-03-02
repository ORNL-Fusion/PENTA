function [Gammas,QoTs,Gamma_U,Gamma_L,gamma_PS]=calc_flux_SN_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
    cmesh,emesh,D11_star_slice,D31_star_slice,D33_star_slice,Smax,flags,Avec,Bsq,U2,Flows,lmat,B_0)
%This function calculates arrays of fluxes using the S-N method.
%The fluxes are defined as Gamma = [Gamma_e Gamma_i1 ...] 
% qoTs=qa/Ta=[q_e/Te q_i1/Ti1 ...]
%Note that the S&N definition of q/T is only the diffusive flux.
% JL 2/2010

%constants
elem_charge = 1.602176487e-19;
coefname='';
%expand flags structure
log_interp=flags.log_interp;
use_quad=flags.use_quad;

%none of the coeffs currently use logopt
logopt=0;

%loop over species and calculate fluxes
for ispec1=1:num_species
    %species parameters
    ma=masses(ispec1);
    qa=charges(ispec1);
    qa2=charges(ispec1)^2;
    na=ns(ispec1);
    vta=v_ths(ispec1);
    Ta=Ts(ispec1);
    
    %all other species parameters
    qb2=charges(1:end~=ispec1).^2;
    nb=ns(1:end~=ispec1);
    vtb=v_ths(1:end~=ispec1);
    
    %L11 L12 L21 and L22
    norm_factor=na*(ma^2*vta^3/(2*qa^2));
    nu_exp=0;
    intfac=max( D11_star_slice - (2/3)*cmesh.'*U2 + D31_star_slice.^2./D33_star_slice/B_0^2,0 );  %Don't allow negative L
    K_exp=1.5; 
    Lmat{ispec1}(1,1)=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,intfac,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);      
    K_exp=1.5+1; 
    Lmat{ispec1}(1,2)=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,intfac,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);
    Lmat{ispec1}(2,1)=Lmat{ispec1}(1,2);  %Use symmetry relation
    K_exp=1.5+2; 
    Lmat{ispec1}(2,2)=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,intfac,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);
    
    %Loop over primary Sonine index (k)
    for kval=0:Smax

        %N_j+1,k+1
        norm_factor=na * ma*vta/qa * (2/3)*Bsq/B_0;
        K_exp=1.5; nu_exp=0;
        intfac=D31_star_slice./D33_star_slice;
        Nvec{ispec1}(1,kval+1)=energy_conv3(0,kval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,intfac,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);        
        K_exp=2.5; nu_exp=0;
        Nvec{ispec1}(2,kval+1)=energy_conv3(0,kval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,intfac,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);        
        
    end %k loop
    
    %calculate PS fluxes
    lX_sum=0;
    for ispec2=1:num_species
        spec2_Aind=(ispec2-1)*3+1;
        lX_sum=lX_sum+(lmat{ispec1,ispec2}(1:2,1:2).*[1 -1;-1 1])*...
            [Avec(spec2_Aind)+(5/2)*Avec(spec2_Aind+1);Avec(spec2_Aind+1)]...
            .*Ts(ispec2)*elem_charge/charges(ispec2);
    end
    
    PS_fluxes=-U2*lX_sum/charges(ispec1);
    gamma_PS(ispec1)=PS_fluxes(1);
    QoT_PS(ispec1)=-PS_fluxes(2)+(5/2)*Ta*elem_charge*gamma_PS(ispec1);   %above gives -q_PS/T
    
    
    %Calculate Gamma_a0 and Gamma_a1
    flow_ind=(ispec1-1)*(Smax+1)+1;
    ind_A=(ispec1-1)*3+1;
    Gamma_U(ispec1)=-sum(Nvec{ispec1}(1,:).*Flows(flow_ind:flow_ind+Smax).');
    Gamma_L(ispec1)=- Lmat{ispec1}(1,1)*Avec(ind_A) - Lmat{ispec1}(1,2)*Avec(ind_A+1);    
    Gammas(ispec1)=   -sum(Nvec{ispec1}(1,:).*Flows(flow_ind:flow_ind+Smax).') - Lmat{ispec1}(1,1)*Avec(ind_A) - Lmat{ispec1}(1,2)*Avec(ind_A+1)+ gamma_PS(ispec1);
    QoTs(ispec1)  =   -sum(Nvec{ispec1}(2,:).*Flows(flow_ind:flow_ind+Smax).') - Lmat{ispec1}(2,1)*Avec(ind_A) - Lmat{ispec1}(2,2)*Avec(ind_A+1) + QoT_PS(ispec1)  ;
end %species loop for flows

