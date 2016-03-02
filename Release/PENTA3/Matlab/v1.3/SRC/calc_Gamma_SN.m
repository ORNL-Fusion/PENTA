function [Gammas,Gamma_U,Gamma_L,Gamma_PS]=calc_Gamma_SN(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
    cmesh,emesh,D11_star_slice,D31_star_slice,D33_star_slice,log_D11_star_slice,Smax,flags,Avec,Bsq,U2,Flows,lmat,dTdrs,dndrs)
%This function calculates arrays of fluxes using the S-N method.
%The fluxes are defined as Gamma = [Gamma_e Gamma_i1 ...] 
% JL 2/2010 - 9/21/2011

%constants
elem_charge = 1.602176487e-19;

%expand flags structure
log_interp=flags.log_interp;
use_quad=flags.use_quad;
flux_cap=flags.flux_cap;

%allocate arrays
Nvec = zeros(Smax,1);
Gamma_U = zeros(num_species,1);
Gamma_L = zeros(num_species,1);
Gamma_PS = zeros(num_species,1);
Gammas = zeros(num_species,1);

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
    
    norm_factor=na*(ma^2*vta^3/(2*qa^2));
    if flux_cap == 1
        % Cap radial coefficient at 0
        intfac=max(D11_star_slice  - (2/3)*U2*cmesh.' + D31_star_slice.^2./D33_star_slice,0);
        L11=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,intfac,log_interp,use_quad,1.5,0,norm_factor,flags);
        L12=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,intfac,log_interp,use_quad,2.5,0,norm_factor,flags);        
    else       
        intfac=log_D11_star_slice;
        L11_1=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,1,cmesh,emesh,intfac,log_interp,use_quad,1.5,0,norm_factor,flags);
        L12_1=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,1,cmesh,emesh,intfac,log_interp,use_quad,2.5,0,norm_factor,flags);
        
        intfac=D31_star_slice.^2./D33_star_slice;        
        L11_2=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,intfac,log_interp,use_quad,1.5,0,norm_factor,flags);        
        L12_2=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,intfac,log_interp,use_quad,2.5,0,norm_factor,flags);                
        
        L11=L11_1+L11_2;
        L12=L12_1+L12_2;
    end
    
    %Loop over primary Sonine index (k)
    norm_factor=na * ma*vta/qa * (2/3)*Bsq;
    for kval=0:Smax
        Nvec(kval+1)=energy_conv3(0,kval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,D31_star_slice./D33_star_slice,log_interp,use_quad,1.5,0,norm_factor,flags);        
    end
            
    %Calculate PS flux

   
    Gamma_PS_std = 0;
    %loop over ALL species
    for ispec2=1:num_species
        
        %loop specific species beta
        n_beta=ns(ispec2);
        T_beta=Ts(ispec2);
        dn_betadr=dndrs(ispec2);
        dT_betadr=dTdrs(ispec2);
        q_beta=charges(ispec2);
        
        C_PS1=(n_beta*dT_betadr + T_beta*dn_betadr)*elem_charge/(n_beta*q_beta);
        C_PS2=dT_betadr*elem_charge/q_beta;
        
        Gamma_PS_std=Gamma_PS_std+C_PS1*lmat{ispec1,ispec2}(1,1)-C_PS2*lmat{ispec1,ispec2}(1,2);  
    end
    
    flow_ind=(ispec1-1)*(Smax+1)+1;
    ind_A=(ispec1-1)*3+1;

    if flux_cap == 1 
        Gamma_PS(ispec1) = (U2/qa)*Gamma_PS_std;
    else
        [I_0]=energy_conv_PS(vta,qa2,qb2,ma,na,nb,loglambda,vtb,1,na,flags);
        [I_1]=energy_conv_PS(vta,qa2,qb2,ma,na,nb,loglambda,vtb,2,na,flags);
        Gamma_PS_flow=(2/3)*(ma*Ta*elem_charge/qa)*( I_0*Avec(ind_A) + I_1*Avec(ind_A+1) );
        Gamma_PS(ispec1) = (U2/qa)*(Gamma_PS_std + Gamma_PS_flow);
    end
    Gamma_U(ispec1)= -sum(Nvec.*Flows(flow_ind:flow_ind+Smax));
    Gamma_L(ispec1)= -L11*Avec(ind_A) - L12*Avec(ind_A+1);
    Gammas(ispec1)=  Gamma_U(ispec1) + Gamma_L(ispec1) + Gamma_PS(ispec1);
end %species loop for flows

