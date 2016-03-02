function [QoTs]=calc_QoT_SN(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
    cmesh,emesh,D11_star_slice,D31_star_slice,D33_star_slice,log_D11_star_slice,Smax,flags,Avec,Bsq,U2,Flows,lmat,dTdrs,dndrs)
% QoTs=Qa/Ta=[Q_e/Te Q_i1/Ti1 ...]
% JL 9/21/11

%constants
elem_charge = 1.602176487e-19;

%expand flags structure
log_interp=flags.log_interp;
use_quad=flags.use_quad;
flux_cap=flags.flux_cap;

%allocate arrays
Nvec = zeros(Smax,1);
QoT_PS = zeros(num_species,1);
QoTs = zeros(num_species,1);

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
    
    norm_factor=na*ma^2*vta^3/(2*qa^2);     
    if flux_cap == 1
        % Cap radial coefficient at 0
        intfac=max(D11_star_slice  - (2/3)*U2*cmesh.' + D31_star_slice.^2./D33_star_slice,0);
        L21=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,intfac,log_interp,use_quad,2.5,0,norm_factor,flags);
        L22=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,intfac,log_interp,use_quad,3.5,0,norm_factor,flags);
    else       
        L21_1=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,1,cmesh,emesh,log_D11_star_slice,log_interp,use_quad,2.5,0,norm_factor,flags);
        L22_1=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,1,cmesh,emesh,log_D11_star_slice,log_interp,use_quad,3.5,0,norm_factor,flags);
        
        intfac=D31_star_slice.^2./D33_star_slice;      
        L21_2=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,intfac,log_interp,use_quad,2.5,0,norm_factor,flags);
        L22_2=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,intfac,log_interp,use_quad,3.5,0,norm_factor,flags);
        
        L21=L21_1+L21_2;
        L22=L22_1+L22_2;
    end    
    
    %Loop over primary Sonine index (k)
    norm_factor=na * ma*vta/qa * (2/3)*Bsq;
    for kval=0:Smax        
        Nvec(kval+1)=energy_conv3(0,kval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,D31_star_slice./D33_star_slice,log_interp,use_quad,2.5,0,norm_factor,flags);            
    end %k loop

    %Calculate PS flux


    QoT_PS_std = 0;
    %loop over ALL species
    for ispec2=1:num_species
        
        %specific species beta
        n_beta=ns(ispec2);
        T_beta=Ts(ispec2);
        dn_betadr=dndrs(ispec2);
        dT_betadr=dTdrs(ispec2);
        q_beta=charges(ispec2);
        
        C_PS1=(n_beta*dT_betadr + T_beta*dn_betadr)*elem_charge/(n_beta*q_beta);
        C_PS2=dT_betadr*elem_charge/q_beta;
        
        QoT_PS_std=QoT_PS_std+C_PS1 * ( (5/2)*lmat{ispec1,ispec2}(1,1) - lmat{ispec1,ispec2}(2,1) )...
            -C_PS2*( (5/2)*lmat{ispec1,ispec2}(1,2) - lmat{ispec1,ispec2}(2,2) );      
    end
    
    flow_ind=(ispec1-1)*(Smax+1)+1;
    ind_A=(ispec1-1)*3+1;

    if flux_cap == 1 
        QoT_PS(ispec1) = (U2/qa)*QoT_PS_std;
    else
        [I_1]=energy_conv_PS(vta,qa2,qb2,ma,na,nb,loglambda,vtb,2,na,flags);
        [I_2]=energy_conv_PS(vta,qa2,qb2,ma,na,nb,loglambda,vtb,3,na,flags);
        QoT_PS_flow=(2/3)*(ma*Ta*elem_charge/qa)*( I_1*Avec(ind_A) + I_2*Avec(ind_A+1) );
        QoT_PS(ispec1) = (U2/qa)*(QoT_PS_std + QoT_PS_flow);
    end
    QoT_U(ispec1)= -sum(Nvec.*Flows(flow_ind:flow_ind+Smax));
    QoT_L(ispec1)= -L21*Avec(ind_A) - L22*Avec(ind_A+1);    
    QoTs(ispec1) = QoT_U(ispec1) + QoT_L(ispec1) + QoT_PS(ispec1)  ;
end %species loop for flows

