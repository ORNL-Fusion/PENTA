function [QoTs,DKES_QoT]=calc_QoT_MBT(num_species,abs_Er,ns,Ts,v_ths,charges,masses,dTdrs,dndrs,loglambda,...
    cmesh,emesh,D11_star_slice,log_D11_star_slice,D31_star_slice,Smax,B_0,flags,Avec,lmat,Flows,U2)
%This function calculates arrays of fluxes using the MBT method 
%The fluxes are defined as Gamma = [Gamma_e Gamma_i1 ...] 
% qoTs=qa/Ta=[q_e/Te q_i1/Ti1 ...]
%
% JL 2/2010-9/222/11

%constants
elem_charge = 1.602176487e-19;

%expand flags structure
log_interp=flags.log_interp;
use_quad=flags.use_quad;
flux_cap=flags.flux_cap;

%Define c coefficients
c_l(1:Smax+1,1)=(3/4)*sqrt(pi)*factorial(0:Smax)./gamma(5/2:Smax+5/2);

% Allocate arrays
QoT_Uflux_fric = zeros(num_species,1);
QoT_Uflux = zeros(num_species,1);
DKES_QoT = zeros(num_species,1);
friction_QoT = zeros(num_species,1);
QoT_PS = zeros(num_species,1);
QoTs = zeros(num_species,1);
QoT_coef_k = zeros(num_species,Smax+1);
QoT_coef_l = zeros(Smax+1,1);

%loop over species
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
    
    %Calculate flux without momentum conservation    
    norm_factor=na*ma^2*vta^3/(2*qa^2);     
    if flux_cap == 1
        % Cap radial coefficient at 0
        intfac=max(D11_star_slice  - (2/3)*U2*cmesh.',0);
        L21=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,intfac,log_interp,use_quad,2.5,0,norm_factor,flags);
        L22=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,intfac,log_interp,use_quad,3.5,0,norm_factor,flags);        
    else
        L21=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,1,cmesh,emesh,log_D11_star_slice,log_interp,use_quad,2.5,0,norm_factor,flags);
        L22=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,1,cmesh,emesh,log_D11_star_slice,log_interp,use_quad,3.5,0,norm_factor,flags);
    end
           
    ind_A=(ispec1-1)*3+1;
    if Avec(ind_A+2)==0
        L23=0;
    else
        norm_factor=ma*vta^2/(2*qa*B_0);
        L23=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,D31_star_slice,log_interp,use_quad,2,0,norm_factor,flags);        
    end
    
    DKES_QoT(ispec1) =-L21*Avec(ind_A) - L22*Avec(ind_A+1) + L23*Avec(ind_A+2);
        
    friction_QoT(1:num_species)=0;

    for kval=0:Smax
        
        %calculate flux due to species A flow (not friction component)
        norm_factor=na*ma/qa;
        QoT_coef_k(ispec1,kval+1)=energy_conv3(kval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,D31_star_slice,log_interp,use_quad,2,1,norm_factor,flags);
        
        %calculate friction contribution
        for lval=0:Smax             
            norm_factor=1/qa;                        
            QoT_coef_l(lval+1,1)= energy_conv3(lval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,D31_star_slice,log_interp,use_quad,2,0,norm_factor,flags);
        end
              
        for ispec2=1:num_species
            QoT_l_sum= sum(c_l(1:Smax+1,1).*QoT_coef_l(1:Smax+1,1) .*lmat{ispec1,ispec2}(1:Smax+1,kval+1));            
            flow_ind2=(ispec2-1)*(Smax+1)+1+kval;
            friction_QoT(ispec2)= friction_QoT(ispec2)  + Flows(flow_ind2)*QoT_l_sum;                    
        end %species 2 loop
    end %k loop
    
    %calculate flux due to flows
    flow_ind1=(ispec1-1)*(Smax+1)+1;
    QoT_Uflux_fric(ispec1)= - sum(friction_QoT(:)); 
    QoT_Uflux(ispec1)=-sum(QoT_coef_k(ispec1,:) .*Flows(flow_ind1:flow_ind1+Smax).');

    %PS flux

    
    QoT_PS_std = 0;
    %loop over ALL species
    for ispec2=1:num_species
        
        %let's call the specific species beta
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
    ind_A=(ispec1-1)*3+1;

    if flux_cap == 1 
        QoT_PS(ispec1) = (U2/qa)*QoT_PS_std;
    else
        [I_1]=energy_conv_PS(vta,qa2,qb2,ma,na,nb,loglambda,vtb,2,na,flags);
        [I_2]=energy_conv_PS(vta,qa2,qb2,ma,na,nb,loglambda,vtb,3,na,flags);
        QoT_PS_flow=(2/3)*(ma*Ta*elem_charge/qa)*( I_1*Avec(ind_A) + I_2*Avec(ind_A+1) );
        QoT_PS(ispec1) = (U2/qa)*(QoT_PS_std + QoT_PS_flow);
    end
    %Sum all flux components
    QoTs(ispec1)  =DKES_QoT(ispec1)  + QoT_Uflux(ispec1) + QoT_Uflux_fric(ispec1) + QoT_PS(ispec1);    
end %species loop for flux

