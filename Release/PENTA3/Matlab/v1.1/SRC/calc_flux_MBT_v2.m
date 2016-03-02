function [Gammas,QoTs,DKES_flux,Uflux,Uflux_fric,PS_flux,DKES_QoT]=calc_flux_MBT_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,dTdrs,dndrs,loglambda,...
    cmesh,emesh,D11_star_slice,D31_star_slice,D13_star_slice,Smax,B_0,flags,Avec,lmat,Flows,G2)
%This function calculates arrays of fluxes using the MBT method 
%The fluxes are defined as Gamma = [Gamma_e Gamma_i1 ...] 
% qoTs=qa/Ta=[q_e/Te q_i1/Ti1 ...]
%
% JL 2/2010

%constants
elem_charge = 1.602176487e-19;

%expand flags structure
log_interp=flags.log_interp;
use_quad=flags.use_quad;

%none of the coeffs currently use logopt
logopt=0;

%Define c coefficients
c_l(1:Smax+1,1)=(3/4)*sqrt(pi)*factorial(0:Smax)./gamma(5/2:Smax+5/2);

%loop over species
for ispec1=1:num_species

    %species parameters
    ma=masses(ispec1);
    qa=charges(ispec1);
    qa2=charges(ispec1)^2;
    na=ns(ispec1);
    vta=v_ths(ispec1);
    Ta=Ts(ispec1);
    dTadr=dTdrs(ispec1);
    dnadr=dndrs(ispec1);
    
    %all other species parameters
    qb2=charges(1:end~=ispec1).^2;
    nb=ns(1:end~=ispec1);
    vtb=v_ths(1:end~=ispec1);
    
    %Calculate the thermal monoenergetic coefficients without MC
    
    %L11
    norm_factor=vta^3*ma^2*0.5/qa^2;
    K_exp=1.5; nu_exp=0;
    coefname='11';
    L11=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D11_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);
    
    %L12
    norm_factor=vta^3*ma^2*0.5/qa^2;
    K_exp=1.5+1; nu_exp=0;
    coefname='11';
    L12=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D11_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);
    
    %L21
    L21=L12;
   
    %L22
    norm_factor=vta^3*ma^2*0.5/qa^2;
    K_exp=1.5+2; nu_exp=0;
    coefname='11';
    L22=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D11_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);
    
    
    ind_A=(ispec1-1)*3+1;
    
    if Avec(ind_A+2)==0
        L13=0;
        L23=0;
    else
    %L13
    norm_factor=vta^2*ma*0.5/qa/B_0;
    K_exp=1.0; nu_exp=0;
    coefname='13';
    L13=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D13_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);
    
    %L23
    norm_factor=vta^2*ma*0.5/qa/B_0;
    K_exp=1+1; nu_exp=0;
    coefname='13';
    L23=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D13_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);
    
    end
    %Calculate flux without momentum conservation
%     ind_A=(ispec1-1)*3+1;
    DKES_flux(ispec1)=-na*(L11*Avec(ind_A)+ L12*Avec(ind_A+1)+ L13*Avec(ind_A+2));
    DKES_QoT(ispec1) =-na*(L21*Avec(ind_A)+ L22*Avec(ind_A+1)+ L23*Avec(ind_A+2));
    
    
    friction_flux(1:num_species)=0;
    friction_QoT(1:num_species)=0;

    for kval=0:Smax
        
        %calculate flux due to species A flow (not friction component)
        %nu_a*D31*
        norm_factor=na*ma/qa;
        K_exp=1; nu_exp=1;
        coefname='31';
        flux_coef_k(ispec1,kval+1)=energy_conv3(kval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D31_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);
        K_exp=1+1; nu_exp=1;
        coefname='31';
        QoT_coef_k(ispec1,kval+1)=energy_conv3(kval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D31_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);
        
        %calculate friction contribution
        for lval=0:Smax             
            %D31*
            norm_factor=na*ma/qa;                        
            K_exp=1; nu_exp=0;
            coefname='31';
            flux_coef_l{ispec1}(lval+1,1)=energy_conv3(lval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D31_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);
            K_exp=1+1; nu_exp=0;
            coefname='31';
            QoT_coef_l{ispec1}(lval+1,1)= energy_conv3(lval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D31_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);
        end
              
        for ispec2=1:num_species

            flux_l_sum=sum(c_l(1:Smax+1,1).*flux_coef_l{ispec1}(1:Smax+1,1).*lmat{ispec1,ispec2}(1:Smax+1,kval+1)/na/ma);
            QoT_l_sum= sum(c_l(1:Smax+1,1).*QoT_coef_l{ispec1}(1:Smax+1,1) .*lmat{ispec1,ispec2}(1:Smax+1,kval+1)/na/ma);
            
            flow_ind2=(ispec2-1)*(Smax+1)+1+kval;

            friction_flux(ispec2)=friction_flux(ispec2) + Flows(flow_ind2)*flux_l_sum;
            friction_QoT(ispec2)= friction_QoT(ispec2)  + Flows(flow_ind2)*QoT_l_sum;
            
            
        end %species 2 loop
    end %k loop
    
    %calculate flux due to flows
    flow_ind1=(ispec1-1)*(Smax+1)+1;
    
    Uflux_fric(ispec1)    = - sum(friction_flux(:));
    QoT_Uflux_fric(ispec1)= - sum(friction_QoT(:));
    
    Uflux(ispec1)    =-sum(flux_coef_k(ispec1,:).*Flows(flow_ind1:flow_ind1+Smax).');    
    QoT_Uflux(ispec1)=-sum(QoT_coef_k(ispec1,:) .*Flows(flow_ind1:flow_ind1+Smax).');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate PS flux
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Define coefficient
    PS_coeff=na*G2*ma/qa/B_0^2;
    
    %integrate I coeffs
    nu_exp=1;
    norm_factor=1;
    
    cmesh_PS=[1e-20 1e10;1e-20 1e10];
    emesh_PS=[1e-20 1e-20;1e10 1e10];
    coef2d_slice_PS=[1 1;1 1];
    
    K_exp=1.5;
    coefname='I';
    [I_0]=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh_PS,emesh_PS,coef2d_slice_PS,0,0,K_exp,nu_exp,norm_factor,coefname);
    
    K_exp=2.5;
    [I_1]=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh_PS,emesh_PS,coef2d_slice_PS,0,0,K_exp,nu_exp,norm_factor,coefname);
    
    K_exp=3.5;
    [I_2]=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh_PS,emesh_PS,coef2d_slice_PS,0,0,K_exp,nu_exp,norm_factor,coefname);
    
    ind_A=(ispec1-1)*3+1;
    term4(1,1)=(2/3)*(Ta*elem_charge/qa)*( I_0*Avec(ind_A) + I_1*Avec(ind_A+1) );
    term4(2,1)=(2/3)*(Ta*elem_charge/qa)*( I_1*Avec(ind_A) + I_2*Avec(ind_A+1) );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate PS flux new method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Calculate effective friction force
    
    big_PS_sum=[0;0];
    %loop over ALL species
    for ispec2=1:num_species
        
        %let's call the specific species beta
        n_beta=ns(ispec2);
        T_beta=Ts(ispec2);
        dn_betadr=dndrs(ispec2);
        dT_betadr=dTdrs(ispec2);
        q_beta=charges(ispec2);
        
        if ispec1==ispec2
            C_PS1=0;
        else
            C_PS1=( (na*dTadr+Ta*dnadr)*elem_charge/(na*qa) - (n_beta*dT_betadr+T_beta*dn_betadr)*elem_charge/(n_beta*q_beta) );
        end
        C_PS2=dT_betadr*elem_charge/(q_beta);
        
        big_PS_sum(1,1)=big_PS_sum(1,1)+C_PS1*lmat{ispec1,ispec2}(1,1)+C_PS2*lmat{ispec1,ispec2}(1,2);
        big_PS_sum(2,1)=big_PS_sum(2,1)+C_PS1 * ( (5/2)*lmat{ispec1,ispec2}(1,1) - lmat{ispec1,ispec2}(2,1) )...
            +C_PS2*( (5/2)*lmat{ispec1,ispec2}(1,2) - lmat{ispec1,ispec2}(2,2) );      
    end
    X_PS2=-1/(na*ma)*big_PS_sum+term4;
    
    PS_flux(ispec1)=PS_coeff*X_PS2(1,1);
    PS_QoT(ispec1)=PS_coeff*X_PS2(2,1);
    
    %Sum all flux components
    Gammas(ispec1)=DKES_flux(ispec1)+Uflux(ispec1)+Uflux_fric(ispec1)+PS_flux(ispec1);
    QoTs(ispec1)=DKES_QoT(ispec1)+QoT_Uflux(ispec1)+QoT_Uflux_fric(ispec1)+PS_QoT(ispec1);
    
end %species loop for flux

