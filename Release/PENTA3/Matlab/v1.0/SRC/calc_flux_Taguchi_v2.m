function [Gammas,QoTs]=calc_flux_Taguchi_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
    cmesh,emesh,D11_star_slice,D31_star_slice,Smax,flags,Avec,U2,Flows,lmat,B_0)
%This function calculates arrays of fluxes using the T method
% JL 2/2010

elem_charge = 1.602176487e-19;

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
    intfac=max( D11_star_slice - (2/3)*cmesh.'*U2,0);  %Don't allow negative L
    K_exp=1.5; 
    Lmat{ispec1}(1,1)=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,intfac,log_interp,use_quad,K_exp,nu_exp,norm_factor,'LT');      
    K_exp=1.5+1; 
    Lmat{ispec1}(1,2)=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,intfac,log_interp,use_quad,K_exp,nu_exp,norm_factor,'LT');
    Lmat{ispec1}(2,1)=Lmat{ispec1}(1,2);  %Use symmetry relation
    K_exp=1.5+2; 
    Lmat{ispec1}(2,2)=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,intfac,log_interp,use_quad,K_exp,nu_exp,norm_factor,'LT');
    
    
    %Loop over primary Sonine index (k)
    for kval=0:Smax

        %N_j+1,k+1
        norm_factor=na * ma/qa;
        K_exp=1; nu_exp=1;
        Nvec{ispec1}(1,kval+1)=energy_conv3(0,kval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D31_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,'31');        
        K_exp=2; nu_exp=1;
        Nvec{ispec1}(2,kval+1)=energy_conv3(0,kval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D31_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,'31');                
        
        %loop over secondary Sonine index (l) and sum
        for ispec2 = 1:num_species
            Mvec{ispec1,ispec2}(1,kval+1)=0;
            Mvec{ispec1,ispec2}(2,kval+1)=0;
            for lval=0:Smax
                
                c_l=(3/4)*sqrt(pi)*factorial(lval)/gamma(lval+5/2);

                norm_factor=na * ma/qa;
                K_exp=1; nu_exp=0;
                Jcoef{ispec1}(1,lval+1)=energy_conv3(0,lval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D31_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,'31');
                K_exp=2; nu_exp=0;
                Jcoef{ispec1}(2,lval+1)=energy_conv3(0,lval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D31_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,'31');                
                
                Mvec{ispec1,ispec2}(1,kval+1) = Mvec{ispec1,ispec2}(1,kval+1)...
                    + c_l*Jcoef{ispec1}(1,lval+1)*lmat{ispec1,ispec2}(lval+1,kval+1)/(na*ma);
                
                Mvec{ispec1,ispec2}(2,kval+1) = Mvec{ispec1,ispec2}(2,kval+1)...
                    + c_l*Jcoef{ispec1}(2,lval+1)*lmat{ispec1,ispec2}(lval+1,kval+1)/(na*ma);
            end
        end
        
    end %k loop

    %Calculate friction flux contributions
    for ispec2=1:num_species
        flow_ind2=(ispec2-1)*(Smax+1)+1;
        Ufric_Gamma(ispec2)=-sum(Mvec{ispec1,ispec2}(1,:).*Flows(flow_ind2:flow_ind2+Smax).');
        Ufric_QoT(ispec2)=  -sum(Mvec{ispec1,ispec2}(2,:).*Flows(flow_ind2:flow_ind2+Smax).');        
    end
    
    Gamma_fric(ispec1)=sum(Ufric_Gamma);
    QoT_fric(ispec1)=sum(Ufric_QoT);
    
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
    
    
    ind_A=(ispec1-1)*3+1;
    mono_flux(ispec1) = - Lmat{ispec1}(1,1)*Avec(ind_A) - Lmat{ispec1}(1,2)*Avec(ind_A+1);
    mono_QoT(ispec1)  = - Lmat{ispec1}(2,1)*Avec(ind_A) - Lmat{ispec1}(2,2)*Avec(ind_A+1);

    %calculate total flux
    flow_ind1=(ispec1-1)*(Smax+1)+1;
    flux_Ua(ispec1)=-sum(Nvec{ispec1}(1,:).*Flows(flow_ind1:flow_ind1+Smax).');
    QoT_Ua(ispec1)=-sum(Nvec{ispec1}(2,:).*Flows(flow_ind1:flow_ind1+Smax).');
    
    Gammas(ispec1)=flux_Ua(ispec1)+Gamma_fric(ispec1)+mono_flux(ispec1) + gamma_PS(ispec1);
    QoTs(ispec1)  =QoT_Ua(ispec1)+ QoT_fric(ispec1)+mono_QoT(ispec1) + QoT_PS(ispec1);
  
end %species loop for flows

