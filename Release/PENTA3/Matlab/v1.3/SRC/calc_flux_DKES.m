function [Gammas,QoTs]=calc_flux_DKES(num_species,abs_Er,ns,v_ths,charges,masses,loglambda,...
    cmesh,emesh,D11_star_slice,D13_star_slice,flags,Avec,B_0)
%This function calculates arrays of fluxes using the DKES method
% JL 2/2010


%expand flags structure
log_interp=flags.log_interp;
use_quad=flags.use_quad;


log_D11_star_slice = log(D11_star_slice);

%loop over species and calculate fluxes
for ispec1=1:num_species
    %species parameters
    ma=masses(ispec1);
    qa=charges(ispec1);
    qa2=charges(ispec1)^2;
    na=ns(ispec1);
    vta=v_ths(ispec1);
    
    %all other species parameters
    qb2=charges(1:end~=ispec1).^2;
    nb=ns(1:end~=ispec1);
    vtb=v_ths(1:end~=ispec1);
    
     %L11
    norm_factor=vta^3*ma^2*0.5/qa^2;
    [L11]=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,1,cmesh,emesh,log_D11_star_slice,log_interp,use_quad,1.5,0,norm_factor,flags);
    [L12]=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,1,cmesh,emesh,log_D11_star_slice,log_interp,use_quad,2.5,0,norm_factor,flags);
    [L22]=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,1,cmesh,emesh,log_D11_star_slice,log_interp,use_quad,3.5,0,norm_factor,flags);
    L21=L12;
    
    ind_A=(ispec1-1)*3+1;
    if Avec(ind_A+2)==0
        L13=0;
        L23=0;
    else
        norm_factor=vta^2*ma*0.5/qa/B_0;
        [L13]=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,D13_star_slice,log_interp,use_quad,1,0,norm_factor,flags);
        [L23]=energy_conv3(0,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,D13_star_slice,log_interp,use_quad,2,0,norm_factor,flags);
    end
    %Calculate flux without momentum conservation    
    Gammas(ispec1)=-na*(L11*Avec(ind_A)+ L12*Avec(ind_A+1)+ L13*Avec(ind_A+2));
    QoTs(ispec1) =-na*(L21*Avec(ind_A)+ L22*Avec(ind_A+1)+ L23*Avec(ind_A+2));

    
end %species loop for flows

