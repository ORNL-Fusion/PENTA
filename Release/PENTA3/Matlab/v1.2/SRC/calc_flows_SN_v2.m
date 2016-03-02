function Flows=calc_flows_SN_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
    cmesh,emesh,D31_star_slice,D33_star_slice,Smax,flags,Avec,lmat,Bsq,B_0)
%This function calculates an array of flows using the S-N method
%The flows are defined as [Ue,0 Ue,1 ... Ue,Smax Ui1,0 ...]
%where Ua,k = <Bu_||ak>/<B**2>
% JL 2/2010

%constants
elem_charge = 1.602176487e-19;
        coefname='';

%define inline function for Kronecker delta 
delta=inline('i==j','i','j');

%expand flags structure
log_interp=flags.log_interp;
use_quad=flags.use_quad;

%none of the coeffs currently use logopt
logopt=0;

%loop over species and calculate flows
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
    
    
    %The j loop defines each equation for the current species in the system
    for jval=0:Smax
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Define the RHS of the equation system
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       
        %B_j+1,1
        norm_factor=na * vta*ma*B_0/qa;
        K_exp=0.5+1; nu_exp=0;
        intfac=D31_star_slice./D33_star_slice;

        Cmat{ispec1}(jval+1,1)=energy_conv3(jval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,intfac,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);

        %B_j+1,2
        norm_factor=na * vta*ma*B_0/qa;
        K_exp=0.5+2; nu_exp=0;
        intfac=D31_star_slice./D33_star_slice;
        Cmat{ispec1}(jval+1,2)=energy_conv3(jval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,intfac,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);                     
        
        %Sum RHS
        ind_A=(ispec1-1)*3+1;
        ind_RHS=(ispec1-1)*(Smax+1)+jval+1;
        RHS(ind_RHS,1) = -Ta*elem_charge*( Cmat{ispec1}(jval+1,1)*Avec(ind_A)   + ...
                                                         Cmat{ispec1}(jval+1,2)*Avec(ind_A+1) - ...
                                                         delta(jval,0)*(3/2)*na*Avec(ind_A+2)...
                                                     );                                
        %Loop over primary Sonine index (k)
        for kval=0:Smax
            
               %first UA term
            norm_factor=na*ma*B_0*vta;
            K_exp=1.5; nu_exp=0;
            intfac=(2*Bsq/3./cmesh.'-D33_star_slice).*cmesh.'./D33_star_slice;
            Dmat{ispec1}(jval+1,kval+1)=energy_conv3(jval,kval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,intfac,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);


        end %k loop
        
    end %j loop
end %species loop for flows

%Form matrix to solve for flows.
for ispec1=1:num_species    %row
    for ispec2=1:num_species  %column
        flow_mat_add=delta(ispec1,ispec2)*Dmat{ispec1}-B_0*(3/2)*lmat{ispec1,ispec2};
        ind1=Smax*(ispec1-1)+ispec1;
        ind2=Smax*(ispec2-1)+ispec2;
        flow_mat(ind1:ind1+Smax,ind2:ind2+Smax)=flow_mat_add;
    end
end

%Solve system to get flows.
Flows=inv(flow_mat)*RHS;