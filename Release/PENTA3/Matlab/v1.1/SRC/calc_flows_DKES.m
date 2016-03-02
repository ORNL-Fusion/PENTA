function Flows=calc_flows_DKES(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
    cmesh,emesh,D31_star_slice,Smax,flags,Avec,B_0)
%This function calculates arrays of flows using the DKES method
% JL 2/2010

%expand flags structure
log_interp=flags.log_interp;
use_quad=flags.use_quad;

%none of the coeffs currently use logopt
logopt=0;

coefname='';

%loop over species and calculate flows

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
    
    %The j loop defines each equation for the current species in the system
    for jval=0:Smax
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Define the RHS of the equation system
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Integrate the three coefficients              
        
        %B_j+1,1 
        norm_factor=vta^2*ma/(2*qa*B_0);
        K_exp=1; nu_exp=0;
        Cmat{ispec1}(jval+1,1)=energy_conv3(jval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D31_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);
        
        %B_j+1,2
        norm_factor=vta^2*ma/(2*qa*B_0);
        K_exp=1+1; nu_exp=0;
        Cmat{ispec1}(jval+1,2)=energy_conv3(jval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D31_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);

        ind_A=(ispec1-1)*3+1;

        if Avec(ind_A+2)==0
            Cmat{ispec1}(jval+1,3)=0;
        else
            %B_j+1,3
            norm_factor=vta/(2*B_0^2);
            K_exp=0.5; nu_exp=0;
            Cmat{ispec1}(jval+1,3)=energy_conv3(jval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D31_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,coefname);
        end
        %Sum RHS
%         ind_A=(ispec1-1)*3+1;
        ind_RHS=(ispec1-1)*(Smax+1)+jval+1;
        RHS(ind_RHS,1) =  -Cmat{ispec1}(jval+1,1)*Avec(ind_A) - Cmat{ispec1}(jval+1,2)*Avec(ind_A+1) - Cmat{ispec1}(jval+1,3)*Avec(ind_A+2);
        
        
    end %j loop
end %species loop for flows
        
%Solve system to get flows.
Flows=RHS;