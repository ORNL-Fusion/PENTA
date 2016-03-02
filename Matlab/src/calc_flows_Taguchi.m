function Flows=calc_flows_Taguchi(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
    cmesh,emesh,D31_star_slice,log_D33_star_slice,Smax,flags,Avec,lmat,Bsq,B_0)
%This function calculates the Sonine expanded flows using the Taguchi
%method.
%1/2010-9/2011 JL

%define inline function for Kronecker delta 
delta=inline('i==j','i','j');

%constants
elem_charge = 1.602176487e-19;

%expand flags structure
log_interp=flags.log_interp;
use_quad=flags.use_quad;

%Define c coefficients
c_l(1:Smax+1,1)=(3/4)*sqrt(pi)*factorial(0:Smax)./gamma(5/2:Smax+5/2);

%loop over species and calculate flows
for ispec1=1:num_species
    %species parameters
    ma=masses(ispec1);
    qa=charges(ispec1);
    qa2=charges(ispec1)^2;
    na=ns(ispec1);
    Ta=Ts(ispec1);
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
        Cmat1=energy_conv3(jval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,D31_star_slice,log_interp,use_quad,1,0,na,flags);
        Cmat2=energy_conv3(jval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,0,cmesh,emesh,D31_star_slice,log_interp,use_quad,2,0,na,flags);                
        ind_A=(ispec1-1)*3+1;
        if Avec(ind_A+2)==0
            Cmat3=0;
        else
            norm_factor=na * qa/(B_0*ma*vta);
            Cmat3=energy_conv3(jval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,1,cmesh,emesh,log_D33_star_slice,log_interp,use_quad,0.5,0,norm_factor,flags);
        end
        
        %Sum RHS
        ind_RHS=(ispec1-1)*(Smax+1)+jval+1;
        RHS(ind_RHS,1) =  -Cmat1*Avec(ind_A) - Cmat2*Avec(ind_A+1) + Cmat3*Avec(ind_A+2);        
        
        %Loop over primary Sonine index (k)
        for kval=0:Smax            
 
            %integrate H_{jk}^a, which uses D33^*
            norm_factor=na*qa/(vta*Ta*elem_charge);
            Hcoeff=energy_conv3(jval,kval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,1,cmesh,emesh,log_D33_star_slice,log_interp,use_quad,0.5,1,norm_factor,flags);
            
            Amat{ispec1}(jval+1,kval+1)=delta(jval,kval)*na*qa*Bsq/(c_l(kval+1)*Ta*elem_charge) - Hcoeff;                     
            
            %loop over secondary Sonine index (l) and sum
            for ispec2 =1:num_species
                Bmat{ispec1,ispec2}(jval+1,kval+1)=0;
                for lval=0:Smax
                   
                    %integrate J_{j,l}^a
                    norm_factor=qa/(Ta*elem_charge*ma*vta);
                    Jcoeff=energy_conv3(jval,lval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,1,cmesh,emesh,log_D33_star_slice,log_interp,use_quad,0.5,0,norm_factor,flags);         
                    
                    Bmat{ispec1,ispec2}(jval+1,kval+1)=Bmat{ispec1,ispec2}(jval+1,kval+1)...
                        - c_l(lval+1)*Jcoeff*lmat{ispec1,ispec2}(lval+1,kval+1);
                    
                end %l loop
            end %species 2 loop            
        end %k loop        
    end %j loop
end %species loop for flows        

%form matrix to solve for flows.
for ispec1=1:num_species    %row
    for ispec2=1:num_species  %column
        flow_mat_add=delta(ispec1,ispec2)*Amat{ispec1}+Bmat{ispec1,ispec2};
            
        ind1=Smax*(ispec1-1)+ispec1;
        ind2=Smax*(ispec2-1)+ispec2;
        flow_mat(ind1:ind1+Smax,ind2:ind2+Smax)=flow_mat_add;
    end
end

%Solve system to get flows.
Flows=flow_mat\RHS;