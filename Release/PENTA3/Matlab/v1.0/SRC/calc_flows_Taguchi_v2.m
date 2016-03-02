function Flows=calc_flows_Taguchi_v2(num_species,abs_Er,ns,Ts,v_ths,charges,masses,loglambda,...
    cmesh,emesh,D31_star_slice,D33_star_slice,Smax,flags,Avec,lmat,Bsq,B_0)
%This function calculates the Sonine expanded flows using the Taguchi
%method.
%1/2010-5/2010 JL

%define inline function for Kronecker delta 
delta=inline('i==j','i','j');

%constants
elem_charge = 1.602176487e-19;

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
        
        %B_j+1,1 
        norm_factor=na * vta^2*ma/(2*qa*B_0);
        K_exp=1; nu_exp=0;       
        Cmat{ispec1}(jval+1,1)=energy_conv3(jval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D31_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,'31');
        
        %B_j+1,2
        norm_factor=na * vta^2*ma/(2*qa*B_0);
        K_exp=1+1; nu_exp=0;
        Cmat{ispec1}(jval+1,2)=energy_conv3(jval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D31_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,'31');

                
        ind_A=(ispec1-1)*3+1;
        if Avec(ind_A+2)==0
            Cmat{ispec1}(jval+1,3)=0;
        else
            norm_factor=na * vta/(2*B_0^2);
            K_exp=0.5; nu_exp=0;
            Cmat{ispec1}(jval+1,3)=energy_conv3(jval,0,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D33_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,'33');
        end
        
        %Sum RHS
%         ind_A=(ispec1-1)*3+1;
        ind_RHS=(ispec1-1)*(Smax+1)+jval+1;
        RHS(ind_RHS,1) =  -Cmat{ispec1}(jval+1,1)*Avec(ind_A) - Cmat{ispec1}(jval+1,2)*Avec(ind_A+1) - Cmat{ispec1}(jval+1,3)*Avec(ind_A+2);        
        
        %Loop over primary Sonine index (k)
        for kval=0:Smax            

            %calculate c coefficient for k
            c_k=(3/4)*sqrt(pi)*factorial(kval)/gamma(kval+5/2);
            
            %integrate H_{jk}^a, which uses D33^*
            norm_factor=na*vta/(2*B_0^2);
            K_exp=0.5; nu_exp=1;
            Hcoeff{ispec1}(jval+1,kval+1)=energy_conv3(jval,kval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D33_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,'33');
            
            Amat{ispec1}(jval+1,kval+1)=delta(jval,kval)*na*Bsq/(c_k*B_0) - B_0*ma/(Ta*elem_charge)*Hcoeff{ispec1}(jval+1,kval+1);                     
            
            %loop over secondary Sonine index (l) and sum
            for ispec2 =1:num_species
                Bmat{ispec1,ispec2}(jval+1,kval+1)=0;
                for lval=0:Smax
                    c_l=(3/4)*sqrt(pi)*factorial(lval)/gamma(lval+5/2);
                    
                    %integrate J_{j,l}^a
                    norm_factor=na*vta/(2*B_0^2);
                    K_exp=0.5; nu_exp=0;
                    Jcoeff{ispec1}(jval+1,lval+1)=energy_conv3(jval,lval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,D33_star_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,'33');
                    
                    Bmat{ispec1,ispec2}(jval+1,kval+1)=Bmat{ispec1,ispec2}(jval+1,kval+1)...
                        - B_0*ma/(Ta*elem_charge)*c_l*Jcoeff{ispec1}(jval+1,lval+1)*lmat{ispec1,ispec2}(lval+1,kval+1)/(na*ma);
                    
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
Flows=inv(flow_mat)*RHS;