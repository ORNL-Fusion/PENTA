function [Lmat,Mmat,Nmat]=define_thermal_mats(Er,ns,Ts,v_ths,charges,masses,loglambda,cmul_mat,efield_mat,coef2d_mat,Smax,log_interp,use_quad,DKES_limit)

elem_charge = 1.602176487e-19;

%define |Er|
abs_Er=abs(Er);

%total number of species
num_species=length(charges);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now integrate coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DKES_limit
    for ispec=1:num_species
        %species parameters
        ma=masses(ispec);
        qa2=charges(ispec)^2;
        na=ns(ispec);
        vta=v_ths(ispec);
        Ta=Ts(ispec);
        
        %all other species parameters
        qb2=charges(1:end~=ispec).^2;
        nb=ns(1:end~=ispec);
        vtb=v_ths(1:end~=ispec);
        
%         jvals=[1 0 -1];
        jvals=[1 2 3];
        logopt=0;mat_ind=1;
        [cmesh,emesh]=meshgrid(cmul_mat{mat_ind},efield_mat{mat_ind});
        logemesh=log(emesh);logcmesh=log(cmesh);
        coef2d_slice=coef2d_mat{mat_ind};
        L_vals{ispec}=energy_conv2(jvals,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh,log_interp,use_quad);
        
        %correct for different normalizing coefficient for L* and D11*
        %and for Z
        L_vals{ispec}=L_vals{ispec}*(1/qa2);
    end        
    
    
    
    
    %assemble L,M,N mats
    for ispec=1:num_species
        Lmat{ispec}=[L_vals{ispec}(1) -L_vals{ispec}(2);L_vals{ispec}(2) -L_vals{ispec}(3)];
        Mmat{ispec}=NaN;
        Nmat{ispec}=NaN;
    end
    
        

else
    
    
    for ispec=1:num_species
        %species parameters
        ma=masses(ispec);
        qa2=charges(ispec)^2;
        na=ns(ispec);
        vta=v_ths(ispec);
        
        %all other species parameters
        qb2=charges(1:end~=ispec).^2;
        nb=ns(1:end~=ispec);
        vtb=v_ths(1:end~=ispec);
        
        %integrate L
        jvals=1:3;
        skip=0;
        for ii=2:Smax
            jstart=jvals(end)+1+skip;
            jadd=jstart:jstart+1;
            jvals=[jvals jadd];
            skip=skip+1;
        end
        logopt=0;mat_ind=1;
        [cmesh,emesh]=meshgrid(cmul_mat{mat_ind},efield_mat{mat_ind});
        logemesh=log(emesh);logcmesh=log(cmesh);
        coef2d_slice=coef2d_mat{mat_ind};
        L_vals{ispec}=energy_conv2(jvals,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh,log_interp,use_quad);
        %correct for Z
        L_vals{ispec}=L_vals{ispec}*(elem_charge^2/qa2);
        
        %integrate M
        jend=3;
        for ii=2:Smax
            jend=jend+(ii+1);
        end
        jvals=1:jend;
        
        logopt=1;mat_ind=2;
        [cmesh,emesh]=meshgrid(cmul_mat{mat_ind},efield_mat{mat_ind});
        logemesh=log(emesh);logcmesh=log(cmesh);
        coef2d_slice=coef2d_mat{mat_ind};
        M_vals{ispec}=energy_conv2(jvals,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh,log_interp,use_quad);
        
        %integrate N
        %same jvals as above
        logopt=0;mat_ind=3;
        [cmesh,emesh]=meshgrid(cmul_mat{mat_ind},efield_mat{mat_ind});
        logemesh=log(emesh);logcmesh=log(cmesh);
        coef2d_slice=coef2d_mat{mat_ind};
        N_vals{ispec}=energy_conv2(jvals,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh,log_interp,use_quad);
        %correct for Z
        N_vals{ispec}=N_vals{ispec}*(elem_charge/sqrt(qa2));
        
        
    end
    
    
    %assemble L,M,N mats
    for ispec=1:num_species
        Lmat{ispec}=[L_vals{ispec}(1) -L_vals{ispec}(2);L_vals{ispec}(2) -L_vals{ispec}(3)];
        for ii=2:Smax
            Lmat{ispec}(ii+1,1:2)=[L_vals{ispec}(ii*2) -L_vals{ispec}(ii*2+1)];
        end
    end
    
    for ispec=1:num_species
        Mmat{ispec}=[M_vals{ispec}(1)  M_vals{ispec}(2);M_vals{ispec}(2)  M_vals{ispec}(3)];
        Nmat{ispec}=[N_vals{ispec}(1)  N_vals{ispec}(2);N_vals{ispec}(2)  N_vals{ispec}(3)];
        ind=4;
        for ii=2:Smax
            for jj=1:ii+1
                Mmat{ispec}(jj,ii+1)=M_vals{ispec}(ind);
                Mmat{ispec}(ii+1,jj)=M_vals{ispec}(ind);
                Nmat{ispec}(jj,ii+1)=N_vals{ispec}(ind);
                Nmat{ispec}(ii+1,jj)=N_vals{ispec}(ind);
                ind=ind+1;
            end
            
        end
    end
    
end