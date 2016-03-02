function [Lmat,Mmat,Nmat]=define_thermal_mats(Er,ns,v_ths,charges,masses,loglambda,cmul_mat,efield_mat,coef2d_mat,Smax,log_interp,use_quad)

%define |Er|
abs_Er=abs(Er);

%total number of species
num_species=length(charges);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now integrate coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    if Smax == 1
        jvals=1:3;
    elseif Smax == 2
        jvals=1:5;
    elseif Smax == 3
        jvals=[1:5 7:8];
    elseif Smax == 4
        jvals=[1:5 7:8 11:12];
    end
    logopt=0;mat_ind=1;
    [cmesh,emesh]=meshgrid(cmul_mat{mat_ind},efield_mat{mat_ind});
    logemesh=log(emesh);logcmesh=log(cmesh);
    coef2d_slice=coef2d_mat{mat_ind};
    L_vals{ispec}=energy_conv2(jvals,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh,log_interp,use_quad);
    
    %integrate M
    if Smax == 1
        jvals=1:3;
    elseif Smax == 2
        jvals=1:6;
    elseif Smax == 3
        jvals=1:10;
    elseif Smax == 4
        jvals=1:15;
    end
    logopt=1;mat_ind=2;
    [cmesh,emesh]=meshgrid(cmul_mat{mat_ind},efield_mat{mat_ind});
    logemesh=log(emesh);logcmesh=log(cmesh);
    coef2d_slice=coef2d_mat{mat_ind};
    M_vals{ispec}=energy_conv2(jvals,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh,log_interp,use_quad);
    
    %integrate N
    if Smax == 1
        jvals=1:3;
    elseif Smax == 2
        jvals=1:6;
    elseif Smax == 3
        jvals=1:10;
    elseif Smax == 4
        jvals=1:15;        
    end
    logopt=0;mat_ind=3;
    [cmesh,emesh]=meshgrid(cmul_mat{mat_ind},efield_mat{mat_ind});
    logemesh=log(emesh);logcmesh=log(cmesh);
    coef2d_slice=coef2d_mat{mat_ind};
    N_vals{ispec}=energy_conv2(jvals,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh,log_interp,use_quad);
    
    
end


%assemble L,M,N mats

%ions
for ispec=1:num_species
    Lmat{ispec}=[L_vals{ispec}(1) -L_vals{ispec}(2);L_vals{ispec}(2) -L_vals{ispec}(3)];
    Mmat{ispec}=[M_vals{ispec}(1)  M_vals{ispec}(2);M_vals{ispec}(2)  M_vals{ispec}(3)];
    Nmat{ispec}=[N_vals{ispec}(1)  N_vals{ispec}(2);N_vals{ispec}(2)  N_vals{ispec}(3)];
end


if Smax >= 2
       
    for ispec=1:num_species
        
        Lmat{ispec}(3,1)= L_vals{ispec}(4);
        Lmat{ispec}(3,2)=-L_vals{ispec}(5);
        
        Nmat{ispec}(1,3)=N_vals{ispec}(4);
        Nmat{ispec}(3,1)=N_vals{ispec}(4);
        Nmat{ispec}(2,3)=N_vals{ispec}(5);
        Nmat{ispec}(3,2)=N_vals{ispec}(5);
        Nmat{ispec}(3,3)=N_vals{ispec}(6);
        
        Mmat{ispec}(1,3)=M_vals{ispec}(4);
        Mmat{ispec}(3,1)=M_vals{ispec}(4);
        Mmat{ispec}(2,3)=M_vals{ispec}(5);
        Mmat{ispec}(3,2)=M_vals{ispec}(5);
        Mmat{ispec}(3,3)=M_vals{ispec}(6);
        
    end
    
end

if Smax >= 3
    for ispec=1:num_species
        
        Lmat{ispec}(4,1)= L_vals{ispec}(6);  %here 6 and 7 correspond to jvals 7:8 above
        Lmat{ispec}(4,2)=-L_vals{ispec}(7);
        
        Nmat{ispec}(1,4)=N_vals{ispec}(7);
        Nmat{ispec}(4,1)=N_vals{ispec}(7);
        Nmat{ispec}(2,4)=N_vals{ispec}(8);
        Nmat{ispec}(4,2)=N_vals{ispec}(8);
        Nmat{ispec}(3,4)=N_vals{ispec}(9);
        Nmat{ispec}(4,3)=N_vals{ispec}(9);
        Nmat{ispec}(4,4)=N_vals{ispec}(10);
        
        Mmat{ispec}(1,4)=M_vals{ispec}(7);
        Mmat{ispec}(4,1)=M_vals{ispec}(7);
        Mmat{ispec}(2,4)=M_vals{ispec}(8);
        Mmat{ispec}(4,2)=M_vals{ispec}(8);
        Mmat{ispec}(3,4)=M_vals{ispec}(9);
        Mmat{ispec}(4,3)=M_vals{ispec}(9);
        Mmat{ispec}(4,4)=M_vals{ispec}(10);
    end
end

if Smax >= 4
    for ispec=1:num_species
        
        Lmat{ispec}(5,1)= L_vals{ispec}(8);  %here 8 and 9 correspond to jvals 11:12 above
        Lmat{ispec}(5,2)=-L_vals{ispec}(9);
        
        Nmat{ispec}(1,5)=N_vals{ispec}(11);
        Nmat{ispec}(5,1)=N_vals{ispec}(11);
        Nmat{ispec}(2,5)=N_vals{ispec}(12);
        Nmat{ispec}(5,2)=N_vals{ispec}(12);
        Nmat{ispec}(3,5)=N_vals{ispec}(13);
        Nmat{ispec}(5,3)=N_vals{ispec}(13);
        Nmat{ispec}(4,5)=N_vals{ispec}(14);
        Nmat{ispec}(5,4)=N_vals{ispec}(14);
        Nmat{ispec}(5,5)=N_vals{ispec}(15);
        
        Mmat{ispec}(1,5)=M_vals{ispec}(11);
        Mmat{ispec}(5,1)=M_vals{ispec}(11);
        Mmat{ispec}(2,5)=M_vals{ispec}(12);
        Mmat{ispec}(5,2)=M_vals{ispec}(12);
        Mmat{ispec}(3,5)=M_vals{ispec}(13);
        Mmat{ispec}(5,3)=M_vals{ispec}(13);
        Mmat{ispec}(4,5)=M_vals{ispec}(14);
        Mmat{ispec}(5,4)=M_vals{ispec}(14);
        Mmat{ispec}(5,5)=M_vals{ispec}(15);
    end
end
