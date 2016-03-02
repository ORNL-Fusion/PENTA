function [FFmat]=define_flux_force_mat(num_ion_species,lmat,Lmat,Mmat,Nmat,Bsq,ne,ni,Z_ion)

%constants
elem_charge = 1.602176487e-19;

num_species=1+num_ion_species;

%define flux matrix
for ispec1=1:num_species    %row
    for ispec2=1:num_species  %column
        if ispec1==ispec2
            flux_mat_add=Mmat{ispec1}*inv(Nmat{ispec1}) - lmat{ispec1,ispec2}*inv(Nmat{ispec1})*Bsq;
        else
            flux_mat_add=-lmat{ispec1,ispec2}*inv(Nmat{ispec2})*Bsq;
        
        end
        flux_mat(1+(ispec1-1)*2:1+(ispec1-1)*2+1,1+(ispec2-1)*2:1+(ispec2-1)*2+1)=flux_mat_add;
    end
end

%define force matrix
for ispec1=1:num_species    %row
    for ispec2=1:num_species  %column
        if ispec1==ispec2
            Xmat_add=Mmat{ispec1}*inv(Nmat{ispec1})*Lmat{ispec1} - lmat{ispec1,ispec2}*inv(Nmat{ispec1})*Lmat{ispec1}*Bsq-Nmat{ispec1};
        else
            Xmat_add=-lmat{ispec1,ispec2}*inv(Nmat{ispec2})*Lmat{ispec2}*Bsq;
        
        end
        Xmat(1+(ispec1-1)*2:1+(ispec1-1)*2+1,1+(ispec2-1)*2:1+(ispec2-1)*2+1)=Xmat_add;
    end
end

%add last column of Xmat
Xmat_end(1:2,1)=[-ne*elem_charge*sqrt(Bsq); 0;];
for ispec1=1:num_species-1 
    Xmat_end(3+(ispec1-1)*2:4+(ispec1-1)*2,1)=[ni(ispec1)*Z_ion(ispec1)*elem_charge*sqrt(Bsq); 0;];
end

Xmat=[Xmat Xmat_end];

FFmat=inv(flux_mat)*Xmat;