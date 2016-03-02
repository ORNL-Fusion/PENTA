function [FFmat]=define_flux_force_mat(num_ion_species,lmat,Lmat,Mmat,Nmat,Bsq,ne,ni,Z_ion,Smax)

%constants
elem_charge = 1.602176487e-19;

num_species=1+num_ion_species;

%define flux matrix
for ispec1=1:num_species    %row
    for ispec2=1:num_species  %column
        if ispec1==ispec2     %diagonals
            flux_mat_add=( Mmat{ispec1} - Bsq*lmat{ispec1,ispec2} ) * inv(Nmat{ispec1});
        else
            flux_mat_add=-Bsq*lmat{ispec1,ispec2}*inv(Nmat{ispec2});
        end
        ind1=1+(ispec1-1)*2+(ispec1-1)*(Smax-1);
        ind2=1+(ispec2-1)*2+(ispec2-1)*(Smax-1);
        flux_mat(ind1:ind1+Smax,ind2:ind2+Smax)=flux_mat_add;
    end
end

%define force matrix
for ispec1=1:num_species    %row
    for ispec2=1:num_species  %column
        if ispec1==ispec2
            Ntmp=[Nmat{ispec1}(1:Smax+1,1), -Nmat{ispec1}(1:Smax+1,2)];
            Xmat_add=( Mmat{ispec1} - Bsq*lmat{ispec1,ispec2} ) * inv(Nmat{ispec1})*Lmat{ispec1} - Ntmp;
        else
            Xmat_add=-Bsq*lmat{ispec1,ispec2}*inv(Nmat{ispec2})*Lmat{ispec2};
        
        end
        ind1=1+(ispec1-1)*2+(ispec1-1)*(Smax-1); 
        ind2=1+(ispec2-1)*2;
        Xmat(ind1:ind1+Smax,ind2:ind2+1)=Xmat_add;
    end
end

%add last column of Xmat
Xmat_end(1:Smax+1,1)=[-ne*elem_charge*sqrt(Bsq); zeros(Smax,1);];
for ispec1=1:num_species-1 
    ind1=Smax+2+(ispec1-1)*2+(ispec1-1)*(Smax-1);
    Xmat_end(ind1:ind1+Smax,1)=[ni(ispec1)*Z_ion(ispec1)*elem_charge*sqrt(Bsq); zeros(Smax,1);];
end

Xmat=[Xmat Xmat_end];

FFmat=inv(flux_mat)*Xmat;