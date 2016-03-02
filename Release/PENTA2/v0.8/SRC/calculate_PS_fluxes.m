function [gamma_PS,q_PS]=calculate_PS_fluxes(Z_ion,num_ion_species,lmat,Xvec,U2)

%constants
elem_charge = 1.602176487e-19;


%calculate PS fluxes
charges=elem_charge*[-1 Z_ion];
for ispec1=1:num_ion_species+1
    lX_sum=0;
    for ispec2=1:num_ion_species+1
        spec2_ind=(ispec2-1)*2+1;
        lX_sum=lX_sum+lmat{ispec1,ispec2}*Xvec(spec2_ind:spec2_ind+1)./charges(ispec2);
    end
    
    PS_fluxes=-U2*lX_sum/charges(ispec1);
    gamma_PS(ispec1)=PS_fluxes(1);
    q_PS(ispec1)=PS_fluxes(2);
end