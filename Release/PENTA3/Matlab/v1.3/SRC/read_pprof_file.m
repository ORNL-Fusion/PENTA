function [pprof_info_allsurf]=...
    read_pprof_file(data_path,pprof_char,num_ion_species,plot_pprof)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function reads the file containing the plasma profile information,
%e.g. 'plasma_profiles_hsx.dat'
%1/26/2010 JL

%read plasma profile info (pprof = 'plasma profile')
pprof_data=dlmread([data_path '/plasma_profiles_' pprof_char '.dat']);

num_pprof_pts=pprof_data(1,1);  
ra_pprof=pprof_data(2:num_pprof_pts,1);
ne_pprof=pprof_data(2:num_pprof_pts,2)*1e18;   %make units m^-3
Te_pprof=pprof_data(2:num_pprof_pts,3);

for ispec=1:num_ion_species

    ni_pprof(ispec,:)=pprof_data(2:num_pprof_pts,4+2*(ispec-1))*1e18;       %make units m^-3
    Ti_pprof(ispec,:)=pprof_data(2:num_pprof_pts,5+2*(ispec-1));
    
end

%smooth cubic spline fit the profiles and derivatives
% pp_ne=csaps(ra_pprof,ne_pprof);
% pp_Te=csaps(ra_pprof,Te_pprof);
% pp_ne_prime=fnder(pp_ne);
% pp_Te_prime=fnder(pp_Te);

pp_ne=spline(ra_pprof,ne_pprof);
pp_Te=spline(ra_pprof,Te_pprof);
[breaks_ne,coefs_ne,l_ne,k_ne,d_ne] = unmkpp(pp_ne);
[breaks_Te,coefs_Te,l_Te,k_Te,d_Te] = unmkpp(pp_Te);
pp_ne_prime = mkpp(breaks_ne,repmat(k_ne-1:-1:1,d_ne*l_ne,1).*coefs_ne(:,1:k_ne-1),d_ne);
pp_Te_prime = mkpp(breaks_Te,repmat(k_Te-1:-1:1,d_Te*l_Te,1).*coefs_Te(:,1:k_Te-1),d_Te);


for ispec=1:num_ion_species
%     pp_Ti(ispec,:)=csaps(ra_pprof,Ti_pprof(ispec,:));
%     pp_Ti_prime(ispec,:)=fnder(pp_Ti(ispec,:));
    pp_Ti(ispec,:)=spline(ra_pprof,Ti_pprof(ispec,:));
    [breaks_Ti,coefs_Ti,l_Ti,k_Ti,d_Ti] = unmkpp(pp_Ti);    
    pp_Ti_prime = mkpp(breaks_Ti,repmat(k_Ti-1:-1:1,d_Ti*l_Ti,1).*coefs_Ti(:,1:k_Ti-1),d_Ti);    
    
%     pp_ni(ispec,:)=csaps(ra_pprof,ni_pprof(ispec,:));
%     pp_ni_prime(ispec,:)=fnder(pp_ni(ispec,:));
    pp_ni(ispec,:)=spline(ra_pprof,ni_pprof(ispec,:));
    [breaks_ni,coefs_ni,l_ni,k_ni,d_ni] = unmkpp(pp_ni);    
    pp_ni_prime = mkpp(breaks_ni,repmat(k_ni-1:-1:1,d_ni*l_ni,1).*coefs_ni(:,1:k_ni-1),d_ni); 
end


%make structure of pprof info for returning
pprof_info_allsurf.pp_ne=pp_ne;
pprof_info_allsurf.pp_Te=pp_Te;
pprof_info_allsurf.pp_ni=pp_ni;
pprof_info_allsurf.pp_Ti=pp_Ti;
pprof_info_allsurf.pp_ne_prime=pp_ne_prime;
pprof_info_allsurf.pp_Te_prime=pp_Te_prime;
pprof_info_allsurf.pp_ni_prime=pp_ni_prime;
pprof_info_allsurf.pp_Ti_prime=pp_Ti_prime;

if plot_pprof
    roa_plot=0:0.001:1;
    
    figure;hold on; box on;
    plot(roa_plot,ppval(pp_ne,roa_plot)./1e18,'r.-')
    for ispec=1:num_ion_species
        plot(roa_plot,ppval(pp_ni(ispec,:),roa_plot)./1e18,'b.-')
    end
    xlabel('\rho');ylabel('n (10^{18} m^{-3})')
    legend('n_e','ions')

    
    figure;hold on; box on;
    plot(roa_plot,ppval(pp_Te,roa_plot),'r.-')
    for ispec=1:num_ion_species
        plot(roa_plot,ppval(pp_Ti(ispec,:),roa_plot),'b.-')
    end
    xlabel('\rho');ylabel('T (eV)')
    legend('T_e','ions')
    axis([0,1,0,max([ppval(pp_Ti(ispec,:),roa_plot) ...
        ppval(pp_Te,roa_plot)])*1.2])
end
