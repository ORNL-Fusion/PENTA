function [Er_roots]=find_Er_roots_PENTA(gamma_e,gamma_i,Er_test_vals,Z_ion,roa_surf)
%this file calculates the radial electric field from the ambipolarity
%contstraint for given ion and electron particle fluxes vs Er.  
%
% 6/2009 JL

plotit=0;

num_efield=length(Er_test_vals);

for is=1:length(Z_ion)
    Zmat(is,:)=Z_ion(is)*ones(1,num_efield);
end

gamma_all_i=sum(Zmat.*gamma_i(:,:),1);

%find where the sign switches
root_vr=zeros(1,num_efield);
for ie=1:(num_efield-1)
    a_test=sign(gamma_e(ie)- gamma_all_i(:,ie) );
    b_test=sign(gamma_e(ie+1)-gamma_all_i(ie+1));
    if a_test~=b_test
        root_vr(ie)=1;
    end
end
I_vr=find(root_vr==1);

%find number of roots at this r/a
num_roots=length(I_vr);

% Catch num_roots ~= 1,3
if num_roots~=1 && num_roots~=3
    figure
    hold on
    plot(Er_test_vals/100,gamma_all_i,'b.--')
    plot(Er_test_vals/100,gamma_e,'r.-')
    legend('Z*ion','electron')
    title(['r/a = ' num2str(roa_surf)])
    xlabel('Er (V/cm)')
%     ylabel('\Gamma m^{-2}s^{-1}')
    ylabel('J/n/e')
    box on
    error('Number of Er roots not equal to 1 or 3')
end


%find Er at the roots from linear interpolation
for iroot=1:length(I_vr)
    ind=I_vr(iroot);
    Er_roots(iroot)=fzero(@diff_flux,Er_test_vals(ind),[],Er_test_vals(ind-1:ind+1),gamma_e(ind-1:ind+1),gamma_all_i(ind-1:ind+1));
end


if plotit
    %plot very rough
    figure; hold on;box on;
    plot(Er_test_vals/100,gamma_all_i/1e20,'b.--')
    plot(Er_test_vals/100,gamma_e/1e20,'r.-')
    xlabel('Er (V/cm)')
%     ylabel('\Gamma (10^{20} m^{-2}s^{-1})')
    ylabel('J/n/e')
    plot(Er_test_vals(I_vr)./100,gamma_all_i(I_vr)/1e20,'ko')
    for iroot=1:length(I_vr)
        ind=I_vr(iroot);
        plot([Er_roots(iroot) Er_roots(iroot)]/100,[gamma_e(ind)*.5 gamma_e(ind)*2]/1e20,'g-')
    end
end


function [dflux]=diff_flux(Er,Er_range,gamma_e,gamma_all_i)
%Define gamma_e - gamma_i using linear interpolation
df=gamma_e-gamma_all_i;
if Er > Er_range(end)
    dflux=df(end);
elseif  Er < Er_range(1)
    dflux=df(1);
else
    dflux=interp1(Er_range,df,Er);
end

