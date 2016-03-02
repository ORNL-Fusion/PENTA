clearvars;


Svals=[1:9];
% 
% for ii=1:length(Svals)
%     Smax=Svals(ii);
%     tic;
%     [roa_vals,Er_ambi{Smax},gamma_e{Smax},q_e{Smax},gamma_i{Smax},q_i{Smax},J_bs{Smax},J_E_e{Smax},J_E_i{Smax},J_E_cl,FFmat{Smax},X_vec{Smax},BUprl_e{Smax},BUprl_i{Smax}]=run_PENTA([],Smax);
%     disp(['Took ' num2str(toc) ' seconds.'])
%     times(ii)=toc;
% end



load convergence_data_s10_qhs1T_8_27_09_with_flows.mat

terms=Svals+1;

%unroll shit
for ii=1:length(Svals)
    ion_roots(ii)=Er_ambi{ii}(1);
    un_roots(ii)=Er_ambi{ii}(2);
    electron_roots(ii)=Er_ambi{ii}(3);
    
    gamma_e_ion(ii)=gamma_e{ii}(1);
    gamma_e_electron(ii)=gamma_e{ii}(3);
    
    q_e_ion(ii)=q_e{ii}(1);
    q_e_electron(ii)=q_e{ii}(3);
    
    q_i_ion(ii)=q_i{ii}(1);
    q_i_electron(ii)=q_i{ii}(3);
    
    J_bs_ion(ii)=J_bs{ii}(1);
    J_bs_electron(ii)=J_bs{ii}(3);
    
    BUprl_e_ion(ii)=BUprl_e{ii}(1);
    BUprl_e_electron(ii)=BUprl_e{ii}(3);
    
    BUprl_i_ion(ii)=BUprl_i{ii}(1);
    BUprl_i_electron(ii)=BUprl_i{ii}(3);    
end

figure;hold on;box on
plot(terms,ion_roots/100,'bo')
% plot(terms,un_roots/100,'go')
% plot(terms,electron_roots/100,'ro')
xlabel('Number of terms')
ylabel('E_r (V/cm)')
title('Ion root solutions')


figure;hold on;box on
plot(terms,q_i_ion,'ro')

figure;hold on;box on
plot(terms,J_bs_ion/1000,'ro')
xlabel('Number of terms')
ylabel('J_{bs} (kA/m^2)')
title('Ion root solutions')

figure;hold on;box on
plot(terms,times,'ro')
xlabel('Number of terms')
ylabel('Run Time (seconds)')


figure;hold on;box on
plot(terms(1:end-1),abs(diff(J_bs_ion)./J_bs_ion(2:end)*100),'r-o')
plot(terms(1:end-1),abs(diff(J_bs_electron)./J_bs_electron(2:end)*100),'r-s')
plot(terms(1:end-1),abs(diff(ion_roots)./ion_roots(2:end)*100),'ko-')
plot(terms(1:end-1),abs(diff(electron_roots)./electron_roots(2:end)*100),'ks-')
plot(terms(1:end-1),abs(diff(q_e_ion)./q_e_ion(2:end)*100),'m-o')
plot(terms(1:end-1),abs(diff(q_e_electron)./q_e_electron(2:end)*100),'m-s')
plot(terms(1:end-1),abs(diff(q_i_ion)./q_i_ion(2:end)*100),'b-o')
plot(terms(1:end-1),abs(diff(q_i_electron)./q_i_electron(2:end)*100),'b-s')
plot(terms(1:end-1),abs(diff(BUprl_i_ion)./BUprl_i_ion(2:end)*100),'g-o')
plot(terms(1:end-1),abs(diff(BUprl_e_ion)./BUprl_e_ion(2:end)*100),'g-s')
xlabel('Number of terms')
ylabel('Relative difference (%)')
set(gca,'yscale','log')
legend('J_{bs} iroot','J_{bs} eroot','E_r iroot','E_r eroot',...
    'q_e iroot','q_e eroot','q_i iroot','q_i eroot','<Bu_i_{||}>','<Bu_e_{||}>')

% plot(terms,J_bs_electron/1000,'ro')
% title('Electron root solutions')