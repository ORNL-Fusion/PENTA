clearvars;

% run_ident='qhs1T_parabola';
num_ion_species=2;
pprof_char='2';

path_fortran='D:\Transport\PENTA\tests\7_2009_PENTA_Matlab_impurity_tests\qhs1T_parabola';
path_matlab='D:\Transport\PENTA\tests\7_2009_PENTA_Matlab_impurity_tests\qhs1T_parabola';




fvEr_data_f=dlmread([path_fortran '/fluxes_vs_Er'],'',2,0);


roa_f=fvEr_data_f(:,1);
Er_f=fvEr_data_f(:,2);
gamma_e_f=fvEr_data_f(:,3);

tmp_ind=4:4+num_ion_species-1;
gamma_i_f=fvEr_data_f(:,tmp_ind);
tmp_ind=4+num_ion_species;
q_e_f=fvEr_data_f(:,tmp_ind);
tmp_ind=5+num_ion_species:5+2*num_ion_species-1;
q_i_f=fvEr_data_f(:,tmp_ind);

figure;hold on;box on;
plot(Er_f,gamma_e_f/1e20,'r.-')
plot(Er_f,gamma_i_f/1e20,'.-')
xlabel('Er')
title(['r/a = ' num2str(roa_f(1))])
ylabel('10^2^0 \Gamma')

figure;hold on;box on;
plot(Er_f,gamma_e_f/1e20,'r.-')
plot(Er_f,gamma_i_f(:,1)/1e20+6*gamma_i_f(:,2)/1e20,'.-')
xlabel('Er')
title(['r/a = ' num2str(roa_f(1))])
ylabel('10^2^0 sum(\Gamma/q)')

figure;hold on;box on;
plot(Er_f,q_e_f,'r.-')
plot(Er_f,q_i_f,'.-')
xlabel('Er')
title(['r/a = ' num2str(roa_f(1))])
ylabel('q')









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PROFILES

% %r/a Te ne dnedr dTedr Ti ni dnidr dTidr
% pprof_data_f=dlmread([path_fortran '/plasma_profiles_check'],'',2,0);
% 
% roa_f=pprof_data_f(:,1);
% Te_f=pprof_data_f(:,2);
% ne_f=pprof_data_f(:,3);
% dnedr_f=pprof_data_f(:,4);
% dTedr_f=pprof_data_f(:,5);
% 
% tmp_ind=6:6+num_ion_species-1;
% Ti_f=pprof_data_f(:,tmp_ind);
% tmp_ind=tmp_ind+num_ion_species;
% ni_f=pprof_data_f(:,tmp_ind);
% tmp_ind=tmp_ind+num_ion_species;
% dnidr_f=pprof_data_f(:,tmp_ind);
% tmp_ind=tmp_ind+num_ion_species;
% dTidr_f=pprof_data_f(:,tmp_ind);
% 
% %load matlab output
% load([path_matlab '/pprof_check.mat'])
% roa_m=roa_vals;
% Te_m=Te;
% ne_m=ne;
% dTedr_m=dTedr;
% dTidr_m=dTidr;
% dnedr_m=dnedr;
% dnidr_m=dnidr;
% Ti_m=Ti;
% ni_m=ni;
% 
% %read actual input file
% pprof_data=dlmread([path_fortran '/plasma_profiles' pprof_char '.dat']);
% num_pprof_pts=pprof_data(1,1);  
% ra_pprof=pprof_data(2:num_pprof_pts,1);
% ne_pprof=pprof_data(2:num_pprof_pts,2)*1e18;   %make units m^-3
% Te_pprof=pprof_data(2:num_pprof_pts,3);
% 
% for ispec=1:num_ion_species
%         Ti_pprof(ispec,:)=pprof_data(2:num_pprof_pts,4+2*(ispec-1));
%         ni_pprof(ispec,:)=pprof_data(2:num_pprof_pts,5+2*(ispec-1))*1e18;   %make units m^-3
% end
% 
% 
% colors=['r','b','g','m'];
% %plot stuff
% figure; hold on; box on;
% plot(roa_m,Te_m,'k.')
% plot(roa_f,Te_f,'ko')
% plot(ra_pprof,Te_pprof,'m.-')
% plot(ra_pprof,Ti_pprof,'m.-')
% for is=1:num_ion_species
%     plot(roa_m,Ti_m(:,is),[colors(is) '.'])
%     plot(roa_f,Ti_f(:,is),[colors(is) 'o'])
% end
% xlabel('r/a'),ylabel('T (eV)')
% 
% %density
% figure; hold on; box on;
% plot(roa_m,ne_m,'k.')
% plot(roa_f,ne_f,'ko')
% plot(ra_pprof,ne_pprof,'m.-')
% plot(ra_pprof,ni_pprof,'m.-')
% for is=1:num_ion_species
%     plot(roa_m,ni_m(:,is),[colors(is) '.'])
%     plot(roa_f,ni_f(:,is),[colors(is) 'o'])
% end
% xlabel('r/a'),ylabel('n')
% 
% %dTdr
% figure; hold on; box on;
% plot(roa_m,dTedr_m,'k.')
% plot(roa_f,dTedr_f,'ko')
% for is=1:num_ion_species
%     plot(roa_m,dTidr_m(:,is),[colors(is) '.'])
%     plot(roa_f,dTidr_f(:,is),[colors(is) 'o'])
% end
% xlabel('r/a'),ylabel('dTdr (eV/m)')
% 
% %dndr
% figure; hold on; box on;
% plot(roa_m,dnedr_m,'k.')
% plot(roa_f,dnedr_f,'ko')
% for is=1:num_ion_species
%     plot(roa_m,dnidr_m(:,is),[colors(is) '.'])
%     plot(roa_f,dnidr_f(:,is),[colors(is) 'o'])
% end
% xlabel('r/a'),ylabel('dndr')
