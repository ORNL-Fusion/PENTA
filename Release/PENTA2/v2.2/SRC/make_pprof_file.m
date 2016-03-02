clear all;


% ts_date='9_24_07'; ts_shots=[27 43 51 52 53 58 59 60 62 63 64 66 67 68]; ts_seven=0; %44kW qhs1T
% ts_plotit=1; ts_gate_width=50; ts_atten=0; ts_num_iter=0; 
% [roa_ts,Te,Te_std_low,Te_std_high,ne,ne_std_low,ne_std_high] = get_te_ne_ml_JL_saved_vars(ts_date,ts_shots,ts_plotit,ts_gate_width,ts_atten,ts_num_iter,ts_seven);

% %44kW QHS B=1T
% roa_ts=[   0.016106504818415   0.098684871626455   0.194538984807418   0.290170400745854   0.384965393926578   0.480236971163569   0.578547082093118   0.679871516191837   0.787234845650106   0.906017548124103];
% Te=[1632        1363         891         394         347         203         121         121          76          60];
% ne=1e12*[   4.863000000000000   4.288000000000000   4.514000000000000   4.447000000000000   4.163000000000000   3.509000000000000   2.575000000000000   2.301000000000000   1.760000000000000   1.471000000000000];
% ni1=1e12*[   4.863000000000000   4.288000000000000   4.514000000000000   4.447000000000000   4.163000000000000   3.509000000000000   2.575000000000000   2.301000000000000   1.760000000000000   1.471000000000000];

%44kW QHS B=1T ne=4
ts_date='8_14_09'; ts_shots=[22 27 29 32 33 36:38 40  74 76 78:80 82:83 85:88];  %44kW qhs1T 8/14/09
[roa_ts,Te,Te_std_low,Te_std_high,ne,ne_std_low,ne_std_high] = make_TS_profile(ts_date,ts_shots);
Ti1=70-30*roa_ts;
ni1=ne;

% % 44kW QHS B=1T ne=5
% ts_date='8_5_09'; ts_shots=[ 44 47 53 54 55 56 58 64 69 70 75  78 80];  %44kW qhs1T 8/5/09 ne=5
% [roa_ts,Te,Te_std_low,Te_std_high,ne,ne_std_low,ne_std_high] = make_TS_profile(ts_date,ts_shots);
% Ti1=70-30*roa_ts;
% ni1=ne;


% %Parabolic - 2 species
% roa_ts=0:.05:1;
% Te=1000-500*roa_ts.^2;
% ne=10e12-1e12*roa_ts.^2;
% Ti1=1000-500*roa_ts.^2;
% Ti2=1000-200*roa_ts.^2;
% ni1=ne*1;
% ni2=ne*.0;


% %Parabolic - 1 species
% roa_ts=0:.05:1;
% Te=3000-1000*roa_ts.^2;
% ne=10e12-1e12*roa_ts.^2;
% Ti1=3000-1000*roa_ts.^2;
% ni1=ne;


% endpath='D:\Transport\PENTA\tests\7_2009_PENTA_Matlab_impurity_tests\qhs1T_parabola';
% endpath='D:\Transport\ideal_tokamak\PENTA_data';
endpath='D:\Transport\PENTA\penta_runs\PENTA_I_runs\qhs1T_44kW_ne4_8_2009';
% endpath='D:\Transport\PENTA\penta_runs\PENTA_I_runs\qhs1T_44kW_ne5_8_2009';


pprof_char='a';

fid=fopen([endpath '\plasma_profiles' pprof_char '.dat'],'w');


numpts=200;

s_for_fit=linspace(0,1,numpts);
rho_for_fit=sqrt(s_for_fit);

% % use doublesided fit to profiles
roa_dbl=[-roa_ts(end:-1:1),roa_ts];

Te_dbl=[Te(end:-1:1),Te];
ppTe=csaps(roa_dbl,Te_dbl);

Ti1_dbl=[Ti1(end:-1:1),Ti1];
ppTi1=csaps(roa_dbl,Ti1_dbl);

% Ti2_dbl=[Ti2(end:-1:1),Ti2];
% ppTi2=csaps(roa_dbl,Ti2_dbl);

ne_dbl=[ne(end:-1:1),ne];
ppne=csaps(roa_dbl,ne_dbl,.99);

ni1_dbl=[ni1(end:-1:1),ni1];
ppni1=csaps(roa_dbl,ni1_dbl,.99);

% ni2_dbl=[ni2(end:-1:1),ni2];
% ppni2=csaps(roa_dbl,ni2_dbl);


%evaluate at input points
Te_at_test=ppval(ppTe,rho_for_fit);
ne_at_test=ppval(ppne,rho_for_fit)/1e12;

Ti1_at_test=ppval(ppTi1,rho_for_fit);
ni1_at_test=ppval(ppni1,rho_for_fit)/1e12;

% Ti2_at_test=ppval(ppTi2,rho_for_fit);
% ni2_at_test=ppval(ppni2,rho_for_fit)/1e12;


figure;hold on;box on;
plot(roa_ts,Te,'ro')
plot(rho_for_fit,ppval(ppTe,rho_for_fit),'.-')
xlabel('\rho');ylabel('T_e (eV)')


rho_veryfine=linspace(0,1,1000);
figure;hold on;box on;
plot(roa_ts,ne,'ro')
plot(roa_ts,ni1,'bo')
% plot(roa_ts,ni2,'ko')
plot(rho_veryfine,ppval(ppne,rho_veryfine),'r.-')
plot(rho_veryfine,ppval(ppni1,rho_veryfine),'b.-')
% plot(rho_veryfine,ppval(ppni2,rho_veryfine),'k.-')
% plot(rho_for_fit,ni2_at_test*1e12,'kx-')
xlabel('\rho');ylabel('n_e /cc')

figure;hold on;box on;
plot(roa_ts,Ti1,'bo-')
% plot(roa_ts,Ti2,'ko-')
plot(rho_veryfine,ppval(ppTi1,rho_veryfine),'b.-')
% plot(rho_veryfine,ppval(ppTi2,rho_veryfine),'k.-')
xlabel('\rho');ylabel('T_i (eV)')

fprintf(fid,'%i\n',numpts);
for ind=1:numpts
%     fprintf(fid,'%15.7e\t%15.7e\t%15.7e\t%15.7e\t%15.7e\t%15.7e\t%15.7e\n',rho_for_fit(ind),ne_at_test(ind),Te_at_test(ind),Ti1_at_test(ind),ni1_at_test(ind),Ti2_at_test(ind),ni2_at_test(ind));
    fprintf(fid,'%15.7e\t%15.7e\t%15.7e\t%15.7e\t%15.7e\n',rho_for_fit(ind),ne_at_test(ind),Te_at_test(ind),Ti1_at_test(ind),ni1_at_test(ind));

end



fclose(fid);