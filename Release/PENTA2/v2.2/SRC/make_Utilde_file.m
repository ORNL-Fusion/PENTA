%1/08 JL
clear all;
inital_dir=pwd;

%surfaces to run this for
% surf_vals=[2:25 30:10:50 70 90 120:20:200];
% surf_vals=[2:25 30:10:200];  %qhs1T
surf_vals=20;


roa_surf=sqrt((surf_vals - 1)/200);

%set some values that will be used later
Bscale_correction=false;  
itok=false;
% DKES_data_type='best';  %num_efield=110;
DKES_data_type='bcorrect'; 
plotit=0;

total_cmul_vals=42;

%file name stuff
% ext='qhs1T';
ext='tok';
% ext='qhs1T_ideal';
% ext='f1410p1T';

dkes_runs_path='D:\Transport\DKES\DKES_runs';
dkes_results_dir=[dkes_runs_path '\' ext];
dkes_results_dir='D:\Transport\ideal_tokamak\DKES_data'


%constants
e_ovr_c=1.602e-19;

%save path
save_path='D:\Transport\PENTA\tests\7_2009_PENTA_Matlab_impurity_tests\qhs1T_parabola';
save_path='D:\Transport\ideal_tokamak\PENTA_data'

%display some stuff
disp(['Working on files with extension: ' ext])


fid_U=fopen([save_path '\Utilde2_profile'],'w');
fprintf(fid_U,'%i\n',length(surf_vals))

    
    disp(['Loading dkes from location: ' dkes_results_dir])

%loop over each surface
for is=1:length(surf_vals)
    
    mysurf=surf_vals(is);
%     disp(['Working on surface ' num2str(is) ' of ' num2str(length(surf_vals))])

    %load dkes file for this surface
    %load either extrapolated dkes file, or CALC_ALL_CMUL file

    try
        cd(dkes_results_dir)
        DKES_filename=['DKES_data_' ext '_s' num2str(mysurf) '_' DKES_data_type];
        load(DKES_filename)
    catch
        disp(['tried to load filename: ' DKES_filename])
        error('DKES folder or file does not exist')
    end
        
    %approximate Utwiddle^2 for this surface ~1.5*D11 as (cmul -> inf)
%     disp('verify u~^2 approximation!')
    utilde_sq(is)=1.5*D11_mean(end,1)/cmul_vals(end); 

%     figure;hold on;box on;
    plot(cmul_vals,D11_mean(:,5))
%     plot(cmul_vals,D11_mean(:,1)./cmul_vals.')
    set(gca,'xscale','log','yscale','log')
%     title(num2str(roa_surf(is))) 
   
    
  
    
    fprintf(fid_U,'%8.5f %8.5f\n', roa_surf(is),utilde_sq(is));
    
    
    
   

    
    
end %surf loop



fclose(fid_U);

figure;hold on; box on;
plot(roa_surf,utilde_sq,'r.-')
xlabel('r/a')
ylabel('<U^2>')


cd(inital_dir);