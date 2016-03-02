%this file loads the Sugama_data files created by
%convert_dkes_to_sugama.m and tabulates them in the file
%format required by PENTA

%-------------------------------------------------
%required variables
%path_penta_input - folder where the files will be created
%penta_input_suffix - ext on dkes and vmec files
%surf_vals - surfaces over which to run this
%DKES_data_type
%dkes_runs_path

%surfaces to run this for
% surf_vals=[2:21 40:10:200];
surf_vals=[10];

%starting path for output directory (initial ui directory)
penta_runs_path='D:\Transport\PENTA\penta_runs';


% penta_input_suffix='qhs1T';
% DKES_data_type='best'; 

penta_input_suffix='tok';
DKES_data_type='bcorrect'; 


%file where DKES_data file is located
% dkes_runs_path='D:\Transport\DKES\DKES_runs\qhs1T';
dkes_runs_path='D:\Transport\ideal_tokamak\DKES_data';

path_penta_input=uigetdir(penta_runs_path, 'Indicate the folder where the penta input files will be created');


%number of collisionalities
% num_cmul=35;   
num_cmul=42;   

%loop over surface
for is=1:length(surf_vals)
    
    %select test surface
    mysurf=surf_vals(is);
    disp(['Working on surface ' num2str(is) ' of ' num2str(length(surf_vals))])
    
    %open input file
    input_file_name=['DKES_data_' penta_input_suffix '_s' num2str(mysurf) '_' DKES_data_type '.mat'];
    load([dkes_runs_path '\' input_file_name]);    
    
    %open output files
    fid_D11star_file=fopen([path_penta_input '\D11_star_' penta_input_suffix '_s' num2str(mysurf)],'w');
    fid_D13star_file=fopen([path_penta_input '\D13_star_' penta_input_suffix '_s' num2str(mysurf)],'w');
    fid_D33star_file=fopen([path_penta_input '\D33_star_' penta_input_suffix '_s' num2str(mysurf)],'w');

    num_efield=length(efield_vals);

    %write num_cmul num_efield
    fprintf(fid_D11star_file,'      %i        %i\n',num_cmul,num_efield);
    fprintf(fid_D13star_file,'      %i        %i\n',num_cmul,num_efield);
    fprintf(fid_D33star_file,'      %i        %i\n',num_cmul,num_efield);
    
    %write cmul values in increasing order
    for ic=1:num_cmul
        fprintf(fid_D11star_file,' %15.7e\n',cmul_vals(ic));
        fprintf(fid_D13star_file,' %15.7e\n',cmul_vals(ic));
        fprintf(fid_D33star_file,' %15.7e\n',cmul_vals(ic));
    end
    
    %write efield values
    for ie=1:num_efield
        fprintf(fid_D11star_file,' %15.7e\n',efield_vals(ie));
        fprintf(fid_D13star_file,' %15.7e\n',efield_vals(ie));
        fprintf(fid_D33star_file,' %15.7e\n',efield_vals(ie));
    end
    
    %write the values to the file
    for ie=1:num_efield
        for ic=1:num_cmul
            fprintf(fid_D11star_file,' %15.7e\n',D11_mean(ic,ie));
            fprintf(fid_D13star_file,' %15.7e\n',D13_mean(ic,ie));
            fprintf(fid_D33star_file,' %15.7e\n',D33_mean(ic,ie));
        end
    end
    
    
end

%close file
fclose(fid_D11star_file);    
fclose(fid_D13star_file);
fclose(fid_D33star_file);
fclose('all')
    
    
