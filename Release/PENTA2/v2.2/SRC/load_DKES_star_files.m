function [cmul_mat,efield_mat,coef2d_mat]=load_DKES_star_files(data_path,run_ident,vmec_jval)
% 9/2009 JL

%list of coefficient characters
file_char_vals = {'11','13','33'};

%loop over coeffcient and make matrix of values
for i_char=1:length(file_char_vals)
    
    %current character
    file_char=char(file_char_vals{i_char});
       
    %open file containing monoenergetic coefficients indicated by fname
    fname=[data_path '/D' file_char '_star_' run_ident '_s' num2str(vmec_jval)];
    try
        fdata=dlmread(fname);
    catch
        error(['Could not find file ' fname])
    end
    
    %read cmul, efield data
    num_cmul=fdata(1,1);
    num_efield=fdata(1,2);
    
    cmul_mat{i_char}=fdata(2:2+num_cmul-1,1);
    efield_mat{i_char}=fdata(2+num_cmul:2+num_cmul+num_efield-1,1);
    
    %loop over efield, cmul and create coefficient matrix
    count=0;
    for j = 1:num_efield
        istart=2+num_cmul+num_efield;
        coef2d(:,j)=fdata(istart+count:istart+count+num_cmul-1,1);
        count=count+num_cmul;
    end
    
    coef2d_mat{i_char}=coef2d;
    
end