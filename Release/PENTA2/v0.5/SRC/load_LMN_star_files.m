function [cmul_mat,efield_mat,coef2d_mat]=load_LMN_star_files(data_path,run_ident,vmec_jval)
%This function loads the L*/e^2,M* and N* files named
%[l,m,n]star_lijs_[run_ident]_s[vmec_jval], located in the director
%data_path.  The output matrices are cells ordered with {1,2,3}
%corresponding to {l,m,n}.
% 6/2009 JL


file_char_vals = ['l','m','n'];

for i_char=1:3
    
    file_char=file_char_vals(i_char);
    
    if file_char == 'm'
        logopt = 1;
    else
        logopt = 0;
    end
    
    fname=[data_path '\' file_char 'star_lijs_' run_ident '_s' num2str(vmec_jval)];
    
    %read file containing monoenergetic coefficients indicated by fname
    fdata=dlmread(fname);
    num_cmul=fdata(2,1);
    num_efield=fdata(2,2);
    
    cmul_mat{i_char}=fdata(3:3+num_cmul-1,1);
    efield_mat{i_char}=fdata(3+num_cmul:3+num_cmul+num_efield-1,1);
    
    
    %loop over efield, cmul and create coefficient matrix
    count=0;
    for j = 1:num_efield
        istart=3+num_cmul+num_efield;
        if logopt
            coef2d(:,j)=log(fdata(istart+count:istart+count+num_cmul-1,1));
        else
            coef2d(:,j)=fdata(istart+count:istart+count+num_cmul-1,1);
        end
        count=count+num_cmul;
    end
    
    coef2d_mat{i_char}=coef2d;
    
end