function [cmul_array,efield_array,D11_star_slice,D31_star_slice,D13_star_slice,D33_star_slice]=...
    load_DKES_star_files(data_path,run_ident,vmec_jval,Add_Spitzer_to_D33)
%This function loads the files containing the D* coefficients.
% 9/2009-2/2010 JL

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

%Try including "Spitzer" portion of D33
if Add_Spitzer_to_D33
    %loop over efield and cmul
    ne_loop=size(efield_mat{3},1);
    nc_loop=size(cmul_mat{3},1);
    
    for ie_ind=1:ne_loop
        for ic_ind=1:nc_loop
            
            %calculate L33(Spitzer)=(2/3)/CMUL
            D33_Spitzer=(2/3)/cmul_mat{3}(ic_ind);
            
            %calculate L33(Physical)=L33(Spitzer)-L33
            coef2d_mat{3}(ic_ind,ie_ind)=D33_Spitzer-coef2d_mat{3}(ic_ind,ie_ind);
        end
    end
end

%make sure all cmul and efield are the same
if any(cmul_mat{1}~=cmul_mat{2}) || any(cmul_mat{1}~=cmul_mat{3}) 
    error('All cmul values must be the same for each surface in the D* files')
end

if any(efield_mat{1}~=efield_mat{2}) || any(efield_mat{1}~=efield_mat{3}) 
    error('All efield values must be the same for each surface in the D* files')
end

%Assign output variables
cmul_array=cmul_mat{1};
efield_array=efield_mat{1};

%Define all four monoenergetic coefficients
D11_star_slice=coef2d_mat{1};
D13_star_slice=coef2d_mat{2};
D31_star_slice=-D13_star_slice;  %Use Onsager symmetry
D33_star_slice=coef2d_mat{3};
