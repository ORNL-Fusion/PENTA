function [VMEC_info,ion_params,pprof_info_allsurf,U2_data]=...
    load_input_files(flags,path_info,arad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function loads the VMEC info file, the ion parameters file and the
%plasma profile input file and compiles the data into structures.
%Optionally (through the flag 'load_U2_file') a U2 file is read.  
%
% Note: If arad=[] it is read from VMEC data
% 1/26/2010 J:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Constants
p_mass=1.672621637e-27;             %proton mass

%expand path_info structure
data_path=path_info.data_path;
pprof_ident=path_info.pprof_ident;
run_ident=path_info.run_ident;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ion parameters from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ion_file=[data_path '/ion_params'];
fid_ion=fopen(ion_file);
ion_data=textscan(fid_ion,'%s');
id1=char(ion_data{1}(2));
id2=char(ion_data{1}(3));
id3=char(ion_data{1}(4));
num_ion_species=str2num(id1(find(id1=='=')+1:end));
Z_ion=str2num(id2(find(id2=='=')+1:end));
miomp=str2num(id3(find(id3=='=')+1:end));
fclose(fid_ion);

ion_mass=miomp*p_mass;

%check if ion parameters defined
if length(ion_mass)~= num_ion_species || length(Z_ion)~= num_ion_species
    error('Variables ion_mass and Z_ion must have length==num_ion_species')
end

%Compile ion parameters
ion_params.ion_mass=ion_mass;
ion_params.Z_ion=Z_ion;
ion_params.num_ion_species=num_ion_species;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read plasma profile file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pprof_info_allsurf]=read_pprof_file(data_path,pprof_ident,num_ion_species,flags.plot_pprof);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read VMEC profile info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data1_vmec=dlmread([data_path '/profile_data_' run_ident],'',[0,0,1,1]);
data2_vmec=dlmread([data_path '/profile_data_' run_ident],'',3,0);
%use predefined <a> if given above
if isempty(arad)
    arad=data1_vmec(2,1);    %<a>
end
VMEC_info.arad=arad;
VMEC_info.Rmajor_vmec=data1_vmec(2,2); %<R0>
VMEC_info.js_vmec=data2_vmec(:,1);     %VMEC surface #
VMEC_info.roa_vmec=data2_vmec(:,3);    %sqrt(norm. tor. flux)
VMEC_info.chip_vmec=data2_vmec(:,4);   %deriv of poloidal flux
VMEC_info.psip_vmec=data2_vmec(:,5);   %deriv of toroidal flux
VMEC_info.Btheta_vmec=data2_vmec(:,6); %Boozer pol field
VMEC_info.Bzeta_vmec=data2_vmec(:,7);  %Boozer tor field
VMEC_info.vp_vmec=data2_vmec(:,8);     %volume derivative
VMEC_info.Bsq_vmec=data2_vmec(:,9);    %<B**2>
VMEC_info.iota_vmec=data2_vmec(:,10);  %iota
VMEC_info.B0_vmec=data2_vmec(:,11);    %B(m=0,n=0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read U^2 file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flags.load_U2_file
    U2_data=dlmread([data_path '/Utilde2_profile.dat'],'',1,0);
else
    U2_data=[];
end