function [Lmat,Mmat,Nmat]=define_thermal_mats(Er,ne,vth_e,Z_ion,ion_mass,ni,vth_i,loglambda,cmul_mat,efield_mat,coef2d_mat)

%constants
elem_charge = 1.602176487e-19;
e_mass = 9.10938215e-31;

%define |Er|
abs_Er=abs(Er);

%ELECTRONS
%   Integrate electron La,Ma,Na coefficients using (Eq. 36) to get 
%    Laj,Maj,Naj thermal coefficients:
%   Here species "a" are electrons and "b" are ions, species' mass,
%	 charge, temperature and thermal velocity are defined below and
%	 used in subroutine energy_conv.

elem2=elem_charge^2;

ma=e_mass;
qa2=elem2;
na=ne;
vta=vth_e;

qb2=Z_ion.^2*elem2;
nb=ni;
vtb=vth_i;



%integrate Le
logopt=0; mat_ind=1;
[cmesh,emesh]=meshgrid(cmul_mat{mat_ind},efield_mat{mat_ind});
logemesh=log(emesh);logcmesh=log(cmesh);
coef2d_slice=coef2d_mat{mat_ind};
Le_vals=energy_conv2(1:3,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh);
Le1=Le_vals(1); Le2=Le_vals(2);  Le3=Le_vals(3);
%integrate Me
logopt=1; mat_ind=2;
[cmesh,emesh]=meshgrid(cmul_mat{mat_ind},efield_mat{mat_ind});
logemesh=log(emesh);logcmesh=log(cmesh);
coef2d_slice=coef2d_mat{mat_ind};
Me_vals=energy_conv2(1:3,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh);
Me1=Me_vals(1); Me2=Me_vals(2);  Me3=Me_vals(3);
%integrate Ne
logopt=0;mat_ind=3;
[cmesh,emesh]=meshgrid(cmul_mat{mat_ind},efield_mat{mat_ind});
logemesh=log(emesh);logcmesh=log(cmesh);
coef2d_slice=coef2d_mat{mat_ind};
Ne_vals=energy_conv2(1:3,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh);
Ne1=Ne_vals(1); Ne2=Ne_vals(2);  Ne3=Ne_vals(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now integrate ion coefficients using (36)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ispec=1:length(Z_ion)
    %this ion species parameters
    ma=ion_mass(ispec);
    qa2=Z_ion(ispec).^2*elem2;
    na=ni(ispec);
    vta=vth_i(ispec);
    
    %electron and all other ion species params
    qb2=elem2*[-1 Z_ion(1:end~=ispec)].^2;
    nb=[ne ni(1:end~=ispec)];
    vtb=[vth_e vth_i(1:end~=ispec)];

    
    %integrate Li
    logopt=0;mat_ind=1;
    [cmesh,emesh]=meshgrid(cmul_mat{mat_ind},efield_mat{mat_ind});
    logemesh=log(emesh);logcmesh=log(cmesh);
    coef2d_slice=coef2d_mat{mat_ind};
    Li_vals=energy_conv2(1:3,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh);
    Li1(ispec)=Li_vals(1); Li2(ispec)=Li_vals(2);  Li3(ispec)=Li_vals(3);
    %integrate Me
    logopt=1;mat_ind=2;
    [cmesh,emesh]=meshgrid(cmul_mat{mat_ind},efield_mat{mat_ind});
    logemesh=log(emesh);logcmesh=log(cmesh);
    coef2d_slice=coef2d_mat{mat_ind};
    Mi_vals=energy_conv2(1:3,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh);
    Mi1(ispec)=Mi_vals(1); Mi2(ispec)=Mi_vals(2);  Mi3(ispec)=Mi_vals(3);
    %integrate Ne
    logopt=0;mat_ind=3;
    [cmesh,emesh]=meshgrid(cmul_mat{mat_ind},efield_mat{mat_ind});
    logemesh=log(emesh);logcmesh=log(cmesh);
    coef2d_slice=coef2d_mat{mat_ind};
    Ni_vals=energy_conv2(1:3,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh);
    Ni1(ispec)=Ni_vals(1); Ni2(ispec)=Ni_vals(2);  Ni3(ispec)=Ni_vals(3);
    
end


%assemble L,M,N mats

%electrons
Lmat{1}=[Le1 Le2;Le2 Le3];
Mmat{1}=[Me1 Me2;Me2 Me3];
Nmat{1}=[Ne1 Ne2;Ne2 Ne3];

%ions
for ispec=1:length(Z_ion)
    
    Lmat{ispec+1}=[Li1(ispec) Li2(ispec);Li2(ispec) Li3(ispec)];
    Mmat{ispec+1}=[Mi1(ispec) Mi2(ispec);Mi2(ispec) Mi3(ispec)];
    Nmat{ispec+1}=[Ni1(ispec) Ni2(ispec);Ni2(ispec) Ni3(ispec)];

    
end
    