function Xvec=form_Xvec(Er_test,Te,dTedr,ne,dnedr,Ti,dTidr,ni,dnidr,Z_ion,B_Eprl,Bsq)
%form thermodynamic forces from eqs. (10),(C6) SN
%   Species charge sign has been included, elem_charge has been
%   factored out of Er term because Ta is in eV.
% JL 2008-2010

elem_charge = 1.602176487e-19;

%electrons
Xe1 = -Te*elem_charge*(dnedr/ne + dTedr/Te + Er_test/Te);
Xe2 = -elem_charge*dTedr;

%ions
Xi1 = -elem_charge*Ti.*(dnidr./ni + dTidr./Ti - Z_ion.*Er_test./Ti);
Xi2 = -elem_charge*dTidr;

%E_||
XE=B_Eprl/sqrt(Bsq);

%Form vector
Xvec=[Xe1 Xe2];

for ispec=1:length(Ti)
    Xvec=[Xvec Xi1(ispec) Xi2(ispec)];
end

Xvec=[Xvec XE].';
    
    