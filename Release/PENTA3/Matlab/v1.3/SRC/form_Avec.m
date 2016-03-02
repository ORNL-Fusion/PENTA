function Avec=form_Avec(Er_test,Ts,dTdrs,ns,dndrs,charges,B_Eprl,Bsq,B0_vmec)
% JL 9/21/2011

Avec = zeros(1,length(Ts)*3);

elem_charge = 1.602176487e-19;

for is=1:length(Ts)
    ind_A=(is-1)*3+1;    
    Avec(ind_A)   = dndrs(is)./ns(is) - 1.5*dTdrs(is)/Ts(is) - charges(is)*Er_test/(elem_charge*Ts(is));
    Avec(ind_A+1) = dTdrs(is)/Ts(is);
    Avec(ind_A+2) =  charges(is)*B0_vmec*B_Eprl/(Ts(is)*elem_charge*Bsq);
end
