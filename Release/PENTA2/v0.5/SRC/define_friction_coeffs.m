function lmat=define_friction_coeffs(ion_mass,Z_ion,vth_e,vth_i,ne,ni,Te,Ti,loglambda)

%constants
elem_charge = 1.602176487e-19;
eps0 = 8.854187817e-12;
e_mass = 9.10938215e-31;


masses=[e_mass ion_mass];
charges=elem_charge*[-1 Z_ion];
v_ths=[vth_e vth_i];
Ts=[Te,Ti];
ns=[ne,ni];

tau_coeff = 3*sqrt(pi)*pi*eps0^2;

num_species=length(Z_ion)+1;
for ispec1=1:num_species
    
        ma=masses(ispec1);
        Ta=Ts(ispec1);
        vth_a=v_ths(ispec1);
        charge_a=charges(ispec1);
        na=ns(ispec1);

        
        %only really need these evaluated at ispec2
        Kak=(v_ths./vth_a).^2;
        Mak00=-(1 + ma./masses).*(1 + Kak).^-1.5;
        Mak10=-1.5*(1+ma./masses).*(1+Kak).^-2.5;
        Mak11=-(13/4 + 4*Kak + 15/2*Kak.^2).*(1+Kak).^-2.5;
        Mak01=Mak10;
    
        Nak11=(27/4)*(Ta./Ts).*Kak.*(1+Kak).^-2.5;
        
        tau_ak=tau_coeff*ma^2*vth_a^3./(ns.*charge_a^2.*charges.^2*loglambda);
        
    for ispec2 = ispec1 : ispec1+(num_species - ispec1)
        

        Nab00=-Mak00(ispec2);
        Nab10=-Mak10(ispec2);
        Nab11=Nak11(ispec2);
        
        tau_ab=tau_ak(ispec2);
        
        Kba=(vth_a./v_ths(ispec2)).^2;
        Mba10=-1.5*(1+masses(ispec2)/ma)*(1+Kba)^-2.5;
        Nab01=-Mba10*(Ta*vth_a)/(Ts(ispec2)*v_ths(ispec2));
        
        if ispec1 == ispec2
            lab11=na*ma*( sum(Mak00./tau_ak) + Nab00./tau_ab );
            lab21=na*ma*( sum(Mak10./tau_ak) + Nab10./tau_ab );
            lab22=na*ma*( sum(Mak11./tau_ak) + Nab11./tau_ab );
            lab12=na*ma*( sum(Mak01./tau_ak) + Nab01./tau_ab );
        else
            lab11=na*ma*( Nab00./tau_ab );
            lab21=na*ma*( Nab10./tau_ab );
            lab22=na*ma*( Nab11./tau_ab );
            lab12=na*ma*( Nab01./tau_ab ); 
        end
        
        lmat{ispec1,ispec2}=[lab11 -lab12;-lab21 lab22];
    end
end


for ispec1=2:num_species
    for ispec2= 1:(ispec1 - 1)
        lmat{ispec1,ispec2}=lmat{ispec2,ispec1}.';
    end
end



% %electrons
% tau_coeff = 3*sqrt(pi)*pi*eps0^2;
% tau_ee =tau_coeff*e_mass^2*vth_e^3/(ne*(elem_charge^4)*loglambda);
% Mee00=-1/sqrt(2);
% Mee10=-3/8*sqrt(2);
% Mee11=-(13/4 + 4 + 15/2)*sqrt(2)/8;
% Nee00=-Mee00;
% Nee10=-Mee10;
% Nee11=(27/4)*sqrt(2)/8;
% 
% tau_ei =tau_coeff*e_mass^2*vth_e^3./(ni.*Z_ion.^2*(elem_charge^4)*loglambda);
% Kei=(vth_i./vth_e).^2;
% Mei00=-(1+e_mass./ion_mass).*(1+Kei).^-1.5;
% Mei10=-1.5*(1+e_mass./ion_mass).*(1+Kei).^-2.5;
% Mei11=-(13/4 + 4*Kei + 15/2*Kei.^2).*(1+Kei).^-2.5;
% 
% lee11=ne*e_mass*(Mee00/tau_ee+sum(Mei00./tau_ei)+Nee00/tau_ee);
% lee21=ne*e_mass*(Mee10/tau_ee+sum(Mei10./tau_ei)+Nee10/tau_ee);
% lee22=ne*e_mass*(Mee11/tau_ee+sum(Mei11./tau_ei)+Nee11/tau_ee);
% lee12=lee21;
