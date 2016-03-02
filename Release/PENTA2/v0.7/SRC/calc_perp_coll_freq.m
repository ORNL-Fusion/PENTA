function [ nu_perp_a ] = calc_perp_coll_freq(va,xb,qa2,qb2,ma,nb,Kb,loglambda,vta,Ka)
%Input: Tb in eV, all other units mks
%
%calculates perpendicular collision frequency of particles of velocity
%va of type a on thermal species b (note that if there is more than one
%species, this result should be added to the collision frequency of the
%test particle on that species)
%
%Accepts vector input for va

%From eq 2.86 Callen's Book or Eq 25 Sugama and Nishimura PoP V9 p4637
%(2002)

%constants
eps0_sqrd=7.839664189871123e-023;
sqrtpi=1.772453850905516;


va3=Ka.*va*vta^2;
%calculate reference collision freq of a on b 
nu_ab=nb*qa2*qb2*loglambda./(ma^2*va3*4*pi*eps0_sqrd);


%calculate H(x_ab) (related to rosenbluth potential)
%using error function and analytic deriv. of erf
H_ab=(1-1./(2*Kb)).*erf(xb) + exp(-Kb)./(xb.*sqrtpi);

nu_perp_a=H_ab.*nu_ab;
