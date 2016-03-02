function [ nu_perp_ab_core ] = calc_perp_coll_freq(xb,qb2,nb)
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

% Updated 6/28/2011 to only calculate field species component

sqrtpi=1.772453850905516;
H_ab=(1-1./(2*xb.^2)).*erf(xb) + exp(-xb.^2)./(xb.*sqrtpi);
nu_perp_ab_core=nb*qb2*H_ab;

