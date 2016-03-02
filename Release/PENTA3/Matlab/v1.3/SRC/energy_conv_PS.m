function [Coeff]=energy_conv_PS(vta,qa2,qb2,ma,na,nb,loglambda,vtb,K_exp,norm_factor,flags)
% 9/21/2011 JL


%Define energy limits of integration
Kmin=flags.Kmin;
Kmax=flags.Kmax;
%number of K values for rectangular integration
num_K_vals=flags.num_K_vals;

% %Define energy limits of integration
% Kmin=1e-5;
% Kmax=20;
% 
% %number of K values for rectangular integration
% num_K_vals=10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Perform energy integral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Using regtangular approx.
Kvals=logspace(log10(Kmin),log10(Kmax),num_K_vals);
test=intfun(Kvals,vta,qa2,qb2,ma,na,nb,loglambda,vtb,K_exp);
Coeff_tmp=sum(test(1:end-1).* (Kvals(2:end)-Kvals(1:end-1)) );

%Unnormalize the coefficient.  2/sqrt(pi) is from the convolution, the
%norm_factor must account for the monoenergetic coefficient.  (and density)
Coeff=2*Coeff_tmp*norm_factor/sqrt(pi);  


function [integrand]=intfun(Kval,vta,qa2,qb2,ma,na,nb,loglambda,vtb,K_exp)

% Test species normalized energy and velcocity and velocity
xa=sqrt(Kval);
va=xa*vta;

% Coefficient in front of nu_ab
eps0_sqrd=7.839664189871123e-023;
nu_a_coeff = qa2*loglambda./(4*pi*eps0_sqrd*ma^2*va.^3);

%Define energy dependent collision frequency
nu_perp_aa_core= calc_perp_coll_freq(xa,qa2,na);

%loop over field species to get collision frequency
for ispec=1:length(vtb)
    % Field species norm. energy and velocity
    xb=xa.*(vta/vtb(ispec));

    %Define energy dependent collision frequency
    nu_perp_ab_core(ispec,:)= calc_perp_coll_freq(xb,qb2(ispec),nb(ispec));
end

nu_perp_a=nu_a_coeff.*(nu_perp_aa_core + sum(nu_perp_ab_core,1));

integrand=exp(-Kval).*Kval.^(K_exp+1/2).*nu_perp_a;
