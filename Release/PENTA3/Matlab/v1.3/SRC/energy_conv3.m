function [Coeff]=energy_conv3(jval,kval,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,log_interp,use_quad,K_exp,nu_exp,norm_factor,flags)
%This function performs the energy integral by
%interpolating the specified file(s).
%
% 1/6/2010 - 9/22/2011 JL

%integrate a given coefficient, multiplied, sometimes, by the collision
%frequency and a sonine polynomial

%Calculate log of x and y interpolation variables if necessary
if log_interp
    emesh(emesh==0)=1e-20;   %Catch -inf
    logcmesh=log(cmesh);
    logemesh=log(emesh);
else
    logcmesh=[];
    logemesh=[];    
end

%determine coefficients of Sonine polynomal of order kval
%coeffs_k is ordered as [K^0 K^1 K^2 ... K^kval]
for mtest=0:kval
    ckm=(-1)^mtest*gamma(kval+5/2)/( factorial(kval-mtest)*gamma(mtest+5/2)*factorial(mtest) );
    coeffs_k(mtest+1)=ckm;
end

%determine coefficients of Sonine polynomal of order jval
%coeffs_j is ordered as [K^0 K^1 K^2 ... K^jval]
for mtest=0:jval
    ckm=(-1)^mtest*gamma(jval+5/2)/( factorial(jval-mtest)*gamma(mtest+5/2)*factorial(mtest) );
    coeffs_j(mtest+1)=ckm;
end

%Define energy limits of integration
Kmin=flags.Kmin;
Kmax=flags.Kmax;
%number of K values for rectangular integration
num_K_vals=flags.num_K_vals;


% Kmin=1e-5;
% Kmax=20;
% num_K_vals=10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Perform energy integral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if use_quad
    %%%%% Using quad
    tol=1e-8;
    Coeff_tmp=quad(@intfun,Kmin,Kmax,tol,[],cmesh,emesh,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,jval,kval,vtb,coef2d_slice,logopt,log_interp,logcmesh,logemesh,coeffs_j,coeffs_k,K_exp,nu_exp);
else
    %%%% Using regtangular approx.
    
%     Kvals=linspace(Kmin,Kmax,num_K_vals);
%     dK=Kvals(2)-Kvals(1);
%     Coeff_tmp=sum(intfun(Kvals,cmesh,emesh,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,jval,kval,vtb,coef2d_slice,logopt,log_interp,logcmesh,logemesh,coeffs_j,coeffs_k,K_exp,nu_exp))*dK;
    
    Kvals=logspace(log10(Kmin),log10(Kmax),num_K_vals);
    test=intfun(Kvals,cmesh,emesh,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,jval,kval,vtb,coef2d_slice,logopt,log_interp,logcmesh,logemesh,coeffs_j,coeffs_k,K_exp,nu_exp);
    Coeff_tmp=sum(test(1:end-1).* (Kvals(2:end)-Kvals(1:end-1)) );
end

%Unnormalize the coefficient.  2/sqrt(pi) is from the convolution, the
%norm_factor must account for the monoenergetic coefficient.  
Coeff=2*Coeff_tmp*norm_factor/sqrt(pi);  


function [integrand]=intfun(Kval,cmesh,emesh,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,juse,kuse,vtb,coef2d,logopt,log_interp,logcmesh,logemesh,coeffs_j,coeffs_k,K_exp,nu_exp)
%This function defines the integrand of eq. 36 except for the numerical
%factor in front and includes the energy dependence of the coefficients.

% Test species normalized energy and velcocity and velocity
xa=sqrt(Kval);
va=xa*vta;

% |Er|/v_a
efield=abs_Er./va;

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
cmul_K = nu_perp_a./va;

%Interpolate 2D coefficient matrix to get D(cmul,efield)

extrapval=NaN;
%check for out of bounds
if any(efield>emesh(end,end))
%     disp('big efield')
    efield(efield>emesh(end,end))=emesh(end,end);
end
if any(cmul_K>cmesh(end,end))
%     disp('big cmul')
    cmul_K(cmul_K>cmesh(end,end))=cmesh(end,end);
end
if any(cmul_K<cmesh(1,1))
%     disp('small cmul')
    cmul_K(cmul_K<cmesh(1,1))=cmesh(1,1);
end

if log_interp
    mono_coeff=interp2(logcmesh,logemesh,coef2d.',log(cmul_K),log(efield),'linear',extrapval);
else
    mono_coeff=interp2(cmesh,emesh,coef2d.',cmul_K,efield,'linear',extrapval);
end

%Calculate the Sonine polynomial product
kfun_tmp1=0;
for mtest=0:kuse
    kfun_tmp1=kfun_tmp1+coeffs_k(mtest+1)*Kval.^mtest;
end

kfun_tmp2=0;
for mtest=0:juse
    kfun_tmp2=kfun_tmp2+coeffs_j(mtest+1)*Kval.^mtest;
end

kfun=kfun_tmp1.*kfun_tmp2;

%if the coefficient is log(D) then perform exp(mono_coeff), also take
%exp(-K) at the same time
if logopt  
    kfun2=exp(mono_coeff-Kval);
else
    kfun2=exp(-Kval).*mono_coeff;
end

%Take into account the energy (velocity) dependence of the monoenergetic
%coefficient.  D_11* --> v**3, D_13* --> v**2, D_33 -->v**1
% That makes the K dependence **(3/2, 1, 1/2), respectively.
% There is also the K**1/2 in the integral, giving **(2,3/2,1)

integrand=kfun.*kfun2.*Kval.^(K_exp+1/2).*nu_perp_a.^(nu_exp);
