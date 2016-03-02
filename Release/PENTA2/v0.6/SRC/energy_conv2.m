function [Coeff]=energy_conv2(jvals,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh)
%This function performs the energy integral defined in (36) by
%interpolating the specified lijs file(s).
% 6/2009 JL

%use logarithmic interpolation in x,y (cmul,efield)
log_interp=1;

%loop over jvals
for jval=jvals
        
    %Define energy limits of integration
    Kmin=1e-4;
    Kmax=10;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Perform energy integral of (36)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     %%%%% Using quad
%     tol=1e-8;
%     Coeff_tmp=quad(@intfun,Kmin,Kmax,tol,[],cmesh,emesh,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,jval,vtb,coef2d_slice,logopt,log_interp,logcmesh,logemesh);
%     
    %%%% Using regtangular approx.
    Kvals=linspace(Kmin,Kmax,1000);
    dK=Kvals(2)-Kvals(1);
    Coeff_tmp=sum(intfun(Kvals,cmesh,emesh,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,jval,vtb,coef2d_slice,logopt,log_interp,logcmesh,logemesh)).*dK;
    
    %Unnormalize the coefficient
    Coeff(jval)=Coeff_tmp*na*2*vta*ma/sqrt(pi);
    
end

function [integrand]=intfun(Kval,cmesh,emesh,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,jval,vtb,coef2d,logopt,log_interp,logcmesh,logemesh)
%This function defines the integrand of eq. 36 except for the numerical
%factor in front and includes the energy dependence of the coefficients.

% Test species normalized energy and velcocity and velocity
Ka=Kval;
xa=sqrt(Kval);
va=xa*vta;

% |Er|/v_a
efield=abs_Er./va;

%Define energy dependent collision frequency
nu_perp_aa= calc_perp_coll_freq(va,xa,qa2,qa2,ma,na,Ka,loglambda,vta,Ka);

%loop over field species to get collision frequency
for ispec=1:length(vtb)
    % Field species norm. energy and velocity
    xb=xa.*(vta/vtb(ispec));
    Kb=xb.^2;

    %Define energy dependent collision frequency
    nu_perp_ab(ispec,:)= calc_perp_coll_freq(va,xb,qa2,qb2(ispec),ma,nb(ispec),Kb,loglambda,vta,Ka);
end

nu_perp_a=nu_perp_aa + sum(nu_perp_ab,1);
cmul_K = nu_perp_a./va;

%Interpolate 2D coefficient matrix to get D(cmul,efield)
extrapval=0;
if log_interp
    mono_coeff=interp2(logcmesh,logemesh,coef2d.',log(cmul_K),log(efield),'linear',extrapval);
else
    mono_coeff=interp2(cmesh,emesh,coef2d.',cmul_K,efield,'linear',extrapval);
end

%integrate the coefficients using equation (36).  Note the K^2 instead of
%sqrt(K) is due to the the K^-3/2 dependence of the normalized
%coefficients.  See definition of M* before (54) for example.  The
%n*2/sqrt(pi) factor is included in the un-normalization of the cofficients
%in the main routine above for speed.
if logopt
    integrand=( Kval.^2.*exp(mono_coeff-Kval).*(Kval-5/2).^(jval-1) );
else
    integrand=( Kval.^2.*exp(-Kval).*(Kval-5/2).^(jval-1) ).*mono_coeff;
end
