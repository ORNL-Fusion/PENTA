function [Coeff]=energy_conv2(jvals,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,vtb,logopt,cmesh,emesh,coef2d_slice,logcmesh,logemesh,log_interp,use_quad)
%This function performs the energy integral defined in (36) by
%interpolating the specified lijs file(s).
%
% Here jvals give the combination of Laguerre polynomials:
%   1: 00
%   2: 01, (10)
%   3: 11
%   4: 02, (20)
%   5: 12, (21)
%   6: 22
%   7: 03, (30)
%   8: 13, (31)
%   9: 23, (32)
%   10: 33
%   11: 04, 40
%   12: 14, 41
%   13: 24, 42
%   14: 34, 43
%   15: 44
%
% 6/2009 JL

%loop over jvals
for jval=jvals
        
    %Define energy limits of integration
    Kmin=1e-4;
    Kmax=10;
    
    %number of K values for rectangular integration
    num_K_vals=1000;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Perform energy integral of (36)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if use_quad
        %%%%% Using quad
        tol=1e-8;
        Coeff_tmp=quad(@intfun,Kmin,Kmax,tol,[],cmesh,emesh,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,jval,vtb,coef2d_slice,logopt,log_interp,logcmesh,logemesh);
    else
        %%%% Using regtangular approx.
        Kvals=linspace(Kmin,Kmax,num_K_vals);
        dK=Kvals(2)-Kvals(1);
        Coeff_tmp=sum(intfun(Kvals,cmesh,emesh,vta,qa2,qb2,ma,na,nb,loglambda,abs_Er,jval,vtb,coef2d_slice,logopt,log_interp,logcmesh,logemesh))*dK;
    end
    
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


%   1: 00
%   2: 01, (10)
%   3: 11
%   4: 02, (20)
%   5: 12, (21)
%   6: 22
%   7: 03, (30)
%   8: 13, (31)
%   9: 23, (32)
%   10: 33
%   11: 04, 40
%   12: 14, 41
%   13: 24, 42
%   14: 34, 43
%   15: 44
if jval==1
    kfun=1;
elseif jval==2
    kfun=2.5-Kval;
elseif jval==3
    kfun=6.25-5*Kval+Kval.^2;
elseif jval==4
    kfun=.5*Kval.^2-3.5*Kval+4.375;
elseif jval==5
    kfun=(.5*Kval.^2-3.5*Kval+4.375).*(2.5-Kval);
elseif jval==6
    kfun=(.5*Kval.^2-3.5*Kval+4.375).^2;
elseif jval==7
    kfun=315/48-(63/8)*Kval+(9/4)*Kval.^2-(1/6)*Kval.^3;
elseif jval==8
    kfun=(315/48-(63/8)*Kval+(9/4)*Kval.^2-(1/6)*Kval.^3).*(2.5-Kval);
elseif jval==9
    kfun=(315/48-(63/8)*Kval+(9/4)*Kval.^2-(1/6)*Kval.^3).*(.5*Kval.^2-3.5*Kval+4.375);
elseif jval==10
    kfun=(315/48-(63/8)*Kval+(9/4)*Kval.^2-(1/6)*Kval.^3).^2;
elseif jval==11
    kfun=(3465/384-(693/48)*Kval+(99/16)*Kval.^2-(11/12)*Kval.^3+(1/24)*Kval.^4);
elseif jval==12
    kfun=(3465/384-(693/48)*Kval+(99/16)*Kval.^2-(11/12)*Kval.^3+(1/24)*Kval.^4).*(2.5-Kval);
elseif jval==13
    kfun=(3465/384-(693/48)*Kval+(99/16)*Kval.^2-(11/12)*Kval.^3+(1/24)*Kval.^4).*(.5*Kval.^2-3.5*Kval+4.375);
elseif jval==14
    kfun=(3465/384-(693/48)*Kval+(99/16)*Kval.^2-(11/12)*Kval.^3+(1/24)*Kval.^4).*(315/48-(63/8)*Kval+(9/4)*Kval.^2-(1/6)*Kval.^3);
elseif jval==15
    kfun=(3465/384-(693/48)*Kval+(99/16)*Kval.^2-(11/12)*Kval.^3+(1/24)*Kval.^4).^2;
else
    error(['bad jval' num2str(jval)])
end
    
if logopt
    kfun2=exp(mono_coeff-Kval);
else
    kfun2=exp(-Kval).*mono_coeff;
end

integrand=Kval.^2.*kfun2.*kfun;
