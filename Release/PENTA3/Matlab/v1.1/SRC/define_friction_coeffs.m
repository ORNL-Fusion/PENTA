function lmat=define_friction_coeffs(masses,charges,v_ths,ns,Ts,loglambda,Smax,use_Ji)
%
% This function defines a matrix of classical friction coefficients.  This
% matrix is returned as a cell array of the form lmat{a,b}(i,j) 
% for species a and b and subscripts i,j. Each submatrix lab=lmat{a,b} 
% is of size Smax+1,Smax+1 where Smax is the number of Sonine polynomials 
% used in expanding the l=1 component of fa1.  
%
% The classical friction coefficients are defined in, for example, Eq. (25)
% of [Sugama and Nishimura, PoP 15, 042502 (2008)].  Low orders are given 
% in several references (Hirshman and Sigmar, Helander and Sigmar), but
% general expressions can be derived from [Ji and Held, PoP 13, 102103
% (2006)], see my talks for expressions.  
%
% The logical variable "use_Ji" uses the expressions derived from Ji and
% Held if use_Ji==1.  
%
% JL 2009-2010

%constants
eps0 = 8.854187817e-12;

%coefficient in front of collision time
tau_coeff = 3*sqrt(pi)*pi*eps0^2;

%total number of species
num_species=length(masses);

%%%%Define upper triangular portion of lmat
%loop over test species (species "a")
for ispec1=1:num_species
    
    %parameters for test species
    ma=masses(ispec1);
    Ta=Ts(ispec1);
    vth_a=v_ths(ispec1);
    charge_a=charges(ispec1);
    na=ns(ispec1);
    nama=na*ma;
    
    %loop over field species (species "b", includes b=a)
%     for ispec2 = ispec1 : num_species
    for ispec2 = 1 : num_species        
        
        %loop over row and column of submatrix lab
        for irow=0:Smax
            for icol=0:Smax
                
                %Define matrix of friction coefficients for specified Sonine order
                lab(irow+1,icol+1)=0;
                
                %diagonal terms (Otherwise the M sum is zero)
                if ispec1 == ispec2
                    M_ovr_tau_sum=0;
                    %And now we have to sum over every species, indicated "k"
                    for spec_k=1:num_species
                        vth_k=v_ths(spec_k);
                        Tk=Ts(spec_k);
                        mk=masses(spec_k);
                        nk=ns(spec_k);
                        charge_k=charges(spec_k);
                        
                        %Calculate the collision time for species a on k for this row, col
                        tau_ak=tau_coeff*ma^2*vth_a^3/(nk*charge_a^2*charge_k^2*loglambda);
                        
                        %Define dimensionless ratios
                        chi_ak=vth_k/vth_a;
                        theta_ak=Tk/Ta;
                        mu_ak=mk/ma;
                        Q_ak=1+chi_ak^2;
                        
                        %calc Mak for this row, col
                        if use_Ji
                            Mak=2*calc_A(irow,icol,theta_ak,mu_ak,chi_ak,Q_ak);
                        else
                            Mak=calc_M_helander(irow,icol,theta_ak,mu_ak,chi_ak,Q_ak);
                        end
                        
                        %Add this Mak to the sum
                        M_ovr_tau_sum=M_ovr_tau_sum+Mak/tau_ak;
                    end
                    
                    %For diagonal elements we have a M term and a N term (added below)
                    lab(irow+1,icol+1)=nama*M_ovr_tau_sum;
                    
                end
                
                %define parameters for species b
                vth_b=v_ths(ispec2);
                Tb=Ts(ispec2);
                mb=masses(ispec2);
                charge_b=charges(ispec2);
                nb=ns(ispec2);
                
                %Define dimensionless ratios
                chi_ab=vth_b/vth_a;
                theta_ab=Tb/Ta;
                mu_ab=mb/ma;
                Q_ab=1+chi_ab^2;
                
                %All elements of lab have the N term
                tau_ab=tau_coeff*ma^2*vth_a^3/(nb*charge_a^2*charge_b^2*loglambda);
                if use_Ji
                    Nab=(2/chi_ab)*calc_B(irow,icol,theta_ab,mu_ab,chi_ab,Q_ab);
                else
                    Nab=calc_N_helander(irow,icol,theta_ab,mu_ab,chi_ab,Q_ab);
                end
                lab(irow+1,icol+1)=lab(irow+1,icol+1)+nama*( Nab/tau_ab );
            end %col
        end %row
        
        %assign the submatrix to the cell array
        lmat{ispec1,ispec2}=lab;
        
    end %spec2 loop
end %spec1 loop

% %define lower triangular portion of lmat (lab = lba.')
% for ispec1=2:num_species
%     for ispec2= 1:(ispec1 - 1)
%         lmat{ispec1,ispec2}=lmat{ispec2,ispec1}.';
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                   SUBROUTINES                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The subroutines below define the friction coefficients M and N using
% either expressions derived from Helander and Sigmar or from Ji and Held.
%
% JL 2009-2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Mpk=calc_M_helander(irow,icol,theta,mu,chi,Q)
%Direct expressions for Mpk_ij derived from Helander and Sigmar
% JL 2009

xab=chi;
yab=1/mu;

if irow == 0 && icol == 0
    Mpk=-(1 + yab)*Q^-1.5;
elseif (irow == 0 && icol == 1) || (irow == 1 && icol == 0)
    Mpk=-1.5*(1 + yab)*Q^-2.5;
elseif irow == 1 && icol == 1
    Mpk=-(13/4 + 4*xab^2 + 15/2*xab^4)*Q^-2.5;
elseif (irow == 0 && icol == 2) || (irow == 2 && icol == 0)
    Mpk=-(15/8)*(1 + yab)*Q^-3.5;
elseif (irow == 1 && icol == 2) || (irow == 2 && icol == 1)
    Mpk=-( 69/16 + 6*xab^2 + 63/4*xab^4 )*Q^-3.5;
elseif irow == 2 && icol == 2
    Mpk=-( 433/64 + 17*xab^2 + 459/8*xab^4 + 28*xab^6 + 175/8*xab^8 )*Q.^-4.5;
else
    error('No non-Ji method defined for this Sonine order')
end

function Npk=calc_N_helander(irow,icol,theta,mu,chi,Q)
%Direct expressions for Npk_ij derived from Helander and Sigmar
% Reqires calc_M_helander
% JL 2009

xab=chi;
TaoTb=1/theta;

chi_ba=1/chi;
mu_ba=1/mu;
theta_ba=1/theta;
Q_ba=1+chi_ba^2;

xba=chi_ba;
TboTa=theta;

if irow == 0 && icol == 0
    Npk=-calc_M_helander(0,0,theta,mu,chi,Q);
elseif (irow == 1 && icol == 0)
    Npk=-calc_M_helander(0,1,theta,mu,chi,Q);
elseif (irow == 0 && icol == 1)
    Nba10=-calc_M_helander(0,1,theta_ba,mu_ba,chi_ba,Q_ba);
    Npk=Nba10*(TaoTb)*chi_ba;
elseif irow == 1 && icol == 1
    Npk=(27/4)*TaoTb*xab^2*Q^-2.5;
elseif (irow == 0 && icol == 2)
    Npk=-calc_M_helander(0,2,theta_ba,mu_ba,chi_ba,Q_ba)/xab;
elseif (irow == 2 && icol == 0)
    Npk=-calc_M_helander(2,0,theta,mu,chi,Q);
elseif (irow == 1 && icol == 2)
    Npk=(225/16)*TaoTb*xab^4*Q^-3.5;
elseif (irow == 2 && icol == 1)
    Nba12=(225/16)*TboTa*xba^4*Q_ba^-3.5;
    Npk=Nba12*(TaoTb)*chi_ba;
elseif irow == 2 && icol == 2
    Npk=(2625/64)*TaoTb*xab^4*Q^-4.5;
else
    error('No non-Ji method defined for this Sonine order')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The functions below are used when calculating the coefficients for the Ji
% and Held method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eab=calc_eab(k,Q)
%Simple function to calculate eab
% JL 2009
Xab=Q^-(k+1/2);
eab=gamma(k+1/2)*Xab/gamma(0.5);

function Eba=calc_Eba(k,chi,Q)
%Simple function to calculate Eab
% JL 2009
coeff=factorial(k)/2/chi;
Etmp=0;
for jtest=0:k
    Etmp=Etmp+chi^(1+2*jtest)*calc_eab(jtest,Q)/factorial(jtest);
end
Eba=coeff*Etmp;


function Apk=calc_A(p,k,theta,mu,chi,Q)
%Function to calculate Apk
%Note that this likely breaks for "large" values of p or k because
%factorials and gamma functions of these values (and derived quantities)
%are taken.
% JL 2009

%calculate A for "arbitrary" p,k
Apk=0;
mvals=0:k;
qvals=0:p;
for qtest=qvals
    for mtest=mvals
        
        %calculate c coefficients
        cpq=(-1)^qtest*gamma(p+5/2)/( factorial(p-qtest)*gamma(qtest+5/2)*factorial(qtest) );
        ckm=(-1)^mtest*gamma(k+5/2)/( factorial(k-mtest)*gamma(mtest+5/2)*factorial(mtest) );
        
        %calculate Astar
        Astar_qm=calc_Astar(qtest,mtest,theta,mu,chi,Q);
        
        %calculate portion of  Apk
        Atmp=cpq*ckm*Astar_qm;
        
        %total Apk
        Apk=Apk+Atmp;
    end
end


function Bpk=calc_B(p,k,theta,mu,chi,Q)
%Function to calculate Bpk
%Note that this likely breaks for "large" values of p or k because
%factorials and gamma functions of these values (and derived quantities)
%are taken.
% JL 2009

%calculate B for "arbitrary" p,k
Bpk=0;
mvals=0:k;
qvals=0:p;
for qtest=qvals
    for mtest=mvals
        
        %calculate c coefficients
        cpq=(-1)^qtest*gamma(p+5/2)/( factorial(p-qtest)*gamma(qtest+5/2)*factorial(qtest) );
        ckm=(-1)^mtest*gamma(k+5/2)/( factorial(k-mtest)*gamma(mtest+5/2)*factorial(mtest) );
        
        %calculate Bstar
        Bstar_qm=calc_Bstar(qtest,mtest,theta,mu,chi,Q);
        
        %calculate portion of Bpk
        Btmp=cpq*ckm*Bstar_qm;
        
        %total Bpk
        Bpk=Bpk+Btmp;
    end
end


function Astar=calc_Astar(q,m,theta,mu,chi,Q)
%Calculates variables for Astar
%Note this actually gives Astar * (tau/3/n)
%JL 2009

%lower case a's
a_E0=-2*(1-theta)/mu;
a_E1=-1+(1+2*m)*(1-2*theta)/mu;
a_E2=2*theta*m^2/mu;
a_e0=2*(1-theta)*(1/theta+1/mu);
a_e1=(1+2*m)*(1-1/mu+2*theta/mu);
a_e2=-2*theta*m^2/mu;

%alphas - actually alpha * (tau/3/n)

alpha_E0=calc_Eba(1+q+m,chi,Q);
alpha_E1=calc_Eba(q+m,chi,Q);
if q==0 && m==0     %alpha_E not defined for n=2, l=1, q=m=0 (aE2^00 = 0)
    alpha_E2=0;
else
    alpha_E2=calc_Eba(q+m-1,chi,Q);
end

% chi^(1+2*k)*calc_eab(k,Q) gives eba

alpha_e0=chi^(1+2*(2+q+m))*calc_eab(2+q+m,Q)/chi;
alpha_e1=chi^(1+2*(1+q+m))*calc_eab(1+q+m,Q)/chi;
alpha_e2=chi^(1+2*(q+m))*calc_eab(q+m,Q)/chi;

%actually Astar * (tau/3/n)
Astar=a_E0*alpha_E0 + a_E1*alpha_E1 + a_E2*alpha_E2 + a_e0*alpha_e0 + a_e1*alpha_e1 + a_e2*alpha_e2;


function Bstar=calc_Bstar(q,m,theta,mu,chi,Q)
%Calculates variables for Bstar
%Note this actually gives Bstar * (tau/3/n)
%JL 2009

%lower case b's
b_e0=2*mu/theta^2;
b_e1=-4/5;
b_p0=4+8/5*m+4/3*mu/theta-8/3/theta;
b_m0=-8/3*mu/theta+4/3/theta;
b_m1=8/5;

%betas  - actually beta * (tau/3/n)
beta_e0=chi^(2*q+5)*calc_eab(2+q+m,Q);
beta_e1=chi^(2*q+5)*calc_eab(3+q+m,Q);

back_p0=0;
for j=0:q
    back_p0=back_p0 + factorial(q)/factorial(j)*chi^(2*j)*calc_eab(2+m+j,Q);
end

back_m0=0;
back_m1=0;
for j=0:m
    back_m0=back_m0 + 1/factorial(j)*calc_eab(2+q+j,Q);
    back_m1=back_m1 + 1/factorial(j)*calc_eab(3+q+j,Q);
end

beta_p0=0.5*chi^3 * back_p0;
beta_m0=0.5*chi^(2*q+5)*factorial(m)*back_m0;
beta_m1=0.5*chi^(2*q+5)*factorial(m)*back_m1;

%actually Bstar * (tau/3/n)
Bstar=b_e0*beta_e0 + b_e1*beta_e1 + b_p0*beta_p0 + b_m0*beta_m0 + b_m1*beta_m1;
