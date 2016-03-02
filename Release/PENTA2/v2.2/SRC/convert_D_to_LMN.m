function [cmul_mat,efield_mat,coef2d_mat_out,U2]=convert_D_to_LMN(cmul_mat,efield_mat,coef2d_mat,U2,Bsq)
%note that cmul=nu/v, efield=|Er|/v

%assumes bcorrect has been done and same efield cmul for each coefficient
%1/08 JL

%Set flags

%L* manipulation flags
USE_OLD_MERGE_FACTOR=0;  %This merges the L* coefficient to the Er=0 result above cmul~=0.1
USE_NEW_MERGE=0;         %Begins merging when L* becomes negative
EXTEND_IT=0;             %uses last non negative L* for all higher cmul values
no_PS_flow=0;            %Turns off PS flow term in L*
ZERO_IT=1;               %sets all negative L* to 1e-12

%These flags should be set for legacy reasons
take_log_mstar=1;        %return log(M*)
include_efactor=1;       %return L*/e^2 and N*/e, where e=elem_charge

%plot coefficients
plotit=0;                

%constants
elem_charge = 1.602176487e-19;

%check flags
if (USE_NEW_MERGE && USE_OLD_MERGE_FACTOR) || (EXTEND_IT && USE_OLD_MERGE_FACTOR) || (USE_NEW_MERGE && EXTEND_IT)
    error('Incompatible flag settings in convert_D_to_LMN')
end

%check that different sized arrays are not used for each coefficient
size1=size(coef2d_mat{1});
size2=size(coef2d_mat{2});
size3=size(coef2d_mat{3});

if any([size1~=size2, size1~=size3, size2 ~= size3])
    error('All D coefficients must have identical cmul, efield arrays')
end

%assign cmul, efield values
num_efield=size(efield_mat{1},1);
cmul_vals=cmul_mat{1};

%assign DKES coefficients
D11_mean=coef2d_mat{1};
D13_mean=coef2d_mat{2};
D33_mean=coef2d_mat{3};

%approximate <U~^2>  -- this is passed as an input argument
if no_PS_flow
    U2=0;
else
    U2=1.5*D11_mean(end,1)/cmul_vals(end);
end

%for new merging ignore highest cmul value
if USE_NEW_MERGE
    cmul_vals=cmul_mat{1}(1:end-1);
    D11_mean=coef2d_mat{1}(1:end-1,:);
    D13_mean=coef2d_mat{2}(1:end-1,:);
    D33_mean=coef2d_mat{3}(1:end-1,:);
end

for ie=1:num_efield  %efield loop

    %define denominator used in all three coefficients
    denom=1-1.5*cmul_vals.*D33_mean(:,ie)./Bsq;
    
    %Handle Merge factor settings              
    if USE_OLD_MERGE_FACTOR==1   %perform merging unless Er=0
        merge_hi_nu=true;
        if ie==1
            merge_hi_nu=false;
        end
    else                        %otherwise perform no merging at high cmul
        merge_hi_nu=false;
    end
    
    %set merge_factor
    merge_factor=1;  %default value for merge_factor (1 has no effect on coefficient)
    if merge_hi_nu
        merge_factor  = 1-( tanh( 4*(log(cmul_vals)+2.0) ) +1 )./2; %goes to zero for cmul > ~.1
    end    
    
    %define L*
    lstar(:,ie) = D11_mean(:,ie) - (2/3)*cmul_vals*U2  + 1.5*cmul_vals.*D13_mean(:,ie).^2./(Bsq*denom);
     
    %save Er=0 L* for merging
    if ie == 1
        lstar_E0 = lstar(:,ie);
    elseif USE_NEW_MERGE
        %find cmul where L* becomes negative
        [merge_index]=find(lstar(:,ie)<0);
        %begin merging when L* becomes negative
        if ~isempty(merge_index)
            if merge_index(1)==1
                lstar(:,ie)=lstar(:,ie-1);
                merge_hi_nu=0;
            else              
                cmul_merge=cmul_vals(merge_index(1)-1);
                merge_factor  = 1-( tanh( 4*(log(cmul_vals)-log(cmul_merge)) ) +1 )./2; %4 is like a steepness parameter, other terms are for normalization
                merge_hi_nu=1;
            end
        end
    elseif EXTEND_IT
        %find cmul where L* becomes negative
        [merge_index]=find(lstar(:,ie)<0);
        %use last positive L*
        if ~isempty(merge_index)
            if merge_index(1)==1 || merge_index(1)==2
                lstar(:,ie)=1e-12;
                
            else                
                lstar(merge_index(1):end,ie)=lstar(merge_index(1)-1,ie);
            end
            merge_hi_nu=0;
        end
    end
    
    if merge_hi_nu
        lstar(:,ie)=lstar(:,ie).*merge_factor + lstar_E0.*( 1 - merge_factor );
    end
           
    %define M*
    mstar(:,ie)=cmul_vals.^2.*D33_mean(:,ie)./denom;
    
    %define N*
    nstar(:,ie)=cmul_vals.*D13_mean(:,ie)./denom;

    
end %efield loop
    
%Handle additional flags
if ZERO_IT
    lstar(lstar<1e-10)=1e-10;
end

%plotting
if plotit
    figure;hold on;box on;
    plot(cmul_vals,lstar,'.-')
    xlabel('\nu/v')
    ylabel('L*')
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    title('L*')
    
    figure;hold on;box on;
    plot(cmul_vals,mstar,'.-')
    xlabel('\nu/v')
    ylabel('M*')
    set(gca,'yscale','log')
    set(gca,'xscale','log')
    title('M*')
    
    figure;hold on;box on;
    plot(cmul_vals,nstar,'.-')
    xlabel('\nu/v')
    ylabel('N*')
    set(gca,'xscale','log')
    title('N*')
    
end

%handle even more flags


%Invoke lower limits on M*
if any(any(mstar <= 0))
    disp('setting lower limit on M*')
    mstar(mstar <= 0)=1e-12;
end

if take_log_mstar==1
    mstar=log(mstar);
end

if include_efactor
    lstar=lstar./elem_charge^2;
    nstar=nstar./elem_charge;
end

%set outputs
coef2d_mat_out{1}=lstar;
coef2d_mat_out{2}=mstar;
coef2d_mat_out{3}=nstar;

cmul_mat{1}=cmul_vals;
cmul_mat{2}=cmul_vals;
cmul_mat{3}=cmul_vals;