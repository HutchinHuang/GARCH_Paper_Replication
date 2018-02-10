function [tab_n, tab_alpha, tab_beta, tab_stat, tab_p] = test_the_last(indices,alphainit,betainit)
nbind=length(indices);

tab_n=zeros(nbind,1)';
tab_alpha=zeros(nbind,1)';
tab_beta=zeros(nbind,1)';
tab_stat=zeros(nbind,1)';
tab_p=zeros(nbind,1)';

for ind = 1:nbind
    data = csvread(char(strcat('../csv/',indices(ind))), 1, 3); % 为了防止日期被误读
    close=flipud(data(:,2)); 
    ntot=length(close);
    rend=zeros(ntot,1);
    rend(2:ntot)=log(close(2:ntot)./close(1:(ntot-1)))*100; 

    ser=rend(2:ntot)-mean(rend(2:ntot));
    if ind == 2
        plot(ser);
        title('FIGURE.1: Log returns (in %) of the MCBF stock series');
    end
    
    omegainit=abs(ser(1));

    [para, ~, stat, pval]=estimate_stationary_test(omegainit,alphainit,betainit,ser);
 
    tab_n(ind) = ntot;
    tab_alpha(ind) = para(2);
    tab_beta(ind) = para(3);
    tab_stat(ind) = round(stat,1);
    tab_p(ind) = pval;
end
end