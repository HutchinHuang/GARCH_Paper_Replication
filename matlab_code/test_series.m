function tab_stat = test_series(indices,alphainit,betainit)
nbind=length(indices);
tab_stat=zeros(nbind,1);

for ind = 1:nbind
    data = csvread(char(strcat('../csv/',indices(ind))), 1, 3); % 为了防止日期被误读
    close=flipud(data(:,2)); 
    ntot=length(close);
    rend=zeros(ntot,1);
    rend(2:ntot)=log(close(2:ntot)./close(1:(ntot-1)))*100; 

    ser=rend(2:ntot)-mean(rend(2:ntot));

    omegainit=abs(ser(1));

    [~, ~, stat, ~]=estimate_stationary_test(omegainit,alphainit,betainit,ser);

    tab_stat(ind) = round(stat,1);
end
end