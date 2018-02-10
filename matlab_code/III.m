% Main Programme to Get Table III.
clear
df=7;
omega0=0.1;
alpha0=0.5;
betavec=[0.61;0.64;0.67;0.70;0.73;0.76;0.79];
beta_star=0.7;
niter=length(betavec);
njter=[500;2000;4000];
for jter = 1:length(njter)
    n=njter(jter);
    nrep=1000;
    pval_matrix=zeros(niter,nrep);
    para_omega=zeros(niter,nrep);
    para_alpha=zeros(niter,nrep);
    para_beta=zeros(niter,nrep);
    for i = 1:niter
        beta0=betavec(i);
        for j = 1:nrep
            ser=garch_sim_switch(n, omega0,  alpha0, beta0,df);
            [para, minimum, stat, pval, ok] = estimgarch11_and_test(abs((1+0.1/2-0.1*unifrnd(0,1))*omega0),abs((1+0.1/2-0.1*unifrnd(0,1))*alpha0),abs((1+0.1/2-0.1*unifrnd(0,1))*beta0), ser,beta_star);
            if(ok~=0)
                disp([i,ok,pval]);
            end
            pval_matrix(i,j)=pval;
            para_omega(i,j)=para(1);
            para_alpha(i,j)=para(2);
            para_beta(i,j)=para(3);
        end
    end
    freq_rej=zeros(niter,1);
    for i = 1:niter 
        freq_rej(i)=100*length(find(pval_matrix(i,:)< 0.05))/nrep;
    end
    disp(freq_rej);
end
