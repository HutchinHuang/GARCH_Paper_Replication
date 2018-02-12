% Main Programme to Get Table V
% Almost the same as that to get Table IV, except for the test function.
clear
omega0=0.1;
beta0=0.8;
df=7;
alphavec=[0.18;0.20;0.22;0.25752;0.28;0.30;0.31];
niter=length(alphavec);
njter=[500;2000;4000];
for jter = 1:length(njter)
    n=njter(jter);
    nrep=1000;
    pval_matrix = zeros(niter,nrep);
    para_omega=zeros(niter,nrep);
    para_alpha=zeros(niter,nrep);
    para_beta=zeros(niter,nrep);
    for i = 1:niter
        alpha0=alphavec(i);
        for j = 1:nrep
            ser=garch_sim_switch(n, omega0,  alpha0, beta0,df);
            [para, minimum, stat, pval] = estimate_non_stationary_test(abs((1+0.1/2-0.1*unifrnd(0,1))*omega0), abs((1+0.1/2-0.1*unifrnd(0,1))*alpha0),abs((1+0.1/2-0.1*unifrnd(0,1))*beta0), ser);
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