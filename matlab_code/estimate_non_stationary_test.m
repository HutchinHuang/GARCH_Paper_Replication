% Estimation of GARCH(1,1) 
% Then test H0: gamma0 < 0 (Strictly Stationary)
function [para, minimum, stat, pval] = estimate_non_stationary_test(omega0,alpha0,beta0,eps0)
options = optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');

b0 = [omega0;alpha0;beta0];

lb = zeros(size(b0));
ub = Inf(size(b0));

[para, minimum] = fmincon(@(x)objf3(x, eps0),b0,[],[],[],[],lb,ub,[],options);

omega=para(1);
alpha=para(2);
beta=para(3);

n = length(eps0);

sigma2 = zeros(n,1);
sigma2(1) = omega;

for t = 2:n 
    sigma2(t)=omega+alpha*eps0(t-1).^2+beta*sigma2(t-1);
end
         
eta2=zeros(n,1);
eta2(2:n)=eps0(2:n).^2./sigma2(2:n);

sigma2u=var(log(alpha*eta2(2:n)+beta));
stat=sqrt(n)*mean(log(alpha*eta2(2:n)+beta))/sqrt(sigma2u);
pval=normcdf(stat);
end

