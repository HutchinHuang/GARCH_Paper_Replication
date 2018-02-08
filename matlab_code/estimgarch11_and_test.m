% Estimation of GARCH(1,1) 
% Then test H0: alpha<=alpha.star 
function [para, minimum, stat, pval, ok] = estimgarch11_and_test(omega0,alpha0,beta0,eps0,beta_star)

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
[kappa_eta, ~, Imat] = MatJ(omega,alpha,beta,eps0);
I = Imat;
ok = 0;
if cond(I) < 6.710886474424970e+07
    Iinv = inv(I); 
else
    ok = 1; 
    Iinv = diag(ones(size(I,1),1));
end
kap = kappa_eta; 
den = sqrt((kap-1)*Iinv(2,2));
stat= sqrt(n) * (beta-beta_star)/den;
pval=1-normcdf(stat);
end