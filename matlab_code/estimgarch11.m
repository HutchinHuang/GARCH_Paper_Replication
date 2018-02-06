function para = estimgarch11 (omega0,alpha0,beta0,eps0)
options = optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');

b0 = [omega0;alpha0;beta0];

lb = zeros(size(b0));
ub = Inf(size(b0));

para = fmincon(@(x)objf3(x, eps0),b0,[],[],[],[],lb,ub,[],options);
% para = para';

% res = nlminb(valinit,objf3, lower=c(0.0000001,0,0),upper=c(Inf,Inf,Inf), eps=eps0)
% omega=res$par[1]
% alpha=res$par[2]
% beta=res$par[3]
% list(coef=c(omega,alpha,beta),minimum=res$objective)
end

% estime un GARCH(1,1) avec nlminb omega0=1 alpha0=0.05 beta=0.7 eps0=sim
function qmle = objf3(x, eps)    
    omega = x(1);
    alpha = x(2);
    beta = x(3);
    n = length(eps);
    sigma2 = zeros(n,1);
    sigma2(1) = omega;
    for t = 2:n
    sigma2(t) = omega+alpha*eps(t-1).^2+beta*sigma2(t-1);
    end
    qmle = mean(log(sigma2(2:n))+eps(2:n).^2 ./ sigma2(2:n));
end