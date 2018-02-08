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
