function para = estimgarch11 (omega0,alpha0,beta0,eps0)
options = optimset('Hessian','bfgs','Algorithm','interior-point','Display','off');

b0 = [omega0;alpha0;beta0];

lb = zeros(size(b0));
ub = Inf(size(b0));

para = fmincon(@(x)objf3(x, eps0),b0,[],[],[],[],lb,ub,[],options);

end
