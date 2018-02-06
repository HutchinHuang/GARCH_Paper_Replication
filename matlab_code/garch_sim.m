% simulation d'un GARCH(1,1)
function eps = garch_sim(n, omega, alpha, beta)
eta = randn(n, 1);
ht=zeros(n,1);
eps=ht;
eps(1) = sqrt(omega)*eta(1);
for t = 2:n 
    ht(t) = omega+alpha*eps(t-1).^2+beta*ht(t-1);
    eps(t) = sqrt(ht(t))*eta(t);
end
end