% simulation d'un GARCH(1,1) avec bruit iid Student normalise
function eps = garch_sim_switch(n, omega, alpha, beta, df)
if df>2
    factor = sqrt((df-2)/df);
else
    factor = 1;
end

if df>0
    eta = factor*trnd(df, n, 1);
else
    eta = randn(n, 1);
end

ht=zeros(n,1);
eps = ht;
ht(1) = omega;
eps(1) = sqrt(ht(1))*eta(1);
for t = 2:n 
    ht(t) = omega+alpha*eps(t-1).^2+beta*ht(t-1);
    eps(t) = sqrt(ht(t))*eta(t);
end
end
