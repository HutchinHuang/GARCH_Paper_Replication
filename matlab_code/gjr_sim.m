% simulation of GJR(1,1)
function eps = gjr_sim(n, omega,  alpha1, alpha2, beta,df)
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

eps(1) = sqrt(omega)*eta(1);
for t = 2:n 
    ht(t) = omega+alpha1*(max(0,eps(t-1))).^2+alpha2*(min(0,eps(t-1))).^2+beta*ht(t-1);
    eps(t) = sqrt(ht(t))*eta(t);
end
end
