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