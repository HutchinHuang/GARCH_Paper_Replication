% Estimate J and kappa, with regard to GARCH(1,1)
function [kappa_eta, Jmat, Imat] = MatJ(omega,alpha,beta,eps)
n = length(eps);
sigma2 = zeros(n,1);
sigma2(1) = omega;
dersigma2 = zeros(3,n);
dersigma2(:,1) = [1;0;0];
for t = 2:n
    sigma2(t)=omega+alpha*eps(t-1).^2+beta*sigma2(t-1);
    dersigma2(1:3,t)=[1;eps(t-1)^2;sigma2(t-1)];
    dersigma2(1:3,t)=dersigma2(1:3,t)+beta*dersigma2(1:3,t-1);
end

Jmat=zeros(3,3);

for t = 2:n
Jmat=Jmat+(dersigma2(1:3,t)/sigma2(t))*(dersigma2(1:3,t)/sigma2(t))';
end

Jmat=Jmat/n;

eta= eps ./ sqrt(sigma2);

kappa_eta=mean(eta(2:n).^4);
Imat=Jmat(2:3,2:3)-Jmat(2:3,1)*(Jmat(2:3,1))'/Jmat(1,1);
end