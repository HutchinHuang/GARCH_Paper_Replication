% Main Programme to Get Table 1.
% ===== Part1:Basic Constants =====
clear;
format short;
omega0=1;
alpha0=0.3;
alpha1=0.5;
alpha2=0.7;
omega00=3;
beta0=0.6;
alpha00=alpha0;
alpha01=alpha1;
alpha02=alpha2;
beta00=beta0;
n1=200;
n2=4000;
niter=1000;

% ===== Part2:Garch Models Estimations =====
% 1.1 n1 = 200, garch_model = 2dst
estGARCH2dst = zeros(niter,3);
for i = 1:niter 
    ser=garch_sim(n1,omega0,alpha0,beta0);
    resestimarch1 = estimgarch11(abs(omega00+0.1*unifrnd(0,1)),abs(alpha00+0.1*unifrnd(0,1)),abs(beta00+0.1*unifrnd(0,1)),ser);
    estGARCH2dst(i,1) = resestimarch1(1);
    estGARCH2dst(i,2) = resestimarch1(2);
    estGARCH2dst(i,3) = resestimarch1(3);
end
eromega0=omega0-estGARCH2dst(:,1);
eralpha0=alpha0-estGARCH2dst(:,2);
erbeta0=beta0-estGARCH2dst(:,3);

% 1.2 n1 = 200, garch_model = stst
estGARCHstst = zeros(niter,3);
for i = 1:niter
    ser=garch_sim(n1,omega0,alpha1,beta0);
    resestimarch1 = estimgarch11(abs(omega00+0.1*unifrnd(0,1)), abs(alpha01+0.1*unifrnd(0,1)),abs(beta00+0.1*unifrnd(0,1)),ser);
    estGARCHstst(i,1) = resestimarch1(1);
    estGARCHstst(i,2) = resestimarch1(2);
    estGARCHstst(i,3) = resestimarch1(3);
end
eromega1=omega0-estGARCHstst(:,1); 
eralpha1=alpha1-estGARCHstst(:,2); 
erbeta1=beta0-estGARCHstst(:,3);

% 1.3 n1 = 200, garch_model = exp
estGARCHexp = zeros(niter,3);
for i = 1:niter
    ser=garch_sim(n1,omega0,alpha2,beta0);
    resestimarch1 = estimgarch11(abs(omega00+0.1*unifrnd(0,1)), abs(alpha02+0.1*unifrnd(0,1)),abs(beta00+0.1*unifrnd(0,1)),ser);
    estGARCHexp(i,1) = resestimarch1(1);
    estGARCHexp(i,2) = resestimarch1(2);
    estGARCHexp(i,3) = resestimarch1(3);
end
eromega2=omega0-estGARCHexp(:,1); 
eralpha2=alpha2-estGARCHexp(:,2);
erbeta2=beta0-estGARCHexp(:,3); 

% 2.1 n2 = 4000, garch_model = 2dst
estGARCH2dst2 = zeros(niter,3);
for i = 1:niter
    ser=garch_sim(n2,omega0,alpha0,beta0);
    resestimarch1 = estimgarch11(abs(omega00+0.1*unifrnd(0,1)), abs(alpha00+0.1*unifrnd(0,1)),abs(beta00+0.1*unifrnd(0,1)),ser);
    estGARCH2dst2(i,1) = resestimarch1(1);
    estGARCH2dst2(i,2) = resestimarch1(2);
    estGARCH2dst2(i,3) = resestimarch1(3);
end
eromega02=omega0-estGARCH2dst2(:,1); 
eralpha02=alpha0-estGARCH2dst2(:,2); 
erbeta02=beta0-estGARCH2dst2(:,3);

% 2.2 n2 = 4000, garch_model = stst
estGARCHstst2 = zeros(niter,3);
for i = 1:niter
    ser=garch_sim(n2,omega0,alpha1,beta0);
    resestimarch1 = estimgarch11(abs(omega00+0.1*unifrnd(0,1)),abs(alpha01+0.1*unifrnd(0,1)),abs(beta00+0.1*unifrnd(0,1)),ser);
    estGARCHstst2(i,1) = resestimarch1(1);
    estGARCHstst2(i,2) = resestimarch1(2);
    estGARCHstst2(i,3) = resestimarch1(3);
end
eromega12=omega0-estGARCHstst2(:,1);
eralpha12=alpha1-estGARCHstst2(:,2);
erbeta12=beta0-estGARCHstst2(:,3);

% 2.3 n2 = 4000, garch_model = exp
estGARCHexp2 = zeros(niter,3);
for i = 1:niter
    ser=garch_sim(n2,omega0,alpha2,beta0);
    resestimarch1 = estimgarch11(abs(omega00+0.1*unifrnd(0,1)), abs(alpha02+0.1*unifrnd(0,1)),abs(beta00+0.1*unifrnd(0,1)),ser);
    estGARCHexp2(i,1) = resestimarch1(1);
    estGARCHexp2(i,2) = resestimarch1(2);
    estGARCHexp2(i,3) = resestimarch1(3);
end
eromega22=omega0-estGARCHexp2(:,1); 
eralpha22=alpha2-estGARCHexp2(:,2);
erbeta22=beta0-estGARCHexp2(:,3);

% ===== Part3: Variance Estimations - Derive Tables =====
% Table I - Part1 (n1=200)
Table1_n1 = [mean(eromega0),mean(eralpha0),mean(eralpha0),mean(eromega1),mean(eralpha1),mean(eralpha1),mean(eromega2),mean(eralpha2),mean(eralpha2);...
    mean(eromega0.^2),mean(eralpha0.^2),mean(eralpha0.^2),mean(eromega1.^2),mean(eralpha1.^2),mean(eralpha1.^2),mean(eromega2.^2),mean(eralpha2.^2),mean(eralpha2.^2)];
disp('Table I - Part1 (n1=200):');
disp(round(Table1_n1,2));

% Table I - Part2 (n2=4000)
Table1_n2 = [mean(eromega02),mean(eralpha02),mean(eralpha02),mean(eromega12),mean(eralpha12),mean(eralpha12),mean(eromega22),mean(eralpha22),mean(eralpha22);...
    mean(eromega02.^2),mean(eralpha02.^2),mean(eralpha02.^2),mean(eromega12.^2),mean(eralpha12.^2),mean(eralpha12.^2),mean(eromega22.^2),mean(eralpha22.^2),mean(eralpha22.^2)];
disp('Table I - Part2 (n2=4000):');
disp(round(Table1_n2,2));








