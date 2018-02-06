set.seed(666)
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

for i = 1:niter 
    ser=garch.sim(n1,omega0,alpha0,beta0);
    resestimarch1 = estimgarch11(abs(omega00+0.1*runif(1)),
        abs(alpha00+0.1*runif(1)),abs(beta00+0.1*runif(1)),ser);
    write(resestimarch1$coef,file="estGARCH2dst.dat",ncolumns=3,append=TRUE);
end

for i = 1:niter
    ser=garch.sim(n1,omega0,alpha1,beta0)
    resestimarch1 = estimgarch11(abs(omega00+0.1*runif(1)),
    abs(alpha01+0.1*runif(1)),abs(beta00+0.1*runif(1)),ser)
    write(resestimarch1$coef,file="estGARCHstst.dat",ncolumns=3,append=TRUE)
end

for i = 1:niter
ser=garch.sim(n1,omega0,alpha2,beta0)
resestimarch1 = estimgarch11(abs(omega00+0.1*runif(1)),
abs(alpha02+0.1*runif(1)),abs(beta00+0.1*runif(1)),ser)
write(resestimarch1$coef,file="estGARCHexp.dat",ncolumns=3,append=TRUE)
end

for i = 1:niter
ser=garch.sim(n2,omega0,alpha0,beta0)
resestimarch1 = estimgarch11(abs(omega00+0.1*runif(1)),
abs(alpha00+0.1*runif(1)),abs(beta00+0.1*runif(1)),ser)
write(resestimarch1$coef,file="estGARCH2dst2.dat",ncolumns=3,append=TRUE)
end

for i = 1:niter
ser=garch.sim(n2,omega0,alpha1,beta0)
resestimarch1 = estimgarch11(abs(omega00+0.1*runif(1)),
abs(alpha01+0.1*runif(1)),abs(beta00+0.1*runif(1)),ser)
write(resestimarch1$coef,file="estGARCHstst2.dat",ncolumns=3,append=TRUE)
end

for i = 1:niter
ser=garch.sim(n2,omega0,alpha2,beta0)
resestimarch1 = estimgarch11(abs(omega00+0.1*runif(1)),
abs(alpha02+0.1*runif(1)),abs(beta00+0.1*runif(1)),ser)
write(resestimarch1$coef,file="estGARCHexp2.dat",ncolumns=3,append=TRUE)
end

#### variance des estimations
toutest = read.table("estGARCH2dst.dat")
colnames(toutest)=c("omega","alpha","beta")
var(toutest)
eromega0=omega0-toutest$omega 
eralpha0=alpha0-toutest$alpha 
erbeta0=beta0-toutest$beta 

toutest = read.table("estGARCHstst.dat")
colnames(toutest)=c("omega","alpha","beta")
var(toutest)
eromega1=omega0-toutest$omega 
eralpha1=alpha1-toutest$alpha 
erbeta1=beta0-toutest$beta 

toutest = read.table("estGARCHexp.dat")
colnames(toutest)=c("omega","alpha","beta")
var(toutest)
eromega2=omega0-toutest$omega 
eralpha2=alpha2-toutest$alpha 
erbeta2=beta0-toutest$beta 

round(c(mean(eromega0),mean(eralpha0),mean(eralpha0),
mean(eromega1),mean(eralpha1),mean(eralpha1),
mean(eromega2),mean(eralpha2),mean(eralpha2)),digits=2)
round(c(mean(eromega0^2),mean(eralpha0^2),mean(eralpha0^2),
mean(eromega1^2),mean(eralpha1^2),mean(eralpha1^2),
mean(eromega2^2),mean(eralpha2^2),mean(eralpha2^2)),digits=2)

toutest = read.table("estGARCH2dst2.dat")
colnames(toutest)=c("omega","alpha","beta")
var(toutest)
eromega02=omega0-toutest$omega 
eralpha02=alpha0-toutest$alpha 
erbeta02=beta0-toutest$beta 

toutest = read.table("estGARCHstst2.dat")
colnames(toutest)=c("omega","alpha","beta")
var(toutest)
eromega12=omega0-toutest$omega 
eralpha12=alpha1-toutest$alpha 
erbeta12=beta0-toutest$beta 

toutest = read.table("estGARCHexp2.dat")
colnames(toutest)=c("omega","alpha","beta")
var(toutest)
mean(toutest)
eromega22=omega0-toutest$omega 
eralpha22=alpha2-toutest$alpha 
erbeta22=beta0-toutest$beta 

round(c(mean(eromega02),mean(eralpha02),mean(eralpha02),
mean(eromega12),mean(eralpha12),mean(eralpha12),
mean(eromega22),mean(eralpha22),mean(eralpha22)),digits=2)
round(c(mean(eromega02^2),mean(eralpha02^2),mean(eralpha02^2),
mean(eromega12^2),mean(eralpha12^2),mean(eralpha12^2),
mean(eromega22^2),mean(eralpha22^2),mean(eralpha22^2)),digits=2)