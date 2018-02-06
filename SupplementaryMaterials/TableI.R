# simulation d'un GARCH(1,1)
garch.sim<- function (n, omega,  alpha, beta) {
eta <- rnorm(n)
ht<-rep(0,n)
eps<-ht
eps[1]<-sqrt(omega)*eta[1]
for (t in 2:n) { 
ht[t] <- omega+alpha*eps[t-1]**2+beta*ht[t-1]
eps[t]<-sqrt(ht[t])*eta[t]}
return(eps[1:n])}

# estime un GARCH(1,1) avec nlminb omega0<-1 alpha0<-0.05 beta<-0.7 eps0<-sim

objf3 <- function(x, eps){     
         omega <- x[1]
         alpha <- x[2]
         beta <- x[3]
         n <- length(eps)
         sigma2<-as.numeric(n)
         sigma2[1]<-omega
for (t in 2:n) sigma2[t]<-omega+alpha*eps[t-1]**2+beta*sigma2[t-1]
         qml <- mean(log(sigma2[2:n])+eps[2:n]**2/sigma2[2:n])
         qml }

estimgarch11<- function(omega0,alpha0,beta0,eps0)
     {
valinit<-c(omega0,alpha0,beta0)
res <- nlminb(valinit,objf3, lower=c(0.0000001,0,0),upper=c(Inf,Inf,Inf), eps=eps0)
omega<-res$par[1]
alpha<-res$par[2]
beta<-res$par[3]
list(coef=c(omega,alpha,beta),minimum=res$objective)
}

##############
file.remove("estGARCH2dst.dat","estGARCHstst.dat","estGARCHexp.dat")
file.remove("estGARCH2dst2.dat","estGARCHstst2.dat","estGARCHexp2.dat")
set.seed(666)
omega0<-1
alpha0<-0.3
alpha1<-0.5
alpha2<-0.7
omega00<-3
beta0<-0.6
alpha00<-alpha0
alpha01<-alpha1
alpha02<-alpha2
beta00<-beta0
n1<-200
n2<-4000
niter<-1000

for (i in 1:niter) { 
ser<-garch.sim(n1,omega0,alpha0,beta0)
resestimarch1 <- estimgarch11(abs(omega00+0.1*runif(1)),
    abs(alpha00+0.1*runif(1)),abs(beta00+0.1*runif(1)),ser)
write(resestimarch1$coef,file="estGARCH2dst.dat",ncolumns=3,append=TRUE)
     }
for (i in 1:niter) { 
ser<-garch.sim(n1,omega0,alpha1,beta0)
resestimarch1 <- estimgarch11(abs(omega00+0.1*runif(1)),
abs(alpha01+0.1*runif(1)),abs(beta00+0.1*runif(1)),ser)
write(resestimarch1$coef,file="estGARCHstst.dat",ncolumns=3,append=TRUE)
     }
for (i in 1:niter) { 
ser<-garch.sim(n1,omega0,alpha2,beta0)
resestimarch1 <- estimgarch11(abs(omega00+0.1*runif(1)),
abs(alpha02+0.1*runif(1)),abs(beta00+0.1*runif(1)),ser)
write(resestimarch1$coef,file="estGARCHexp.dat",ncolumns=3,append=TRUE)
     }
for (i in 1:niter) { 
ser<-garch.sim(n2,omega0,alpha0,beta0)
resestimarch1 <- estimgarch11(abs(omega00+0.1*runif(1)),
abs(alpha00+0.1*runif(1)),abs(beta00+0.1*runif(1)),ser)
write(resestimarch1$coef,file="estGARCH2dst2.dat",ncolumns=3,append=TRUE)
     }
for (i in 1:niter) { 
ser<-garch.sim(n2,omega0,alpha1,beta0)
resestimarch1 <- estimgarch11(abs(omega00+0.1*runif(1)),
abs(alpha01+0.1*runif(1)),abs(beta00+0.1*runif(1)),ser)
write(resestimarch1$coef,file="estGARCHstst2.dat",ncolumns=3,append=TRUE)
     }
for (i in 1:niter) { 
ser<-garch.sim(n2,omega0,alpha2,beta0)
resestimarch1 <- estimgarch11(abs(omega00+0.1*runif(1)),
abs(alpha02+0.1*runif(1)),abs(beta00+0.1*runif(1)),ser)
write(resestimarch1$coef,file="estGARCHexp2.dat",ncolumns=3,append=TRUE)
     }
#### variance des estimations
toutest <- read.table("estGARCH2dst.dat")
colnames(toutest)<-c("omega","alpha","beta")
var(toutest)
eromega0<-omega0-toutest$omega 
eralpha0<-alpha0-toutest$alpha 
erbeta0<-beta0-toutest$beta 

toutest <- read.table("estGARCHstst.dat")
colnames(toutest)<-c("omega","alpha","beta")
var(toutest)
eromega1<-omega0-toutest$omega 
eralpha1<-alpha1-toutest$alpha 
erbeta1<-beta0-toutest$beta 

toutest <- read.table("estGARCHexp.dat")
colnames(toutest)<-c("omega","alpha","beta")
var(toutest)
eromega2<-omega0-toutest$omega 
eralpha2<-alpha2-toutest$alpha 
erbeta2<-beta0-toutest$beta 

round(c(mean(eromega0),mean(eralpha0),mean(eralpha0),
mean(eromega1),mean(eralpha1),mean(eralpha1),
mean(eromega2),mean(eralpha2),mean(eralpha2)),digits=2)
round(c(mean(eromega0^2),mean(eralpha0^2),mean(eralpha0^2),
mean(eromega1^2),mean(eralpha1^2),mean(eralpha1^2),
mean(eromega2^2),mean(eralpha2^2),mean(eralpha2^2)),digits=2)

toutest <- read.table("estGARCH2dst2.dat")
colnames(toutest)<-c("omega","alpha","beta")
var(toutest)
eromega02<-omega0-toutest$omega 
eralpha02<-alpha0-toutest$alpha 
erbeta02<-beta0-toutest$beta 

toutest <- read.table("estGARCHstst2.dat")
colnames(toutest)<-c("omega","alpha","beta")
var(toutest)
eromega12<-omega0-toutest$omega 
eralpha12<-alpha1-toutest$alpha 
erbeta12<-beta0-toutest$beta 

toutest <- read.table("estGARCHexp2.dat")
colnames(toutest)<-c("omega","alpha","beta")
var(toutest)
mean(toutest)
eromega22<-omega0-toutest$omega 
eralpha22<-alpha2-toutest$alpha 
erbeta22<-beta0-toutest$beta 

round(c(mean(eromega02),mean(eralpha02),mean(eralpha02),
mean(eromega12),mean(eralpha12),mean(eralpha12),
mean(eromega22),mean(eralpha22),mean(eralpha22)),digits=2)
round(c(mean(eromega02^2),mean(eralpha02^2),mean(eralpha02^2),
mean(eromega12^2),mean(eralpha12^2),mean(eralpha12^2),
mean(eromega22^2),mean(eralpha22^2),mean(eralpha22^2)),digits=2)

##### graphiques
op <- par(mfrow = c(3, 2), cex=0.8) # 3 x 2 pictures on one plot
min<-min(eromega0,eralpha0,erbeta0,eromega1,eralpha1,erbeta1,eralpha2,erbeta2)
max<-max(eromega0,eralpha0,erbeta0,eromega1,eralpha1,erbeta1,eralpha2,erbeta2)
min2<-min(eromega02,eralpha02,erbeta02,eromega12,eralpha12,
erbeta12,eralpha22,erbeta22)
max2<-max(eromega02,eralpha02,erbeta02,eromega12,eralpha12,
erbeta12,eralpha22,erbeta22)

plot(x=c(0.5,1,2,3.5),y=c(min,min,max,max), type="n",ann=FALSE,axes=FALSE)
boxplot(list(eromega0,eralpha0,erbeta0), 
names=c(expression(omega-hat(omega)),expression(alpha-hat(alpha)),
expression(beta-hat(beta))),
xlab = expression(paste('estimation errors when ',
omega,'=1, ',alpha,'=0.3,',beta,'=0.6')),
add=TRUE)
title(expression(paste('Second-order stationarity, n=200')))

plot(x=c(0.5,1,2,3.5),y=c(min2,min2,max2,max2), type="n",ann=FALSE,axes=FALSE)
boxplot(list(eromega02,eralpha02,erbeta02), 
names=c(expression(omega-hat(omega)),expression(alpha-hat(alpha)),
expression(beta-hat(beta))),
xlab = expression(paste('estimation errors when ',
omega,'=1, ',alpha,'=0.3,',beta,'=0.6')),
add=TRUE)
title(expression(paste('Second-order stationarity, n=4000')))

plot(x=c(0.5,1,2,3.5),y=c(min,min,max,max), type="n",ann=FALSE,axes=FALSE)
boxplot(list(eromega1,eralpha1,erbeta1), 
names=c(expression(omega-hat(omega)),expression(alpha-hat(alpha)),
expression(beta-hat(beta))),
xlab = expression(paste('estimation errors when ',
omega,'=1, ',alpha,'=0.5,',beta,'=0.6')),
add=TRUE)
title(expression(paste('Strict stationarity, n=200')))

plot(x=c(0.5,1,2,3.5),y=c(min2,min2,max2,max2), type="n",ann=FALSE,axes=FALSE)
boxplot(list(eromega12,eralpha12,erbeta12), 
names=c(expression(omega-hat(omega)),expression(alpha-hat(alpha)),
expression(beta-hat(beta))),
xlab = expression(paste('estimation errors when ',
omega,'=1, ',alpha,'=0.5,')),
add=TRUE)
title(expression(paste('Strict stationarity, n=4000')))

plot(x=c(0.5,1,2,3.5),y=c(min,min,max,max), type="n",ann=FALSE,axes=FALSE)
boxplot(list(eromega2,eralpha2,erbeta2), 
names=c(expression(omega-hat(omega)),expression(alpha-hat(alpha)),
expression(beta-hat(beta))),
xlab = expression(paste('estimation errors when ',
omega,'=1, ',alpha,'=0.7,',beta,'=0.6')),
add=TRUE)
title(expression(paste('Non stationarity, n=200')))

plot(x=c(0.5,1,2,3.5),y=c(min2,min2,max2,max2), type="n",ann=FALSE,axes=FALSE)
boxplot(list(eromega22,eralpha22,erbeta22), 
names=c(expression(omega-hat(omega)),expression(alpha-hat(alpha)),
expression(beta-hat(beta))),
xlab = expression(paste('estimation errors when ',
omega,'=1, ',alpha,'=0.7,',beta,'=0.6')),
add=TRUE)
title(expression(paste('Non stationarity, n=4000')))
par(op)



