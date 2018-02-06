# simulation d'un GARCH(1,1) avec bruit iid Student normalise
garch.sim<- function (n, omega,  alpha, beta,df) {
ifelse(df>2,factor<-sqrt((df-2)/df),factor<-1)
ifelse(df>0,eta <- factor*rt(n,df=df),eta<-rnorm(n))
ht<-rep(0,n)
eps<-ht
ht[1]<-omega
eps[1]<-sqrt(ht[1])*eta[1]
for (t in 2:n) { 
ht[t] <- omega+alpha*eps[t-1]**2+beta*ht[t-1]
eps[t]<-sqrt(ht[t])*eta[t]}
return(eps[1:n])}

# estime J et kappa, pour un GARCH(1,1)
MatJ <- function(omega,alpha,beta,eps){
n <- length(eps)
sigma2<-as.numeric(n)
dersigma2<-matrix(nrow=3,ncol=n)
sigma2[1]<-omega
dersigma2[2:3,1]<-0
dersigma2[1,1]<-1
for (t in 2:n){
sigma2[t]<-omega+alpha*eps[t-1]**2+beta*sigma2[t-1]
dersigma2[1:3,t]<-c(1,eps[t-1]^2,sigma2[t-1])
dersigma2[1:3,t]<-dersigma2[1:3,t]+beta*dersigma2[1:3,t-1]
}
Jmat<-matrix(0,nrow=3,ncol=3)
for (t in 2:n){
Jmat<-Jmat+(dersigma2[1:3,t]/sigma2[t])%*%t(dersigma2[1:3,t]/sigma2[t])}
Jmat<-Jmat/n
eta<-rep(0,n)
eta[1:n]<-eps[1:n]/sqrt(sigma2[1:n])
#mean(eta[2:n]); mean(eta[2:n]^2)
kappa.eta<-mean(eta[2:n]^4)
I<-Jmat[2:3,2:3]-Jmat[2:3,1]%*%t(Jmat[2:3,1])/Jmat[1,1]
list(kappa=kappa.eta,Jmat=Jmat,Imat=I)
}

# estime un GARCH(1,1) 
# calcule la statistique du test H0: alpha<=alpha.star 

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

estimgarch11ettest<- function(omega0,alpha0,beta0,eps0,beta.star,
 grand=1/sqrt(.Machine$double.eps))
     {
valinit<-c(omega0,alpha0,beta0)
res <- nlminb(valinit,objf3, lower=c(0.0000001,0,0),upper=c(Inf,Inf,Inf), eps=eps0)
omega<-res$par[1]
alpha<-res$par[2]
beta<-res$par[3]
n<-length(eps0)
         sigma2<-vector(length=n)
         eta2<-vector(length=n)
        sigma2[1]<-omega
for (t in 2:n) sigma2[t]<-omega+alpha*eps0[t-1]**2+beta*sigma2[t-1]
         eta2[2:n]<-eps0[2:n]**2/sigma2[2:n]
#mean(eta2[2:n])
Jres<-MatJ(omega,alpha,beta,eps0)
I<-Jres$Imat
ok<-0
if(kappa(I)<grand)Iinv<-solve(I) else {ok<-1; Iinv<-diag(rep(1,nrow(I)))}
kap<-Jres$kappa 
den<-sqrt((kap-1)*Iinv[2,2])
stat<-sqrt(n)*(beta-beta.star)/den
pval<-1-pnorm(stat)
list(coef=c(omega,alpha,beta),minimum=res$objective,stat=stat,pval=pval,ok=ok)
}



##############
set.seed(7)
df<-7
omega0<-0.1
alpha0<-0.5
betavec<-c(0.61,0.64,0.67,0.70,0.73,0.76,0.79)
beta.star<-0.7
niter<-length(betavec)
njter<-c(500,2000,4000)
for(jter in 1:length(njter)){
n<-njter[jter]
nrep<-1000
pval<-matrix(nrow=niter,ncol=nrep)
para.omega<-matrix(nrow=niter,ncol=nrep)
para.alpha<-matrix(nrow=niter,ncol=nrep)
para.beta<-matrix(nrow=niter,ncol=nrep)
for (i in 1:niter) {
beta0<-betavec[i]
for (j in 1:nrep) {
ser<-garch.sim(n, omega0,  alpha0, beta0,df)
resestimgarch11 <- estimgarch11ettest(abs((1+0.1/2-0.1*runif(1))*omega0),
abs((1+0.1/2-0.1*runif(1))*alpha0),abs((1+0.1/2-0.1*runif(1))*beta0),
ser,beta.star)
if(resestimgarch11$ok!=0)print(c(i,resestimgarch11$ok,resestimgarch11$pval
))
pval[i,j]<-resestimgarch11$pval
para.omega[i,j]<-resestimgarch11$coef[1]
para.alpha[i,j]<-resestimgarch11$coef[2]
para.beta[i,j]<-resestimgarch11$coef[3]
     }}
freq.rej<-rep(0,niter)
for (i in 1:niter) {
freq.rej[i]<-100*length(which(pval[i,]< 0.05))/nrep}
print(freq.rej)}

