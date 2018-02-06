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


# estime un GARCH(1,1) 
# calcule la statistique du test de stationnarité stricte

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

estimgarch11ettest<- function(omega0,alpha0,beta0,eps0)
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
sigma2u<-var(log(alpha*eta2[2:n]+beta))
stat<-sqrt(n)*mean(log(alpha*eta2[2:n]+beta))/sqrt(sigma2u)
pval<-1-pnorm(stat)
list(coef=c(omega,alpha,beta),minimum=res$objective,stat=stat,pval=pval)
}


##############
set.seed(7)
omega0<-0.1
beta0<-0.8
df<-7
alphavec<-c(0.18,0.20,0.22,0.25752,0.28,0.30,0.31)
niter<-length(alphavec)
njter<-c(500,2000,4000)
for(jter in 1:length(njter)){
n<-njter[jter]
nrep<-1000
pval<-matrix(nrow=niter,ncol=nrep)
para.omega<-matrix(nrow=niter,ncol=nrep)
para.alpha<-matrix(nrow=niter,ncol=nrep)
para.beta<-matrix(nrow=niter,ncol=nrep)
for (i in 1:niter) {
alpha0<-alphavec[i]
for (j in 1:nrep) {
ser<-garch.sim(n, omega0,  alpha0, beta0,df)
resestimgarch11 <- estimgarch11ettest(abs((1+0.1/2-0.1*runif(1))*omega0),
abs((1+0.1/2-0.1*runif(1))*alpha0),abs((1+0.1/2-0.1*runif(1))*beta0),
ser)
pval[i,j]<-resestimgarch11$pval
para.omega[i,j]<-resestimgarch11$coef[1]
para.alpha[i,j]<-resestimgarch11$coef[2]
para.beta[i,j]<-resestimgarch11$coef[3]
     }}
freq.rej<-rep(0,niter)
for (i in 1:niter) {
freq.rej[i]<-100*length(which(pval[i,]< 0.05))/nrep}

print(freq.rej)}
