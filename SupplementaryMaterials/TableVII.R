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
sigma2u<-var(log(alpha*eta2[2:n]+beta))
stat<-sqrt(n)*mean(log(alpha*eta2[2:n]+beta))/sqrt(sigma2u)
pval<-1-pnorm(stat)
list(coef=c(omega,alpha,beta),minimum=res$objective,stat=stat,pval=pval)
}


############################### Tests de Stationnarité stricte sur series reelles
tests.seriesrelles<-function(indices,name.indices,alphainit,betainit,type='rend'){
nbind<-length(indices)
tab.stat<-rep('&',2*nbind) 

for(ind in (1:nbind)){
### lecture des donnees 
data <- read.table(indices[ind],header=TRUE,sep=",")
close<-rev(data$Close) 
date<-rev(data$Date) 
ntot<-length(close)
rend<-rep(0,ntot); rend[2:ntot]<-log(close[2:ntot]/close[1:(ntot-1)])*100 

ser<-rend[2:ntot]-mean(rend[2:ntot])

omegainit<-abs(ser[1])

test<-estimgarch11ettest(omegainit,alphainit,betainit,ser)

tab.stat[2*(ind-1)+1]<-round(test$stat,1)
}
tab.stat
}
#####################################################
alphainit<-0.05
betainit<-0.9

indices<-c("cac.csv","dax.csv","dja.csv","dowjones.csv","djt.csv",
    "dju.csv","ftse.csv","nasdaq.csv","nikkei.csv","smi.csv","sp500.csv")
name.indices<-c("CAC","DAX","DJA","DJI","DJT","DJU","FTSE","Nasdaq","Nikkei","SMI","SP500")


res<-tests.seriesrelles(indices,name.indices,alphainit,betainit)
res

