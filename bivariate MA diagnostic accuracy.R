
###
### This R syntax performs meta-analysis of the diagnostic accuracy of a binary diagnostic test
###      averaged over a collection of independent studies
###
### The data that each study is supposed to contribute is in the form of the two-by-two table of
###      the diagnostic test versus the gold-standard and summarized with the numbers of patients
###      for whom the diagnostic test is (1) correctly positive (true positives, TN), (2) incorrectly 
###      negative (false negatives, FN), correctly negative (true negatives, TN), and incorrectly
###      positive (false positives, FP).
###
### Example data are in the R-package mada and taking as an example mada's-AuditC, the data that must be
###      provided are in the form of a data.frame/matrix and each row represents a separate study:
###
###                TP FN   FP   TN
###            1   47  9  101  738
###            2  126 51  272 1543
###            3   19 10   12  192
###            4   36  3   78  276
###            .   .   .   .   .
###
### The analysis uses an extension of the "bivariate random-effects model" first developed by Reitsma et al (2005).
###     The numbers of TP and TN are modeled using binomial distributions with study-specific
###     sensitivity and specificity, and the logit-transformed study-specific sens and spec are assumed to be sampled
###     from a bivariate normal distribution. Its mean and covariance-matrix are the structural parameters of the model
###     that are estimated.
###
### We use a Bayesian algorithm with Gibbs sampling to estimate the mean and covariance parameters. The current code 
###     depends on the JAGS-program (Just Another Gibbs Sampler), which needs to be installed (mcmc-jags.sourceforge.net),
###     and on the R-package rjags, which needs to be attached.
###     There are no other dependencies, but if example data are wanted, the R-package mada is recommended. This 
###     package can analyze the same kind of data and is useful for comparison sake anyway.
###
### Results of the analysis are estimates of the mean and covariance and 95% credibility intervals. In addition transformations
###     of these estimates are calculated, in particular estimates of the parameters of the Rutter-Gatsonis sROC-model, the ROC-curve 
###     itself and the area-under-curve. Results are illustrated with the common ROC-plot for this type of data including 
###     confidence and prediction ellipses for the mean.
###


###############################################################################################################################################################
## some analysis- and plot-constants
alfa=0.05
symbool=16
maxsize=1.5
minsize=0.5
x=seq(0.01,0.99,0.01)
P=length(x)
zalfa=qnorm((1-alfa/2))
delta=0.5


###############################################################################################################################################################
# open connection to jags, load other useful packages, and define functions
library(rjags)
logit=function(x){log(x/(1-x))}
antilogit=function(x){1/(1+exp(-x))}
confellips=function(mu,sigma,alfa,npoints) {
   es <- eigen(sigma)
   e1 <- es$vec %*% diag(sqrt(es$val))
   r1 <- sqrt(qchisq(1 - alfa, 2))
   theta <- seq(0, 2 * pi, len = npoints)
   v1 <- cbind(r1 * cos(theta), r1 * sin(theta))
   pts = t(mu - (e1 %*% t(v1)))
}


###############################################################################################################################################################
## define the binomial-bivariate normal model in a jags-model statement
modelstring1 <- "
   model {
      # likelihood of the data
      for (i in 1:N) {
         logitprobs[i,1:2] ~ dmnorm(mu[1:2],invSigma[1:2,1:2])
         sens[i] <- 1/(1+exp(-logitprobs[i,1]))
         spec[i] <- 1/(1+exp(-logitprobs[i,2]))
         TP[i] ~ dbinom(sens[i],ncases[i])
         TN[i] ~ dbinom(spec[i],ncontrols[i])
      }
      # priors
      mu[1:2] ~ dmnorm(zeroos[1:2],precision[1:2,1:2])
      invSigma[1:2,1:2] ~ dwish(R[1:2,1:2],dfx)
      # functions of model-parameters
      Sigma[1:2,1:2] <- inverse(invSigma[1:2,1:2])
      
      # Rutter and Gatsonis parameters: Harbord et al. Biostatistics, 2007, 8, 2, 239-251.
      beta <- 0.5*log(Sigma[2,2]/Sigma[1,1])
      LAMBDA <- sqrt(sqrt(Sigma[2,2])/sqrt(Sigma[1,1]))*mu[1] + sqrt(sqrt(Sigma[1,1])/sqrt(Sigma[2,2]))*mu[2]
      THETA  <- 0.5*(sqrt(sqrt(Sigma[2,2])/sqrt(Sigma[1,1]))*mu[1] - sqrt(sqrt(Sigma[1,1])/sqrt(Sigma[2,2]))*mu[2])
      skwtheta <- 0.5*(sqrt(Sigma[1,1]*Sigma[2,2])-Sigma[1,2])
      skwalfa  <- 2.0*(sqrt(Sigma[1,1]*Sigma[2,2])+Sigma[1,2])

      # ROC given 1-specificity (auc1) or given sensitivity (auc2)
      b1 <- -Sigma[1,2]/Sigma[2,2]
      a1 <- mu[1] + b1 * mu[2]
      b2 <- -Sigma[1,2]/Sigma[1,1]
      a2 <- mu[2] + b2 * mu[1]      
      y1[1] <- 1/(1+exp(-(a1 + b1 * log(x[1]/(1-x[1])))))
      y2[1] <- 1/(1+exp(-(a2 + b2 * log(x[1]/(1-x[1])))))
      y3[1] <- 1/(1+exp(-(LAMBDA*exp(beta/2) + exp(beta)*log(x[1]/(1-x[1])))))
      opp1[1] <- (x[1]-0)*(y1[1]-0)/2
      opp2[1] <- (x[1]-0)*(y2[1]-0)/2
      opp3[1] <- (x[1]-0)*(y3[1]-0)/2
      for (j in 2:P) {
         y1[j] <- 1/(1+exp(-(a1 + b1 * log(x[j]/(1-x[j])))))
         opp1[j] <- (x[j]-x[(j-1)])*y1[(j-1)]+(x[j]-x[(j-1)])*(y1[j]-y1[(j-1)])/2
         y2[j] <- 1/(1+exp(-(a2 + b2 * log(x[j]/(1-x[j])))))
         opp2[j] <- (x[j]-x[(j-1)])*y2[(j-1)]+(x[j]-x[(j-1)])*(y2[j]-y2[(j-1)])/2
         y3[j] <- 1/(1+exp(-(LAMBDA*exp(beta/2) + exp(beta)*log(x[j]/(1-x[j])))))
         opp3[j] <- (x[j]-x[(j-1)])*y3[(j-1)]+(x[j]-x[(j-1)])*(y3[j]-y3[(j-1)])/2
      }
      opp1[(P+1)] <- (1-x[P])*y1[P]+(1-x[P])*(1-y1[P])/2      
      auc1 <- sum(opp1)
      opp2[(P+1)] <- (1-x[P])*y2[P]+(1-x[P])*(1-y2[P])/2      
      auc2 <- sum(opp2)      
      opp3[(P+1)] <- (1-x[P])*y3[P]+(1-x[P])*(1-y3[P])/2      
      aucRG <- sum(opp3) 
   }"


###############################################################################################################################################################
## Import example data from the mada package
library(mada)
data("AuditC")
data("Dementia")
data("IAQ")
data("SAQ")
data("smoking")
d=AuditC
#d=Dementia
#d=IAQ
#d=SAQ
#d=smoking

## or READ DATA from a local file, foir instance
#d = read.csv("arfi_af.csv")


###############################################################################################################################################################
## calculate number of cases and number of true positives
ncases=d$TP+d$FN
TP=d$TP
## number of controls and number of true negatives
ncontrols=d$TN+d$FP
TN=d$TN


###############################################################################################################################################################
##  1. overview of data and statistics (add delta=0.5 to all TP, FN, TN, FP numbers)
sens=(TP+delta)/(ncases+2*delta)
sesens=sqrt(sens*(1-sens)/ncases)
spec=(TN+delta)/(ncontrols+2*delta)
sespec=sqrt(spec*(1-spec)/ncontrols)
data.frame(TP,ncases,sens,sesens,TN,ncontrols,spec,sespec)

## plot of the logit-statistics in ROC-space
#c(mean(logit(sens)),mean(logit(spec)))
#c(var(logit(sens)),var(logit(spec)))
#cov(logit(sens),logit(spec))
#cor((spec),(sens),method="spearman")
#summary(lm(logit(sens) ~ logit(spec),weights=1/(ncontrols*spec*(1-spec))))
symbolsize=minsize+(maxsize-minsize)*((ncases+ncontrols)-min((ncases+ncontrols)))/(max((ncases+ncontrols))-min((ncases+ncontrols)))
symbolsize=rep(1,length(sens))
dev.new()
plot(logit(1-spec),logit(sens),type="n",pch=symbool,xlab="logit(1 - specificity)",ylab="logit(sensitivity)",main="",xlim=c(-8,+8),ylim=c(-8,+8))
for (i in 1:length(ncases)) {
   points(logit(1-spec)[i],logit(sens)[i],pch=symbool,cex=symbolsize[i])
   lines(logit(1-c(spec[i],spec[i])),
      c((logit(sens[i])-zalfa/(ncases[i]*sens[i]*(1-sens[i]))),(logit(sens[i])+zalfa/(ncases[i]*sens[i]*(1-sens[i])))),
      lty=3,lwd=1)
   lines(c((logit(1-spec[i])-zalfa/(ncontrols[i]*spec[i]*(1-spec[i]))),(logit(1-spec[i])+zalfa/(ncontrols[i]*spec[i]*(1-spec[i])))),
      logit(c(sens[i],sens[i])),lty=3,lwd=1)
}


###############################################################################################################################################################
## 2. do the analysis
## supply the data in a list-format  (----> what should R look like? <----)
dd1=list(N=length(ncases), ncases=ncases,TP=TP, ncontrols=ncontrols,TN=TN,
         zeroos=c(0,0), precision=diag(0.001,2), R=diag(0.1,2), dfx=2, x=x, P=P)
## check and compile the model-syntax
model1=jags.model(textConnection(modelstring1), data=dd1, n.chains = 2, n.adapt=0)
## burn-in
update(model1, n.iter=50000)
## posterior sampling
output.raw1 = coda.samples(model=model1, variable.names=c("mu","Sigma"), n.iter=50000, thin=10)
output.raw1a = coda.samples(model=model1, variable.names=c("sens","spec"), n.iter=50000, thin=10)
output.raw1b = coda.samples(model=model1, variable.names=c("auc1","auc2","aucRG","b1"), n.iter=50000, thin=10)
output.raw1c = coda.samples(model=model1, variable.names=c("beta","LAMBDA","THETA","skwtheta","skwalfa"), n.iter=50000, thin=10)
## check convergence
dev.new()
#plot(output.raw1,ask=T)
plot(output.raw1b)
## posterior statistics
summary(output.raw1)              # parameters of the bivariate model
summary(output.raw1a)             # study-specific sensitivity and specificity
summary(output.raw1b)             # auc-values: given specificity or given sensitivity or according to the Rutter-Gatsonis SROC-model
plot(output.raw1c,ask=T)
summary(output.raw1c)             # parameters of the Rutter-Gatsonis SROC-model

# derive some estimates to construct the ROC-plot
namen=row.names(summary(output.raw1)$statistics)
##statistiek=sapply(namen,function(x){strsplit(x,"[",fixed=TRUE)[[1]][1]})
mu1=summary(output.raw1)$statistics[which(namen=="mu[1]"),1]
sekw1=summary(output.raw1)$statistics[which(namen=="mu[1]"),2]^2        
mu2=summary(output.raw1)$statistics[which(namen=="mu[2]"),1]
sekw2=summary(output.raw1)$statistics[which(namen=="mu[2]"),2]^2        
s1kw=summary(output.raw1)$statistics[which(namen=="Sigma[1,1]"),1]
s2kw=summary(output.raw1)$statistics[which(namen=="Sigma[2,2]"),1]
s12= summary(output.raw1)$statistics[which(namen=="Sigma[2,1]"),1]
covmatrix=matrix(c(s2kw,-s12,-s12,s1kw),nrow=2,ncol=2)
se12=cov(output.raw1[[1]][,5],output.raw1[[1]][,6])#/sqrt(length(output.raw1[[1]][,6]))       
sematrix=matrix(c(sekw2, -se12, -se12,sekw1),nrow=2,ncol=2)
symbolsize=minsize+(maxsize-minsize)*((ncases+ncontrols)-min((ncases+ncontrols)))/(max((ncases+ncontrols))-min((ncases+ncontrols)))

## ROC plot with 95% confidence ellipse for the mean and 95% prediction ellipse
dev.new()
plot((1-spec),sens,type="n",pch=symbool,xlab="specificity",ylab="sensitivity",main="",xlim=c(0,1),ylim=c(0,1),xaxt="n")
axis(1,at=seq(0,1,0.2),labels=c("1.0","0.8","0.6","0.4","0.2","0.0"))
for (i in 1:length(ncases)) {
   points(1-spec[i],sens[i],pch=symbool,cex=symbolsize[i])
   lines(1-c(spec[i],spec[i]),
      antilogit(c((logit(sens[i])-zalfa/(ncases[i]*sens[i]*(1-sens[i]))),(logit(sens[i])+zalfa/(ncases[i]*sens[i]*(1-sens[i]))))),
      lty=3,lwd=1)
   lines(antilogit(c((logit(1-spec[i])-zalfa/(ncontrols[i]*spec[i]*(1-spec[i]))),(logit(1-spec[i])+zalfa/(ncontrols[i]*spec[i]*(1-spec[i]))))),
      c(sens[i],sens[i]),lty=3,lwd=1)
}
b=(-s12/s2kw)
a=mu1+b*mu2
help=seq(min(1-spec),max(1-spec),0.01)
help1=antilogit(a+b*logit(help))
help1[help1>max(sens)]=NA
help1[help1<min(sens)]=NA
lines(help,help1,col=2)
points(antilogit(-mu2),antilogit(mu1),pch=symbool,col=2,cex=0.75)
punten1=confellips(mu=c(-mu2,mu1),sigma=sematrix,alfa=alfa,npoints=1000)
lines(antilogit(punten1[,1]),antilogit(punten1[,2]),col=2)
punten2=confellips(mu=c(-mu2,mu1),sigma=covmatrix,alfa=alfa,npoints=1000)
lines(antilogit(punten2[,1]),antilogit(punten2[,2]),col=2,lty=3)
## Rutter & Gatsonis ROC curve
beta = 0.5*log(s1kw/s2kw)
LAMBDA = sqrt(sqrt(s2kw)/sqrt(s1kw))*mu1 + sqrt(sqrt(s1kw)/sqrt(s2kw))*mu2
yyy = antilogit(LAMBDA*exp(beta/2) + exp(beta)*logit(help))
lines(help,yyy,col=4)
legend(x=0.5,y=0.2,legend=c("naive","Rutter-Gatsonis SROC-model"),col=c(2,4),lty=1,bty="n",adj=0)



###############################################################################################################################################################
## 3. do a similar analysis with mada 
#library(mada)
xx=data.frame(TP=TP,FN=ncases-TP,FP=ncontrols-TN,TN=TN)
madares=reitsma(xx, method = "reml",predict=T,sroc.type="naive")
summary(madares)
AUC(madares,sroc.type="naive")    # 0.688  vgl de mijne is 0.713
dev.new()
plot(madares,predict=T,type="naive")
for (i in 1:length(ncases)) {
   points(1-spec[i],sens[i],pch=symbool,cex=symbolsize[i])
   lines(1-c(spec[i],spec[i]),
      antilogit(c((logit(sens[i])-zalfa/(ncases[i]*sens[i]*(1-sens[i]))),(logit(sens[i])+zalfa/(ncases[i]*sens[i]*(1-sens[i]))))),
      lty=3,lwd=1)
   lines(antilogit(c((logit(1-spec[i])-zalfa/(ncontrols[i]*spec[i]*(1-spec[i]))),(logit(1-spec[i])+zalfa/(ncontrols[i]*spec[i]*(1-spec[i]))))),
      c(sens[i],sens[i]),lty=3,lwd=1)
}
madares=reitsma(xx, method = "reml",predict=T,sroc.type="ruttergatsonis")
summary(madares)
AUC(madares,sroc.type="ruttergatsonis")   #0.91
dev.new()
plot(madares,predict=T,type="ruttergatsonis")
for (i in 1:length(ncases)) {
   points(1-spec[i],sens[i],pch=symbool,cex=symbolsize[i])
   lines(1-c(spec[i],spec[i]),
      antilogit(c((logit(sens[i])-zalfa/(ncases[i]*sens[i]*(1-sens[i]))),(logit(sens[i])+zalfa/(ncases[i]*sens[i]*(1-sens[i]))))),
      lty=3,lwd=1)
   lines(antilogit(c((logit(1-spec[i])-zalfa/(ncontrols[i]*spec[i]*(1-spec[i]))),(logit(1-spec[i])+zalfa/(ncontrols[i]*spec[i]*(1-spec[i]))))),
      c(sens[i],sens[i]),lty=3,lwd=1)
}



