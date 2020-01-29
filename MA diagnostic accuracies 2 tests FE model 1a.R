
########## FIXED-EFFECTS VARIANT


###
### meta-analysis of the comparative diagnostic accuracies of two diagnostic tests. 
###
### In this version we meta-analyze the results of a set of independent studies that all compared the same two
###     diagnostic tests (A and B) to a reference standard and that reported the results with the numbers of
###     true positives (TP), true negatives (TN), false positives (FP) and false negatives (FN) for both diagnostic tests. 
###     In this version we do not use the correlations of the outcomes of both tests in the same patients, because such 
###     information is unfortunately rarely reported in scientific reports of this kind of diagnostic accuracy studies. 
###     So in effect the outcomes of tests A and B are handled as if they are independent and the analysis is done as if
###     the patients tested with test A are different than the patients tested with test B.
###
### Input data is in the form of a data.frame/matrix where every study is a separate row and there are 8 columns with the
###     TP, FN, TN, FP for test A and the TP, FN, TN, FP for test B. An example is given below.
###

d = as.data.frame(matrix(
     c(14, 26,  90, 10, 28, 12,  68, 32,
       18, 37,  50,  0, 40, 15,  26, 24,
        2,  4, 185,  8,  4,  2, 159, 34,
       24, 30,  55,  4, 42, 15, 101, 38,
       38, 33,  11,  3, 64,  7,  13,  4),
       nrow=5,ncol=8,byrow=TRUE)
    )
names(d)=c("TPsA","FNsA","TNsA","FPsA","TPsB","FNsB","TNsB","FPsB")

###
### The meta-analysis model is based on binomial distributions of the TP and TN of tests A and B given the overall -non varying- 
###     sensitivities and specificities of tests A and B. 
###     So the model is a fixed-effects model.
###
### The current code depends on the JAGS-program (Just Another Gibbs Sampler), which needs to be installed (mcmc-jags.sourceforge.net),
###     and on the R-package rjags, which needs to be attached. There are no other dependencies.
###
### Results of the analysis are estimates of the logit-transformed sensitivities and specificities and 95% credibility intervals. 
###     In addition transformations of these estimates are calculated, in particular estimates of odds ratio's of the sensitivity of 
###     test A vs test B and of the specificity of test A vs test B. Also, the difference of the two mean sensitivities and of the two specificities
###     are calculated together with their 95% credibility intervals.
###


### attach the rjags-package
library(rjags)

### get the data from the data-frame d
N=nrow(d) 
TPsA   =d$TPsA          #c(14,18,2,24,38)
nAcases=d$TPsA+d$FNsA   #c(40,55,6,54,71)
TPsB   =d$TPsB          #c(28,40,4,42,64)
nBcases=d$TPsB+d$FNsB   #c(40,55,6,57,71)
TNsA      =d$TNsA          #c( 90,50,185,55,11)
nAcontrols=d$TNsA+d$FPsA   #c(100,50,193,59,14)
TNsB      =d$TNsB          #c(68,26,159,101,13)
nBcontrols=d$TNsB+d$FPsB   #c(100,50,193,139,17)

### define non-informative hyper prior for 
# the logit-transformed parameters (multivariate normal)
zero4=c(0,0,0,0)   
prec4=diag(0.001,nrow=4,ncol=4)  

##### put data and hyper priors in a list-object to be transferred to JAGS
jagsdatalist=list(
  N=N, TPsA=TPsA,TPsB=TPsB,TNsA=TNsA,TNsB=TNsB, 
  nAcases=nAcases,nBcases=nBcases,nAcontrols=nAcontrols,nBcontrols=nBcontrols,
  zero4=zero4,prec4=prec4)

### model in JAGS code
modelstring = "
   model {
      sensB <- 1/(1+exp(-1*(meanalpha[1])))
      sensA <- 1/(1+exp(-1*(meanalpha[1]+meanalpha[2])))
      specB <- 1/(1+exp(-1*(meanalpha[3])))
      specA <- 1/(1+exp(-1*(meanalpha[3]+meanalpha[4])))

      for (i in 1:N){ 
         TPsA[i] ~ dbin(sensA,nAcases[i])
         TPsB[i] ~ dbin(sensB,nBcases[i])
         TNsA[i] ~ dbin(specA,nAcontrols[i])
         TNsB[i] ~ dbin(specB,nBcontrols[i])
       }

      meanalpha[1:4] ~ dmnorm(zero4[1:4],prec4[1:4,1:4])

      meansensitivityB <- 1/(1+exp(-1*(meanalpha[1])))
      meansensitivityA <- 1/(1+exp(-1*(meanalpha[1]+meanalpha[2])))
      meanspecificityB <- 1/(1+exp(-1*(meanalpha[3])))
      meanspecificityA <- 1/(1+exp(-1*(meanalpha[3]+meanalpha[4])))

      sensdifferenceBvsA <- 1/(1+exp(-1*(meanalpha[1])))  -  1/(1+exp(-1*(meanalpha[1]+meanalpha[2])))
      specdifferenceBvsA <- 1/(1+exp(-1*(meanalpha[3])))  -  1/(1+exp(-1*(meanalpha[3]+meanalpha[4])))
      sensoddsratioAvsB <- exp(meanalpha[2])
      specoddsratioAvsB <- exp(meanalpha[4])
      
   }  
"

### useful functions
logit = function(x){log(x/(1-x))}
invlogit = function(x) {exp(x)/(1+exp(x))}
confellips=function(mu,sigma,alfa,npoints) {
   es <- eigen(sigma)
   e1 <- es$vec %*% diag(sqrt(es$val))
   r1 <- sqrt(qchisq(1 - alfa, 2))
   theta <- seq(0, 2 * pi, len = npoints)
   v1 <- cbind(r1 * cos(theta), r1 * sin(theta))
   pts = t(mu - (e1 %*% t(v1)))
}


##### do the JAGS analysis
inits = list(meanalpha=c(0,0,0,0))                                                                     # starting values
m = jags.model(textConnection(modelstring), data=jagsdatalist, inits, n.chains=2)                      # compile thee analysis 
                                                                                                       #    (i.e. check whether model and data fit together)
update(m, 50000)                                                                                       # do 50000 iterations 
x = coda.samples(m, c("meanalpha","sensdifferenceBvsA","specdifferenceBvsA",
       "sensoddsratioAvsB","specoddsratioAvsB"), n.iter=50000)                                         # draw 50000 parameters from the posterior-distributions
plot(x,ask=TRUE)                                                                                       # check for convergence 
           



                                                                                            #     trace-plots should display stable chaos
# summary of the structural and derived parameters
summary(x)                                                                                             # calculate statistics

dev.new()
par(mfrow=c(2,2))
####### since the convergence was OK for the datasets tried out so far, the results below are based only on sample drawn from the first chain
hist(x[[1]][,which(colnames(x[[1]])=="sensoddsratioAvsB")],xlab="odds ratio A vs B",main="sensitivity")   # a histogram of the samples for one parameter
quantile(x[[1]][,which(colnames(x[[1]])=="sensoddsratioAvsB")],probs=c(0.025,0.25,0.50,0.75,0.975))       # a 95% credibility interval for this parameter
hist(x[[1]][,which(colnames(x[[1]])=="specoddsratioAvsB")],xlab="odds ratio A vs B",main="specificity")   # a histogram of the samples for one parameter
quantile(x[[1]][,which(colnames(x[[1]])=="specoddsratioAvsB")],probs=c(0.025,0.25,0.50,0.75,0.975))       # a 95% credibility interval for this parameter
hist(x[[1]][,which(colnames(x[[1]])=="sensdifferenceBvsA")],xlab="difference B vs A",main="sensitivity")  # a histogram of the samples for another parameter
quantile(x[[1]][,which(colnames(x[[1]])=="sensdifferenceBvsA")],probs=c(0.025,0.25,0.50,0.75,0.975))      # a 95% credibility interval for this other parameter
hist(x[[1]][,which(colnames(x[[1]])=="specdifferenceBvsA")],xlab="difference B vs A",main="specificity")  # a histogram of the samples for another parameter
quantile(x[[1]][,which(colnames(x[[1]])=="specdifferenceBvsA")],probs=c(0.025,0.25,0.50,0.75,0.975))      # a 95% credibility interval for this other parameter
par(mfrow=c(1,1))


x1 = coda.samples(m,c("meansensitivityA","meansensitivityB","meanspecificityA","meanspecificityB"), n.iter=50000)
#help1=t(apply(x1[[1]],2,function(y){quantile(y,probs=c(0.025,0.5,0.975))}))
help1=summary(x1)$quantiles[,c(1,3,5)]


dev.new()
#pdf("fig1.pdf")
layout(matrix(c(1,2,3,4),nrow=2,ncol=2))
plot(1,1,xlab="sensitivity",ylab="study number",main="",xlim=c(0,1),ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[1,2],(-N-2),pch=16,cex=1)
lines(c(help1[1,1],help1[1,3]),c(-N-2,-N-2))
pp=(TPsA+0.5)/(nAcases+1)
lowxx=invlogit(logit(pp) - 1.96*sqrt(1/(nAcases*pp*(1-pp))))
higxx=invlogit(logit(pp) + 1.96*sqrt(1/(nAcases*pp*(1-pp))))
points(pp,(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   lines( c(lowxx[i],higxx[i]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[1,2],lty=3,lwd=0.75)
plot(1,1,xlab="sensitivity",ylab="study number",main="",xlim=c(0,1),ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[2,2],(-N-2),pch=16,cex=1)
lines(c(help1[2,1],help1[2,3]),c(-N-2,-N-2))
pp=(TPsB+0.5)/(nBcases+1)
lowxx=invlogit(logit(pp) - 1.96*sqrt(1/(nBcases*pp*(1-pp))))
higxx=invlogit(logit(pp) + 1.96*sqrt(1/(nBcases*pp*(1-pp))))
points(pp,(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   lines( c(lowxx[i],higxx[i]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[2,2],lty=3,lwd=0.75)
plot(1,1,xlab="specificity",ylab="study number",main="",xlim=c(0,1),ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[3,2],(-N-2),pch=16,cex=1)
lines(c(help1[3,1],help1[3,3]),c(-N-2,-N-2))
pp=(TNsA+0.5)/(nAcontrols+1)
lowxx=invlogit(logit(pp) - 1.96*sqrt(1/(nAcontrols*pp*(1-pp))))
higxx=invlogit(logit(pp) + 1.96*sqrt(1/(nAcontrols*pp*(1-pp))))
points(pp,(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   lines( c(lowxx[i],higxx[i]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[3,2],lty=3,lwd=0.75)
plot(1,1,xlab="specificity",ylab="study number",main="",xlim=c(0,1),ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[4,2],(-N-2),pch=16,cex=1)
lines(c(help1[4,1],help1[4,3]),c(-N-2,-N-2))
pp=(TNsB+0.5)/(nBcontrols+1)
lowxx=invlogit(logit(pp) - 1.96*sqrt(1/(nBcontrols*pp*(1-pp))))
higxx=invlogit(logit(pp) + 1.96*sqrt(1/(nBcontrols*pp*(1-pp))))
points(pp,(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   lines( c(lowxx[i],higxx[i]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[4,2],lty=3,lwd=0.75)
mtext("test A",side=3,outer=TRUE,line=-2)
mtext("test B",side=3,outer=TRUE,line=-23)
#dev.off()

dev.new()
#pdf("fig2.pdf")
#help0=t(apply(x[[1]],2,function(y){quantile(y,probs=c(0.025,0.5,0.975))}))
help0=summary(x)$quantiles[,c(1,3,5)]
layout(matrix(c(1,2,3,4),nrow=2,ncol=2))
plot(1,1,xlab="difference of sensitivities B vs A",ylab="study number",xlim=c(-1,1),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help0[which(row.names(help0)=="sensdifferenceBvsA"),2],(-N-2),pch=16,cex=1)
lines(c(help0[which(row.names(help0)=="sensdifferenceBvsA"),1],help0[which(row.names(help0)=="sensdifferenceBvsA"),3]),c(-N-2,-N-2))
pp2=(TPsA+0.5)/(nAcases+1)
pp1=(TPsB+0.5)/(nBcases+1)
points(pp1-pp2,(-N:-1),pch="+",cex=1)
lowxx=(pp1-pp2)-1.96*sqrt((pp2*(1-pp2)/nAcases) + (pp1*(1-pp1)/nBcases))
higxx=(pp1-pp2)+1.96*sqrt((pp2*(1-pp2)/nAcases) + (pp1*(1-pp1)/nBcases))
for (i in 1:N) {
   lines( c(lowxx[i],higxx[i]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help0[which(row.names(help0)=="sensdifferenceBvsA"),2],lty=3,lwd=0.75)



plot(1,1,xlab="difference of specificities B vs A",ylab="study number",xlim=c(-1,1),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help0[which(row.names(help0)=="specdifferenceBvsA"),2],(-N-2),pch=16,cex=1)
lines(c(help0[which(row.names(help0)=="specdifferenceBvsA"),1],help0[which(row.names(help0)=="specdifferenceBvsA"),3]),c(-N-2,-N-2))
pp2=(TNsA+0.5)/(nAcontrols+1)
pp1=(TNsB+0.5)/(nBcontrols+1)
points(pp1-pp2,(-N:-1),pch="+",cex=1)
lowxx=(pp1-pp2)-1.96*sqrt((pp2*(1-pp2)/nAcontrols) + (pp1*(1-pp1)/nBcontrols))
higxx=(pp1-pp2)+1.96*sqrt((pp2*(1-pp2)/nAcontrols) + (pp1*(1-pp1)/nBcontrols))
for (i in 1:N) {
   lines( c(lowxx[i],higxx[i]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help0[which(row.names(help0)=="specdifferenceBvsA"),2],lty=3,lwd=0.75)


pp2=(TPsA+0.5)/(nAcases+1)
pp1=(TPsB+0.5)/(nBcases+1)
ORs=(pp2/(1-pp2)) / (pp1/(1-pp1))
lowxx=exp(log(ORs)-1.96*sqrt( (1/(TPsA+0.5))+(1/(d$FNsA+0.5))+(1/(TPsB+0.5))+(1/(d$FNsB+0.5)) ))
higxx=exp(log(ORs)+1.96*sqrt( (1/(TPsA+0.5))+(1/(d$FNsA+0.5))+(1/(TPsB+0.5))+(1/(d$FNsB+0.5)) ))
minx=0 #min(lowxx,na.rm=TRUE)
maxx=1.1 #max(higxx,na.rm=TRUE)
plot(1,1,xlab="OR of sensitivities A vs B",ylab="study number",xlim=c(minx,maxx),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help0[which(row.names(help0)=="sensoddsratioAvsB"),2],(-N-2),pch=16,cex=1)
lines(c(help0[which(row.names(help0)=="sensoddsratioAvsB"),1],help0[which(row.names(help0)=="sensoddsratioAvsB"),3]),c(-N-2,-N-2))
points(ORs,(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   lines( c(lowxx[i],higxx[i]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help0[which(row.names(help0)=="sensoddsratioAvsB"),2],lty=3,lwd=0.75)
pp2=(TNsA+0.5)/(nAcontrols+1)
pp1=(TNsB+0.5)/(nBcontrols+1)
ORs=(pp2/(1-pp2)) / (pp1/(1-pp1))
lowxx=exp(log(ORs)-1.96*sqrt( (1/(TNsA+0.5))+(1/(d$FPsA+0.5))+(1/(TNsB+0.5))+(1/(d$FPsB+0.5)) ))
higxx=exp(log(ORs)+1.96*sqrt( (1/(TNsA+0.5))+(1/(d$FPsA+0.5))+(1/(TNsB+0.5))+(1/(d$FPsB+0.5)) ))
minxx=min(lowxx,na.rm=TRUE)
maxxx=10#max(higxx,na.rm=TRUE)
plot(1,1,xlab="OR of specificities A vs B",ylab="study number",xlim=c(minxx,maxxx),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help0[which(row.names(help0)=="specoddsratioAvsB"),2],(-N-2),pch=16,cex=1)
lines(c(help0[which(row.names(help0)=="specoddsratioAvsB"),1],help0[which(row.names(help0)=="specoddsratioAvsB"),3]),c(-N-2,-N-2))
points(ORs,(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   lines( c(lowxx[i],higxx[i]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help0[which(row.names(help0)=="specoddsratioAvsB"),2],lty=3,lwd=0.75)
#dev.off()

dev.new()
#pdf("fig3.pdf")
plot(1,1,xlab="1-specificity",ylab="sensitivity",xlim=c(0,1),ylim=c(0,1),type="n",main="means and 95% confidence ellipses in ROC-space for both tests")
pp1=(TPsA+0.5)/(nAcases+1)
lowxx1=invlogit(logit(pp1) - 1.96*sqrt(1/(nAcases*pp1*(1-pp1))))
higxx1=invlogit(logit(pp1) + 1.96*sqrt(1/(nAcases*pp1*(1-pp1))))
pp2=(TNsA+0.5)/(nAcontrols+1)
lowxx2=invlogit(logit(pp2) - 1.96*sqrt(1/(nAcontrols*pp2*(1-pp2))))
higxx2=invlogit(logit(pp2) + 1.96*sqrt(1/(nAcontrols*pp2*(1-pp2))))
points(1-pp2,pp1,pch="+",col=2)
for (i in 1:N) {
   lines(1-c(pp2[i],pp2[i]),c(lowxx1[i],higxx1[i]),lty=2,lwd=0.75,col=2)
   lines(1-c(lowxx2[i],higxx2[i]),c(pp1[i],pp1[i]),lty=2,lwd=0.75,col=2)
}
pp1s=(TPsB+0.5)/(nBcases+1)
lowxx1s=invlogit(logit(pp1s) - 1.96*sqrt(1/(nBcases*pp1s*(1-pp1s))))
higxx1s=invlogit(logit(pp1s) + 1.96*sqrt(1/(nBcases*pp1s*(1-pp1s))))
pp2s=(TNsB+0.5)/(nBcontrols+1)
lowxx2s=invlogit(logit(pp2s) - 1.96*sqrt(1/(nBcontrols*pp2s*(1-pp2s))))
higxx2s=invlogit(logit(pp2s) + 1.96*sqrt(1/(nBcontrols*pp2s*(1-pp2s))))
points(1-pp2s,pp1s,pch="+",col=3)
for (i in 1:N) {
   lines(1-c(pp2s[i],pp2s[i]),c(lowxx1s[i],higxx1s[i]),lty=2,lwd=0.75,col=3)
   lines(1-c(lowxx2s[i],higxx2s[i]),c(pp1s[i],pp1s[i]),lty=2,lwd=0.75,col=3)
}
legend(x=0.8,y=0.1,legend=c("test A","test B"),pch=16,col=c(2,3),bty="n")
stat=summary(x)$statistics
help0=summary(x)$quantiles[,c(1,3,5)]
points(1-invlogit(help0[which(row.names(help0)=="meanalpha[3]"),2]),invlogit(help0[which(row.names(help0)=="meanalpha[1]"),2]),col=3,pch=16,cex=1.5)
mu1=help0[which(row.names(help0)=="meanalpha[3]"),2]
mu2=help0[which(row.names(help0)=="meanalpha[1]"),2]
sematrix=matrix(c(stat[which(row.names(stat)=="meanalpha[3]"),2]^2,-cov(x[[1]][,1],x[[1]][,3]),-cov(x[[1]][,1],x[[1]][,3]),stat[which(row.names(stat)=="meanalpha[1]"),2]^2),ncol=2,nrow=2)
punten1=confellips(mu=c(-mu1,mu2),sigma=sematrix,alfa=0.05,npoints=1000)
lines(invlogit(punten1[,1]),invlogit(punten1[,2]),col=3)
points(1-invlogit(help0[which(row.names(help0)=="meanalpha[3]"),2]+help0[which(row.names(help0)=="meanalpha[4]"),2]),
         invlogit(help0[which(row.names(help0)=="meanalpha[1]"),2]+help0[which(row.names(help0)=="meanalpha[2]"),2]),col=2,pch=16,cex=1.5)
mu1=help0[which(row.names(help0)=="meanalpha[3]"),2] + help0[which(row.names(help0)=="meanalpha[4]"),2]
mu2=help0[which(row.names(help0)=="meanalpha[1]"),2] + help0[which(row.names(help0)=="meanalpha[2]"),2]
sematrix=matrix(c(var(x[[1]][,1]+x[[1]][,2]),-cov(x[[1]][,1]+x[[1]][,2],x[[1]][,3]+x[[1]][,4]),-cov(x[[1]][,1],x[[1]][,3]),var(x[[1]][,3]+x[[1]][,4])),ncol=2,nrow=2)
punten1=confellips(mu=c(-mu1,mu2),sigma=sematrix,alfa=0.05,npoints=1000)
lines(invlogit(punten1[,1]),invlogit(punten1[,2]),col=2)
#dev.off()

