
########## RANDOM-EFFECTS VARIANT


###
### meta-analysis of the comparative diagnostic accuracies of two diagnostic tests. 
###
### In this version we meta-analyze the results of a set of independent studies that all compared the same two
###     diagnostic tests (A and B) to a reference standard and that reported the results with the numbers of
###     true positives (TP), true negatives (TN), false positives (FP) and false negatives (FN) for both diagnostic tests. 
###     In this version we do not (cannot) use the correlations of the outcomes of both tests in the same patients, because such 
###     information is unfortunately rarely reported in scientific reports of this kind of diagnostic accuracy studies. 
###     So in effect the outcomes of tests A and B are handled as if they are independent and the analysis is done as if
###     the patients tested with test A are different than the patients tested with test B.
###
### Input data is in the form of a data.frame where every study is represemted by a separate row and there are 8 columns with the
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
### The meta-analysis model is based on binomial distributions of the TP and TN of tests A and B given the study-specific
###     sensitivities and specificities of tests A and B. Logit-transforms of these four study-specific parameters are assumed 
###     to have been samped from a mutivariate normal distribution which mean and covariance-matrix is estimated.
###     So the model is a random-effects model.
###
### The current code depends on the JAGS-program (Just Another Gibbs Sampler), which needs to be installed (mcmc-jags.sourceforge.net),
###     and on the R-package rjags, which needs to be attached. There are no other dependencies.
###
### Results of the analysis are estimates of the mean and covariance and 95% credibility intervals. In addition transformations
###     of these estimates are calculated, in particular estimates of odds ratio's of the mean sensitivity of test A vs test B and
###     of the mean specificity of test A vs test B. Also, the difference of the two mean sensitivities and of the two specificities
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

### define non-informative hyper priors for 
# I. the mean of the logit-transformed parameters (multivariate normal)
zero4=c(0,0,0,0)   
prec4=diag(0.001,nrow=4,ncol=4)  
# and II. the inverse of the covariance matrix of the study-specific parameters
InvTau=diag(0.1,nrow=4,ncol=4)

##### put data and hyper priors in a list-object to be transferred to JAGS
jagsdatalist=list(
  N=N, TPsA=TPsA,TPsB=TPsB,TNsA=TNsA,TNsB=TNsB, 
  nAcases=nAcases,nBcases=nBcases,nAcontrols=nAcontrols,nBcontrols=nBcontrols,
  zero4=zero4,prec4=prec4,InvTau=InvTau)

### model in JAGS code
modelstring = "
   model {
      for (i in 1:N){ 
         alpha[i,1:4] ~ dmnorm(meanalpha[1:4],InvSigma[1:4,1:4])
         sensB[i] <- 1/(1+exp(-1*(alpha[i,1])))
         sensA[i] <- 1/(1+exp(-1*(alpha[i,1]+alpha[i,2])))
         specB[i] <- 1/(1+exp(-1*(alpha[i,3])))
         specA[i] <- 1/(1+exp(-1*(alpha[i,3]+alpha[i,4])))
         sensORAvsB[i] <- (sensA[i]/(1-sensA[i])) / (sensB[i]/(1-sensB[i]))
         deltasensBvsA[i] <- sensB[i] - sensA[i]
         specORAvsB[i] <- (specA[i]/(1-specA[i])) / (specB[i]/(1-specB[i]))
         deltaspecBvsA[i] <- specB[i] - specA[i]
      }

      for (i in 1:N){ 
         TPsA[i] ~ dbin(sensA[i],nAcases[i])
         TPsB[i] ~ dbin(sensB[i],nBcases[i])
         TNsA[i] ~ dbin(specA[i],nAcontrols[i])
         TNsB[i] ~ dbin(specB[i],nBcontrols[i])
       }

      meanalpha[1:4] ~ dmnorm(zero4[1:4],prec4[1:4,1:4])
      InvSigma[1:4,1:4] ~ dwish(InvTau[1:4,1:4],4)

      meansensitivityB <- 1/(1+exp(-1*(meanalpha[1])))
      meansensitivityA <- 1/(1+exp(-1*(meanalpha[1]+meanalpha[2])))
      meanspecificityB <- 1/(1+exp(-1*(meanalpha[3])))
      meanspecificityA <- 1/(1+exp(-1*(meanalpha[3]+meanalpha[4])))

      sensdifferenceBvsA <- 1/(1+exp(-1*(meanalpha[1])))  -  1/(1+exp(-1*(meanalpha[1]+meanalpha[2])))
      specdifferenceBvsA <- 1/(1+exp(-1*(meanalpha[3])))  -  1/(1+exp(-1*(meanalpha[3]+meanalpha[4])))
      sensoddsratioAvsB <- exp(meanalpha[2])
      specoddsratioAvsB <- exp(meanalpha[4])
            
      Sigma[1:4,1:4] <- inverse(InvSigma[1:4,1:4])
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
inits = list(meanalpha=c(0,0,0,0), InvSigma=diag(1,nrow=4,ncol=4))                                     # starting values
m = jags.model(textConnection(modelstring), data=jagsdatalist, inits, n.chains=2)                      # compile the analysis 
                                                                                                       #    (i.e. check whether model and data fit together)
update(m, 50000)                                                                                       # do 50000 iterations 
x = coda.samples(m, c("meanalpha","Sigma","sensdifferenceBvsA","specdifferenceBvsA",
       "sensoddsratioAvsB","specoddsratioAvsB"), n.iter=50000)                                         # draw 50000 parameters from the posterior-distributions
#plot(x,ask=TRUE)                                                                                       # check for convergende 
                                                                                                       #     trace-plots should display stable chaos
# summary of the structural and derived parameters
summary(x)[[1]]                                                                                        # calculate statistics

dev.new()
par(mfrow=c(2,2))
hist(x[[1]][,which(colnames(x[[1]])=="sensoddsratioAvsB")],xlab="odds ratio A vs B",main="sensitivity")   # a histogram of the samples for one parameter
quantile(x[[1]][,which(colnames(x[[1]])=="sensoddsratioAvsB")],probs=c(0.025,0.25,0.50,0.75,0.975))       # a 95% credibility interval for this parameter
hist(x[[1]][,which(colnames(x[[1]])=="specoddsratioAvsB")],xlab="odds ratio A vs B",main="specificity")   # a histogram of the samples for a second parameter
quantile(x[[1]][,which(colnames(x[[1]])=="specoddsratioAvsB")],probs=c(0.025,0.25,0.50,0.75,0.975))       # a 95% credibility interval for this parameter
#par(mfrow=c(1,1))

#par(mfrow=c(1,2))
hist(x[[1]][,which(colnames(x[[1]])=="sensdifferenceBvsA")],xlab="difference B vs A",main="sensitivity")  # a histogram of the samples for another parameter
quantile(x[[1]][,which(colnames(x[[1]])=="sensdifferenceBvsA")],probs=c(0.025,0.25,0.50,0.75,0.975))      # a 95% credibility interval for this other parameter
hist(x[[1]][,which(colnames(x[[1]])=="specdifferenceBvsA")],xlab="difference B vs A",main="specificity")  # a histogram of the samples for yet another parameter
quantile(x[[1]][,which(colnames(x[[1]])=="specdifferenceBvsA")],probs=c(0.025,0.25,0.50,0.75,0.975))      # a 95% credibility interval for this other parameter
par(mfrow=c(1,1))


x1 = coda.samples(m,c("meansensitivityA","meansensitivityB","meanspecificityA","meanspecificityB"), n.iter=50000)
help1=t(apply(x1[[1]],2,function(y){quantile(y,probs=c(0.025,0.5,0.975))}))
x2 = coda.samples(m, c("sensA","sensB","specA","specB"), n.iter=50000) 
help2=t(apply(x2[[1]],2,function(y){quantile(y,probs=c(0.025,0.5,0.975))}))

dev.new()
#pdf("fig1.pdf")
layout(matrix(c(1,2,3,4),nrow=2,ncol=2))
plot(1,1,xlab="sensitivity",ylab="study number",main="",xlim=c(0,1),ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[1,2],(-N-2),pch=16,cex=1)
lines(c(help1[1,1],help1[1,3]),c(-N-2,-N-2))
points(help2[1:N,2],(-N:-1),pch=16,cex=0.75)
points(TPsA/(nAcases),(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   lines( c(help2[i,1],help2[i,3]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[1,2],lty=3,lwd=0.75)
plot(1,1,xlab="sensitivity",ylab="study number",main="",xlim=c(0,1),ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[2,2],(-N-2),pch=16,cex=1)
lines(c(help1[2,1],help1[2,3]),c(-N-2,-N-2))
points(help2[(N+1):(2*N),2],(-N:-1),pch=16,cex=0.75)
points(TPsB/(nBcases),(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   lines( c(help2[N+i,1],help2[N+i,3]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[2,2],lty=3,lwd=0.75)
plot(1,1,xlab="specificity",ylab="study number",main="",xlim=c(0,1),ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[3,2],(-N-2),pch=16,cex=1)
lines(c(help1[3,1],help1[3,3]),c(-N-2,-N-2))
points(help2[(2*N+1):(3*N),2],(-N:-1),pch=16,cex=0.75)
points(TNsA/(nAcontrols),(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   lines( c(help2[2*N+i,1],help2[2*N+i,3]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[3,2],lty=3,lwd=0.75)
plot(1,1,xlab="specificity",ylab="study number",main="",xlim=c(0,1),ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[4,2],(-N-2),pch=16,cex=1)
lines(c(help1[4,1],help1[4,3]),c(-N-2,-N-2))
points(help2[(3*N+1):(4*N),2],(-N:-1),pch=16,cex=0.75)
points(TNsB/(nBcontrols),(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   lines( c(help2[3*N+i,1],help2[3*N+i,3]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[4,2],lty=3,lwd=0.75)
mtext("test A",side=3,outer=TRUE,line=-2)
mtext("test B",side=3,outer=TRUE,line=-23)
#dev.off()

dev.new()
#pdf("fig2.pdf")
x3 = coda.samples(m,c("sensORAvsB","specORAvsB","deltasensBvsA","deltaspecBvsA"), n.iter=50000)
help3=t(apply(x3[[1]],2,function(y){quantile(y,probs=c(0.025,0.5,0.975))}))
help0=t(apply(x[[1]],2,function(y){quantile(y,probs=c(0.025,0.5,0.975))}))
layout(matrix(c(1,2,3,4),nrow=2,ncol=2))
plot(1,1,xlab="difference of sensitivities B vs A",ylab="study number",xlim=c(-1,1),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help0[which(row.names(help0)=="sensdifferenceBvsA"),2],(-N-2),pch=16,cex=1)
lines(c(help0[which(row.names(help0)=="sensdifferenceBvsA"),1],help0[which(row.names(help0)=="sensdifferenceBvsA"),3]),c(-N-2,-N-2))
points((TPsB/nBcases - TPsA/nAcases),(-N:-1),pch="+",cex=0.75)
points(help3[1:N,2],(-N:-1),pch=16,cex=0.75)
for (i in 1:N) {
   lines( c(help3[i,1],help3[i,3]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help0[which(row.names(help0)=="sensdifferenceBvsA"),2],lty=3,lwd=0.75)
plot(1,1,xlab="difference of specificities B vs A",ylab="study number",xlim=c(-1,1),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help0[which(row.names(help0)=="specdifferenceBvsA"),2],(-N-2),pch=16,cex=1)
lines(c(help0[which(row.names(help0)=="specdifferenceBvsA"),1],help0[which(row.names(help0)=="specdifferenceBvsA"),3]),c(-N-2,-N-2))
points((TNsB/nBcontrols - TNsA/nAcontrols),(-N:-1),pch="+",cex=0.75)
points(help3[(N+1):(2*N),2],(-N:-1),pch=16,cex=0.75)
for (i in 1:N) {
   lines( c(help3[N+i,1],help3[N+i,3]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help0[which(row.names(help0)=="specdifferenceBvsA"),2],lty=3,lwd=0.75)
minx=0 #min(help3[(2*N+1):(3*N),1],na.rm=TRUE)
maxx=1.1 #max(help3[(2*N+1):(3*N),3],na.rm=TRUE)
plot(1,1,xlab="OR of sensitivities A vs B",ylab="study number",xlim=c(minx,maxx),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help0[which(row.names(help0)=="sensoddsratioAvsB"),2],(-N-2),pch=16,cex=1)
lines(c(help0[which(row.names(help0)=="sensoddsratioAvsB"),1],help0[which(row.names(help0)=="sensoddsratioAvsB"),3]),c(-N-2,-N-2))
points(exp(logit(TPsA/nAcases) - logit(TPsB/nBcases)),(-N:-1),pch="+",cex=0.75)
points(help3[(2*N+1):(3*N),2],(-N:-1),pch=16,cex=0.75)
for (i in 1:N) {
   lines( c(help3[2*N+i,1],help3[2*N+i,3]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help0[which(row.names(help0)=="sensoddsratioAvsB"),2],lty=3,lwd=0.75)
minxx=0 #min(help3[(3*N+1):(4*N),1],na.rm=TRUE)
maxxx=10 #max(help3[(3*N+1):(4*N),3],na.rm=TRUE)
plot(1,1,xlab="OR of specificities A vs B",ylab="study number",xlim=c(minxx,maxxx),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help0[which(row.names(help0)=="specoddsratioAvsB"),2],(-N-2),pch=16,cex=1)
lines(c(help0[which(row.names(help0)=="specoddsratioAvsB"),1],help0[which(row.names(help0)=="specoddsratioAvsB"),3]),c(-N-2,-N-2))
points(exp(logit((TNsA+0.5)/(nAcontrols+1)) - logit((TNsB+0.5)/(nBcontrols+1))),(-N:-1),pch="+",cex=0.75)
points(help3[(3*N+1):(4*N),2],(-N:-1),pch=16,cex=0.75)
for (i in 1:N) {
   lines( c(help3[3*N+i,1],help3[3*N+i,3]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help0[which(row.names(help0)=="specoddsratioAvsB"),2],lty=3,lwd=0.75)
#dev.off()


dev.new()
#pdf("fig3.pdf")
plot(1,1,xlab="1-specificity",ylab="sensitivity",xlim=c(0,1),ylim=c(0,1),type="n",main="ROC curves for both tests")
points(1-help2[(2*N+1):(3*N),2],help2[1:N,2],pch=16,col=2)
for (i in 1:N) {
   lines(1-c(help2[2*N+i,2],help2[2*N+i,2]),c(help2[i,1],help2[i,3]),lty=2,lwd=0.75)
   lines(1-c(help2[2*N+i,1],help2[2*N+i,3]),c(help2[i,2],help2[i,2]),lty=2,lwd=0.75)
}
points(1-help2[(3*N+1):(4*N),2],help2[(N+1):(2*N),2],pch=16,col=3)
for (i in 1:N) {
   lines(1-c(help2[3*N+i,2],help2[3*N+i,2]),c(help2[N+i,1],help2[N+i,3]),lty=2,lwd=0.75)
   lines(1-c(help2[3*N+i,1],help2[3*N+i,3]),c(help2[N+i,2],help2[N+i,2]),lty=2,lwd=0.75)
}
legend(x=0.8,y=0.1,legend=c("test A","test B"),pch=16,col=c(2,3),bty="n")
help0=t(apply(x[[1]],2,function(y){quantile(y,probs=c(0.025,0.5,0.975))}))
b2=help0[which(row.names(help0)=="Sigma[1,3]"),2]/help0[which(row.names(help0)=="Sigma[3,3]"),2]
a2= help0[which(row.names(help0)=="meanalpha[1]"),2] - b2 * help0[which(row.names(help0)=="meanalpha[3]"),2]
b1=(help0[which(row.names(help0)=="Sigma[1,3]"),2]+help0[which(row.names(help0)=="Sigma[1,4]"),2]+help0[which(row.names(help0)=="Sigma[2,3]"),2]+help0[which(row.names(help0)=="Sigma[2,4]"),2])/(help0[which(row.names(help0)=="Sigma[3,3]"),2]+help0[which(row.names(help0)=="Sigma[4,4]"),2]+2*help0[which(row.names(help0)=="Sigma[3,4]"),2])
a1=(help0[which(row.names(help0)=="meanalpha[1]"),2] + help0[which(row.names(help0)=="meanalpha[2]"),2]) - b1*(help0[which(row.names(help0)=="meanalpha[3]"),2] + help0[which(row.names(help0)=="meanalpha[4]"),2])
yyy2=c()
yyy1=c()
lowestyyy2=min(1-help2[(3*N+1):(4*N),3],na.rm=TRUE)
highestyyy2=max(1-help2[(3*N+1):(4*N),1],na.rm=TRUE)
lowestyyy1=min(1-help2[(2*N+1):(3*N),3],na.rm=TRUE)
highestyyy1=max(1-help2[(2*N+1):(3*N),1],na.rm=TRUE)
xxx=(1:99)/100
yyy2=exp(a2-b2*log(xxx/(1-xxx)))/(1+exp(a2-b2*log(xxx/(1-xxx))))
yyy1=exp(a1-b1*log(xxx/(1-xxx)))/(1+exp(a1-b1*log(xxx/(1-xxx))))
lines(xxx[xxx>=lowestyyy2 & xxx<=highestyyy2],yyy2[xxx>=lowestyyy2 & xxx<=highestyyy2],col=3)
lines(xxx[xxx>=lowestyyy1 & xxx<=highestyyy1],yyy1[xxx>=lowestyyy1 & xxx<=highestyyy1],col=2)
#dev.off()


