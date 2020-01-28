
########## RANDOM-EFFECTS VARIANT


###
### meta-analysis of the comparative diagnostic accuracies of two diagnostic tests. 
###
### In this version we meta-analyze the results of a set of independent studies that all compared the same two
###     diagnostic tests (A and B) to a reference standard. 
###     Some studies only reported results with the numbers of true positives (TP), true negatives (TN), 
###     false positives (FP) and false negatives (FN) for both diagnostic tests. These are type-1 studies.
###     Other studies reported the complete 2-by-2 table of the outcomes of tests A and B separately for cases 
###     and controls. These are type-2 studies. Using these data we can also estimate the correlations of the 
###     outcomes of both tests in the same patients, separately for cases and controls. 
###     The outcomes of tests A and B in the studies that only reported TP, FN, TN and FP results are handled as if 
###     they are independent and the analysis is done as if the patients tested with test A are different than the 
###     patients tested with test B.
###
### Input data for type-1 studies is in the form of a data.frame where every study is a separate row and there 
###     are 8 columns with the TP, FN, TN, FP for test A and the TP, FN, TN, FP for test B. An example is given below.
###     Input for type-2 studies is a second data.frame where every study is a row and there are also 8 colums with
###     the numbers of case-patients that are positive on both tests (Capp), that are positive for A and negative for B
###     Capm, that are negative for A but positive for B (Camp) and that are negative for both (Camm). Similar numbers
##3     are given for the control patients (Copp, Copm, Comp, Comm). An example is given below.
###

d1 = as.data.frame(matrix(
     c(14, 26,  90, 10, 28, 12,  68, 32,
       18, 37,  50,  0, 40, 15,  26, 24,
        2,  4, 185,  8,  4,  2, 159, 34,
       24, 30,  55,  4, 42, 15, 101, 38,
       38, 33,  11,  3, 64,  7,  13,  4),
       nrow=5,ncol=8,byrow=TRUE)
    )
names(d1)=c("TPsA","FNsA","TNsA","FPsA","TPsB","FNsB","TNsB","FPsB")

d2 = as.data.frame(matrix(
     c( 3,0,0,0,1,3,15,78,
       51,5,1,6,2,1, 5,56,
       22,5,1,1,1,1, 3,25),
      nrow=3,ncol=8,byrow=TRUE)
   )
names(d2)=c("Capp","Capm","Camp","Camm","Copp","Copm","Comp","Comm")

###
### For type-1 studies the meta-analysis model is based on binomial distributions of the TP and TN of tests A and B 
###     given the overall -non varying- sensitivities and specificities of tests A and B. 
### For type-2 studies the meta-analysis model is also based on the binomial distributions of the TP and TN of tests A and B
###     whereas the numbers of case-patients that are positive on both tests and the numbers of patients that are negative
###     on both tests are based on the hypergeometric distributions given the marginals of the 2-by-2 tables. Logit-transforms 
###     of the six study-specific parameters are assumed to have been samped from a mutivariate normal distribution which 
###     mean and covariance-matrix is estimated. So the model is a random-effects model.
###
### The current code depends on the JAGS-program (Just Another Gibbs Sampler), which needs to be installed (mcmc-jags.sourceforge.net),
###     and on the R-package rjags, which needs to be attached. There are no other dependencies.
###
### Results of the analysis are estimates of the logit-transformed sensitivities and specificities and 95% credibility intervals,
###     and log-odds ratios of tests A vs B for cases and controls. In total there are therefore six parameters.
###     In addition transformations of these estimates are calculated, in particular estimates of odds ratio's of the sensitivity of 
###     test A vs test B and of the specificity of test A vs test B. Also, the difference of the two mean sensitivities and of the two specificities
###     are calculated together with their 95% credibility intervals.
###


###  some settings
delta=0.5                 # add 0.5 to all TN, FP, TP, FN?
iternr = 50000            # maximum number of iterations
xas=seq(0.01,0.99,0.01)   # steps on the x-axis for calculation of roc-curves and AUCs
P=length(xas)


### translate d2 into d1, combine into d
d2$TPsA=d2$Capp+d2$Capm
d2$FNsA=d2$Camp+d2$Camm
d2$TNsA=d2$Comm+d2$Comp
d2$FPsA=d2$Copm+d2$Copp
d2$TPsB=d2$Capp+d2$Camp
d2$FNsB=d2$Capm+d2$Camm
d2$TNsB=d2$Comm+d2$Copm
d2$FPsB=d2$Comp+d2$Copp
d=rbind(d1,d2[,names(d1)])


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
N2=nrow(d2)

### define non-informative hyper priors for 
# the logit-transformed parameters and the log-ORs (multivariate normal)
zero6=c(0,0,0,0,0,0)   
prec6=diag(0.001,nrow=6,ncol=6)  
# and the inverse of the variance-covariance matrix of the study-specific parameters
InvTau=diag(0.1,nrow=6,ncol=6)

##### put data and hyper priors in a list-object to be transferred to JAGS
jagsdatalist=list(
  N=N, TPsA=TPsA,TPsB=TPsB,TNsA=TNsA,TNsB=TNsB, 
  nAcases=nAcases,nBcases=nBcases,nAcontrols=nAcontrols,nBcontrols=nBcontrols,
  N2=N2,
  Capp=d2$Capp,CapA=d2$Capp+d2$Capm,CapB=d2$Capp+d2$Camp,Ncases   =d2$Capp+d2$Capm+d2$Camp+d2$Camm,
  Comm=d2$Comm,ComA=d2$Comm+d2$Comp,ComB=d2$Comm+d2$Copm,Ncontrols=d2$Comm+d2$Comp+d2$Copm+d2$Copp,
  zero6=zero6,prec6=prec6,InvTau=InvTau,x=xas,P=P)

### model in JAGS code
modelstring = "
   model {
      # sample parameters from the population distribution and transform
      for (i in 1:N) {
         alpha[i,1:6] ~ dmnorm(meanalpha[1:6],InvSigma[1:6,1:6])
         sensB[i] <- 1/(1+exp(-1*(alpha[i,1])))
         sensA[i] <- 1/(1+exp(-1*(alpha[i,1]+alpha[i,2])))
         specB[i] <- 1/(1+exp(-1*(alpha[i,3])))
         specA[i] <- 1/(1+exp(-1*(alpha[i,3]+alpha[i,4])))
         ORcases[i]    <- exp(alpha[i,5])
         ORcontrols[i] <- exp(alpha[i,6])
      }
      # likelihoods of the observed data: marginals
      # for all type-1 and type-2 studies
      for (i in 1:N){ 
         TPsA[i] ~ dbin(sensA[i],nAcases[i])
         TPsB[i] ~ dbin(sensB[i],nBcases[i])
         TNsA[i] ~ dbin(specA[i],nAcontrols[i])
         TNsB[i] ~ dbin(specB[i],nBcontrols[i])
      }
      # likelihoods of the 2-by-2 tables for cases and controls
      # only for type-2 studies
      for ( i in 1:N2) {
         Capp[i] ~ dhyper(CapA[i],(Ncases[i]-CapA[i]),   CapB[i],ORcases[i])
         Comm[i] ~ dhyper(ComA[i],(Ncontrols[i]-ComA[i]),ComB[i],ORcontrols[i])
      }
      # noninformative priors
      meanalpha[1:6] ~ dmnorm(zero6[1:6],prec6[1:6,1:6])
      InvSigma[1:6,1:6] ~ dwish(InvTau[1:6,1:6],6)

      # transformations of population parameters
      Sigma[1:6,1:6] <- inverse(InvSigma[1:6,1:6])

      meansensitivityB <- 1/(1+exp(-1*(meanalpha[1])))
      meansensitivityA <- 1/(1+exp(-1*(meanalpha[1]+meanalpha[2])))
      meanspecificityB <- 1/(1+exp(-1*(meanalpha[3])))
      meanspecificityA <- 1/(1+exp(-1*(meanalpha[3]+meanalpha[4])))

      sensdifferenceBvsA <- 1/(1+exp(-1*(meanalpha[1])))  -  1/(1+exp(-1*(meanalpha[1]+meanalpha[2])))
      specdifferenceBvsA <- 1/(1+exp(-1*(meanalpha[3])))  -  1/(1+exp(-1*(meanalpha[3]+meanalpha[4])))
      sensoddsratioAvsB <- exp(meanalpha[2])
      specoddsratioAvsB <- exp(meanalpha[4])
      
      # parameters of the ROC-curve: named naive by mada
      slopeA     <- -1*(Sigma[1,3]+Sigma[1,4]+Sigma[2,3]+Sigma[2,4])/(Sigma[3,3]+Sigma[4,4]+2*Sigma[3,4])
      slopeB     <- -1*Sigma[3,1]/Sigma[3,3]
      interceptB <-  meanalpha[1] + slopeB*meanalpha[3]
      interceptA <- (meanalpha[1]+meanalpha[2]) + slopeA*(meanalpha[3]+meanalpha[4])

      # Rutter and Gatsonis parameters: Harbord et al. Biostatistics, 2007, 8, 2, 239-251.
      beta_A     <- 0.5*log((Sigma[3,3]+Sigma[4,4]+2*Sigma[3,4])/(Sigma[1,1]+Sigma[2,2]+2*Sigma[1,2]))
      LAMBDA_A   <- sqrt(sqrt(Sigma[3,3]+Sigma[4,4]+2*Sigma[3,4])/sqrt(Sigma[1,1]+Sigma[2,2]+2*Sigma[1,2])) * (meanalpha[1]+meanalpha[2]) +
                    sqrt(sqrt(Sigma[1,1]+Sigma[2,2]+2*Sigma[1,2])/sqrt(Sigma[3,3]+Sigma[4,4]+2*Sigma[3,4])) * (meanalpha[3]+meanalpha[4])
      THETA_A    <- 0.5*(sqrt(sqrt(Sigma[3,3]+Sigma[4,4]+2*Sigma[3,4])/sqrt(Sigma[1,1]+Sigma[2,2]+2*Sigma[1,2])) * (meanalpha[1]+meanalpha[2]) -
                         sqrt(sqrt(Sigma[1,1]+Sigma[2,2]+2*Sigma[1,2])/sqrt(Sigma[3,3]+Sigma[4,4]+2*Sigma[3,4])) * (meanalpha[3]+meanalpha[4]))
      skwtheta_A <- 0.5*(sqrt((Sigma[3,3]+Sigma[4,4]+2*Sigma[3,4])*(Sigma[1,1]+Sigma[2,2]+2*Sigma[1,2]))-(Sigma[1,3]+Sigma[1,4]+Sigma[2,3]+Sigma[2,4]))
      skwalfa_A  <- 2*(sqrt((Sigma[3,3]+Sigma[4,4]+2*Sigma[3,4])*(Sigma[1,1]+Sigma[2,2]+2*Sigma[1,2]))+(Sigma[1,3]+Sigma[1,4]+Sigma[2,3]+Sigma[2,4]))

      beta_B     <- 0.5*log(Sigma[3,3]/Sigma[1,1])
      LAMBDA_B   <- sqrt(sqrt(Sigma[3,3])/sqrt(Sigma[1,1]))*meanalpha[1] + sqrt(sqrt(Sigma[1,1])/sqrt(Sigma[3,3]))*meanalpha[3]
      THETA_B    <- 0.5*(sqrt(sqrt(Sigma[3,3])/sqrt(Sigma[1,1]))*meanalpha[1] - sqrt(sqrt(Sigma[1,1])/sqrt(Sigma[3,3]))*meanalpha[3])
      skwtheta_B <- 0.5*(sqrt(Sigma[1,1]*Sigma[3,3])-Sigma[1,3])
      skwalfa_B  <- 2.0*(sqrt(Sigma[1,1]*Sigma[3,3])+Sigma[1,3])

      # AUCs
      y1[1]   <- 1/(1+exp(-(interceptA + slopeA * log(x[1]/(1-x[1])))))
      y2[1]   <- 1/(1+exp(-(interceptB + slopeB * log(x[1]/(1-x[1])))))
      y3A[1] <- 1/(1+exp(-(LAMBDA_A*exp(beta_A/2) + exp(beta_A)*log(x[1]/(1-x[1])))))
      y3B[1] <- 1/(1+exp(-(LAMBDA_B*exp(beta_B/2) + exp(beta_B)*log(x[1]/(1-x[1])))))
      opp1[1] <- (x[1]-0)*(y1[1]-0)/2
      opp2[1] <- (x[1]-0)*(y2[1]-0)/2
      opp3A[1] <- (x[1]-0)*(y3A[1]-0)/2
      opp3B[1] <- (x[1]-0)*(y3B[1]-0)/2
      for (j in 2:P) {
         y1[j]   <- 1/(1+exp(-(interceptA + slopeA * log(x[j]/(1-x[j])))))
         opp1[j] <- (x[j]-x[(j-1)])*y1[(j-1)]+(x[j]-x[(j-1)])*(y1[j]-y1[(j-1)])/2
         y2[j]   <- 1/(1+exp(-(interceptB + slopeB * log(x[j]/(1-x[j])))))
         opp2[j] <- (x[j]-x[(j-1)])*y2[(j-1)]+(x[j]-x[(j-1)])*(y2[j]-y2[(j-1)])/2
         y3A[j] <- 1/(1+exp(-(LAMBDA_A*exp(beta_A/2) + exp(beta_A)*log(x[j]/(1-x[j])))))
         opp3A[j] <- (x[j]-x[(j-1)])*y3A[(j-1)]+(x[j]-x[(j-1)])*(y3A[j]-y3A[(j-1)])/2
         y3B[j] <- 1/(1+exp(-(LAMBDA_B*exp(beta_B/2) + exp(beta_B)*log(x[j]/(1-x[j])))))
         opp3B[j] <- (x[j]-x[(j-1)])*y3B[(j-1)]+(x[j]-x[(j-1)])*(y3B[j]-y3B[(j-1)])/2
      }
      opp1[(P+1)] <- (1-x[P])*y1[P]+(1-x[P])*(1-y1[P])/2      
      aucA <- sum(opp1)
      opp2[(P+1)] <- (1-x[P])*y2[P]+(1-x[P])*(1-y2[P])/2      
      aucB <- sum(opp2)
      delta_AUC <- aucA-aucB
      opp3A[(P+1)] <- (1-x[P])*y3A[P]+(1-x[P])*(1-y3A[P])/2      
      aucRG_A <- sum(opp3A) 
      opp3B[(P+1)] <- (1-x[P])*y3B[P]+(1-x[P])*(1-y3B[P])/2      
      aucRG_B <- sum(opp3B) 
      delta_AUC_RG <- aucRG_A - aucRG_B
   }  
"

### some useful functions
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
inits = list(meanalpha=c(0,0,0,0,0,0))                                                                 # starting values
m = jags.model(textConnection(modelstring), data=jagsdatalist, inits, n.chains=2)                      # compile the analysis 
                                                                                                       #    (i.e. check whether model and data fit together)
update(m, iternr)                                                                                      # do iternr iterations 
x = coda.samples(m, c("meanalpha","Sigma",
       "meansensitivityA","meansensitivityB","meanspecificityA","meanspecificityB",
       "sensdifferenceBvsA","specdifferenceBvsA","sensoddsratioAvsB","specoddsratioAvsB",
        "slopeA","slopeB","interceptA","interceptB","aucA","aucB","delta_AUC",
        "beta_A","LAMBDA_A","THETA_A","skwtheta_A","skwalfa_A",
        "beta_B","LAMBDA_B","THETA_B","skwtheta_B","skwalfa_B",
        "aucRG_B","aucRG_A","delta_AUC_RG"), n.iter=iternr)                                             # draw iternr parameters from the posterior-distributions

#dev.new()
#pdf("traceplots.pdf")
#plot(x)   #,ask=TRUE)   
#dev.off()                                                                                               # check for convergence: trace-plots should display stable chaos



####### since the convergence was OK for the datasets tried out so far, the results below are based only on sample drawn from the first chain


# summary of the structural and derived parameters
summary(x)[[1]] 
help1 = t(apply(x[[1]],2,function(y){quantile(y,probs=c(0.025,0.5,0.975))}))                            # calculate some statistics

dev.new()
#pdf("fig0.pdf")
layout(matrix(c(1,3,2,4),nrow=2,ncol=2))
hist(log(x[[1]][,which(colnames(x[[1]])=="sensoddsratioAvsB")]),xlab="log odds ratio of sensitivity A vs B",main="")   # a histogram of the samples for one parameter
quantile(x[[1]][,which(colnames(x[[1]])=="sensoddsratioAvsB")],probs=c(0.025,0.25,0.50,0.75,0.975))                    # a 95% credibility interval for this parameter
hist(log(x[[1]][,which(colnames(x[[1]])=="specoddsratioAvsB")]),xlab="log odds ratio of specificity A vs B",main="")   # a histogram of the samples for one parameter
quantile(x[[1]][,which(colnames(x[[1]])=="specoddsratioAvsB")],probs=c(0.025,0.25,0.50,0.75,0.975))                    # a 95% credibility interval for this parameter
hist(x[[1]][,which(colnames(x[[1]])=="sensdifferenceBvsA")],xlab="difference of sensitivity B vs A",main="")      # a histogram of the samples for another parameter
quantile(x[[1]][,which(colnames(x[[1]])=="sensdifferenceBvsA")],probs=c(0.025,0.25,0.50,0.75,0.975))              # a 95% credibility interval for this other parameter
hist(x[[1]][,which(colnames(x[[1]])=="specdifferenceBvsA")],xlab="difference of specificity B vs A",main="")      # a histogram of the samples for another parameter
quantile(x[[1]][,which(colnames(x[[1]])=="specdifferenceBvsA")],probs=c(0.025,0.25,0.50,0.75,0.975))              # a 95% credibility interval for this other parameter
title("posterior distributions of population parameters",line=-2,outer=TRUE)   
#dev.off()

x2 = coda.samples(m,c("sensA","sensB","specA","specB","ORcases","ORcontrols","alpha"), n.iter=iternr)
help2=t(apply(x2[[1]],2,function(y){quantile(y,probs=c(0.025,0.5,0.975))}))

dev.new()
#pdf("fig1.pdf")
layout(matrix(c(1,2,3,4),nrow=2,ncol=2))
par(mar=c(5, 4, 4, 2) + 0.1)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(1,1,xlab="sensitivity of test A",ylab="study number",main="",xlim=c(0,1),ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[which(row.names(help1)=="meansensitivityA"),2],(-N-2),pch=16,cex=1)
lines(c(help1[which(row.names(help1)=="meansensitivityA"),1],help1[which(row.names(help1)=="meansensitivityA"),3]),c(-N-2,-N-2))
pp=(TPsA+delta)/(nAcases+2*delta)
points(pp,(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   points(help2[which(row.names(help2)==paste("sensA[",i,"]",sep="")),2],-N+(i-1),pch=16)
   lines(c(help2[which(row.names(help2)==paste("sensA[",i,"]",sep="")),1],help2[which(row.names(help2)==paste("sensA[",i,"]",sep="")),3]) ,c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[which(row.names(help1)=="meansensitivityA"),2],lty=3,lwd=0.75)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(1,1,xlab="sensitivity of test B",ylab="study number",main="",xlim=c(0,1),ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[which(row.names(help1)=="meansensitivityB"),2],(-N-2),pch=16,cex=1)
lines(c(help1[which(row.names(help1)=="meansensitivityB"),1],help1[which(row.names(help1)=="meansensitivityB"),3]),c(-N-2,-N-2))
pp=(TPsB+delta)/(nBcases+2*delta)
points(pp,(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   points(help2[which(row.names(help2)==paste("sensB[",i,"]",sep="")),2],-N+(i-1),pch=16)
   lines(c(help2[which(row.names(help2)==paste("sensB[",i,"]",sep="")),1],help2[which(row.names(help2)==paste("sensB[",i,"]",sep="")),3]) ,c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[which(row.names(help1)=="meansensitivityB"),2],lty=3,lwd=0.75)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(1,1,xlab="specificity of test A",ylab="study number",main="",xlim=c(0,1),ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[which(row.names(help1)=="meanspecificityA"),2],(-N-2),pch=16,cex=1)
lines(c(help1[which(row.names(help1)=="meanspecificityA"),1],help1[which(row.names(help1)=="meanspecificityA"),3]),c(-N-2,-N-2))
pp=(TNsA+delta)/(nAcontrols+2*delta)
points(pp,(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   points(help2[which(row.names(help2)==paste("specA[",i,"]",sep="")),2],-N+(i-1),pch=16)
   lines(c(help2[which(row.names(help2)==paste("specA[",i,"]",sep="")),1],help2[which(row.names(help2)==paste("specA[",i,"]",sep="")),3]) ,c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[which(row.names(help1)=="meanspecificityA"),2],lty=3,lwd=0.75)
par(mar=c(5, 4, 4, 2) + 0.1)
plot(1,1,xlab="specificity of test B",ylab="study number",main="",xlim=c(0,1),ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[which(row.names(help1)=="meanspecificityB"),2],(-N-2),pch=16,cex=1)
lines(c(help1[which(row.names(help1)=="meanspecificityB"),1],help1[which(row.names(help1)=="meanspecificityB"),3]),c(-N-2,-N-2))
pp=(TNsB+delta)/(nBcontrols+2*delta)
points(pp,(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   points(help2[which(row.names(help2)==paste("specB[",i,"]",sep="")),2],-N+(i-1),pch=16)
   lines(c(help2[which(row.names(help2)==paste("specB[",i,"]",sep="")),1],help2[which(row.names(help2)==paste("specB[",i,"]",sep="")),3]) ,c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[which(row.names(help1)=="meanspecificityB"),2],lty=3,lwd=0.75)
#mtext("test A",side=3,outer=TRUE,line=-2)
#mtext("test B",side=3,outer=TRUE,line=-23)
par(mar=c(5, 4, 4, 2) + 0.1)
title("posterior medians and 95% credibility intervals",line=-2,outer=TRUE)
mtext("(observed values are indicated by '+')",line=-3,outer=TRUE,cex=0.75)
#dev.off()

dev.new()
#pdf("fig2.pdf")
layout(matrix(c(1,2,3,4),nrow=2,ncol=2))
plot(1,1,xlab="difference of sensitivities B vs A",ylab="study number",xlim=c(-1,1),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[which(row.names(help1)=="sensdifferenceBvsA"),2],(-N-2),pch=16,cex=1)
lines(c(help1[which(row.names(help1)=="sensdifferenceBvsA"),1],help1[which(row.names(help1)=="sensdifferenceBvsA"),3]),c(-N-2,-N-2))
pp2=(TPsA+delta)/(nAcases+2*delta)
pp1=(TPsB+delta)/(nBcases+2*delta)
points(pp1-pp2,(-N:-1),pch="+",cex=1)
for (i in 1:N) {
   temp1 = quantile(as.numeric(x2[[1]][,which(colnames(x2[[1]])==paste("sensB[",i,"]",sep=""))]) - 
                    as.numeric(x2[[1]][,which(colnames(x2[[1]])==paste("sensA[",i,"]",sep=""))]),probs=c(0.025,0.50,0.975))
   points(temp1[2],-N+(i-1),pch=16)
   lines( c(temp1[1],temp1[3]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[which(row.names(help1)=="sensdifferenceBvsA"),2],lty=3,lwd=0.75)
plot(1,1,xlab="difference of specificities B vs A",ylab="study number",xlim=c(-1,1),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(help1[which(row.names(help1)=="specdifferenceBvsA"),2],(-N-2),pch=16,cex=1)
lines(c(help1[which(row.names(help1)=="specdifferenceBvsA"),1],help1[which(row.names(help1)=="specdifferenceBvsA"),3]),c(-N-2,-N-2))
pp2=(TNsA+delta)/(nAcontrols+2*delta)
pp1=(TNsB+delta)/(nBcontrols+2*delta)
points(pp1-pp2,(-N:-1),pch="+",cex=1)
for (i in 1:N) {
   temp1 = quantile(as.numeric(x2[[1]][,which(colnames(x2[[1]])==paste("specB[",i,"]",sep=""))]) - 
                    as.numeric(x2[[1]][,which(colnames(x2[[1]])==paste("specA[",i,"]",sep=""))]),probs=c(0.025,0.50,0.975))
   points(temp1[2],-N+(i-1),pch=16)
   lines( c(temp1[1],temp1[3]),c(-N+(i-1),-N+(i-1)))
}
abline(v=help1[which(row.names(help1)=="specdifferenceBvsA"),2],lty=3,lwd=0.75)
pp2=(TPsA+delta)/(nAcases+2*delta)
pp1=(TPsB+delta)/(nBcases+2*delta)
ORs=(pp2/(1-pp2)) / (pp1/(1-pp1))
SEs=sqrt( (1/(TPsA+delta))+(1/(d$FNsA+delta))+(1/(TPsB+delta))+(1/(d$FNsB+delta)) )
pp1x=pp1[(nrow(d1)+1):N]
pp2x=pp2[(nrow(d1)+1):N]
d2n=(d2$Capp+d2$Camp+d2$Capm+d2$Camm)
ses = sqrt((1/((d2n+2*delta)*pp1x*(1-pp1x))) + (1/((d2n+2*delta)*pp2x*(1-pp2x))) - 2*((d2$Capp+0.5*delta)*(d2$Camm+0.5*delta)-(d2$Capm+0.5*delta)*(d2$Camp+0.5*delta))/(((d2n+2*delta)^3)*pp1x*(1-pp1x)*pp2x*(1-pp2x)))
SEs[(nrow(d1)+1):N]=ses
lowxx=(log(ORs)-1.96*SEs)
higxx=(log(ORs)+1.96*SEs)
minx=min(lowxx,na.rm=TRUE)
maxx=max(higxx,na.rm=TRUE)
plot(1,1,xlab="log OR of sensitivities A vs B",ylab="study number",xlim=c(minx,maxx),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(log(help1[which(row.names(help1)=="sensoddsratioAvsB"),2]),(-N-2),pch=16,cex=1)
lines(log(c(help1[which(row.names(help1)=="sensoddsratioAvsB"),1],help1[which(row.names(help1)=="sensoddsratioAvsB"),3])),c(-N-2,-N-2))
points(log(ORs),(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   temp1= as.numeric(x2[[1]][,which(colnames(x2[[1]])==paste("sensA[",i,"]",sep=""))])
   temp2= as.numeric(x2[[1]][,which(colnames(x2[[1]])==paste("sensB[",i,"]",sep=""))])
   temp3 = log((temp1*(1-temp2))/((1-temp1)*temp2))
   temp4 = quantile(temp3,probs=c(0.025,0.50,0.975))
   points(temp4[2],-N+(i-1),pch=16)
   lines( c(temp4[1],temp4[3]),c(-N+(i-1),-N+(i-1)))
}
abline(v=log(help1[which(row.names(help1)=="sensoddsratioAvsB"),2]),lty=3,lwd=0.75)
pp2=(TNsA+delta)/(nAcontrols+2*delta)
pp1=(TNsB+delta)/(nBcontrols+2*delta)
ORs=(pp2/(1-pp2)) / (pp1/(1-pp1))
SEs=sqrt( (1/(TNsA+delta))+(1/(d$FPsA+delta))+(1/(TNsB+delta))+(1/(d$FPsB+delta)) )
pp1x=pp1[(nrow(d1)+1):N]
pp2x=pp2[(nrow(d1)+1):N]
d2n=(d2$Copp+d2$Comp+d2$Copm+d2$Comm)
ses = sqrt((1/((d2n+2*delta)*pp1x*(1-pp1x))) + (1/((d2n+2*delta)*pp2x*(1-pp2x))) - 2*((d2$Copp+0.5*delta)*(d2$Comm+0.5*delta)-(d2$Copm+0.5*delta)*(d2$Comp+0.5*delta))/(((d2n+2*delta)^3)*pp1x*(1-pp1x)*pp2x*(1-pp2x)))
SEs[(nrow(d1)+1):N]=ses
lowxx=(log(ORs)-1.96*SEs)
higxx=(log(ORs)+1.96*SEs)
minx=min(lowxx,na.rm=TRUE)
maxx=max(higxx,na.rm=TRUE)
plot(1,1,xlab="log OR of specificities A vs B",ylab="study number",xlim=c(minx,maxx),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2),-N:-1),labels=c("overall",1:N))
points(log(help1[which(row.names(help1)=="specoddsratioAvsB"),2]),(-N-2),pch=16,cex=1)
lines(log(c(help1[which(row.names(help1)=="specoddsratioAvsB"),1],help1[which(row.names(help1)=="specoddsratioAvsB"),3])),c(-N-2,-N-2))
points(log(ORs),(-N:-1),pch="+",cex=0.75)
for (i in 1:N) {
   temp1= as.numeric(x2[[1]][,which(colnames(x2[[1]])==paste("specA[",i,"]",sep=""))])
   temp2= as.numeric(x2[[1]][,which(colnames(x2[[1]])==paste("specB[",i,"]",sep=""))])
   temp3 = log((temp1*(1-temp2))/((1-temp1)*temp2))
   temp4 = quantile(temp3,probs=c(0.025,0.50,0.975))
   points(temp4[2],-N+(i-1),pch=16)
   lines( c(temp4[1],temp4[3]),c(-N+(i-1),-N+(i-1)))
}
abline(v=log(help1[which(row.names(help1)=="specoddsratioAvsB"),2]),lty=3,lwd=0.75)
title("posterior medians and 95% credibility intervals",line=-2,outer=TRUE)
mtext("(observed values are indicated by '+')",line=-3,outer=TRUE,cex=0.75)
#dev.off()

stats=summary(x[[1]])$stat
dev.new()
#pdf("fig3.pdf")
plot(1,1,xlab="1-specificity",ylab="sensitivity",xlim=c(0,1),ylim=c(0,1),type="n",main="means and 95% confidence ellipses in ROC-space for both tests")
for (i in 1:N) {
   points(1-help2[which(row.names(help2)==paste("specA[",i,"]",sep="")),2],help2[which(row.names(help2)==paste("sensA[",i,"]",sep="")),2],pch=16,col=2,cex=0.5)
   lines(c(1-help2[which(row.names(help2)==paste("specA[",i,"]",sep="")),1],1-help2[which(row.names(help2)==paste("specA[",i,"]",sep="")),3]) ,
         c(help2[which(row.names(help2)==paste("sensA[",i,"]",sep="")),2],help2[which(row.names(help2)==paste("sensA[",i,"]",sep="")),2]),lty=3,col=2,lwd=0.75)
   lines(c(1-help2[which(row.names(help2)==paste("specA[",i,"]",sep="")),2],1-help2[which(row.names(help2)==paste("specA[",i,"]",sep="")),2]) ,
         c(help2[which(row.names(help2)==paste("sensA[",i,"]",sep="")),1],help2[which(row.names(help2)==paste("sensA[",i,"]",sep="")),3]),lty=3,col=2,lwd=0.75)
   points(1-help2[which(row.names(help2)==paste("specB[",i,"]",sep="")),2],help2[which(row.names(help2)==paste("sensB[",i,"]",sep="")),2],pch=16,col=3,cex=0.5)
   lines(c(1-help2[which(row.names(help2)==paste("specB[",i,"]",sep="")),1],1-help2[which(row.names(help2)==paste("specB[",i,"]",sep="")),3]) ,
         c(help2[which(row.names(help2)==paste("sensB[",i,"]",sep="")),2],help2[which(row.names(help2)==paste("sensB[",i,"]",sep="")),2]),lty=3,col=3,lwd=0.75)
   lines(c(1-help2[which(row.names(help2)==paste("specB[",i,"]",sep="")),2],1-help2[which(row.names(help2)==paste("specB[",i,"]",sep="")),2]) ,
         c(help2[which(row.names(help2)==paste("sensB[",i,"]",sep="")),1],help2[which(row.names(help2)==paste("sensB[",i,"]",sep="")),3]),lty=3,col=3,lwd=0.75)
}
legend(x=0.8,y=0.1,legend=c("test A","test B"),pch=16,col=c(2,3),bty="n")
points(1-help1[which(row.names(help1)=="meanspecificityA"),2],help1[which(row.names(help1)=="meansensitivityA"),2],col=2,pch=16,cex=1.5)
points(1-help1[which(row.names(help1)=="meanspecificityB"),2],help1[which(row.names(help1)=="meansensitivityB"),2],col=3,pch=16,cex=1.5)
mu1= stats[which(row.names(stats)=="meanalpha[3]"),1]
mu2= stats[which(row.names(stats)=="meanalpha[1]"),1]
var1= stats[which(row.names(stats)=="meanalpha[3]"),2]^2
var2= stats[which(row.names(stats)=="meanalpha[1]"),2]^2
covv=cov(x[[1]][,which(colnames(x[[1]])=="meanalpha[1]")],x[[1]][,which(colnames(x[[1]])=="meanalpha[3]")])
sematrix=matrix(c(var1,-covv,-covv,var2),ncol=2,nrow=2)
punten1=confellips(mu=c(-mu1,mu2),sigma=sematrix,alfa=0.05,npoints=1000)
lines(invlogit(punten1[,1]),invlogit(punten1[,2]),col=3,lty=2)
var1= stats[which(row.names(stats)=="Sigma[3,3]"),1]
var2= stats[which(row.names(stats)=="Sigma[1,1]"),1]
covv= stats[which(row.names(stats)=="Sigma[1,3]"),1]
sematrix=matrix(c(var1,-covv,-covv,var2),ncol=2,nrow=2)
punten1=confellips(mu=c(-mu1,mu2),sigma=sematrix,alfa=0.05,npoints=1000)
lines(invlogit(punten1[,1]),invlogit(punten1[,2]),col=3,lty=1)
minxxx=invlogit(min(punten1[,1]))
maxxxx=invlogit(max(punten1[,1]))
minyyy=invlogit(min(punten1[,2]))
maxyyy=invlogit(max(punten1[,2]))
xxx=(1-seq(0.01,0.99,0.01))
yyy=invlogit(stats[which(row.names(stats)=="interceptB"),1] + stats[which(row.names(stats)=="slopeB"),1]*logit(1-seq(0.01,0.99,0.01)))
lines(xxx[xxx>minxxx & xxx<maxxxx & yyy>minyyy & yyy<maxyyy],yyy[xxx>minxxx & xxx<maxxxx & yyy>minyyy & yyy<maxyyy],lty=1,col=3)
mu1= stats[which(row.names(stats)=="meanalpha[3]"),1] + stats[which(row.names(stats)=="meanalpha[4]"),1]
mu2= stats[which(row.names(stats)=="meanalpha[1]"),1] + stats[which(row.names(stats)=="meanalpha[2]"),1]
var1= stats[which(row.names(stats)=="meanalpha[3]"),2]^2 + stats[which(row.names(stats)=="meanalpha[4]"),2]^2
var2= stats[which(row.names(stats)=="meanalpha[1]"),2]^2 + stats[which(row.names(stats)=="meanalpha[2]"),2]^2
covv=cov(x[[1]][,which(colnames(x[[1]])=="meanalpha[1]")]+x[[1]][,which(colnames(x[[1]])=="meanalpha[2]")],
         x[[1]][,which(colnames(x[[1]])=="meanalpha[3]")]+x[[1]][,which(colnames(x[[1]])=="meanalpha[4]")])
sematrix=matrix(c(var1,-covv,-covv,var2),ncol=2,nrow=2)
punten1=confellips(mu=c(-mu1,mu2),sigma=sematrix,alfa=0.05,npoints=1000)
lines(invlogit(punten1[,1]),invlogit(punten1[,2]),col=2,lty=2)
var1= stats[which(row.names(stats)=="Sigma[3,3]"),1] + stats[which(row.names(stats)=="Sigma[4,4]"),1] + 2*stats[which(row.names(stats)=="Sigma[3,4]"),1]
var2= stats[which(row.names(stats)=="Sigma[1,1]"),1] + stats[which(row.names(stats)=="Sigma[2,2]"),1] + 2*stats[which(row.names(stats)=="Sigma[1,2]"),1]
covv= stats[which(row.names(stats)=="Sigma[1,3]"),1] + stats[which(row.names(stats)=="Sigma[1,4]"),1] + 
      stats[which(row.names(stats)=="Sigma[2,3]"),1] + stats[which(row.names(stats)=="Sigma[2,4]"),1]
sematrix=matrix(c(var1,-covv,-covv,var2),ncol=2,nrow=2)
punten1=confellips(mu=c(-mu1,mu2),sigma=sematrix,alfa=0.05,npoints=1000)
lines(invlogit(punten1[,1]),invlogit(punten1[,2]),col=2,lty=1)
minxxx=invlogit(min(punten1[,1]))
maxxxx=invlogit(max(punten1[,1]))
minyyy=invlogit(min(punten1[,2]))
maxyyy=invlogit(max(punten1[,2]))
xxx=(1-seq(0.01,0.99,0.01))
yyy=invlogit(stats[which(row.names(stats)=="interceptA"),1] + stats[which(row.names(stats)=="slopeA"),1]*logit(1-seq(0.01,0.99,0.01)))
lines(xxx[xxx>minxxx & xxx<maxxxx & yyy>minyyy & yyy<maxyyy],yyy[xxx>minxxx & xxx<maxxxx & yyy>minyyy & yyy<maxyyy],lty=1,col=2)
#dev.off()



dev.new()
#pdf("fig4.pdf)
layout(matrix(c(1,2),nrow=1,ncol=2))
logORs=log(((d2$Capp+delta)*(d2$Camm+delta)) / ((d2$Capm+delta)*(d2$Camp+delta)))
selogORs = sqrt((1/(d2$Capp+delta)) + (1/(d2$Camp+delta)) + (1/(d2$Capm+delta)) + (1/(d2$Camm+delta)))
minxx=min((logORs-1.96*selogORs))
#minxx=-6
maxxx=max((logORs+1.96*selogORs))
#maxxx=6
plot(1,1,xlab="log Odds Ratio in Cases",ylab="study number",xlim=c(minxx,maxxx),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2):(-1)),labels=c("overall","",1:N))
for (i in 1:N) {
   points(log(help2[which(row.names(help2)==paste("ORcases[",i,"]",sep="")),2]),-N+(i-1),pch=16)
   lines(log(c(help2[which(row.names(help2)==paste("ORcases[",i,"]",sep="")),1],help2[which(row.names(help2)==paste("ORcases[",i,"]",sep="")),3])),
         c(-N+(i-1),-N+(i-1)))
}
points(help1[which(row.names(help1)=="meanalpha[5]"),2],-N-2,pch=16)
lines(c(help1[which(row.names(help1)=="meanalpha[5]"),1],help1[which(row.names(help1)=="meanalpha[5]"),3]),c(-N-2,-N-2))
abline(v=help1[which(row.names(help1)=="meanalpha[5]"),2],lty=3)
points(logORs,(-(N-nrow(d1)):(-1)),pch="+",col=1)
logORs=log(((d2$Comm+delta)*(d2$Copp+delta)) / ((d2$Copm+delta)*(d2$Comp+delta)))
selogORs = sqrt((1/(d2$Copp+delta)) + (1/(d2$Comp+delta)) + (1/(d2$Copm+delta)) + (1/(d2$Comm+delta)))
minxx=min((logORs-1.96*selogORs))
maxxx=max((logORs+1.96*selogORs))
plot(1,1,xlab="log Odds Ratio in Controls",ylab="study number",xlim=c(minxx,maxxx),main="",ylim=c(-(N+2),-1),type="n",yaxt="n")
axis(2,at=c(-(N+2):(-1)),labels=c("overall","",1:N))
for (i in 1:N) {
   points(log(help2[which(row.names(help2)==paste("ORcontrols[",i,"]",sep="")),2]),-N+(i-1),pch=16)
   lines(log(c(help2[which(row.names(help2)==paste("ORcontrols[",i,"]",sep="")),1],help2[which(row.names(help2)==paste("ORcontrols[",i,"]",sep="")),3])),
         c(-N+(i-1),-N+(i-1)))
}
points(help1[which(row.names(help1)=="meanalpha[6]"),2],-N-2,pch=16)
lines(c(help1[which(row.names(help1)=="meanalpha[6]"),1],help1[which(row.names(help1)=="meanalpha[6]"),3]),c(-N-2,-N-2))
abline(v=help1[which(row.names(help1)=="meanalpha[6]"),2],lty=3)
points(logORs,(-(N-nrow(d1)):(-1)),pch="+",col=1)
title("posterior medians and 95% credibility intervals",outer=TRUE, line=-2)
mtext("(observed values are indicated by '+')",line=-3,outer=TRUE,cex=0.75)
#dev.off()


### do separate analyses with the mada package, just for fun and comparison sake
library(mada)
xx=data.frame(TP=d$TPsA,FN=d$FNsA,FP=d$FPsA,TN=d$TNsA)
madares=reitsma(xx, method = "reml",predict=T,sroc.type="naive")
summary(madares)
AUC(madares,sroc.type="naive")
dev.new()
par(mfrow=c(1,2))
plot(madares,predict=T,type="naive",main="test A: naive")
plot(madares,predict=T,type="ruttergatsonis",main="test A: Rutter-Gatsonis")

xx=data.frame(TP=d$TPsB,FN=d$FNsB,FP=d$FPsB,TN=d$TNsB)
madares=reitsma(xx, method = "reml",predict=T,sroc.type="naive")
summary(madares)
AUC(madares,sroc.type="naive")
dev.new()
par(mfrow=c(1,2))
plot(madares,predict=T,type="naive",main="test B: naive")
plot(madares,predict=T,type="ruttergatsonis",main="test B: Rutter-Gatsonis")














