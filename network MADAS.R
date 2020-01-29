

### some settings
analyze_single_studies = T      # incluse single test studies in the analysis (TRUE) or not (FALSE) ?
delta = 0.5                     # add 0.5 to cell-numbers? (only used for plotting 95% CIs of the observed sens & spec)
iternr = 50000                  # maximum number of iterations
xas=seq(0.01,0.99,0.01)         # steps on the x-axis for calculation of roc-curves and AUCs
P=length(xas)


### read the data
###    read from different files (1) the data of studies that evaluated one or multiple tests but do not report the 2-by-2 tables of test-results against each other
###                          and (2) the data of studies that evaluated two tests and do provide the 2-by-2 tables of test-results against each other
setwd("/users/ahzwinderman/desktop/talen/diagnostic-accuracy-meta-analysis/netwerk MADAS")
d1=read.csv("book1.csv",sep=",",header=TRUE)
d2=read.csv("book2.csv",sep=",",header=TRUE)

### combine into 1 data.frame
d=rbind(d1,
      data.frame(studynr=d2$studynr,testnr=d2$testnr1,TP=d2$Capp+d2$Capm,cases=d2$Capp+d2$Capm+d2$Camp+d2$Camm,TN=d2$Comm+d2$Comp,controls=d2$Copp+d2$Copm+d2$Comp+d2$Comm,studyname=d2$studyname,testname=d2$testname1),
      data.frame(studynr=d2$studynr,testnr=d2$testnr2,TP=d2$Capp+d2$Camp,cases=d2$Capp+d2$Capm+d2$Camp+d2$Camm,TN=d2$Comm+d2$Copm,controls=d2$Copp+d2$Copm+d2$Comp+d2$Comm,studyname=d2$studyname,testname=d2$testname2)
      )
rangorde=order(d$studynr,d$testnr)
d=d[rangorde,]
d                              # view data

### some summarizing numbers
N=length(unique(d$studynr))    # number of studies
N
K=length(unique(d$testnr))     # number of tests
K
table((d$studyname))           # how many tests per study
table((d$testname))            # how many studies per test


k1=length(which(table(d$studyname)==1))      # number of studies evaluating 1 test
k1
k2=length(which(table(d$studyname)> 1))      # number of studies evaluatiung multiple tests
k2

### exclude single test studies????
nrssinglestudies=which(d$studyname %in% names(which(table(d$studyname)==1)))
dx=d[-nrssinglestudies,]
dy=d
#table(d$testname[nrssinglestudies])
if (analyze_single_studies == FALSE) {d=dx}
nrow(d)                                      # number of test results


### re-organize: just to be sure, and define a new study-number
rangorde=order(d$studyname,d$testname)
d=d[rangorde,]
d$nummer=1
for (j in 2:nrow(d)) {
   if (d$studyname[j] == d$studyname[(j-1)]) {d$nummer[j]=d$nummer[j-1]}
   if (d$studyname[j] != d$studyname[(j-1)]) {d$nummer[j]=d$nummer[j-1]+1}
}
d

### useful functions
logit=function(x) {log(x/(1-x))}
invlogit=function(x) {exp(x)/(1+exp(x))}
confellips=function(mu,sigma,alfa,npoints) {
   es <- eigen(sigma)
   e1 <- es$vec %*% diag(sqrt(es$val))
   r1 <- sqrt(qchisq(1 - alfa, 2))
   theta <- seq(0, 2 * pi, len = npoints)
   v1 <- cbind(r1 * cos(theta), r1 * sin(theta))
   pts = t(mu - (e1 %*% t(v1)))
}
library(rjags)

### model in JAGS code
modelstring = "
   model {
      ### sample parameters for each study from the population distribution 
      ###     and transform these to sensitivity, specificity and odds ratio's
      for (i in 1:Nstudy) {
         alpha[i,1:npar] ~ dmnorm(meanalpha[1:npar],InvSigma[1:npar,1:npar])
         for (j in 1:K) {
            sens[i,j]     <- 1/(1+exp(-1*(alpha[i,j])))
            spec[i,j]     <- 1/(1+exp(-1*(alpha[i,(K+j)])))
         }
         ORcases[i]    <- exp(alpha[i,(2*K+1)])
         ORcontrols[i] <- exp(alpha[i,(2*K+2)])
      }
      ### likelihoods of the observed data of each row in the data-matrix: 
      ###    these are the marginals for all (type-1 and type-2) studies
      for (i in 1:Ntotaal){ 
         TP[i] ~ dbin(sens[studynr[i],testnr[i]],cases[i])
         TN[i] ~ dbin(spec[studynr[i],testnr[i]],controls[i])
      }
      ### likelihoods of the 2-by-2 tables for cases and controls
      ### only for type-2 studies
      for ( i in 1:N2) {
         Capp[i] ~ dhyper(Cap1[i],Cam1[i],Ca1p[i],ORcases[studynr2[i]])
         Comm[i] ~ dhyper(Com1[i],Cop1[i],Co1m[i],ORcontrols[studynr2[i]])
      }
      ### noninformative priors
      meanalpha[1:npar] ~ dmnorm(zeroos[1:npar],prec[1:npar,1:npar])
      InvSigma[1:npar,1:npar] ~ dwish(InvTau[1:npar,1:npar],npar)
      ### transformations of population parameters
      Sigma[1:npar,1:npar] <- inverse(InvSigma[1:npar,1:npar])
      ### mean sensitivities and specificities of the K tests
      ### slopes and intercepts of the ROC-curves of the K test
      ### Rutter and Gatsonis parameters (Harbord et al. Biostatistics, 2007, 8, 2, 239-251): beta, Lambda, Theta, skwtheta, skwalfa
      ### Naive and Rutter-Gatsonis AUCs of the K tests
      for (j in 1:K) {
         meansensitivity[j] <- 1/(1+exp(-1*(meanalpha[j])))
         meanspecificity[j] <- 1/(1+exp(-1*(meanalpha[(K+j)])))
         slope[j]     <- -1*Sigma[j,(K+j)]/Sigma[(K+j),(K+j)]
         intercept[j] <-  meanalpha[j] + slope[j]*meanalpha[(K+j)]
         # Rutter en Gatsonis
         beta[j]      <- 0.5*log(Sigma[(K+j),(K+j)]/Sigma[j,j])
         Lambda[j]    <- sqrt(sqrt(Sigma[(K+j),(K+j)])/sqrt(Sigma[j,j]))*meanalpha[j] + sqrt(sqrt(Sigma[j,j])/sqrt(Sigma[(K+j),(K+j)]))*meanalpha[(K+j)]
         Theta[j]     <- 0.5*(sqrt(sqrt(Sigma[(K+j),(K+j)])/sqrt(Sigma[j,j]))*meanalpha[j] - sqrt(sqrt(Sigma[j,j])/sqrt(Sigma[(K+j),(K+j)]))*meanalpha[(K+j)])
         skwtheta[j]  <- 0.5*(sqrt(Sigma[j,j]*Sigma[(K+j),(K+j)])-Sigma[j,(K+j)])
         skwalfa[j]   <- 2.0*(sqrt(Sigma[j,j]*Sigma[(K+j),(K+j)])+Sigma[j,(K+j)])
         #AUCs
         yna[j,1]  <- 1/(1+exp(-(intercept[j] + slope[j] * log(x[1]/(1-x[1])))))
         yrg[j,1]  <- 1/(1+exp(-(Lambda[j]*exp(beta[j]/2) + exp(beta[j])*log(x[1]/(1-x[1])))))
         oppna[j,1] <- (x[1]-0)*(yna[j,1]-0)/2
         opprg[j,1] <- (x[1]-0)*(yrg[j,1]-0)/2
         for (jj in 2:P) {
            yna[j,jj]   <- 1/(1+exp(-(intercept[j] + slope[j] * log(x[jj]/(1-x[jj])))))
            oppna[j,jj] <- (x[jj]-x[(jj-1)])*yna[j,(jj-1)]+(x[jj]-x[(jj-1)])*(yna[j,jj]-yna[j,(jj-1)])/2
            yrg[j,jj]   <- 1/(1+exp(-(Lambda[j]*exp(beta[j]/2) + exp(beta[j])*log(x[jj]/(1-x[jj])))))
            opprg[j,jj] <- (x[jj]-x[(jj-1)])*yrg[j,(jj-1)]+(x[jj]-x[(jj-1)])*(yrg[j,jj]-yrg[j,(jj-1)])/2
         }
         oppna[j,(P+1)] <- (1-x[P])*yna[j,P]+(1-x[P])*(1-yna[j,P])/2      
         auc_na[j] <- sum(oppna[j,1:(P+1)])
         opprg[j,(P+1)] <- (1-x[P])*yrg[j,P]+(1-x[P])*(1-yrg[j,P])/2      
         auc_rg[j] <- sum(opprg[j,1:(P+1)])
      }
   }  
"

### define non-informative hyper priors for 
# the logit-transformed parameters and the log-ORs (multivariate normal)
zeroos=rep(0,2*K+2)  
prec=diag(0.001,nrow=2*K+2,ncol=2*K+2)  
# and the inverse of the variance-covariance matrix of the study-specific parameters
InvTau=diag(0.1,nrow=2*K+2,ncol=2*K+2)

### put data and hyper priors in a list-object to be transferred to JAGS
jagsdatalist=list(npar=2*K+2, K=K, Nstudy=N, Ntotaal=nrow(d),
  TP=d$TP, cases=d$cases, TN=d$TN, controls=d$controls, testnr=d$testnr,
  studynr=d$nummer, N2=nrow(d2), studynr2=unique(d[(nrow(d)-nrow(d2)*2+1):nrow(d),"nummer"]),
  Capp=d2$Capp, Cap1=(d2$Capp+d2$Capm), Cam1=(d2$Camm+d2$Camp), Ca1p=(d2$Capp+d2$Camp),
  Comm=d2$Comm, Com1=(d2$Comm+d2$Comp), Cop1=(d2$Copp+d2$Copm), Co1m=(d2$Comm+d2$Copm),
  zeroos=zeroos, prec=prec, InvTau=InvTau, 
  x=xas, P=P)  

##### do the JAGS analysis
inits = list(meanalpha=rep(0,2*K+2))                                                                   # starting values
m = jags.model(textConnection(modelstring), data=jagsdatalist, inits, n.chains=2)                      # compile the analysis 
update(m, iternr)                                                                                      # do iternr iterations 
x = coda.samples(m, c("meanalpha","Sigma", "meansensitivity","meanspecificity",
        "slope","intercept","auc_na","auc_rg", "beta","Lambda","Theta","skwtheta","skwalfa"
         ), n.iter=iternr)                                                                             # draw iternr parameters from the posterior-distributions

#dev.new()
#pdf("traceplots.pdf")
#plot(x)   #,ask=TRUE)   
#dev.off()                                                                                               # check for convergence: trace-plots should display stable chaos


# summary of the structural and derived parameters
helpm1=summary(x)
helpm1
help0=helpm1$statistics
help1=helpm1$quantiles[,c(1,3,5)]


# plot de posterior medians van sens/spec (en 95% cred.intervals) instead of the observed sens/spec ???

#### plot in ROC-space
sens=(d$TP+delta)/(d$cases+2*delta)
senslow=invlogit(logit(sens)-1.96/(d$cases*sens*(1-sens)))
senshigh=invlogit(logit(sens)+1.96/(d$cases*sens*(1-sens)))
spec=(d$TN+delta)/(d$controls+2*delta)
speclow=invlogit(logit(spec)-1.96/(d$controls*spec*(1-spec)))
spechigh=invlogit(logit(spec)+1.96/(d$controls*spec*(1-spec)))
plot(1,1,type="n",xlab="1-specifivity",ylab="sensitivity",xlim=c(0,1),ylim=c(0,1))
for (i in 1:nrow(d)) {
   kleur=d$testnr[i]+1  #as.numeric(d$testname[i])+1
   points(1-spec[i],sens[i],pch="+",cex=0.75,col=kleur)
   lines(c(1-spec[i],1-spec[i]),c(senslow[i],senshigh[i]),lty=3,lwd=0.75,col=kleur)
   lines(c(1-spechigh[i],1-speclow[i]),c(sens[i],sens[i]),lty=3,lwd=0.75,col=kleur)
}
for (j in 1:K) {
   points(1-help1[which(row.names(help1)==paste("meanspecificity[",j,"]",sep="")),2],help1[which(row.names(help1)==paste("meansensitivity[",j,"]",sep="")),2],cex=1.5,pch=16,col=j+1)
   mu2= help0[which(row.names(help0)==paste("meanalpha[",j,"]",sep="")),1]
   mu1= help0[which(row.names(help0)==paste("meanalpha[",(K+j),"]",sep="")),1]
   var2= help0[which(row.names(help0)==paste("meanalpha[",j,"]",sep="")),2]^2
   var1= help0[which(row.names(help0)==paste("meanalpha[",(K+j),"]",sep="")),2]^2
   covv=cov(x[[1]][,which(colnames(x[[1]])==paste("meanalpha[",j,"]",sep=""))],x[[1]][,which(colnames(x[[1]])==paste("meanalpha[",(K+j),"]",sep=""))])             # covariance based only on samples in chain 1
   sematrix=matrix(c(var1,-covv,-covv,var2),ncol=2,nrow=2)
   punten1=confellips(mu=c(-mu1,mu2),sigma=sematrix,alfa=0.05,npoints=1000)
   lines(invlogit(punten1[,1]),invlogit(punten1[,2]),col=j+1,lty=2)
   var1= help0[which(row.names(help0)==paste("Sigma[",(K+j),",",(K+j),"]",sep="")),1]
   var2= help0[which(row.names(help0)==paste("Sigma[",j,",",j,"]",sep="")),1]
   covv= cov(x[[1]][,which(colnames(x[[1]])==paste("Sigma[",j,",",j,"]",sep=""))],x[[1]][,which(colnames(x[[1]])==paste("Sigma[",(K+j),",",(K+j),"]",sep=""))])    # covariance based only on samples in chain 1
   sematrix=matrix(c(var1,-covv,-covv,var2),ncol=2,nrow=2)
   punten1=confellips(mu=c(-mu1,mu2),sigma=sematrix,alfa=0.05,npoints=1000)
   lines(invlogit(punten1[,1]),invlogit(punten1[,2]),col=j+1,lty=1)
   minxxx=invlogit(min(punten1[,1]))
   maxxxx=invlogit(max(punten1[,1]))
   minyyy=invlogit(min(punten1[,2]))
   maxyyy=invlogit(max(punten1[,2]))
   xxx=(1-seq(0.01,0.99,0.01))
   yyy=invlogit(help0[which(row.names(help0)==paste("intercept[",j,"]",sep="")),1] + help0[which(row.names(help0)==paste("slope[",j,"]",sep="")),1]*logit(1-seq(0.01,0.99,0.01)))
   lines(xxx[xxx>minxxx & xxx<maxxxx & yyy>minyyy & yyy<maxyyy],yyy[xxx>minxxx & xxx<maxxxx & yyy>minyyy & yyy<maxyyy],lty=1,col=j+1)
}
legend(x=0.7,y=0.4,legend=(names(table(d$testname))),col=as.numeric(names(table(d$testnr)))+1,pch=16,lty=1,bty="n")

# differences between AUCs of all pairs of tests
#    median differences and the 2.5th and 97.5th percentiles
AUCs=matrix(NA,nrow=K*(K-1)/2,ncol=9)
delta_auc_na=delta_auc_rg=rep(NA,3)
tel=0
for (j in 1:(K-1))  {
   for (jj in (j+1):K) {
      tel=tel+1
      temp1=temp2=NA
      temp1=x[[1]][,which(colnames(x[[1]])==paste("auc_na[",j,"]",sep=""))] - x[[1]][,which(colnames(x[[1]])==paste("auc_na[",jj,"]",sep=""))]
      temp2=x[[2]][,which(colnames(x[[1]])==paste("auc_na[",j,"]",sep=""))] - x[[1]][,which(colnames(x[[2]])==paste("auc_na[",jj,"]",sep=""))]
      delta_auc_na[1:3] = quantile(c(temp1,temp2),probs=c(0.025,0.5,0.975))
      temp1=temp2=NA
      temp1=x[[1]][,which(colnames(x[[1]])==paste("auc_rg[",j,"]",sep=""))] - x[[1]][,which(colnames(x[[1]])==paste("auc_rg[",jj,"]",sep=""))]
      temp2=x[[2]][,which(colnames(x[[1]])==paste("auc_rg[",j,"]",sep=""))] - x[[1]][,which(colnames(x[[2]])==paste("auc_rg[",jj,"]",sep=""))]
      delta_auc_rg[1:3] = quantile(c(temp1,temp2),probs=c(0.025,0.5,0.975))
      AUCs[tel,1:9]=c(round(c(tel,j,jj),0),c(delta_auc_na[1:3],delta_auc_rg[1:3]))
   }
}
colnames(AUCs)=c("pairnr","test1","test2","delta_na_low","delta_na_median","delta_na_high","delta_rg_low","delta_rg_median","delta_nrg_high")
AUCs

# differences between senss of all pairs of tests
#    median differences and the 2.5th and 97.5th percentiles

# differences between specs of all pairs of tests
#    median differences and the 2.5th and 97.5th percentiles



