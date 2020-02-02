

########
######## Meta-analysis of a set of studies of the diagnostic accuracy a quantitative diagnostic test, that report the accuracies-results 
########           at a number of different thresholds. 
########
######## Per study with k different thresholds theta_1, theta_2, ..., theta_k the observed numbers of cases and controls in the ordered
########           categories defined by theta_j < Y < theta_(j+1) is reconstructed. The vector of numbers of cases (and controls) 
########           is assumed to follow a multinomial distribution with parameters pi_1 to pi_(k+1). These category-probabilities may
########           vary between studies. So the meta-analysis model is a random-effects model. 
######## The vectors of "logit"-transformed category-probabilities (i.e. softmax transforms) for the cases and controls of a study are 
########           assumed to be drawn from multivariate normal distributions (for cases and controls).
########
######## The parameters are estimated using an MCMC-algorithm in a Bayesian framework with uninformed prior distributions.
########

######## The data are in the form of a data.frame with studyid, cutoff-value and the number of true positive (TP),
########          false nagetive (FN), true negative (TN) and false positive (FP) per cutoof-values. An example is given below.
########
########   study_id    n cutoff   TP   FN  FP  TN
########          1   56   5.96   15    2   1  38
########          1   56   5.94   16    1   3  36
########          2  158   7.70   28   10  38  82
########          2  158   9.80    7   31   0 120
########          3   82   9.80   13    2   5  62
########          4  192  10.17   35    9  15 133
########          4  192  -4.25   44    0 130  18
########          4  192   8.97   43    1  86  62
########          4  192   9.30   42    2  64  84
########          4  192   9.65   40    4  37 111
########          4  192  10.69   27   17   7 141
########          4  192  28.35   13   31   1 147
########          4  192  35.59    7   37   0 148
########          5 3183   9.80 1672  588 249 674
########          5 3183  11.30  452 1808  18 905
########          6  417   8.64  150   17 144 106
########          6  417  10.00   79   88  26 224
########          7   31   9.00    6    1   4  20
########

######## some constants
delta=0.0       ### add 0.5 to sens and spec cell-numbers (for plotting only)

setwd("/users/ahzwinderman/desktop/jenny")

######## read the data
qq1=read.csv("2.LITMUS-DEF-Advanced fibrosis -Converted to Siemens.csv",header=TRUE)

####################################################################
######## The R-package diagmeta can analyse this type of data, but it assumes that the results in the same study at different thresholds are independent
library(diagmeta)
#library(lme4)
Steinhauserresults = diagmeta(TP=qq1$TP,FP=qq1$FP,TN=qq1$TN,FN=qq1$FN,cutoff=qq1$cutoff, studlab=qq1$study, model="DIDS", distr="logistic",
    equalvar=FALSE, log.cutoff=FALSE, method.weights="invvar",incr=0.5,lambda=0.5)
print(Steinhauserresults)
diagstats(Steinhauserresults,cutoff=sort(unique(qq1$cutoff)))
par(mfrow=c(2,2))
plot(Steinhauserresults, which = "cdf",lines = TRUE, line.optcut = TRUE, ci = TRUE)
plot(Steinhauserresults, which = "roc",lines = TRUE, line.optcut = TRUE, ci = TRUE)
plot(Steinhauserresults, which = "sroc",lines = TRUE, line.optcut = TRUE, ci = TRUE)
plot(Steinhauserresults, which = "youden",lines = TRUE, line.optcut = TRUE, ci = TRUE)
#plot(Steinhauserresults, which = "survival",lines = TRUE, line.optcut = TRUE, ci = TRUE)
####################################################################


qq1=qq1[order(qq1$study_id,qq1$cutoff),]                                              ## order with respect to study-number
qq1$sens=(qq1$TP+delta)/(1+qq1$TP+qq1$FN)                                             ## add sensitivity and specificity
qq1$spec=(qq1$TN+delta)/(1+qq1$TN+qq1$FP)
studyid=unique(qq1$study_id)                                                          ## determine unique study-numbers
N=length(studyid)                                                                     ## determine how many studies
ncases=aggregate(round(qq1$n*qq1$prev,0),by=list(qq1$study_id),FUN=mean)[,2]          ## calculate numbers of cases and controls per study
ncontrols=aggregate(qq1$n,by=list(qq1$study_id),FUN=mean)[,2]-ncases
######## reorganize the data
datamat=list()                
k_i=c()       # calculate the number of cutpoints per study
for (i in 1:length(studyid)) 
{
	vv1=subset(qq1,study_id==studyid[i]) 
	k_i[i]=nrow(vv1)
	datamat[[i]]=vv1[,c("cutoff","spec","sens")]
}
######## the reorganized data
datamat


###################################################################
######## derive the observed numbers of cases (y1) and controls (y0) in the study-specific categories
kp1=k_i+1   ## (kp1[i] is the number of categories in study i)
y1 = y0 = array(0,dim=c(N,max(kp1)))
for (i in 1:N) {
   temp1=c(round(datamat[[i]][,2]*ncontrols[i]),ncontrols[i])
   temp1[2:(k_i[i]+1)] = temp1[2:(k_i[i]+1)] -temp1[1:(k_i[i])]
   y0[i,1:length(temp1)]=temp1
   temp2=c(ncases[i],round(datamat[[i]][,3]*ncases[i]))
   temp2[1:(k_i[i])] = temp2[1:(k_i[i])]-temp2[2:(k_i[i]+1)]
   y1[i,1:length(temp2)] = temp2
}
######## determine the total number of unique cutoff-values
######## and add an extra value that is larger than the maximum cutoff-value used in the N studies
xtra=1
theta=c()
for (i in 1:N) {theta=c(theta,datamat[[i]][,1])}
theta=c(sort(unique(theta)),max(theta)+xtra)
m=length(theta)
######## determine the relationship between tilde{pi} and pi
designmat=array(0,dim=c(N,m,max(kp1)))
for (i in 1:N) {
   for (j in 1:m) {
      if (theta[j] <= datamat[[i]][1,1]) {designmat[i,j,1]=1}
      if (k_i[i]>1) {
         for (el in 2:k_i[i]) {
            if ( theta[j]<= datamat[[i]][el,1]) {designmat[i,j,el]=1}
            if ( theta[j]<=datamat[[i]][(el-1),1]) {designmat[i,j,el]=0}
         }
      }
      if (theta[j] > datamat[[i]][k_i[i],1]) {designmat[i,j,kp1[i]]=1}
   }
}
dimnames(designmat)[[2]]=theta
######## check study i
##i=1
#designmat[i,,]


###################################################################
######## JAGS
library(rjags)

modelstring0 <- "
   model {
      for (i in 1:N) {
         ######## draw the study-specific parameters from the multivariate normal distribution
         a0[i, 1:(m-1)]  ~ dmnorm(mu0[1:(m-1)], invSigma0[1:(m-1), 1:(m-1)])
         a1[i, 1:(m-1)]  ~ dmnorm(mu1[1:(m-1)], invSigma1[1:(m-1), 1:(m-1)])
         ######## calculate the category-probabilities: pi
         teller0[i,1] <- 1
         for (j in 2:m) {teller0[i,j] <- exp(a0[i,(j-1)])}
         somteller0[i] <- sum(teller0[i,1:m])
         for (j in 1:m) {pi0[i,j] <- teller0[i,j]/somteller0[i]}         
         teller1[i,1] <- 1
         for (j in 2:m) {teller1[i,j] <- exp(a1[i,(j-1)])}
         somteller1[i] <- sum(teller1[i,1:m])
         for (j in 1:m) {pi1[i,j] <- teller1[i,j]/somteller1[i]}  
         ######## derive the study-specific category-probabilities: pitilde
         for (el in 1:(kp1[i])) {
            for (j in 1:m) {help1[i,j,el] <- designmat[i,j,el] * pi1[i,j]}
            pitilde1[i,el] <- sum(help1[i,1:m,el])
            for (j in 1:m) {help0[i,j,el] <- designmat[i,j,el] * pi0[i,j]}
            pitilde0[i,el] <- sum(help0[i,1:m,el])
         }
         y0[i,1:kp1[i]] ~ dmulti(pitilde0[i,1:kp1[i]],ncontrols[i])
         y1[i,1:kp1[i]] ~ dmulti(pitilde1[i,1:kp1[i]],ncases[i])
      }
      
      ######## hyperpriors
      mu0[1:(m-1)] ~ dmnorm(zeroos[1:(m-1)],prec[1:(m-1), 1:(m-1)])
      invSigma0[1:(m-1), 1:(m-1)] ~ dwish(R[1:(m-1), 1:(m-1)], df1)
      mu1[1:(m-1)] ~ dmnorm(zeroos[1:(m-1)],prec[1:(m-1), 1:(m-1)])
      invSigma1[1:(m-1), 1:(m-1)] ~ dwish(R[1:(m-1), 1:(m-1)], df1)
      
      ######## average category-probabilities (averaged over the N studies)
      teller0x[1] <- 1
      for (j in 2:m) {teller0x[j] <- exp(mu0[(j-1)])}
      somteller0x <- sum(teller0x[1:m])
      for (j in 1:m) {pi0x[j] <- teller0x[j]/somteller0x}         
      teller1x[1] <- 1
      for (j in 2:m) {teller1x[j] <- exp(mu1[(j-1)])}
      somteller1x <- sum(teller1x[1:m])
      for (j in 1:m) {pi1x[j] <- teller1x[j]/somteller1x}  
      sumpi1x[1] <- pi1x[1]
      sumpi0x[1] <- pi0x[1]
      for (j in 2:m) {sumpi1x[j] <- sum(pi1x[1:j])}  
      for (j in 2:m) {sumpi0x[j] <- sum(pi0x[1:j])}  
      ######## averaged spec, sens, Youden
      for (i in 1:(m-1)) {
         spec[i] <- sum(pi0x[1:i])
         emspec[i] <- 1-spec[i]
         sens[i] <- sum(pi1x[(i+1):m])
         Youden[i] <- 2 * ((1-lambda)*spec[i] + lambda*sens[i]) - 1
         emsp[(i+1)] <- emspec[i]
         se[(i+1)] <- sens[i]
      }
      se[1] <- 1
      emsp[1] <- 1
      se[(m+1)] <- 0
      emsp[(m+1)] <-0 
      ######## AUC/c-stat    
      w[1] <- 0
      for (j in 2:(m-1)) {w[j] <- (emsp[j]-emsp[(j-1)]) * se[(j-1)] + (emsp[j]-emsp[(j-1)]) * (se[j]-se[(j-1)])/2}
      cstat <- -1*(sum(w[1:(m-1)]))
}
"

######## construct a list-object to transport data from R to JAGS
dd0 = list(N=N, kp1=kp1, ncases=ncases, ncontrols=ncontrols, m=m, designmat=designmat, y0=y0, y1=y1, zeroos=rep(0,(m-1)), prec=diag(0.001,(m-1)), 
             R=diag(0.1,(m-1)), df1=(m-1),lambda=0.5)
######## check model-syntax, and compile the model-syntax
model0=jags.model(textConnection(modelstring0), data=dd0, n.chains = 2, n.adapt=0)
######## burn-in
update(model0, n.iter=50000)             # burn-in
######## posterior sampling
output.raw0 = coda.samples(model=model0,variable.names=c("pi1x","pi0x","sumpi1x","sumpi0x","spec","sens","Youden","cstat"), n.iter=50000, thin=10)
######## check convergence
#par(ask=T)
#plot(output.raw0)

######## show estimates
abcd=summary(output.raw0)
abcd

#save.image("jenny1b.Rdata")
#load("jenny1b.Rdata")

######## some preparation to display further results
namen=row.names(abcd$statistics)
statistiek=sapply(namen,function(x){strsplit(x,"[",fixed=TRUE)[[1]][1]})
nummer=sapply(sapply(namen,function(x){strsplit(x,"[",fixed=TRUE)[[1]][2]}),function(x){strsplit(x,"]",fixed=TRUE)[[1]][1]})

######## statistics: sens
cbind(theta[1:(m-1)], 
      abcd$statistics[which(statistiek=="sens"),],
      abcd$quantiles[which(statistiek=="sens"),])

######## statistics: spec
cbind(theta[1:(m-1)], 
      abcd$statistics[which(statistiek=="spec"),],
      abcd$quantiles[which(statistiek=="spec"),])

######## statistics: Youden
hh1=cbind(theta[1:(m-1)], 
      abcd$statistics[which(statistiek=="Youden"),],
      abcd$quantiles[which(statistiek=="Youden"),])
hh1
######## cutoff with largest (weighted) Youden
theta[which(hh1[,2] == max(hh1[,2]))]

######## statistics: c-statistic
c(abcd$statistics[which(statistiek=="cstat"),],
      abcd$quantiles[which(statistiek=="cstat"),]
     )


######## plaatje 1: cumulative distributions
dev.new()
#pdf("plaatje1.pdf")
par(mfrow=c(2,1))
par(mar=c(2,4,4,2)+0.1)
plot(1,1,type="n",xaxt="n",xlim=c(min(theta)-1.1,max(theta[1:(m-1)])+1),ylim=c(0,1),xlab="biomarker",ylab="cumulative distribution",main="controls")
axis(1,at=theta[1:(m-1)],labels=theta[1:(m-1)])
for (i in 1:N) {
   lines(c(min(theta)-1,datamat[[i]][,1],max(theta[1:(m-1)])+1),c(cumsum(y0[i,1:kp1[i]])/ncontrols[i],1),type="s",col=i+1,pcabcdh=NULL)
}
legend(x=15,y=0.8,legend=unique(qq1$study),col=1+(1:N),lty=1,bty="n",cex=0.6)
lines(c(min(theta)-1,theta[1:(m-1)]) , abcd$statistics[which(statistiek=="sumpi0x"),1], type="s", col=1, lwd=3)
lines(c(min(theta)-1,theta[1:(m-1)]) , abcd$quantiles[which(statistiek=="sumpi0x"),1], type="s", col=1, lwd=2,lty=2)
lines(c(min(theta)-1,theta[1:(m-1)]) , abcd$quantiles[which(statistiek=="sumpi0x"),5], type="s", col=1, lwd=2,lty=2)
par(mar=c(5,4,2,2)+0.1)
plot(1,1,type="n",xaxt="n",xlim=c(min(theta)-1.1,max(theta[1:(m-1)])+1),ylim=c(0,1),xlab="biomarker",ylab="cumulative distribution",main="cases")
axis(1,at=theta[1:(m-1)],labels=theta[1:(m-1)])
for (i in 1:N) {
   lines(c(min(theta)-1,datamat[[i]][,1],max(theta[1:(m-1)])+1),c(cumsum(y1[i,1:kp1[i]])/ncases[i],1), type="s", col=(i+1), pch=NULL,lwd=0.5+2*(ncases/sum(ncases)))
}
lines(c(min(theta)-1,theta[1:(m-1)]) , abcd$statistics[which(statistiek=="sumpi1x"),1], type="s", col=1, lwd=3)
lines(c(min(theta)-1,theta[1:(m-1)]) , abcd$quantiles[which(statistiek=="sumpi1x"),1], type="s", col=1, lwd=2,lty=2)
lines(c(min(theta)-1,theta[1:(m-1)]) , abcd$quantiles[which(statistiek=="sumpi1x"),5], type="s", col=1, lwd=2,lty=2)
par(mfrow=c(1,1))
par(mar=c(5,4,4,2)+0.1)
#dev.off()

######## plaatje 2
dev.new()
#pdf("plaatje2.pdf")
plot(1,1,xlab="1-specificity",ylab="sensitivity",xlim=c(0,1),ylim=c(0,1),type="n")
for (i in 1:N) {
	datamatt=datamat[[i]]
	datamatt=datamatt[order(datamatt$sens,1-datamatt$spec,datamatt$cutoff),]
	lines(c(0,1-datamatt[,2],1),c(0,datamatt[,3],1),col=i+1,pch=2+i)
}
sens=abcd$statistics[which(statistiek=="sens"),1]
sdsens=abcd$statistics[which(statistiek=="sens"),2]
lower=abcd$quantiles[which(statistiek=="sens"),1]
upper=abcd$quantiles[which(statistiek=="sens"),5]
emspec=1-abcd$statistics[which(statistiek=="spec"),1]
rang=order(emspec)
sens=sens[rang]
emspec=emspec[rang]
lower=lower[rang]
upper=upper[rang]
lines(c(0,emspec,1),c(0,sens,1),col=1,lwd=3)
lines(c(0,emspec,1),c(0,lower,1),col=1,lwd=1.5,lty=2)
lines(c(0,emspec,1),c(0,upper,1),col=1,lwd=1.5,lty=2)
legend(x=0.7,y=0.5,legend=unique(qq1$study),col=1+(1:N),lty=1,bty="n",cex=0.6,adj=0)
#dev.off()

######## plaatje 3: Youden index
dev.new()
#pdf("plaatje3.pdf")
plot(1,1,type="n",xlim=c(min(theta),max(theta[(m-1)])),ylim=c(0,1),xlab="cut point",ylab="Youden index")
for (i in 1:N) {
   lines(datamat[[i]][,1],datamat[[i]][,2]+datamat[[i]][,3]-1,type="b",pch=i,col=i+1)
}
lines(theta[1:(m-1)],abcd$statistics[which(statistiek=="Youden"),1],lwd=3,col=1,lty=1)
lines(theta[1:(m-1)],abcd$quantiles[which(statistiek=="Youden"),1],lwd=2,col=1,lty=3)
lines(theta[1:(m-1)],abcd$quantiles[which(statistiek=="Youden"),5],lwd=2,col=1,lty=3)
legend(x=20,y=1,legend=unique(qq1$study),col=1+(1:N),lty=1,bty="n",cex=0.6,adj=0)
#dev.off()

######## plaatje 4: densities
dev.new()
#pdf("plaatje4.pdf")
labs=c()
labs[1]=paste("<=",theta[1],sep="")
labs[2:(m)]=paste("(",theta[1:(m-1)],";",theta[2:m],"]",sep="")
labs[m]=paste(">",theta[m-1],sep="")
layout(matrix(c(1,2),nrow=2,ncol=1),heights=c(1,0.3))
par(mar=c(0,4,4,2)+0.1)
plot(c(0:m,m+0.01),c(0,abcd$statistics[which(statistiek=="pi0x"),1],0),type="S",lty=1,col=3,lwd=2,ylim=c(0,0.4),xaxt="n",ylab="density function",xlab="")
lines(c(0:m,m+0.01),c(0,abcd$quantiles[which(statistiek=="pi0x"),5],0),type="S",lty=2,col=3,lwd=0.8)
abline(h=0)
lines(c(0:m,m+0.01),c(0,abcd$statistics[which(statistiek=="pi1x"),1],0),type="S",lty=1,col=2,lwd=2)
lines(c(0:m,m+0.01),c(0,abcd$quantiles[which(statistiek=="pi1x"),5],0),type="S",lty=2,col=2,lwd=0.8)
axis(1,at=(1:m)-0.5,labels=FALSE)
legend(x=1,y=0.4,legend=c("cases","controls"),lty=1,col=c(2,3),bty="n")
par(mar=c(5,4,0,2)+0.1)
plot(0,0,xlab="",ylab="",xaxt="n",yaxt="n",bty="n",type="n",ylim=c(-10,1),xlim=c(0,m+0.01))
text(x=(1:m)-0.5,y=rep(0,m),labels=labs,cex=0.6,adj=c(1,1),srt=45)
mtext("biomarker",1,line=+0)
par(mar=c(5,4,4,2)+0.1)
#dev.off()






