
###
###  Bayesian Accuracy Data Analysis (same as in the R-syntax file "bivariate MA diagnostic accuracy.R", but models extended with covariables)
###

###
### data are to be provided in data.frame- or matrix-form.
### Diagnostic accuracy data of a binary test is expected to be available of a number of independent studies, at study-level.
### Several formats can be read but most efficient case is when for each study the number of true and false positives and negatives are given, preferable with names TN, FP, TP, FN.
### Covariable data are also at study-level and most efficient is the case to provide covariate information in separate matrices. Covariates for controls and cases may be the same or different.
### Covariates can be of various type (numeric, factor, character), but are relatively simply translated to numeric format.
###
### it is important to set the names or colnames attributes of the data.frames or matrices.
### The R mada-package provides an example data-set: data of the first five rows of the smoking-data in mada with the accuracy-data columns TN, FP, TP, FN and 2 covariates called test and population
###
###    TN  FP  TP FN type population
### 1 324  28  21 15  SAQ          S
### 2 969 120  90 10  SAQ          S
### 3 232  26 104  8  SAQ          G
### 4 673  92 332 18  SAQ          G
### 5  77   2   3  0  SAQ          S
###

## Import example data from the mada package
library(mada)
data("smoking")
d0=smoking

# put the covariates in separate matrices (x1_0 for sensitivity/theta and x2_0 for specificity/alpha)
x1_0=x2_0=cbind(as.numeric(d0$type)-1,as.numeric(d0$population)-1)
colnames(x1_0)=colnames(x2_0)=c("type", "population")
rm(smoking)

## or use my own example data of studies for bladder cancer diagnosis
d0=data.frame(TP=c(c(14,18,2,24,38),c(28,40,4,42,64)),        ncases=c(c(40,55,6,54,71),c(40,55,6,57,71)),
              TN=c(c( 90,50,185,55,11),c( 68,26,159,101,13)), ncontrols=c(c(100,50,193,59,14),c(100,50,193,139,17)),
              test=c(rep("A",5),rep("B",5)))
x1_0=x2_0=data.frame(test=d0$test)

## or use the data of Schuetz GM, Zacharopoulou NM, Schlattmann P, Dewey M. Annals of Internal Medicine. 2010;152(3):167-177.
d0=read.csv("schuetz.csv")
d0[1:10,]         # 108 rijen met tp, fp, fn, tn
table(d0$Test )   # 89 CT, 19 M~RI
x1_0=x2_0=data.frame(Test=d0$Test)

## main function


   library(metafor)
   library(mvmeta)
   library(lme4)
   library(rjags)
   source("bada2 functions.R")

   descplot=1    #(=0 is default)
   delta=0.5     #(=0.5 is default, correction factor for the case that sens/spec=0 or 1)
   SEyes=1       #(=0 is default)
   beta=0        #(=0 is default, scaling factor of the RG SROC model)
   
   ff1=checkdata1(d0, delta)
   if (ff1[[2]]<0) {
      if (ff1[[2]]== -1) {print("Dataset with accuracy-statistics is neither a data.frame nor a matrix.")}
      if (ff1[[2]]== -2) {print("Some of the accuracy-statistics columns are not numeric.")}
      if (ff1[[2]]== -3) {print("Some of the accuracy-statistics columns are less than zero.")}
      if (ff1[[2]]== -4) {print("Some of the sensitivity or specificity-statistics columns are below zero or larger than 1.")}
      return(d0)
   }
   d=ff1[[1]]
   rm(ff1)
   summarizingsensspec(d, beta=0)
   optimalbeta=initial_estimate_beta(d, -3, 3)
   dev.new()
   if (descplot==1) {descriptiveplot(d, delta=0.5, varyingsymbolsize=0, deltasymbolsize, minsymbolsize, SEyes=1, titel="Observed Sensitivity versus Specificity values")}
   ff2=summarizingcovariates(nrow(d),x1_0,x2_0)
   if (ff2[[3]]<0) {
      print("Covariate data is not usable.")
      return(d0)
   }
   ff3a=separate_sensspec_approximate(d, ff2[[1]], ff2[[2]])        # lm and metafor
      mu=Sigma=SE2matrix=covxy=NA
      mu=c(ff3a[[1]][1,1],ff3a[[3]][1,1])
      covxy = ff3a[[2]][1]*ff3a[[4]][1]*ff3a[[5]][1]
      Sigma=matrix(c(ff3a[[2]][1]^2, covxy, covxy, ff3a[[4]][1]^2), nrow=2, ncol=2)
      SE2matrix=matrix(c(ff3a[[1]][1,2]^2,-0.0,-0.0,ff3a[[3]][1,2]^2), ncol=2, nrow=2)
      dev.new()
      par(mfrow=c(1,2))
      bivmodelplot(mu, Sigma, SE2matrix, d, alfa=0.05, plottype=1, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      bivmodelplot(mu, Sigma, SE2matrix, d, alfa=0.05, plottype=2, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      title("separate analysis of lsens and lspec using univariate models with the unweighted normal distribution", line=-2.5, outer=T, cex=0.8)
      mu=Sigma=SE2matrix=covxy=NA
      mu=c(ff3a[[1]][1,3],ff3a[[3]][1,3])
      covxy = ff3a[[2]][2]*ff3a[[4]][2]*ff3a[[5]][2]
      Sigma=matrix(c(ff3a[[2]][2]^2, covxy, covxy, ff3a[[4]][2]^2), nrow=2, ncol=2)
      SE2matrix=matrix(c(ff3a[[1]][1,4]^2,-0.0,-0.0,ff3a[[3]][1,4]^2), ncol=2, nrow=2)
      dev.new()
      par(mfrow=c(1,2))
      bivmodelplot(mu, Sigma, SE2matrix, d, alfa=0.05, plottype=1, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      bivmodelplot(mu, Sigma, SE2matrix, d, alfa=0.05, plottype=2, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      title("separate analysis of lsens and lspec using univariate models with the normal-normal distribution", line=-2.5, outer=T, cex=0.8)
   ff3b=separate_gem_s_approximate(d, ff2[[1]], ff2[[2]], b=beta)        # lm and metafor and b=beta or optimalbeta
      mu=Sigma=SE2matrix=covxy=NA
      mu=c(ff3b[[1]][1,1],ff3b[[3]][1,1])
      covxy = 0
      Sigma=matrix(c(ff3b[[2]][1]^2, covxy, covxy, ff3b[[4]][1]^2), nrow=2, ncol=2)
      SE2matrix=matrix(c(ff3b[[1]][1,2]^2,-0.0,-0.0, ff3b[[3]][1,2]^2), ncol=2, nrow=2)
      dev.new()
      par(mfrow=c(2,2))
      RGmodelplot(mu, Sigma, SE2matrix, beta, d, alfa=0.05, plottype=3, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      RGmodelplot(mu, Sigma, SE2matrix, beta, d, alfa=0.05, plottype=1, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      RGmodelplot(mu, Sigma, SE2matrix, beta, d, alfa=0.05, plottype=2, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      title("separate analysis of (lsens+lspec)/2 and lsens-lspec using the RG SROC model with the normal distribution", line=-2.5, outer=T, cex=0.8)
      mu=Sigma=SE2matrix=covxy=NA
      mu=c(ff3b[[1]][1,3],ff3b[[3]][1,3])
      covxy = 0
      Sigma=matrix(c(ff3b[[2]][2]^2, covxy, covxy, ff3b[[4]][2]^2), nrow=2, ncol=2)
      SE2matrix=matrix(c(ff3b[[1]][1,4]^2,-0.0,-0.0, ff3b[[3]][1,4]^2), ncol=2, nrow=2)
      dev.new()
      par(mfrow=c(2,2))
      RGmodelplot(mu, Sigma, SE2matrix, beta, d, alfa=0.05, plottype=3, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      RGmodelplot(mu, Sigma, SE2matrix, beta, d, alfa=0.05, plottype=1, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      RGmodelplot(mu, Sigma, SE2matrix, beta, d, alfa=0.05, plottype=2, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      title("separate analysis of (lsens+lspec)/2 and lsens-lspec using the RG SROC model with the normal-normal distribution", line=-2.5, outer=T, cex=0.8)
   ff3c=joint_analysis_approximate(d,ff2[[1]],ff2[[2]],0)                   # mvmeta: same set of covs for sens/theta and spec/alpha (for the time being)
   #madares=reitsma(d, formula = cbind(tsens, tfpr) ~ type + population)    # same as ff3c, but results without normalization of covs
   #madares=reitsma(d, formula = cbind(tsens, tfpr) ~ test)                 # same as ff3c, maar met niet genormaliseerde covariaten
   #summary(madares)
      mu=Sigma=SE2matrix=covxy=NA
      mu=c(ff3c[[1]][1,1], ff3c[[1]][1,3])
      Sigma=ff3c[[2]]
      covxy=summary(ff3c[[5]])$vcov[1,(2+ncol(ff2[[1]]))]
      SE2matrix=matrix(c(ff3c[[1]][1,2]^2, covxy, covxy, ff3c[[1]][1,4]^2), ncol=2, nrow=2)
      dev.new()
      par(mfrow=c(1,2))
      bivmodelplot(mu, Sigma, SE2matrix, d, alfa=0.05, plottype=1, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      bivmodelplot(mu, Sigma, SE2matrix, d, alfa=0.05, plottype=2, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      title("joint analysis of lsens and lspec using the bivariate model with the normal-normal distribution", line=-2.5, outer=T, cex=0.8)
      mu=Sigma=SE2matrix=covxy=NA
      mu=c(ff3c[[3]][1,1], ff3c[[3]][1,3])
      Sigma=ff3c[[4]]           ### cov thetaen alfa hoeft niet perse nul te zijn, nu wel maar is misschien beter om niet te doen omdat beta op nul wordt gezet
      covxy=summary(ff3c[[6]])$vcov[1,(2+ncol(ff2[[1]]))]
      SE2matrix=matrix(c(ff3c[[3]][1,2]^2, covxy, covxy, ff3c[[3]][1,4]^2), ncol=2, nrow=2)
      dev.new()
      par(mfrow=c(2,2))
      RGmodelplot(mu, Sigma, SE2matrix, beta, d, alfa=0.05, plottype=3, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      RGmodelplot(mu, Sigma, SE2matrix, beta, d, alfa=0.05, plottype=1, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      RGmodelplot(mu, Sigma, SE2matrix, beta, d, alfa=0.05, plottype=2, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      title("joint analysis of (lsens+lspec)/2 and lsens-lspec using the RG SROC model with the normal-normal distribution", line=-2.5, outer=T, cex=0.8)

   ff4=sens_spec_binomial(d,ff2[[1]], ff2[[2]])                    # glmer: separate and joint analyses (with same set of covs for sens and spec at joint analysis (for the time being))
      mu=Sigma=SE2matrix=covxy=NA
      mu=c(ff4[[1]][1,1], ff4[[3]][1,1])
      covxy=0                                                                                              # what is the covariance between sens and spec at separate analysis ?????
      covxy=cor(predict(ff4[[7]],type="link"), predict(ff4[[8]],type="link"))*sqrt(ff4[[2]]*ff4[[4]])      # not likely correct
      covxy=cor(resid(ff4[[7]],type="response"), resid(ff4[[8]],type="response"))*sqrt(ff4[[2]]*ff4[[4]])  # not likely correct
      Sigma=matrix(c(ff4[[2]], covxy, covxy, ff4[[4]]), nrow=2, ncol=2)
      SE2matrix=matrix(c(ff4[[1]][1,2]^2, 0, 0, ff4[[3]][1,2]^2), ncol=2, nrow=2)
      dev.new()
      par(mfrow=c(1,2))
      bivmodelplot(mu, Sigma, SE2matrix, d, alfa=0.05, plottype=1, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      bivmodelplot(mu, Sigma, SE2matrix, d, alfa=0.05, plottype=2, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      title("separate analysis of lsens and lspec using the bivariate model with the binomial-normal distribution", line=-2.5, outer=T, cex=0.8)
      mu=Sigma=SE2matrix=covxy=NA
      mu=c(ff4[[5]][1,1], ff4[[5]][1,3])
      Sigma=ff4[[6]]
      covxy=summary(ff4[[9]])$vcov[1,2]
      SE2matrix=matrix(c(ff4[[5]][1,2]^2, covxy, covxy, ff4[[5]][1,4]^2), ncol=2, nrow=2)
      dev.new()
      par(mfrow=c(1,2))
      bivmodelplot(mu, Sigma, SE2matrix, d, alfa=0.05, plottype=1, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      bivmodelplot(mu, Sigma, SE2matrix, d, alfa=0.05, plottype=2, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes, titel="")
      title("joint analysis of lsens and lspec using the bivariate model with the binomial-normal distribution", line=-2.5, outer=T, cex=0.8)
      
   ff5=JAGSbivmodel(d, ff2[[1]], ff2[[2]], niter=50000, thin=10)
      dev.new()
      par(mfrow=c(1,2))
      Bivmodelplot_JAGS(output=ff5, d=d, plottype=1, maxaantal=100, alfa=0.05, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, titel="") 
      auceesBiv=Bivmodelplot_JAGS(output=ff5, d=d, plottype=2, maxaantal=100, alfa=0.05, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, titel="") 
      title("joint analysis of lsens and lspec using the bivariate model with the Bayesian binomial-normal distribution", line=-2.5, outer=T, cex=0.8)
   ff6=JAGS_RGmodel(d, ff2[[1]], ff2[[2]], niter=50000, thin=10)
      dev.new()
      par(mfrow=c(2,2))
      RGmodelplot_JAGS(output=ff6, d=d, plottype=3, maxaantal=100, alfa=0.05, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, titel="") 
      RGmodelplot_JAGS(output=ff6, d=d, plottype=1, maxaantal=100, alfa=0.05, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, titel="") 
      auceesRG=RGmodelplot_JAGS(output=ff6, d=d, plottype=2, maxaantal=100, alfa=0.05, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, titel="Rutter-Gatsonis SROC model: JAGS") 
      title("joint analysis of (lsens+lspec)/2 and lsens-lspec using the RG SROC model with the Bayesian binomial-normal distribution", line=-2.5, outer=T, cex=0.8)



##########################################################################################################################################################################################
## By the way, RG model: see Harbord et al. Biostatistics, 2007
##
## The RG model is a development of the Moses-Littenberg method, which weighing was incorrect and R&G fixed this.
##
##   for study i
##      logit(sens)[i] = (theta[i]+0.5*alpha[i]) * exp(-0.5*beta)
##      logit(spec)[i] = (theta[i]-0.5*alpha[i]) * exp( 0.5*beta)
##
##   Wtf are these parameters theta and alpha? 
##
##   Suppose beta=0, then exp(0.5*beta)=1. Now calculate both the sum and the difference of logit(sens) and logit(spec)
##
##      logit(sens[i])+logit(spec)[i] = log(DOR)[i] = (theta[i]+0.5*alpha[i]) + (theta[i]-0.5*alpha[i]) = 2*theta[i]
##      logit(sens[i])-logit(spec)[i] = (theta[i]+0.5*alpha[i]) - (theta[i]-0.5*alpha[i]) = alpha[i]
##
##   Apparently, theta has the interpretation of being half the logDOR; interpretation of alpha is not straightforward. But perhaps it helps to
##   point to the Bland-Altman plot/procedure that also is based on associating the difference versus the sum of two repeated measures of the same variable.
##   (Don't know how the parallel should be taken further, in the BA-plot the correlation between the difference and the sum is used to test whether the
##    variances of the 2 measures are different, but here theta and alpha are fixed to be independent of each other....)
##
##   Anyway, when theta is modelled as a function of covariates, the interpretation is pretty clear, namely to see how the logDOR depends on (study-)characteristics.
##   I am less clear on what it means to let alpha depend on covariates.... It is possible and therefore we do (?).
##
##   The parameter beta has a kind of covariance function in RG's model, because theta and alpha are fixed to be independent and by introducing beta makes RG's model similar 
##      to the bivariate model (when there are no covariates). It functions as a shape parameter for the form of the ROC-curve. If beta=0 then 'true' logDOR is said to be constant.
##########################################################################################################################################################################################


