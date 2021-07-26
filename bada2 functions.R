
#controlsettings=function(d) {
   alfa=0.05
   talfa=qnorm((1-alfa/2))
   delta=0.5
   varyingsymbolsize=1
   minsymbolsize=0.5
   deltasymbolsize=2.5
   SEyes=0
   npoints=1000
   beta=0
   descplot=0
   #return(list(alfa,talfa,delta,varyingsymbolsize,minsymbolsize,deltasymbolsize,SEyes,npoints,beta, descplot))
#}


Bivmodelplot_JAGS=function(output=ff6, d=d, plottype=3, maxaantal=100, alfa=0.05, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, titel="") {
   npoints=1000
   symbolsize=1
   if (varyingsymbolsize==1) {symbolsize=minsymbolsize+deltasymbolsize*(d$ncases+d$ncontrols)/max(d$ncases+d$ncontrols, na.rm=T)}
   j=1
   z=output[[j]]
   for (j in 2:length(output)) {z=rbind(z,output[[j]])}
   z=as.data.frame(z)
   zz=z[sample(1:nrow(z),min(nrow(z),maxaantal),replace=F),]
   zz=as.data.frame(zz)
   mu=c(mean(z$'mu[1]', na.rm=T), mean(z$'mu[2]'))
   Sigma=matrix(c(mean(z$'Sigma[1,1]', na.rm=T), mean(z$'Sigma[1,2]', na.rm=T), mean(z$'Sigma[2,1]', na.rm=T), mean(z$'Sigma[2,2]', na.rm=T)),nrow=2,ncol=2)
   SE2=matrix(c(sd(z$'mu[1]', na.rm=T),cov(z$'mu[1]',z$'mu[2]'),cov(z$'mu[1]',z$'mu[2]'),sd(z$'mu[2]', na.rm=T)),nrow=2,ncol=2)
   hh1=confellips(mu,Sigma,alfa,npoints)
   hh2=confellips(mu,SE2,alfa,npoints)
   if (plottype==1) {
      minx=min(c(hh1[,2],logit(d$spec)),na.rm=T)
      maxx=max(c(hh1[,2],logit(d$spec)),na.rm=T)
      miny=min(c(hh1[,1],logit(d$sens)),na.rm=T)
      maxy=max(c(hh1[,1],logit(d$sens)),na.rm=T)
      plot(0,0, xlab="logit(specificity)", ylab="logit(sensitivity)", xlim=c(minx, maxx), ylim=c(miny,maxy), type="n", main=titel)
      slope = zz$'Sigma[1,2]' / zz$'Sigma[2,2]'
      intercept = zz$'mu[1]' - slope * zz$'mu[2]'
      for (i in 1:length(slope)) {abline(a=intercept[i], b=slope[i], col="lightgrey", lwd=0.7)}
      points(logit(d$spec), logit(d$sens), pch=1, col=1, cex=symbolsize)
      points(mu[2], mu[1], col=2, cex=2, pch=16)
      lines(hh1[,2], hh1[,1], lty=2, col=2, lwd=1.25)
      lines(hh2[,2], hh2[,1], lty=3, col=2, lwd=1.25)
      gemslope=Sigma[1,2]/Sigma[2,2]
      gemintercept=mu[1] - gemslope*mu[2]
      abline(a=gemintercept, b=gemslope, col=3, lwd=2)
   }
   if (plottype != 1) {
      minx=min(c(hh1[,2],logit(d$spec)),na.rm=T)
      maxx=max(c(hh1[,2],logit(d$spec)),na.rm=T)
      xxx=seq(minx, maxx, 0.1)
      plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Specificity",ylab="Sensitivity",xaxt="n", main=titel)
      axis(1,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("1.0","0.8","0.6","0.4","0.2","0.0"))
      abline(a=0,b=1,lty=2,col=1,lwd=1)
      slope=zz$'Sigma[1,2]' / zz$'Sigma[2,2]'
      intercept=zz$'mu[1]' - slope * zz$'mu[2]'
      aucees=c()
      for (i in 1:length(slope)) {
         lines(antilogit(-xxx), antilogit(intercept[i]+slope[i]*xxx), col="lightgrey", lwd=0.7)
         aucees[i]=calculateauc(antilogit(-xxx), antilogit(intercept[i]+slope[i]*xxx), 0, 1)
      }
      points(1-d$spec,d$sens, pch=1, col=1, cex=symbolsize)
      points(antilogit(-mu[2]), antilogit(mu[1]), col=2, cex=2, pch=16)
      lines(antilogit(-hh1[,2]), antilogit(hh1[,1]), lty=2, col=2, lwd=1.25)
      lines(antilogit(-hh2[,2]), antilogit(hh2[,1]), lty=3, col=2, lwd=1.25)
      gemslope=Sigma[1,2]/Sigma[2,2]
      gemintercept=mu[1] - gemslope*mu[2]
      lines(antilogit(-xxx),antilogit(gemintercept+gemslope*xxx), col=3, lwd=2)
      meanauc=calculateauc(antilogit(-xxx),antilogit(gemintercept+gemslope*xxx), 0, 1)
      gg=quantile(aucees,probs=c(0.025,0.975),na.rm=T)
      text(x=0.4,y=0.3,paste("AUC = ",round(meanauc,4)," (95% CI: ", round(gg[1],4)," - ",round(gg[2],4),")",sep=""),adj=0)
      return(aucees)
   }
}

RGmodelplot_JAGS=function(output=ff6, d=d, plottype=3, maxaantal=100, alfa=0.05, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, titel="") {
   npoints=1000
   symbolsize=1
   if (varyingsymbolsize==1) {symbolsize=minsymbolsize+deltasymbolsize*(d$ncases+d$ncontrols)/max(d$ncases+d$ncontrols, na.rm=T)}
   j=1
   z=output[[j]]
   for (j in 2:length(output)) {z=rbind(z,output[[j]])}
   z=as.data.frame(z)
   b=mean(z$beta, na.rm=T)
   zz=z[sample(1:nrow(z),min(nrow(z),maxaantal),replace=F),]
   zz=as.data.frame(zz)
   thet=(logit(d$sens)*exp(b/2) + logit(d$spec)*exp(-b/2))/2
   alph=(logit(d$sens)*exp(b/2) - logit(d$spec)*exp(-b/2))
   mu=c(mean(z$'mu[1]', na.rm=T), mean(z$'mu[2]'))
   Sigma=matrix(c(mean(z$'Sigma[1]', na.rm=T),0,0,mean(z$'Sigma[2]', na.rm=T)),nrow=2,ncol=2)
   SE2=matrix(c(sd(z$'mu[1]', na.rm=T),cov(z$'mu[1]',z$'mu[2]'),cov(z$'mu[1]',z$'mu[2]'),sd(z$'mu[2]', na.rm=T)),nrow=2,ncol=2)
   hh1=confellips(mu,Sigma,alfa,npoints)
   hh2=confellips(mu,SE2,alfa,npoints)
   if (plottype==3) {
      minx=min(c(hh1[,1],thet),na.rm=T)
      maxx=max(c(hh1[,1],thet),na.rm=T)
      miny=min(c(hh1[,2],alph),na.rm=T)
      maxy=max(c(hh1[,2],alph),na.rm=T)
      plot(0,0, xlab="Theta", ylab="Alpha", xlim=c(minx, maxx), ylim=c(miny,maxy), type="n", main=titel)
      for (i in 1:nrow(zz)) {abline(h=zz$'mu[2]'[i], col="lightgrey", lwd=0.7)}
      points(thet, alph, pch=1, cex=symbolsize)
      points(mu[1], mu[2], col=2, cex=2, pch=16)
      lines(hh1[,1],hh1[,2], lty=2, col=2, lwd=1.25)
      lines(hh2[,1],hh2[,2], lty=3, col=2, lwd=1.25)
      gemslope=Sigma[1,2]/Sigma[1,1]
      gemintercept=mu[2]-gemslope*mu[1]
      abline(a=gemintercept, b=gemslope, col=3, lwd=2)
   }
   if (plottype==1) {
      minx=min(c((hh1[,1]-hh1[,2]/2)*exp(b/2),(thet-alph/2)*exp(b/2)),na.rm=T)
      maxx=max(c((hh1[,1]-hh1[,2]/2)*exp(b/2),(thet-alph/2)*exp(b/2)),na.rm=T)
      miny=min(c((hh1[,1]+hh1[,2]/2)*exp(-b/2),(thet+alph/2)*exp(-b/2)),na.rm=T)
      maxy=max(c((hh1[,1]+hh1[,2]/2)*exp(-b/2),(thet+alph/2)*exp(-b/2)),na.rm=T)
      plot(0,0, xlab="logit(specificity)", ylab="logit(sensitivity)", xlim=c(minx, maxx), ylim=c(miny,maxy), type="n", main=titel)
      slope=(zz$'Sigma[1]' - zz$'Sigma[2]'/4)/ (exp(zz$'beta')*(zz$'Sigma[1]'+zz$'Sigma[2]'/4))
      intercept=(zz$'mu[1]'+zz$'mu[2]'/2)*exp(-zz$'beta'/2) - slope *  (zz$'mu[1]'-zz$'mu[2]'/2)*exp(zz$'beta'/2)
      for (i in 1:length(slope)) {abline(a=intercept[i], b=slope[i], col="lightgrey", lwd=0.7)}
      points((thet-alph/2)*exp(b/2), (thet+alph/2)*exp(-b/2), pch=1, col=1, cex=symbolsize)
      points((mu[1]-mu[2]/2)*exp(b/2), (mu[1]+mu[2]/2)*exp(-b/2), col=2, cex=2, pch=16)
      lines((hh1[,1]-hh1[,2]/2)*exp(b/2),(hh1[,1]+hh1[,2]/2)*exp(-b/2), lty=2, col=2, lwd=1.25)
      lines((hh2[,1]-hh2[,2]/2)*exp(b/2),(hh2[,1]+hh2[,2]/2)*exp(-b/2), lty=3, col=2, lwd=1.25)
      gemslope=(Sigma[1,1]-Sigma[2,2]/4)/(exp(b)*(Sigma[1,1]+Sigma[2,2]/4))
      gemintercept=(mu[1]+mu[2]/2)*exp(-b/2) - gemslope*(mu[1]-mu[2]/2)*exp(b/2)
      abline(a=gemintercept, b=gemslope, col=3, lwd=2)
   }
   if (plottype != 1 & plottype !=3) {
      minx=min(c((hh1[,1]-hh1[,2]/2)*exp(b/2),(thet-alph/2)*exp(b/2)),na.rm=T)
      maxx=max(c((hh1[,1]-hh1[,2]/2)*exp(b/2),(thet-alph/2)*exp(b/2)),na.rm=T)
      xxx=seq(minx, maxx, 0.1)
      plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),xlab="Specificity",ylab="Sensitivity",xaxt="n", main=titel)
      axis(1,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("1.0","0.8","0.6","0.4","0.2","0.0"))
      abline(a=0,b=1,lty=2,col=1,lwd=1)
      slope=(zz$'Sigma[1]' - zz$'Sigma[2]'/4)/ (exp(zz$'beta')*(zz$'Sigma[1]'+zz$'Sigma[2]'/4))
      intercept=(zz$'mu[1]'+zz$'mu[2]'/2)*exp(-zz$'beta'/2) - slope *  (zz$'mu[1]'-zz$'mu[2]'/2)*exp(zz$'beta'/2)
      aucees=c()
      for (i in 1:length(slope)) {
         lines(antilogit(-xxx), antilogit(intercept[i]+slope[i]*xxx), col="lightgrey", lwd=0.7)
         aucees[i]=calculateauc(antilogit(-xxx), antilogit(intercept[i]+slope[i]*xxx), 0, 1)
      }
      points(1-d$spec,d$sens, pch=1, col=1, cex=symbolsize)
      points(antilogit(-(mu[1]-mu[2]/2)*exp(b/2)), antilogit((mu[1]+mu[2]/2)*exp(-b/2)), col=2, cex=2, pch=16)
      lines(antilogit(-(hh1[,1]-hh1[,2]/2)*exp(b/2)), antilogit((hh1[,1]+hh1[,2]/2)*exp(-b/2)), lty=2, col=2, lwd=1.25)
      lines(antilogit(-(hh2[,1]-hh2[,2]/2)*exp(b/2)), antilogit((hh2[,1]+hh2[,2]/2)*exp(-b/2)), lty=3, col=2, lwd=1.25)
      gemslope=(Sigma[1,1]-Sigma[2,2]/4)/(exp(b)*(Sigma[1,1]+Sigma[2,2]/4))
      gemintercept=(mu[1]+mu[2]/2)*exp(-b/2) - gemslope*(mu[1]-mu[2]/2)*exp(b/2)
      lines(antilogit(-xxx),antilogit(gemintercept+gemslope*xxx), col=3, lwd=2)
      meanauc=calculateauc(antilogit(-xxx),antilogit(gemintercept+gemslope*xxx), 0, 1)
      gg=quantile(aucees,probs=c(0.025,0.975),na.rm=T)
      text(x=0.4,y=0.3,paste("AUC = ",round(meanauc,4)," (95% CI: ", round(gg[1],4)," - ",round(gg[2],4),")",sep=""),adj=0)
      return(aucees)
   }
}

JAGSbivmodel = function(d,x1,x2,niter=50000, thin=10) {
   write.table("Bayesian analysis of the bivariate sens-spec model:", quote=F, row.names=F, col.names=F)
   P1=ncol(x1)
   P2=ncol(x2)
   dd2=list(N=length(d$ncases), ncases=d$ncases, TP=d$TP, ncontrols=d$ncontrols, TN=d$TN, x1=x1, x2=x2, P1=P1, P2=P2,
            zeroos=rep(0.001,(P1+P2+2)), precision=diag(0.001,(P1+P2+2)), R=diag(0.1,2), dfx=2)
   write.table("compiling", quote=F, row.names=F, col.names=F)
   model2=jags.model(textConnection(modelstring2), data=dd2, n.chains = 2, n.adapt=0)
   write.table("burn-in", quote=F, row.names=F, col.names=F)
   update(model2, n.iter=niter)
   write.table("sampling from posteriors", quote=F, row.names=F, col.names=F)
   output.raw2 = coda.samples(model=model2, variable.names=c("mu","Sigma", "beta1", "beta2"), n.iter=niter, thin=thin)
   print(summary(output.raw2))
   return(output.raw2)
}

JAGS_RGmodel = function(d,x1,x2,niter=50000, thin=10) {
   write.table("Bayesian analysis of the bivariate Rutter-Gatsonis SROC model:", quote=F, row.names=F, col.names=F)
   P1=ncol(x1)
   P2=ncol(x2)
   dd3=list(N=length(d$ncases), ncases=d$ncases, TP=d$TP, ncontrols=d$ncontrols, TN=d$TN, x1=x1, x2=x2, P1=P1, P2=P2,
         zeroos=rep(0.001,(P1+P2+2)), precision=diag(0.001,(P1+P2+2)))
   write.table("compiling", quote=F, row.names=F, col.names=F)
   model3=jags.model(textConnection(modelstring3), data=dd3, n.chains = 2, n.adapt=0)
   write.table("burn-in", quote=F, row.names=F, col.names=F)
   update(model3, n.iter=50000)
   write.table("sampling from posteriors", quote=F, row.names=F, col.names=F)
   output.raw3 = coda.samples(model=model3, variable.names=c("mu","Sigma", "beta1", "beta2", "beta"), n.iter=50000, thin=10)
   print(summary(output.raw3))
   return(output.raw3)
}

modelstring2 <- "
   model {
      # likelihood of the data
      mu[1:2] <- modelparms[1:2]
      beta1[1:P1] <- modelparms[3:(P1+2)]
      beta2[1:P2] <- modelparms[(P1+3):(P1+P2+2)]
      for (i in 1:N) {
         intercepts[i,1:2] ~ dmnorm(mu[1:2],invSigma[1:2,1:2])

         lp[i,1] <- inprod(x1[i,1:P1], beta1[1:P1])               # x1 is a N-by-P1 matrix of numeric covariate values of P1 covariates for the N studies
         lp[i,2] <- inprod(x2[i,1:P2], beta2[1:P2])               # x2 is a N-by-P2 matrix of numeric covariate values of P2 covariates for the N studies 
                                                                  # x1 and x2 are covariates for sens and spec respectively, they are not necessarily the same
         sens[i] <- 1/(1+exp(-(intercepts[i,1]+lp[i,1])))
         spec[i] <- 1/(1+exp(-(intercepts[i,2]+lp[i,2])))
         TP[i] ~ dbinom(sens[i],ncases[i])
         TN[i] ~ dbinom(spec[i],ncontrols[i])
      }

      # priors
      modelparms[1:(P1+P2+2)] ~ dmnorm(zeroos[1:(P1+P2+2)],precision[1:(P1+P2+2),1:(P1+P2+2)])
      invSigma[1:2,1:2] ~ dwish(R[1:2,1:2],dfx)

      # functions of model-parameters
      Sigma[1:2,1:2] <- inverse(invSigma[1:2,1:2])

   }"

modelstring3 <- "
   model {
      # likelihood of the data
      mu[1:2] <- modelparms[1:2]
      beta1[1:P1] <- modelparms[3:(P1+2)]
      beta2[1:P2] <- modelparms[(P1+3):(P1+P2+2)]
      for (i in 1:N) {
         intercepts[i,1] ~ dnorm(mu[1],invSigma[1])               # independent intercepts for the theta- and alpha-parameters of the RG SROC model
         intercepts[i,2] ~ dnorm(mu[2],invSigma[2])

         lp[i,1] <- inprod(x1[i,1:P1], beta1[1:P1])               # x1 is a N-by-P1 matrix of numeric covariate values of P1 covariates for the N studies
         lp[i,2] <- inprod(x2[i,1:P2], beta2[1:P2])               # x2 is a N-by-P2 matrix of numeric covariate values of P2 covariates for the N studies 
                                                                  # x1 and x2 are covariates for sens and spec respectively, they are not necessarily the same
         theta[i] <- intercepts[i,1] +lp[i,1]
         alpha[i] <- intercepts[i,2] +lp[i,2]

         sens[i] <- 1/(1+exp(-((theta[i]+0.5*alpha[i])*exp(-beta/2))))      # beta is the fixed shape-parameter of the RG SROC model (may vary between studies too)
         spec[i] <- 1/(1+exp(-((theta[i]-0.5*alpha[i])*exp( beta/2))))
         TP[i] ~ dbinom(sens[i],ncases[i])
         TN[i] ~ dbinom(spec[i],ncontrols[i])
      }

      # priors
      modelparms[1:(P1+P2+2)] ~ dmnorm(zeroos[1:(P1+P2+2)],precision[1:(P1+P2+2),1:(P1+P2+2)])
      beta ~ dnorm(0,0.001)
      invSigma[1] ~ dgamma(0.01,0.1)
      invSigma[2] ~ dgamma(0.01,0.1)

      # functions of model-parameters
      Sigma[1] <- 1/invSigma[1]
      Sigma[2] <- 1/invSigma[2]

   }"

sens_spec_binomial=function(d, x1, x2) {
   write.table(" ",quote=F,row.names=F,col.names=F)
   write.table("Separate analysis logit(Sensitivity) and logit(Specificity) using binomial-normal model:",quote=F,row.names=F,col.names=F)
   m1=glmer(cbind(d$TP,d$FN)~x1+(1|d$studienr), family="binomial")
   write.table("logit(Sensitivity):",quote=F,row.names=F,col.names=F)
   h11=summary(m1)$coefficients[,1:2]
   row.names(h11)=c("Intercept",colnames(x1))
   colnames(h11)=c("beta", "SE")
   print(h11, quote=F)
   write.table(" ",quote=F,row.names=F,col.names=F)
   tau2=as.numeric(summary(m1)$varcor)
   write.table(paste("between-study variance = ",format(tau2,digits=4,nsmall=4),sep=""), quote=F, row.names=F, col.names=F)
   m2=glmer(cbind(d$TN,d$FP)~x1+(1|d$studienr), family="binomial")
   write.table(" ",quote=F,row.names=F,col.names=F)
   write.table("logit(Specificity):",quote=F,row.names=F,col.names=F)
   h12=summary(m2)$coefficients[,1:2]
   row.names(h12)=c("Intercept",colnames(x1))
   colnames(h12)=c("beta", "SE")
   print(h12, quote=F)
   write.table(" ",quote=F,row.names=F,col.names=F)
   tau3=as.numeric(summary(m2)$varcor)
   write.table(paste("between-study variance = ",format(tau3,digits=4,nsmall=4),sep=""), quote=F, row.names=F, col.names=F)
   write.table(" ",quote=F,row.names=F,col.names=F)
   write.table("Joint analysis logit(Sensitivity) and logit(Specificity) using binomial-bivariate normal model:",quote=F,row.names=F,col.names=F)
   xx=x1
   if (length(which(colnames(x2) %in% colnames(x1))) < ncol(x2)) {xx=cbind(xx,x2[,-which(colnames(x2) %in% colnames(x1))])}
   d1=data.frame(d[,c("TP","FN","studienr")], sens=1, spec=0)
   d2=data.frame(d[,c("TN","FP","studienr")], sens=0, spec=1)
   names(d1)=names(d2)=c("pos","neg","studienr","sens","spec")
   d3=rbind(d1,d2)
   xx3=rbind(xx,xx)
#   xx3=xx3[order(d3$studienr,d3$sens),]
#   d3=d3[order(d3$studienr,d3$sens),]    
   m3=glmer(cbind(pos,neg)~0+sens+sens:xx3+spec+spec:xx3+(0+sens+spec|studienr), family="binomial", data=d3)
   h9=summary(m3)$coefficients[,1:2]
   h13=cbind(rbind(h9[1,1:2],h9[3:(ncol(xx3)+2),1:2]),rbind(h9[2,1:2],h9[(ncol(xx3)+3):(2+2*(ncol(xx3))),1:2]))
   row.names(h13)=c("Intercept",colnames(xx3))
   colnames(h13)=rep(c("beta","SE"),2)
   write.table("            logit(Sensitivity)   logit(Specificity):",quote=F,row.names=F,col.names=F)
   print(h13,quote=F)
   h14=as.matrix(summary(m3)$varcor[[1]][1:2,1:2])
   row.names(h14)=colnames(h14)=c("lsens","lspec")
   write.table(" ",quote=F,row.names=F,col.names=F)
   write.table("between-study covariance matrix:", quote=F, row.names=F, col.names=F)
   print(h14,quote=F)
   return(list(h11,tau2,h12,tau3,h13,h14,m1,m2,m3))
   #return(list(h11,tau2,h12,tau3,m1,m2))
}

RGmodelplot=function(mu, Sigma, SE2matrix, beta=0, d, alfa=0.05, plottype=2, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes=0, titel="") {
   npoints=1000
   talfa=qnorm(1-alfa/2)
   symbolsize=1
   if (varyingsymbolsize==1) {symbolsize=minsymbolsize+ deltasymbolsize*(d$ncases+d$ncontrols)/max(d$ncases+d$ncontrols, na.rm=T)}
   hh1=hh2=matrix(c(NA,NA),nrow=1,ncol=2)
   if (class(Sigma)[1]=="matrix") {hh1=confellips(mu,Sigma,alfa,npoints)}
   if (class(SE2matrix)[1]=="matrix") {hh2=confellips(mu,SE2matrix,alfa,npoints)}
   minxx = min(c(hh1[,1],(logit(d$sens)*exp(beta/2)+logit(d$spec)*exp(-beta/2))/2),na.rm=T)
   maxxx = max(c(hh1[,1],(logit(d$sens)*exp(beta/2)+logit(d$spec)*exp(-beta/2))/2),na.rm=T)
   minyy = min(c(hh1[,2],(logit(d$sens)*exp(beta/2)-logit(d$spec)*exp(-beta/2))),na.rm=T)
   maxyy = max(c(hh1[,2],(logit(d$sens)*exp(beta/2)-logit(d$spec)*exp(-beta/2))),na.rm=T)

   if (plottype==1) {
      minxx1 = min(c(logit(d$spec),(hh1[,1]-hh1[,2]/2)*exp( beta/2)), na.rm=T)
      maxxx1 = max(c(logit(d$spec),(hh1[,1]-hh1[,2]/2)*exp( beta/2)), na.rm=T)
      minyy1 = min(c(logit(d$sens),(hh1[,1]+hh1[,2]/2)*exp(-beta/2)), na.rm=T)
      maxyy1 = max(c(logit(d$sens),(hh1[,1]+hh1[,2]/2)*exp(-beta/2)), na.rm=T)
      plot(1,1,xlab="logit(specificity)", ylab="logit(sensitivity)", type="n", main=titel, 
         xlim=c(minxx1,maxxx1), ylim=c(minyy1, maxyy1))
      points(logit(d$spec), logit(d$sens), pch=1, col=1, cex=symbolsize)
      if (SEyes==1) {
         for (j in 1:nrow(d)) {
            lines(c(logit(d$spec[j]),logit(d$spec[j])), 
                  c(logit(d$sens[j])-talfa/sqrt(d$ncases[j]*d$sens[j]*(1-d$sens[j])), logit(d$sens[j])+talfa/sqrt(d$ncases[j]*d$sens[j]*(1-d$sens[j]))), lty=3, col="grey")
            lines(c(logit(d$spec[j])-talfa/sqrt(d$ncontrols[j]*d$spec[j]*(1-d$spec[j])), logit(d$spec[j])+talfa/sqrt(d$ncontrols[j]*d$spec[j]*(1-d$spec[j]))),
                  c(logit(d$sens[j]),logit(d$sens[j])), lty=3, col="grey")
         }
      }
      points((mu[1]-mu[2]/2)*exp(beta/2),(mu[1]+mu[2]/2)*exp(-beta/2), col=2, pch=16, cex=1.5)
      if (class(Sigma)[1]=="matrix") {
         lines((hh1[,1]-hh1[,2]/2)*exp(beta/2), (hh1[,1]+hh1[,2]/2)*exp(-beta/2), lty=3, col=2, lwd=1.25)
         slope=(Sigma[1,1]-Sigma[2,2]/4)/(exp(beta)*(Sigma[1,1]+Sigma[2,2]/4))
         intercept=exp(-beta/2)*(mu[1]+mu[2]/2)-slope*exp(beta/2)*(mu[1]-mu[2]/2)
         xas=seq(minxx1,maxxx1, 0.1)
         lines(xas,xas*slope+intercept, col=3, lwd=2)
      }
      if (class(SE2matrix)[1]=="matrix") {lines((hh2[,1]-hh2[,2]/2)*exp(beta/2), (hh2[,1]+hh2[,2]/2)*exp(-beta/2), lty=2, col=2, lwd=1.25)}
   }
   if (plottype!=1 & plottype!=3) {
      plot(1,1,xlab="specificity",ylab="sensitivity",xlim=c(0,1),ylim=c(0,1),type="n",xaxt="n",main=titel)
      axis(1,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("1.0","0.8","0.6","0.4","0.2","0.0"))
      abline(a=0,b=1,lty=2,col=1)
      points(1-d$spec, d$sens, pch=1, cex=symbolsize)
      if (SEyes==1) {
         for (j in 1:nrow(d)) {
            lines(c(1-d$spec[j],1-d$spec[j]), antilogit(c(logit(d$sens[j])-talfa/sqrt(d$ncases[j]*d$sens[j]*(1-d$sens[j])),logit(d$sens[j])+talfa/sqrt(d$ncases[j]*d$sens[j]*(1-d$sens[j])))), lty=3, col="grey")
            lines(antilogit(c(logit(1-d$spec[j])-talfa/sqrt(d$ncontrols[j]*d$spec[j]*(1-d$spec[j])), logit(1-d$spec[j])+talfa/sqrt(d$ncontrols[j]*d$spec[j]*(1-d$spec[j])))), c(d$sens[j],d$sens[j]), lty=3, col="grey")
         }
      }
      points(antilogit(-(mu[1]-mu[2]/2)*exp(beta/2)), antilogit((mu[1]+mu[2]/2)*exp(-beta/2)), col=2, cex=1.5, lwd=1.5, pch=16)
      if(class(Sigma)[1]=="matrix") {
         lines(antilogit(-(hh1[,1]-hh1[,2]/2)*exp(beta/2)),antilogit((hh1[,1]+hh1[,2]/2)*exp(-beta/2)), lty=3, col=2, lwd=1.25)
         slope=Sigma[1,2]/Sigma[1,1]
         intercept=mu[2]-slope*mu[1]
         xas=seq(minxx, maxxx, 0.1)
         yas=xas*slope+intercept
         xasx=sort(antilogit(-(xas-yas/2)*exp(beta/2)))
         yasy=sort(antilogit((xas+yas/2)*exp(-beta/2)))
         lines(xasx, yasy, col=3, lty=1, lwd=2)
         auc=-1
         auc=calculateauc(xasx, yasy, 0, 1)
         if (auc != -1) {text(x=0.5, y=0.3, paste("AUC = ", round(auc,4)), adj=0)}
      }
      if(class(SE2matrix)[1]=="matrix") {
         hh2=confellips(mu,SE2matrix,alfa,npoints)
         lines(antilogit(-(hh2[,1]-hh2[,2]/2)*exp(beta/2)),antilogit((hh2[,1]+hh2[,2]/2)*exp(-beta/2)), lty=2, col=2, lwd=1.25)
      }
   }
   if (plottype==3) {
      plot(1,1,xlab="theta=(logit(sensitivity)+logit(specificity))/2", ylab="alpha=logit(sensitivity)-logit(specificity)", type="n", main=titel, 
         xlim=c(minxx,maxxx), ylim=c(minyy, maxyy))
      points((logit(d$spec)*exp(-beta/2)+logit(d$sens)*exp(beta/2))/2,-logit(d$spec)*exp(-beta/2)+logit(d$sens)*exp(beta/2), pch=1, col=1, cex=symbolsize)
      if (SEyes==1) {
         for (j in 1:nrow(d)) {
            lines(c((logit(d$sens[j])*exp(beta/2)+logit(d$spec[j])*exp(-beta/2))/2,(logit(d$sens[j])*exp(beta/2)+logit(d$spec[j])*exp(-beta/2))/2), 
                  c((logit(d$sens[j])*exp(beta/2)-logit(d$spec[j])*exp(-beta/2))-talfa*sqrt(exp(beta)/(d$ncases[j]*d$sens[j]*(1-d$sens[j]))+exp(-beta)/(d$ncontrols[j]*d$spec[j]*(1-d$spec[j]))),
                    (logit(d$sens[j])*exp(beta/2)-logit(d$spec[j])*exp(-beta/2))+talfa*sqrt(exp(beta)/(d$ncases[j]*d$sens[j]*(1-d$sens[j]))+exp(-beta)/(d$ncontrols[j]*d$spec[j]*(1-d$spec[j])))) 
                  , lty=3, col="grey")
            lines(c((logit(d$sens[j])*exp(beta/2)+logit(d$spec[j])*exp(-beta/2))/2-talfa*sqrt(exp(beta)/(4*d$ncases[j]*d$sens[j]*(1-d$sens[j]))+exp(-beta)/(4*d$ncontrols[j]*d$spec[j]*(1-d$spec[j]))),
                    (logit(d$sens[j])*exp(beta/2)+logit(d$spec[j])*exp(-beta/2))/2+talfa*sqrt(exp(beta)/(4*d$ncases[j]*d$sens[j]*(1-d$sens[j]))+exp(-beta)/(4*d$ncontrols[j]*d$spec[j]*(1-d$spec[j])))),
                   c(logit(d$sens[j])*exp(beta/2)-logit(d$spec[j])*exp(-beta/2),logit(d$sens[j])*exp(beta/2)-logit(d$spec[j])*exp(-beta/2)), lty=3, col="grey")
         }
      }
      points(mu[1],mu[2],col=2,pch=16,cex=1.5)
      if (class(Sigma)[1]=="matrix") {
         lines(hh1[,1],hh1[,2], lty=3, col=2, lwd=1.25)
         slope=Sigma[1,2]/Sigma[1,1]
         intercept=mu[2]-slope*mu[1]
         xas=seq(minxx,maxxx, 0.1)
         lines(xas,xas*slope+intercept, col=3, lwd=2)
      }
      if (class(SE2matrix)[1]=="matrix") {lines(hh2[,1],hh2[,2], lty=2, col=2, lwd=1.25)
      }
   }
}

bivmodelplot=function(mu, Sigma, SE2matrix, d, alfa=0.05, plottype=2, varyingsymbolsize=1, minsymbolsize=0.5, deltasymbolsize=2.5, SEyes=0, titel="") {
   npoints=1000
   talfa=qnorm(1-alfa/2)
   if (plottype==1) {
      symbolsize=1
      if (varyingsymbolsize==1) {symbolsize=minsymbolsize+ deltasymbolsize*(d$ncases+d$ncontrols)/max(d$ncases+d$ncontrols, na.rm=T)}
      hh1=matrix(c(NA,NA),nrow=1,ncol=2)
      if (class(Sigma)[1]=="matrix") {hh1=confellips(mu,Sigma,alfa,npoints)}
      plot(1,1,xlab="logit(specificity)", ylab="logit(sensitivity)", type="n", main=titel, 
         xlim=c(min(c(hh1[,1],logit(d$spec)),na.rm=T), max(c(hh1[,1],logit(d$spec)),na.rm=T)), ylim=c(min(c(hh1[,2],logit(d$sens)),na.rm=T), max(c(hh1[,2],logit(d$sens)),na.rm=T)))
      points(logit(d$spec),logit(d$sens), pch=1, col=1, cex=symbolsize)
      if (SEyes==1) {
         for (j in 1:nrow(d)) {
            lines(logit(c(d$spec[j],d$spec[j])), c(logit(d$sens[j])-talfa/sqrt(d$ncases[j]*d$sens[j]*(1-d$sens[j])),logit(d$sens[j])+talfa/sqrt(d$ncases[j]*d$sens[j]*(1-d$sens[j]))), lty=3, col="grey")
            lines((c(logit(d$spec[j])-talfa/sqrt(d$ncontrols[j]*d$spec[j]*(1-d$spec[j])), logit(d$spec[j])+talfa/sqrt(d$ncontrols[j]*d$spec[j]*(1-d$spec[j])))), logit(c(d$sens[j],d$sens[j])), lty=3, col="grey")
         }
      }
      points(mu[2],mu[1],col=2,pch=16,cex=1.5)
      if (class(Sigma)[1]=="matrix") {
         lines(hh1[,2],hh1[,1], lty=3, col=2, lwd=1.25)
         slope=Sigma[1,2]/Sigma[2,2]
         intercept=mu[1]-slope*mu[2]
         xas=seq(min(c(hh1[,2],logit(d$spec)),na.rm=T), max(c(hh1[,2],logit(d$spec)),na.rm=T), 0.1)
         lines(xas,xas*slope+intercept, col=3, lwd=2)
      }
      if (class(SE2matrix)[1]=="matrix") {
         hh2=confellips(mu,SE2matrix,alfa,npoints)
         lines(hh2[,2],hh2[,1], lty=2, col=2, lwd=1.25)
      }
   }
   if (plottype!=1) {
      symbolsize=1
      if (varyingsymbolsize==1) {symbolsize=minsymbolsize+ deltasymbolsize*(d$ncases+d$ncontrols)/max(d$ncases+d$ncontrols, na.rm=T)}
      plot(1,1,xlab="specificity",ylab="sensitivity",xlim=c(0,1),ylim=c(0,1),type="n",xaxt="n",main=titel)
      axis(1,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("1.0","0.8","0.6","0.4","0.2","0.0"))
      abline(a=0,b=1,lty=2,col=1)
      points(1-d$spec, d$sens, pch=1, cex=symbolsize)
      if (SEyes==1) {
         for (j in 1:nrow(d)) {
            lines(c(1-d$spec[j],1-d$spec[j]), antilogit(c(logit(d$sens[j])-talfa/sqrt(d$ncases[j]*d$sens[j]*(1-d$sens[j])),logit(d$sens[j])+talfa/sqrt(d$ncases[j]*d$sens[j]*(1-d$sens[j])))), lty=3, col="grey")
            lines(antilogit(c(logit(1-d$spec[j])-talfa/sqrt(d$ncontrols[j]*d$spec[j]*(1-d$spec[j])), logit(1-d$spec[j])+talfa/sqrt(d$ncontrols[j]*d$spec[j]*(1-d$spec[j])))), c(d$sens[j],d$sens[j]), lty=3, col="grey")
         }
      }
      points(antilogit(-mu[2]),antilogit(mu[1]), pch=16, col=2, lwd=1.5, cex=1.5)
      if(class(Sigma)[1]=="matrix") {
         hh1=confellips(mu,Sigma,alfa,npoints)
         lines(antilogit(-hh1[,2]),antilogit(hh1[,1]), lty=3, col=2, lwd=1.25)
         slope=Sigma[1,2]/Sigma[2,2]
         intercept=mu[1]-slope*mu[2]
         xas=seq(min(c(logit(d$spec),hh1[,2])),max(c(logit(d$spec),hh1[,2])),0.1)
         lines(antilogit(-xas), antilogit(xas*slope+intercept), col=3, lty=1, lwd=2)
         auc=-1
         auc=calculateauc(antilogit(-xas),antilogit(xas*slope+intercept), 0, 1)
         if (auc != -1) {text(x=0.5, y=0.3, paste("AUC = ", round(auc,4)), adj=0)}
      }
      if(class(SE2matrix)[1]=="matrix") {
         hh2=confellips(mu,SE2matrix,alfa,npoints)
         lines(antilogit(-hh2[,2]),antilogit(hh2[,1]), lty=2, col=2, lwd=1.25)
      }
   }
}

joint_analysis_approximate=function(d,x1,x2,b)  {
   xx=x1
   if (length(which(colnames(x2) %in% colnames(x1))) < ncol(x2)) {xx=cbind(xx,x2[,-which(colnames(x2) %in% colnames(x1))])}
   i=1
   S=list()
   S[[i]] = matrix(c(1/(d$ncases[i]*d$sens[i]*(1-d$sens[i])),0,0,1/(d$ncontrols[i]*d$spec[i]*(1-d$spec[i]))),nrow=2,ncol=2)
   for (i in 2:nrow(d)) {S[[i]] = matrix(c(1/(d$ncases[i]*d$sens[i]*(1-d$sens[i])),0,0,1/(d$ncontrols[i]*d$spec[i]*(1-d$spec[i]))),nrow=2,ncol=2)}
   m1=mvmeta(cbind(logit(d$sens),logit(d$spec)) ~ xx, S=S)
   write.table(" ",quote=F, row.names=F,col.names=F)
   write.table("joint approximate analysis of logit(sensitivity) and logit(specificity) ",quote=F, row.names=F,col.names=F)
   write.table(" ",quote=F, row.names=F,col.names=F)
   write.table(format(t(c("","lsens  ","  lspec  ")),width=13), quote=F, row.names=F,col.names=F)
   hh11=cbind(m1$coefficients[,1],sqrt(diag(m1$vcov))[1:(ncol(xx)+1)],m1$coefficients[,2],sqrt(diag(m1$vcov))[(ncol(xx)+2):(2*(ncol(xx)+1))])
   colnames(hh11)=rep(c(" beta"," SE"),2)
   row.names(hh11)=c("Intercept",colnames(xx))
   print(format(hh11, digits=4,nsmall=4, width=7), quote=F)
   write.table(" ",quote=F, row.names=F,col.names=F)
   write.table("random-effects covariance-matrix:", quote=F, row.names=F, col.names=F)
   hh12=m1$Psi
   row.names(hh12)=colnames(hh12)=c("lsens","lspec")
   print(hh12)
   write.table(" ",quote=F, row.names=F,col.names=F)
   i=1
   S=list()
   covxy=exp(b)/(2*d$ncases[i]*d$sens[i]*(1-d$sens[i])) - exp(-b)/(2*d$ncontrols[i]*d$spec[i]*(1-d$spec[i]))
   S[[i]] = matrix(c(exp(b)/(4*d$ncases[i]*d$sens[i]*(1-d$sens[i])) + exp(-b)/(4*d$ncontrols[i]*d$spec[i]*(1-d$spec[i])), covxy, covxy,
                     exp(b)/(d$ncases[i]*d$sens[i]*(1-d$sens[i])) + exp(-b)/(d$ncontrols[i]*d$spec[i]*(1-d$spec[i]))),nrow=2,ncol=2)
   for (i in 2:nrow(d)) {
      covxy=exp(b)/(2*d$ncases[i]*d$sens[i]*(1-d$sens[i])) - exp(-b)/(2*d$ncontrols[i]*d$spec[i]*(1-d$spec[i]))
      S[[i]] = matrix(c(exp(b)/(4*d$ncases[i]*d$sens[i]*(1-d$sens[i])) + exp(-b)/(4*d$ncontrols[i]*d$spec[i]*(1-d$spec[i])), covxy,covxy,
                        exp(b)/(d$ncases[i]*d$sens[i]*(1-d$sens[i])) + exp(-b)/(d$ncontrols[i]*d$spec[i]*(1-d$spec[i]))),nrow=2,ncol=2)
   }
   m2=mvmeta(cbind((logit(d$sens)*exp(b/2)+logit(d$spec)*exp(-b/2))/2,logit(d$sens)*exp(b/2)-logit(d$spec)*exp(-b/2)) ~ xx, S=S, bscov="diag")
   write.table(" ",quote=F, row.names=F,col.names=F)
   write.table("joint approximate analysis of (logit(sensitivity)+logit(specificity))/2 and logit(sensitivity)-logit(specificity)",quote=F, row.names=F,col.names=F)
   write.table(" ",quote=F, row.names=F,col.names=F)
   write.table(format(t(c("","(lsens+lspec)/2","lsens-lspec")),width=15), quote=F, row.names=F,col.names=F)
   hh13=cbind(m2$coefficients[,1],sqrt(diag(m2$vcov))[1:(ncol(xx)+1)],m2$coefficients[,2],sqrt(diag(m2$vcov))[(ncol(xx)+2):(2*(ncol(xx)+1))])
   colnames(hh13)=rep(c(" beta"," SE"),2)
   row.names(hh13)=c("Intercept",colnames(xx))
   print(format(hh13, digits=4,nsmall=4, width=7), quote=F)
   write.table(" ",quote=F, row.names=F,col.names=F)
   write.table("random-effects covariance-matrix:", quote=F, row.names=F, col.names=F)
   hh14=m2$Psi
   row.names(hh14)=colnames(hh14)=c("(lsens+lspec)/2","lsens-lspec")
   print(hh14)
   write.table(" ",quote=F, row.names=F,col.names=F)
   return(list(hh11,hh12,hh13,hh14,m1,m2))
}

separate_gem_s_approximate = function(d,x1,x2,b=0) {
   write.table(" ", row.names=F, col.names=F,quote=F)
   write.table("Separate analysis of theta=(logit(Sensitivity)+logit(Specificity))/2 and alpha=logit(Sensitivity)-logit(Specificity)", row.names=F, col.names=F,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   m1a=lm( (logit(d$sens)*exp(b/2)+logit(d$spec)*exp(-b/2))/2 ~ x1)
   m1d=rma(yi=(logit(d$sens)*exp(b/2)+logit(d$spec)*exp(-b/2))/2, vi=exp(b)/(4*d$ncases*d$sens*(1-d$sens))+exp(-b)/(4*d$ncontrols*d$spec*(1-d$spec)), mods=~x1)
   ff3_1=cbind(summary(m1a)$coefficients[,1], summary(m1a)$coefficients[,2],summary(m1d)$beta, summary(m1d)$se)
   row.names(ff3_1)=c("Intercept",colnames(x2))
   colnames(ff3_1)=rep(c("beta","SE"),2)
   write.table(" ", row.names=F, col.names=F,quote=F)
   write.table(paste("(logit(Sensitivity)+logit(Specificity))/2:"), row.names=F, col.names=F,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   write.table(format(t(c("     ","OLS", "Metafor")),width=15), quote=F, col.names=F, row.names=F)
   print(format(ff3_1,digits=4, nsmall=4, width=7, justify="right"),quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   ff4_1=rbind(summary(m1a)$sigma,sqrt(m1d$tau2))#, m1d$se.tau2)
   row.names(ff4_1)=c("OLS","Metafor")
   colnames(ff4_1)=c("residual SD")
   print(ff4_1,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   m2a=lm( logit(d$sens)*exp(b/2)-logit(d$spec)*exp(-b/2) ~ x2)
   m2d=rma(yi=logit(d$sens)*exp(b/2)-logit(d$spec)*exp(-b/2), vi=exp(b)/(d$ncases*d$sens*(1-d$sens))+exp(-b)/(d$ncontrols*d$spec*(1-d$spec)), mods=~x2)
   ff3_2=cbind(summary(m2a)$coefficients[,1], summary(m2a)$coefficients[,2],summary(m2d)$beta, summary(m2d)$se) 
   row.names(ff3_2)=c("Intercept",colnames(x2))
   colnames(ff3_2)=rep(c("beta","SE"),2)
   write.table("logit(Sensitivity)-logit(Specificity):", row.names=F, col.names=F,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   write.table(format(t(c("     ","OLS", "Metafor")),width=15), quote=F, col.names=F, row.names=F)
   print(format(ff3_2,digits=4, nsmall=4, width=7, justify="right"),quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   ff4_2=rbind(summary(m2a)$sigma,sqrt(m2d$tau2))#, m2d$se.tau2)
   row.names(ff4_2)=c("OLS","Metafor")
   colnames(ff4_2)=c("residual SD")
   print(ff4_2,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   ff5=rbind(cor.test(resid(m1a),resid(m2a))$estimate,
              cor.test(resid(m1d),resid(m2d))$estimate)
   row.names(ff5)=c("OLS","Metafor")
   colnames(ff5)="correlation (lsens+lspec)/2, lsens-lspec"
   print(format(ff5,justify="centre"), quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   return(list(ff3_1, ff4_1, ff3_2, ff4_2, ff5, m1a, m1d, m2a, m2d))
}

separate_sensspec_approximate = function(d,x1,x2) {
   write.table(" ", row.names=F, col.names=F,quote=F)
   write.table("Separate analysis of logit(Sensitivity) and logit(Specificity)", row.names=F, col.names=F,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   m1a=lm( logit(d$sens) ~ x1)
   m1d=rma(yi=logit(d$sens), vi=1/(d$ncases*d$sens*(1-d$sens)), mods=~x1)
   ff3_1=cbind(summary(m1a)$coefficients[,1], summary(m1a)$coefficients[,2],summary(m1d)$beta, summary(m1d)$se)
   row.names(ff3_1)=c("Intercept",colnames(x2))
   colnames(ff3_1)=rep(c("beta","SE"),2)
   write.table(" ", row.names=F, col.names=F,quote=F)
   write.table(paste("logit(Sensitivity):"), row.names=F, col.names=F,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   write.table(format(t(c("     ","OLS", "Metafor")),width=15), quote=F, col.names=F, row.names=F)
   print(format(ff3_1,digits=4, nsmall=4, width=7, justify="right"),quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   ff4_1=rbind(summary(m1a)$sigma,sqrt(m1d$tau2))#, m1d$se.tau2)
   row.names(ff4_1)=c("OLS","Metafor")
   colnames(ff4_1)=c("residual SD")
   print(ff4_1,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   m2a=lm( logit(d$spec) ~ x2)
   m2d=rma(yi=logit(d$spec), vi=1/(d$ncontrols*d$spec*(1-d$spec)), mods=~x2)
   ff3_2=cbind(summary(m2a)$coefficients[,1], summary(m2a)$coefficients[,2],summary(m2d)$beta, summary(m2d)$se)
   row.names(ff3_2)=c("Intercept",colnames(x2))
   colnames(ff3_2)=rep(c("beta","SE"),2)
   write.table("logit(Specificity):", row.names=F, col.names=F,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   write.table(format(t(c("     ","OLS", "Metafor")),width=15), quote=F, col.names=F, row.names=F)
   print(format(ff3_2,digits=4, nsmall=4, width=7, justify="right"),quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   ff4_2=rbind(summary(m2a)$sigma,sqrt(m2d$tau2))#, m2d$se.tau2)
   row.names(ff4_2)=c("OLS","Metafor")
   colnames(ff4_2)=c("residual SD")
   print(ff4_2,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   ff5=rbind(cor.test(resid(m1a),resid(m2a))$estimate,
              cor.test(resid(m1d),resid(m2d))$estimate)
   row.names(ff5)=c("OLS","Metafor")
   colnames(ff5)="correlation logit-sensitivity, logit-specificity"
   print(format(ff5,justify="centre"), quote=F)
   write.table(" ", row.names=F, col.names=F,quote=F)
   return(list(ff3_1, ff4_1, ff3_2, ff4_2, ff5, m1a, m1d, m2a, m2d))
}

descriptiveplot=function(d, delta=0.5, varyingsymbolsize=1, deltasymbolsize, minsymbolsize, SEyes=0, titel="") {
   if (varyingsymbolsize==1) {SEyes=0}
   symbolsize=1
   if (varyingsymbolsize==1) {symbolsize=minsymbolsize+ deltasymbolsize*(d$ncases+d$ncontrols)/max(d$ncases+d$ncontrols, na.rm=T)}
   plot(1,1, xlim=c(0,1), ylim=c(0,1), xlab="Specificity", ylab="Sensitivity", xaxt="n", type="n")
   axis(1,at=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("1.0","0.8","0.6","0.4","0.2","0.0"))
   abline(a=0, b=1, col=1, lty=3, lwd=1)
   points(1-d$spec, d$sens, pch=1, col=1, cex=symbolsize)
   if (SEyes==1) {
      for (j in 1:nrow(d)) {
         lines(c(1-d$spec[j],1-d$spec[j]), antilogit(c(logit(d$sens[j])-talfa/sqrt(d$ncases[j]*d$sens[j]*(1-d$sens[j])),logit(d$sens[j])+talfa/sqrt(d$ncases[j]*d$sens[j]*(1-d$sens[j])))), lty=3, col="grey")
         lines(antilogit(c(logit(1-d$spec[j])-talfa/sqrt(d$ncontrols[j]*d$spec[j]*(1-d$spec[j])), logit(1-d$spec[j])+talfa/sqrt(d$ncontrols[j]*d$spec[j]*(1-d$spec[j])))), c(d$sens[j],d$sens[j]), lty=3, col="grey")
      }
   }
}

summarizingcovariates=function(n,x1,x2) {
   check=0
   if(class(x1)[1]!="data.frame" & class(x1)[1]!="matrix") {
      check=-1
      return (list(x1,x2,check))
   }
   if(class(x2)[1]!="data.frame" & class(x2)[1]!="matrix") {
      check=-1
      return (list(x1,x2,check))
   }
   if (nrow(x1)!=n | nrow(x2)!=n) {
      check=-2
      return(list(x1,x2,check))
   }
   for (j in 1:ncol(x1)) {
      if(class(x1[,j])=="factor") {x1[,j]=as.numeric(x1[,j])-1}
      if(class(x1[,j])=="character") {x1[,j]=as.numeric(as.factor(x1[,j]))-1}
   }
   for (j in 1:ncol(x2)) {
      if(class(x2[,j])=="factor") {x2[,j]=as.numeric(x2[,j])-1}
      if(class(x2[,j])=="character") {x2[,j]=as.numeric(as.factor(x2[,j]))-1}
   }
   write.table("",quote=F,row.names=F,col.names=F)
   write.table("  Covariates for logit(sensitivity):", quote=F, row.names=F, col.names=F)
   write.table("",quote=F,row.names=F,col.names=F)
   ffmin1=cbind(apply(x1, 2,mean,na.rm=T), apply(x1,2, sd,na.rm=T), apply(x1,2,min,na.rm=T), apply(x1,2,max,na.rm=T))
   colnames(ffmin1)=c("mean","sd","min","max")
   print(ffmin1, quote=F)
   write.table("",quote=F,row.names=F,col.names=F)
   write.table("  Covariates for logit(specificity):", quote=F, row.names=F, col.names=F)
   write.table("",quote=F,row.names=F,col.names=F)
   ffmin2=cbind(apply(x2, 2, mean, na.rm=T), apply(x2, 2, sd, na.rm=T), apply(x2, 2, min, na.rm=T), apply(x2, 2, max, na.rm=T))
   colnames(ffmin2)=c("mean","sd","min","max")
   print(ffmin2, quote=F)
   write.table("",quote=F,row.names=F,col.names=F)
   x1=scale(x1)
   x2=scale(x2)
   return(list(x1,x2,check))
}

summarizingsensspec=function(d, beta) {
   write.table("",quote=F,row.names=F,col.names=F)
   write.table(paste("  Number of studies = ", format(nrow(d),digits=0,nsmall=0)), quote=F, row.names=F, col.names=F)
   write.table("",quote=F,row.names=F,col.names=F)
   write.table(paste("  Mean number of cases = ", format(mean(d$ncases,na.rm=T),digits=4,nsmall=4)," varying from ", format(min(d$ncases,na.rm=T),digits=0,nsmall=0)," to ",
       format(max(d$ncases,na.rm=T),digits=0,nsmall=0),sep=""), quote=F, row.names=F, col.names=F)
   write.table(paste("  Mean sensitivity = ", format(mean(d$sens,na.rm=T),digits=4,nsmall=4)," varying from ", format(min(d$sens,na.rm=T),digits=4,nsmall=4)," to ",
       format(max(d$sens,na.rm=T),digits=0,nsmall=0),sep=""), quote=F, row.names=F, col.names=F)
   write.table("",quote=F,row.names=F,col.names=F)
   write.table(paste("  Mean number of controls = ", format(mean(d$ncontrols,na.rm=T),digits=4,nsmall=4)," varying from ", format(min(d$ncontrols,na.rm=T),digits=0,nsmall=0)," to ",
       format(max(d$ncontrols,na.rm=T),digits=0,nsmall=0),sep=""), quote=F, row.names=F, col.names=F)
   write.table(paste("  Mean specificity = ", format(mean(d$spec,na.rm=T),digits=4,nsmall=4)," varying from ", format(min(d$spec,na.rm=T),digits=4,nsmall=4)," to ",
       format(max(d$spec,na.rm=T),digits=0,nsmall=0),sep=""), quote=F, row.names=F, col.names=F)
   write.table(paste("  Spearman correlation sensitivity with specificity = ", format(cor(d$sens,d$spec,method="spearman"),digits=4,nsmall=4)), quote=F, row.names=F, col.names=F)
   write.table("",quote=F,row.names=F,col.names=F)
   write.table("",quote=F,row.names=F,col.names=F)
   write.table(paste("  Mean logit(sensitivity) = ", format(mean(logit(d$sens),na.rm=T),digits=4,nsmall=4)," (SE ",format(sd(logit(d$sens),na.rm=T)/sqrt(nrow(d)),digits=4,nsmall=4),")",
      " SD = ", format(sd(logit(d$sens),na.rm=T), digits=4,nsmall=4), sep=""), quote=F, row.names=F, col.names=F)
   write.table(paste("  Mean logit(specificity) = ", format(mean(logit(d$spec),na.rm=T),digits=4,nsmall=4)," (SE ",format(sd(logit(d$spec),na.rm=T)/sqrt(nrow(d)),digits=4,nsmall=4),")",
      " SD = ", format(sd(logit(d$spec),na.rm=T), digits=4,nsmall=4), sep=""), quote=F, row.names=F, col.names=F)
   write.table(paste("  Pearson correlation logit(sensitivity) with logit(specificity) = ", format(cor(logit(d$sens),logit(d$spec),method="pearson"),digits=4,nsmall=4)), quote=F, row.names=F, col.names=F)
   write.table("",quote=F,row.names=F,col.names=F)
   write.table("",quote=F,row.names=F,col.names=F)
   xxxx=(logit(d$sens)*exp(beta/2) + logit(d$spec)*exp(-beta/2))/2
   yyyy= logit(d$sens)*exp(beta/2) - logit(d$spec)*exp(-beta/2)
   write.table(paste("  Mean theta=(logit(sensitivity)+logit(specificity))/2 = ", format(mean(xxxx,na.rm=T),digits=4,nsmall=4)," (SE ",format(sd(xxxx,na.rm=T)/sqrt(nrow(d)),digits=4,nsmall=4),")",
      " SD = ", format(sd(xxxx,na.rm=T), digits=4,nsmall=4), sep=""), quote=F, row.names=F, col.names=F)
   write.table(paste("  Mean alpha=logit(sensitivity)-logit(specificity) = ", format(mean(yyyy,na.rm=T),digits=4,nsmall=4)," (SE ",format(sd(yyyy,na.rm=T)/sqrt(nrow(d)),digits=4,nsmall=4),")",
      " SD = ", format(sd(yyyy,na.rm=T), digits=4,nsmall=4), sep=""), quote=F, row.names=F, col.names=F)
   write.table(paste("  Pearson correlation theta with alpha = ", format(cor(xxxx,yyyy, method="pearson"),digits=4,nsmall=4)), quote=F, row.names=F, col.names=F)
   write.table("",quote=F,row.names=F,col.names=F)
}

checkdata1 = function(d,delta) {
   #print(dim(d))
   check=0
   if (class(d)[1] != "data.frame" & class(d)[1] != "matrix") {
      check=-1
      return(list(d,check))
   }
   d=as.data.frame(d)
   if (("TN" %in% names(d))==F & "tn" %in% names(d)) {d$TN=d$tn}
   if (("TN" %in% names(d))==F & "Tn" %in% names(d)) {d$TN=d$Tn}
   if (("FN" %in% names(d))==F & "fn" %in% names(d)) {d$FN=d$fn}
   if (("FN" %in% names(d))==F & "Fn" %in% names(d)) {d$FN=d$Fn}
   if (("TP" %in% names(d))==F & "tp" %in% names(d)) {d$TP=d$tp}
   if (("TP" %in% names(d))==F & "Tp" %in% names(d)) {d$TP=d$Tp}
   if (("FP" %in% names(d))==F & "fp" %in% names(d)) {d$FP=d$fp}
   if (("FP" %in% names(d))==F & "Fp" %in% names(d)) {d$FP=d$Fp}
   if (("sens" %in% names(d))==F & "Sens" %in% names(d)) {d$sens=d$Sens}
   if (("sens" %in% names(d))==F & "SENS" %in% names(d)) {d$sens=d$SENS}
   if (("sens" %in% names(d))==F & "Sensitivity" %in% names(d)) {d$sens=d$Sensitivity}
   if (("sens" %in% names(d))==F & "sensitivity" %in% names(d)) {d$sens=d$sensitivity}
   if (("spec" %in% names(d))==F & "Spec" %in% names(d)) {d$spec=d$Spec}
   if (("spec" %in% names(d))==F & "SPEC" %in% names(d)) {d$spec=d$SPEC}
   if (("spec" %in% names(d))==F & "Specificity" %in% names(d)) {d$spec=d$Specificity}
   if (("spec" %in% names(d))==F & "specificity" %in% names(d)) {d$spec=d$specificity}
   if (("ncases" %in% names(d))==F & "NCASES" %in% names(d)) {d$ncases=d$NCASES}
   if (("ncases" %in% names(d))==F & "Ncases" %in% names(d)) {d$ncases=d$Ncases}
   if (("ncases" %in% names(d))==F & "NofCases" %in% names(d)) {d$ncases=d$NofCases}
   if (("ncases" %in% names(d))==F & "Cases" %in% names(d)) {d$ncases=d$Cases}
   if (("ncases" %in% names(d))==F & "cases" %in% names(d)) {d$ncases=d$cases}
   if (("ncontrols" %in% names(d))==F & "NCONTROLS" %in% names(d)) {d$ncontrols=d$NCONTROLS}
   if (("ncontrols" %in% names(d))==F & "Ncontrols" %in% names(d)) {d$ncontrols=d$Ncontrols}
   if (("ncontrols" %in% names(d))==F & "NofControls" %in% names(d)) {d$ncontrols=d$NofControls}
   if (("ncontrols" %in% names(d))==F & "Controls" %in% names(d)) {d$ncontrols=d$Controls}
   if (("ncontrols" %in% names(d))==F & "controls" %in% names(d)) {d$ncontrols=d$controls}
   if ((sum(c("TN","ncontrols","TP","ncases") %in% names(d))==4)) {
      if (class(d$TN*1)!="numeric"|class(d$ncontrols*1)!="numeric"|class(d$TP*1)!="numeric"|class(d$ncases*1)!="numeric") {
         check=-2
         return(list(d,check))
      }
      if (sum(d$TN<0|d$ncontrols<0|d$TP<0|d$ncases<0)!=0) {
         check=-3
         return(list(d,check))
      }   
      d$FP=d$ncontrols - d$TN
      d$FN=d$ncases - d$TP
   }
   if ((sum(c("spec","ncontrols","sens","ncases") %in% names(d))==4)) {
      if (class(d$spec*1)!="numeric"|class(d$ncontrols*1)!="numeric"|class(d$sens*1)!="numeric"|class(d$ncases*1)!="numeric") {
         check=-2
         return(list(d,check))
      }
      if (sum(d$spec<0|d$spec>1|d$ncontrols<0|d$sens<0|d$sens>1|d$ncases<0)!=0) {
         check=-4
         return(list(d,check))
      }
      d$TN=round(d$spec*d$ncontrols,0)
      d$TN=round(d$sens*d$ncases,0)
      d$FP=d$ncontrols - d$TN
      d$FN=d$ncases - d$TP
   }
   if ((sum(c("TN","TP","FN","FP") %in% names(d))==4)) {
      if (class(d$TN*1)!="numeric"|class(d$FP*1)!="numeric"|class(d$FN*1)!="numeric"|class(d$TP*1)!="numeric") {
         check=-2
         return(list(d,check))
      }
      if (sum(d$TN<0|d$FP<0|d$TP<0|d$FN<0)!=0) {
         check=-3
         return(list(d,check))
      }
      d$ncontrols = d$TN + d$FP
      d$ncases = d$TP + d$FN
      d$sens=(d$TP+delta)/(d$TP+d$FN+2*delta)
      d$spec=(d$TN+delta)/(d$TN+d$FP+2*delta)
      d$studienr=1:nrow(d)
      check=1
   }
   return(list(d,check))
}

initial_estimate_beta=function(d,laagste, hoogste) {
   b=seq(laagste, hoogste, 0.001)
   correlaties=sapply(b, correl, xx=logit(d$sens), yy=logit(d$spec))
   plot(b, correlaties, xlab="RG SROC beta parameter", ylab="correlation between theta and alpha", type="l", col=2)
   abline(h=0)
   abline(v=0)
   optimalb=b[which((abs(correlaties)-0) == min((abs(correlaties)-0)))]
   if (length(optimalb)>1) {optimalb=optimalb[which((abs(optimalb)-0)==min((abs(optimalb)-0)))]}
   abline(v=optimalb, col=3)
   text(x=0,y=0.1, paste("optimal initial beta =",round(optimalb,4)))
   return(optimalb)
}

calculateauc=function(x,y, onder=0, boven=1) {
   if (onder<0 | onder>1) {onder=0}
   if (boven<0 | boven>1) {boven=1}
   x=x[x>=0 & x<=1]
   y=y[y>=0 & y<=1]
   auc=-1
   if (length(x)>1 & length(x)==length(y) & class(x)=="numeric" & class(y)=="numeric") {
      x=c(onder,x,boven)
      y=c(onder,y,boven)
      x=x[x>=onder & x<=boven]
      y=y[y>=onder & y<=boven]
      x=sort(x)
      y=sort(y)
      auc=sum(diff(x)*(y[2:length(y)]+y[1:(length(y)-1)])/2)
   }
   return(auc)
}

correl=function(b,xx,yy){cor((xx*exp(b/2) + yy*exp(-b/2))/2, xx*exp(b/2) - yy*exp(-b/2))}
logit=function(x) {log(x/(1-x))}
antilogit=function(x) {1/(1+exp(-x))}

confellips=function(mu,sigma,alfa,npoints) {
   es <- eigen(sigma)
   e1 <- es$vec %*% diag(sqrt(es$val))
   r1 <- sqrt(qchisq(1 - alfa, 2))
   theta <- seq(0, 2 * pi, len = npoints)
   v1 <- cbind(r1 * cos(theta), r1 * sin(theta))
   pts = t(mu - (e1 %*% t(v1)))
}


