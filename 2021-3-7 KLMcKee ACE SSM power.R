
# -------------------------------------------------------------------------
# Simulations for "Hierarchical Biometrical Genetic Analysis of Longitudinal Dynamics"
# Author: Kevin L. McKee
# Contact: kmckee90@gmail.com

# This script runs comparisons for a simple mixed-effects AR1 between two-step analysis with maximum likelihood and hierarchical estimation in Stan. 
# Note: this will take a very long time if not run in parallel chunks.
# -------------------------------------------------------------------------


# Parameters --------------------------------------------------------------
NCores<-4 #Used for multiple MCMC chains
nTrial<-100 #100 iterations is expected to take multiple days to finish (CPU: intel i7 3770)

#Vectors of sample conditions
t.v<-c(50,200)
N.v<-c(50,200)

#Random AR parameter distribution:
VA1<-.8; VC1<-0; VE1<-1-(VA1+VC1)
meanBeta<-.5; varBeta<-.1^2



# Packages and functions --------------------------------------------------
library(xtable)
library(MASS)
library(psych)
library(OpenMx)
library(cmdstanr)
library(foreach)
library(doSNOW)
library(Compositional)


hssmAR1.cmd<-cmdstan_model("ar_ACE.stan")

getMAP<-function(x){
  dens<-mkde(x)
  x.opt<-x[which(dens==max(dens)),]
  return(x.opt)
}

#Simulate AR(1) data by convolution:
rAR<-function(N, b) convolve(exp(log(b)*0:(N-1) ), rnorm(N))



# Simulation --------------------------------------------------------------
#Generate
rBeta<-matrix(c(
  VA1+VC1+VE1, VA1+VC1, 0, 0,
  VA1+VC1, VA1+VC1+VE1, 0, 0,
  0, 0, VA1+VC1+VE1, .5*VA1+VC1,
  0, 0, .5*VA1+VC1, VA1+VC1+VE1
),4,4,byrow=T)

#Set up OpenMx model for constrained estimation of heritability in 2-step approach.
datNames<-c("t1mz","t2mz","t1dz","t2dz")

#Output objects
output.mu<-array(NA, dim=c(3, nTrial, length(t.v),length(N.v), 2))
output.sd<-array(NA, dim=c(3, nTrial, length(t.v),length(N.v), 2))
dimnames(output.mu)<-dimnames(output.sd)<-list(c("a","c","e"), NULL, t.v,  N.v, c("ML","MCMC"))
for(N in N.v)for(t in t.v){
  for(i in 1:nTrial){
    print(c(t,N,i))
    #Generate time series from the random AR coefs:
    betas<-mvrnorm(N, mu=rep(meanBeta, 4), Sigma = rBeta*varBeta, empirical=F)
    t1mz<-vapply(betas[,1], function(x) rAR(t, x), rep(0,t))
    t2mz<-vapply(betas[,2], function(x) rAR(t, x), rep(0,t))
    t1dz<-vapply(betas[,3], function(x) rAR(t, x), rep(0,t))
    t2dz<-vapply(betas[,4], function(x) rAR(t, x), rep(0,t))
    dat<-list("N"=N, "T"=t, "MZ1"=t1mz, "MZ2"=t2mz, "DZ1"=t1dz, "DZ2"=t2dz)
    
    
    #Estimate AR1 coefficient by ML:
    betas.est<-
      cbind(apply(t1mz, 2, function(x)ar(x,order.max=1,aic=F, method="ml" )$ar),
            apply(t2mz, 2, function(x)ar(x,order.max=1,aic=F, method="ml" )$ar),
            apply(t1dz, 2, function(x)ar(x,order.max=1,aic=F, method="ml" )$ar),
            apply(t2dz, 2, function(x)ar(x,order.max=1,aic=F, method="ml" )$ar))
    colnames(betas.est)<-datNames
    
    #Estimate heritabilities in OpenMx
    model<-mxModel("ACE",
                   mxMatrix("Full",nrow=1,ncol=3,labels=c("a","c","e"),values=c(.01,.01,.1), free=T, lbound=1e-8, name="pars"),
                   mxAlgebra(rbind(cbind(a+c+e,a+c,0,0),
                                   cbind(a+c,a+c+e,0,0),
                                   cbind(0,0, a+c+e,.5*a+c),
                                   cbind(0,0, .5*a+c,a+c+e)), name="R"),
                   mxMatrix("Full",nrow=1,ncol=4,labels="m_beta",values=.5,free=T, name="M"),
                   mxData(betas.est, type="raw"),
                   mxExpectationNormal(covariance="R",means="M", dimnames=c(datNames)),
                   mxFitFunctionML())
    
    model<-mxTryHard(model, silent=T)
    output.mu[,i,paste0(t),paste0(N),"ML"]<-coef(model)[1:3]/sum(coef(model)[1:3])
    output.sd[,i,paste0(t),paste0(N),"ML"]<-model$output$standardErrors[1:3]/sum(coef(model)[1:3]) #Unused, will be inaccurate
    
    
    #Hierarchical modeling in cmdstan:
    #Run a single warmup that is shared across that condition for efficiency:
    if(i==1){
      fit <- hssmAR1.cmd$sample(data = dat,
                                chains = 1,
                                parallel_chains = 1,
                                seed=123,
                                adapt_delta= .9,
                                max_treedepth=6,
                                iter_warmup = 500,
                                iter_sampling = 50,
                                save_warmup = TRUE
      )
      d<-fit$sampler_diagnostics()
      ss<-d[1,1,2][[1]]
      im<-fit$inv_metric()
      ws<-list("stepsize"=ss, "invmet"=diag(im[[1]]) )
    }
    #Actual estimation run
    fit2 <- hssmAR1.cmd$sample(data = dat, 
                               chains = 4,
                               parallel_chains = 4,
                               adapt_engaged = TRUE,
                               adapt_delta = .9,
                               max_treedepth=6,
                               iter_warmup = 25,
                               iter_sampling = 500,
                               refresh = 100,
                               metric="diag_e",
                               inv_metric = ws$invmet,
                               step_size = ws$stepsize
    )
    dr.raw<-fit2$draws()
    dr<-list("A"=c(dr.raw[,,c("A")]),
             "C"=c(dr.raw[,,c("C")]),
             "E"=c(dr.raw[,,c("E")]))
    
    output.mu[,i,paste0(t),paste0(N),"MCMC"]<-getMAP(simplify2array(dr[c("A","C","E")]))
    output.sd[,i,paste0(t),paste0(N),"MCMC"]<-sapply(dr, sd)
    saveRDS(list("mu"=output.mu, "sd"=output.sd), "output.RDS")
    
    #Monitor progress:
    pairs.panels(as.data.frame(dr), pch=".", ellipses = FALSE, smooth=FALSE, breaks=24, hist.col="gray", cex.cor = .5, rug=FALSE )
    par(mfrow=c(2,2), mai=c(.75,.75,.25,.1))
    for(j in 1:2)for(n in 1:2) {
      plot(output.mu[1,,j,n,], xlim=c(0,1), ylim=c(0,1), main=paste0("T=",t.v[i],", N=",N.v[n]), xlab="2-Step ML Estimates", ylab="Hierarchical Estimates", pch=16, cex.lab=1.5 )
      abline(0,1, lty=1);abline(h=.8, v=.8, lty=c(2,2))
    } 
    
    gc()
  }
  
}


#Presenting results:
o<-readRDS("output.RDS")
output.mu<-o$mu
output.sd<-o$sd

#Subset outputs 
output.mu.summary<-apply(output.mu, (1:5)[-2], mean)[1,,,]
output.sd.summary<-apply(output.mu, (1:5)[-2], sd)[1,,,]
output.power.summary<- 1-pnorm(qnorm(.975, 0, output.sd.summary), output.mu.summary, output.sd.summary)

#Table:
res.mu<-round(output.mu.summary, 2)
res.ci<-round(1.96*output.sd.summary, 2)
res.f<-array(paste0(res.mu, " (",res.mu - res.ci, ", ",res.mu + res.ci,")" ), dim=c(2,2,2), dimnames=dimnames(output.mu.summary)[2:4])
res.f.combined<-cbind(rbind(res.f[,1,],res.f[,2,]), rbind(round(output.power.summary[,1,], 2), round(output.power.summary[,2,], 2)))
xtable(res.f.combined)
apply(output.mu[1,,2,1,],2,sd)


#Plot:
par(mfrow=c(2,2), mai=c(.75,.75,.25,.1))
for(t in 1:2)for(n in 1:2) {
  plot(output.mu[1,,t,n,], xlim=c(0,1), ylim=c(0,1), cex=.75, main=paste0("T=",t.v[t],", N=",N.v[n]), xlab="2-Step ML Estimates", ylab="Hierarchical Estimates", pch=16, cex.lab=1.5 )
  abline(0,1, lty=1);abline(h=.8, v=.8, lty=c(2,2))
}


