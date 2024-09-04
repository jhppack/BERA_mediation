rm(list=ls())
## This R code demonstrates how to fit BERA.med to data generated from Simulation Setting 1. At the end of the code, 
## BERA.med.summary is used to summarize the posterior samples, providing the posterior mean and a 95% credible interval for the parameters of interest.

library(mvtnorm)
library(tmvtnorm)
library(MCMCpack)

set.seed(101)

# sample size
ns<-100

# number of iterations for the MCMC algorithm
nitr<-30000

### Specification of the true parameters for the simulation study
### please refer to Equation (5) in the manuscript.
tr.W<-matrix(c(0.9, 0.5, rep(0, 5), rep(0, 5), 0.5, 0.5), 7, 2)
# intercept term for M[i,]'
tr.bt02<-c(1, 1)									
# regression coefficients of F[i,]'(=W'X[i,]') for mediator vector M[i,]'
tr.A2<-matrix(c(2, 0, 0, 2), 2, 2)						
# residual covariance matrix for M[i,]'
tr.Sig2<-matrix(c(1,0.3,0.3,1), 2,2)
		
# intercept term for Z[i,]'
tr.bt03<-c(1, 0, 2, -1)								
# regression coefficients of F[i,]'(=W'X[i,]') for latent outcome 
# vector Z[i,]'
tr.A3<-matrix(c(1, 0, 2, 1.5, 1, 1.5, 2, 0), 2, 4)			
# regression coefficients of M[i,]' for Z[i,]
tr.A4<-matrix(c(0, 2, 1, 1, 0.5, 1.5, 0, 2), 2, 4)			
# residual covariance matrix for Z[i,]
tr.Sig3<-matrix(0.3, 4,4)							
diag(tr.Sig3)<-c(0.8, 0.8, 1, 1)
# covariance matrix for covariates X[i,]
tr.Sig.X<-matrix(0.1, 7, 7)							
tr.Sig.X[1:5, 1:5]<-0.3
tr.Sig.X[6:7, 6:7]<-0.3
diag(tr.Sig.X)<-1

## Generate the covariate matrix X
X<-rmvnorm(ns, mean=rep(0, nrow(tr.Sig.X)), sigma=tr.Sig.X)

## Degrees of freedom for outcomes Z and mediators M
# for continuous outcomes
nu.t1<-5			# df for the t-distribution
# for ordinal outcomes
nu.t2<-7.3			# df for the t-distribution
# for mediators
nu.t3<-3


sig2.1<-1
sig2.2<-pi^2*(nu.t2-2)/(3*nu.t2)

## Generate outcomes and mediators conditional on covariate matrix X from multivariate t distribution
z<-matrix(0, ns, ncol(tr.A3))
m<-matrix(0, ns, ncol(tr.A2))
for (i in 1:ns){
m[i,]<-rmvt(1, delta=tr.bt02+X[i,]%*%tr.W%*%tr.A2, sigma=tr.Sig2, df=nu.t3)

phi.i1<-rgamma(1, 5/2, 5/2)
phi.i2<-rgamma(1, nu.t2/2, nu.t2/2)

# using a correlation matrix R (refer to page 6 of Dunson and O'brien 
# paper)
tmp.z<-rmvnorm(1, mean=rep(0, ncol(tr.A4)), sigma=cov2cor(tr.Sig3))	
z[i,1]<-(tr.bt03+X[i,]%*%tr.W%*%tr.A3+m[i,]%*%tr.A4)[1]+tmp.z[1]*sqrt(tr.Sig3[1,1]/phi.i1)
z[i,2]<-(tr.bt03+X[i,]%*%tr.W%*%tr.A3+m[i,]%*%tr.A4)[2]+tmp.z[2]*sqrt(tr.Sig3[2,2]/phi.i1)
z[i,3]<-(tr.bt03+X[i,]%*%tr.W%*%tr.A3+m[i,]%*%tr.A4)[3]+log(pt(tmp.z[3]/sqrt(phi.i2), df=nu.t2)/(1-pt(tmp.z[3]/sqrt(phi.i2), df=nu.t2)))
z[i,4]<-(tr.bt03+X[i,]%*%tr.W%*%tr.A3+m[i,]%*%tr.A4)[4]+log(pt(tmp.z[4]/sqrt(phi.i2), df=nu.t2)/(1-pt(tmp.z[4]/sqrt(phi.i2), df=nu.t2)))
}

## Note that the latent outcomes z[,1:2] are identical to the observed 
## continuous outcomes ys[,1:2], 
y<-z
## Ordinal variables for the third and fourth outcomes
tr.gam1<-c(0, 3, 6)		# true cutpoints for the first ordinal variable
tr.gam2<-c(0, 2, 4)		# true cutpoints for the second ordinal variable
y[,3]<-ifelse(z[,3]>tr.gam1[3], 4, ifelse(z[,3]>tr.gam1[2],3,ifelse(z[,3]>tr.gam1[1], 2, 1)))
y[,4]<-ifelse(z[,4]>tr.gam2[3], 4, ifelse(z[,4]>tr.gam2[2],3,ifelse(z[,4]>tr.gam2[1], 2, 1)))

## Matrix indicating which predictors contribute to latent components
## the first latent component F1 is defined as a liner combination of the first 5 independent variables.
## the second latent component F2 is defined as a liner combination of the last 2 independent variables.
indx.W<-matrix(c(rep(1, 5), rep(0, 2), rep(0, 5), rep(1, 2)), ncol(X), 2)


## Fit BERA-medation to the simulated data.
source("D:/Project/Minjung Kyung/Bayesian ERA/mediation/function for BERA_mediation.R")
rest.out<-BERA.med(y, X, m, indx.W, 2, nitr=30000)

## Summarize posterior samples of parameters of interest.
BERA.med.summary(rest.out, indx.W, burnin=2000, thinning=10)
