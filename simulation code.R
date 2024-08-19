rm(list=ls())
library(mvtnorm)
library(tmvtnorm)
library(MCMCpack)

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
ys<-z
## Ordinal variables for the third and fourth outcomes
tr.gam1<-c(0, 3, 6)		# true cutpoints for the first ordinal variable
tr.gam2<-c(0, 2, 4)		# true cutpoints for the second ordinal variable
ys[,3]<-ifelse(z[,3]>tr.gam1[3], 4, ifelse(z[,3]>tr.gam1[2],3,ifelse(z[,3]>tr.gam1[1], 2, 1)))
ys[,4]<-ifelse(z[,4]>tr.gam2[3], 4, ifelse(z[,4]>tr.gam2[2],3,ifelse(z[,4]>tr.gam2[1], 2, 1)))

### Fit the proposed model.
N<-ns						# sample size
K<-nrow(tr.A3)				# no. of latent variables
Q<-ncol(tr.A3)				# no. of outcomes
Q1<-2						# no. of continuous outcomes
Q2<-Q-Q1					# no. of ordinal outcomes
S<-ncol(tr.A2)				# no. of mediators
P<-ncol(X)					# no. of predictors

## Matrix indicating which predictors contribute to latent components
indx.Ws<-matrix(c(rep(1, 5), rep(0, 2), rep(0, 5), rep(1, 2)),P, K)

### Hyperparameters
# multiplicative constant for the prior variance of matrix W
cW<-10				
# multiplicative constant for the prior variance of matrix A
cA<-100				
# prior mean vector for regression coefficents A2 including the intercept 
# term
E.A2<-rep(0, S+1)			
# prior variance matrix for regression coefficents A2 including the 
# intercept term
V.A2<-cA*diag(S+1)		
# prior mean vector for regression coefficents A1 and A3 including the 
# intercept terms
E.A<-rep(0, K+S+1)		
# prior variance matrix for regression coefficents A1 and A3 including the 
# intercept terms
V.A<-cA*diag(K+S+1)		
# prior mean vector for component weights W
E.Ws<-rep(0, P)			
# prior variance matrix for component weights W
V.Ws<-cW*diag(P)			

# degrees of freedom for a Wishart distribution for the inverse of residual 
# covariance matrices
df.C<-0.01				
# scale matrix for a Wishart distribution for the inverse of # a residual 
# covariance matrices
V.R<-1/df.C*diag(Q)		
# width for a candidate value for a M-H algorithm
alp.MH<-4				


## Initial values
# phi for continuous outcomes
phis1<-rgamma(N, shape=nu.t2/2, rate=nu.t2/2)		
# phi for ordinal outcomes
phis2<-rgamma(N, shape=nu.t2/2, rate=nu.t2/2)		
# phi for mediators
phis3<-rgamma(N, shape=nu.t2/2, rate=nu.t2/2)		
# residual covariance matrix for mediators
Rs2<-diag(S)*1000										
Rs2.inv<-solve(Rs2)
# residual covariance matrix for continuous and ordinal outcomes
Rs3<-diag(Q)*1000										
# scaled for the ordinal variables
diag(Rs3)[3:4]<-1							
Rs3.inv<-solve(Rs3)
# regression coefficents A1 and A3 including the intercept terms
A<-t(rmvnorm(Q, mean=E.A, sigma=V.A)*0.01)		
# regression coefficents A2 including the intercept term
A2<-t(rmvnorm(S, mean=E.A2, sigma=V.A2)*0.01)		
# component weights W
Ws<-t(rmvnorm(K, mean=E.Ws, sigma=V.Ws))			
Ws[indx.Ws==0]<-0
# intercept vector for outcomes
bt03s<-A[1,]							
bt03s.mat<-matrix(bt03s, N, Q, byrow=T)
A34<-A[-1,]
# intercept vector for mediators
bt02s<-A2[1,]							
bt02s.mat<-matrix(bt02s, N, S, byrow=T)
A2<-A2[-1,]
# cutpoints for ordinal outcomes
gam<-matrix(c(-Inf, 0, 0.2, 0.5, Inf), 5, 2)		
# unique values for ordinal outcomes
lvl.ys<-list(sort(unique(ys[,3])), sort(unique(ys[,4])))	

## Grid points for degrees of freedom for multivariate t distributions
can.nu.t<-seq(0.5, 30, by=0.5)
can.nu.t.mat<-matrix(can.nu.t, ns, length(can.nu.t), byrow=T)

## Define matrices and vectors for posterior samples for parameters of 
## interest
post.Ws<-array(0, dim=c(P,K,nitr))
post.A2<-array(0, dim=c(K,S,nitr))
post.A34<-array(0, dim=c(K+S,Q,nitr))
post.Bs2<-array(0, dim=c(P,S,nitr))
post.Bs3<-array(0, dim=c(P,Q,nitr))
post.bt02<-matrix(0, nitr, S)
post.bt03<-matrix(0, nitr, Q)
post.R2<-matrix(0, nitr, S*(S+1)/2)
post.R3<-matrix(0, nitr, Q*(Q+1)/2)
post.gam<-array(0, dim=c(nrow(gam)-2, 2, nitr))
post.nu1<-rep(0, nitr)
post.nu3<-rep(0, nitr)
post.phis1<-matrix(0, nitr, N)
post.phis2<-matrix(0, nitr, N)
post.phis3<-matrix(0, nitr, N)


### MCMC algorithm begins here.
for (ni in 1:nitr){

## Update the latent variables Z for ordinal outcomes.
if(ni>1){tmp<-Z}
Z<-t(sapply(1:N, function(x){
mu.x<-as.numeric(bt03s+X[x,]%*%(Ws%*%A34[1:K,])+m[x,]%*%A34[K+(1:S),])
Sig.x<-diag(c(rep(sqrt(sig2.1/phis1[x]), Q1), rep(sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3%*%diag(c(rep(sqrt(sig2.1/phis1[x]), Q1), rep(sqrt(sig2.2/phis2[x]), Q2)))

## conditional mean and covariance
mu.i<-as.numeric(mu.x[Q1+1:Q2]+Sig.x[Q1+(1:Q2),1:Q1]%*%solve(Sig.x[1:Q1,1:Q1])%*%(ys[x,1:Q1]-mu.x[1:Q1]))
Sig.i<-Sig.x[Q1+(1:Q2),Q1+(1:Q2)]-Sig.x[Q1+(1:Q2),1:Q1]%*%solve(Sig.x[1:Q1,1:Q1])%*%Sig.x[1:Q1, Q1+(1:Q2)]
rtmvnorm(1, lower=gam[ys[x,Q1+1:Q2]+c(0,nrow(gam))], upper=gam[ys[x,Q1+1:Q2]+c(0,nrow(gam))+1], 
mean=mu.i, sigma=Sig.i, algorithm="gibbs", burn.in.samples=1000)}))

## additional step to avoid NaN due to a bad combination of parameters
if(sum(is.na(Z))>0){
Z[apply(is.na(Z),1,any),]<-tmp[apply(is.na(Z),1,any),]
}

## Define all latent variables for both continuous and ordinal outcomes.
Z.all<-cbind(ys[, 1:Q1], Z)

## Update the component weight Wk.
for (k in 1:K){
if (sum(indx.Ws[,k]==1)>1){
post.var.Wk<-solve(
matrix(apply(sapply(1:N, function(x){outer(X[x, indx.Ws[,k]==1], X[x, indx.Ws[,k]==1])*(as.numeric(A2[k,]%*%diag(rep(sqrt(phis3[x]), S))%*%Rs2.inv%*%diag(rep(sqrt(phis3[x]), S))%*%A2[k,])
+as.numeric(A34[k,]%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3.inv%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%A34[k,]))
}), 1, sum), sum(indx.Ws[,k]), sum(indx.Ws[,k]))+solve(V.Ws[indx.Ws[,k]==1,indx.Ws[,k]==1]))

if (K>2){
post.mean.Wk<-post.var.Wk%*%(
apply(sapply(1:N, function(x){X[x, indx.Ws[,k]==1]*(as.numeric(as.numeric(A2[k,]%*%diag(rep(sqrt(phis1[x]), S))%*%Rs2.inv%*%diag(rep(sqrt(phis1[x]), S)))%*%
	as.numeric(m[x,]-bt02s-X[x,]%*%(Ws[,-k]%*%A2[-k,])))
+as.numeric(as.numeric(A34[k,]%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3.inv%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2))))%*%
	as.numeric(Z.all[x,]-bt03s-X[x,]%*%(Ws[,-k]%*%A34[1:K,][-k,])-m[x,]%*%A34[K+(1:S),])))
}), 1, sum))
}else{
post.mean.Wk<-post.var.Wk%*%(
apply(sapply(1:N, function(x){X[x, indx.Ws[,k]==1]*(as.numeric(as.numeric(A2[k,]%*%diag(rep(sqrt(phis1[x]), S))%*%Rs2.inv%*%diag(rep(sqrt(phis1[x]), S)))%*%
	as.numeric(m[x,]-bt02s-X[x,]%*%(Ws[,-k]%o%A2[-k,])))
+as.numeric(as.numeric(A34[k,]%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3.inv%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2))))%*%
	as.numeric(Z.all[x,]-bt03s-X[x,]%*%(Ws[,-k]%o%A34[1:K,][-k,])-m[x,]%*%A34[K+(1:S),])))
}), 1, sum))
}
Ws[indx.Ws[,k]==1,k]<-t(rmvnorm(1, mean=as.vector(post.mean.Wk), sigma=post.var.Wk))
# normalization to satisfy the constrant that t(X%*%Ws[,k])%*%(X%*%Ws[,k])=N.
Ws[,k]<-Ws[,k]*as.numeric(sqrt(N)/sqrt(t(X%*%Ws[,k])%*%(X%*%Ws[,k])))
}else{Ws[indx.Ws[,k]==1,k]<-1}			# if there is only one predictor for a latent variable.
}


## Update the regression coefficients A2 including the intercept for mediators m.
Fs2<-cbind(1, X%*%Ws)
tmp.XX2<-matrix(apply(sapply(1:ns, function(x){t(kronecker(diag(S),matrix(Fs2[x,],1)))%*%diag(rep(sqrt(phis3[x]),S))%*%Rs2.inv%*%diag(rep(sqrt(phis3[x]),S))%*%(kronecker(diag(S),matrix(Fs2[x,],1)))}), 1, sum), dim(Fs2)[2]*S) 
tmp.XZ2<-apply(sapply(1:ns, function(x){t(kronecker(diag(S),matrix(Fs2[x,],1)))%*%diag(rep(sqrt(phis3[x]),S))%*%Rs2.inv%*%diag(rep(sqrt(phis3[x]),S))%*%m[x,]}),1,sum)
A2<-matrix(rmvnorm(1, mean=as.numeric(solve(tmp.XX2+solve(kronecker(diag(S),V.A2)))%*%(tmp.XZ2+solve(kronecker(diag(S),V.A2))%*%rep(E.A2, S))), sigma=solve(tmp.XX2+solve(kronecker(diag(S),V.A2)))), K+1)

## Update the residual covariance matrix  for mediators m.
tmp.ee2<-matrix(apply(sapply(1:ns, function(x){t((m[x,]-Fs2[x,]%*%A2)%*%diag(rep(sqrt(phis3[x]), S)))%*%(m[x,]-Fs2[x,]%*%A2)%*%diag(rep(sqrt(phis3[x]), S))}), 1, sum), S) 
Rs2.inv<-rWishart(1, df.C+ns, solve(solve(V.R[1:S, 1:S])+tmp.ee2))[,,1]
Rs2<-chol2inv(chol(Rs2.inv))

## Update the regression coefficients A including the intercept for latent outcomes Z.all.
Fs3<-cbind(1, X%*%Ws, m)
tmp.XX3<-matrix(apply(sapply(1:ns, function(x){t(kronecker(diag(Q),matrix(Fs3[x,],1)))%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3.inv%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%(kronecker(diag(Q),matrix(Fs3[x,],1)))}), 1, sum), dim(Fs3)[2]*Q) 
tmp.XZ3<-apply(sapply(1:ns, function(x){t(kronecker(diag(Q),matrix(Fs3[x,],1)))%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3.inv%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%Z.all[x,]}),1,sum)
A<-matrix(rmvnorm(1, mean=as.numeric(solve(tmp.XX3+solve(kronecker(diag(Q),V.A)))%*%(tmp.XZ3+solve(kronecker(diag(Q),V.A))%*%rep(E.A, Q))), sigma=solve(tmp.XX3+solve(kronecker(diag(Q),V.A)))), K+S+1)

## Update the residual covariance matrix for latent outcomes Z.all.
tmp.ee3<-matrix(apply(sapply(1:ns, function(x){t((Z.all[x,]-Fs3[x,]%*%A)%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2))))%*%(Z.all[x,]-Fs3[x,]%*%A)%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))}), 1, sum), Q) 
Rs3.inv<-rWishart(1, df.C+ns, solve(solve(V.R)+tmp.ee3))[,,1]
Rs3<-chol2inv(chol(Rs3.inv))

## M-H step for subject-specific precision parameter phis1 for continuous outcomes
lower.bnd1<-apply(cbind(phis1-alp.MH, 0), 1, max)
upper.bnd1<-phis1+alp.MH
can.phis1<-runif(N, lower.bnd1, upper.bnd1) 

target.ratio1<-sapply(1:N, function(x){exp(dmvnorm(Z.all[x,], mean=(Fs3%*%A)[x,], sigma=diag(c(rep(sqrt(sig2.1/can.phis1[x]),Q1), rep(sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3%*%diag(c(rep(sqrt(sig2.1/can.phis1[x]),Q1), rep(sqrt(sig2.2/phis2[x]), Q2))), log=T)+
dgamma(can.phis1[x], shape=nu.t1/2, rate=nu.t1/2, log=T)-
dmvnorm(Z.all[x,], mean=(Fs3%*%A)[x,], sigma=diag(c(rep(sqrt(sig2.1/phis1[x]),Q1), rep(sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3%*%diag(c(rep(sqrt(sig2.1/phis1[x]),Q1), rep(sqrt(sig2.2/phis2[x]), Q2))), log=T)-
dgamma(phis1[x], shape=nu.t1/2, rate=nu.t1/2, log=T))})
phis1<-ifelse(runif(N)>target.ratio1, phis1, can.phis1)
phis1.mat<-matrix(phis1, N, length(can.nu.t))

## M-H step for subject-specific precision parameter phis2 for latent variables for ordinal outcomes
lower.bnd<-apply(cbind(phis2-alp.MH, 0), 1, max)
upper.bnd<-phis2+alp.MH
can.phis2<-runif(N, lower.bnd, upper.bnd) 

target.ratio<-sapply(1:N, function(x){exp(dmvnorm(Z.all[x,], mean=(Fs3%*%A)[x,], sigma=diag(c(rep(sqrt(sig2.1/phis1[x]),Q1), rep(sqrt(sig2.2/can.phis2[x]), Q2)))%*%Rs3%*%diag(c(rep(sqrt(sig2.1/phis1[x]),Q1), rep(sqrt(sig2.2/can.phis2[x]), Q2))), log=T)+
dgamma(can.phis2[x], shape=nu.t2/2, rate=nu.t2/2, log=T)-
dmvnorm(Z.all[x,], mean=(Fs3%*%A)[x,], sigma=diag(c(rep(sqrt(sig2.1/phis1[x]),Q1), rep(sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3%*%diag(c(rep(sqrt(sig2.1/phis1[x]),Q1), rep(sqrt(sig2.2/phis2[x]), Q2))), log=T)-
dgamma(phis2[x], shape=nu.t2/2, rate=nu.t2/2, log=T))})
phis2<-ifelse(runif(N)>target.ratio, phis2, can.phis2)
phis2.mat<-matrix(phis2, N, length(can.nu.t))

## Update subject-specific precision parameter phis3 for mediators m.
phis3<-rgamma(ns, shape=(nu.t3+S)/2, rate=nu.t3/2+1/2*sapply(1:ns, function(x){(m[x,]-Fs2[x,]%*%A2)%*%Rs2.inv%*%t(m[x,]-Fs2[x,]%*%A2)}))
phis3.mat<-matrix(phis3, ns, length(can.nu.t))


## Scale the residual covariance matrix Rs3 and regression coefficients A, so the variances for latent variables corresponding to ordinal outcomes are set to be one. 
scale.factor<-sqrt(c(rep(1, Q1), diag(Rs3)[Q1+1:Q2]))
Rs3<-diag(1/scale.factor)%*%Rs3%*%diag(1/scale.factor)
Rs3.inv<-chol2inv(chol(Rs3))
A<-A*matrix(rep(1/scale.factor, each=S+K+1), S+K+1)		


## Divide A into two parts: bt02s and A2
bt02s<-A2[1,]
bt02s.mat<-matrix(bt02s, N, S, byrow=T)
A2<-A2[-1,]

## Divide A into two parts: bt03s and A34
bt03s<-A[1,]
bt03s.mat<-matrix(bt03s, N, Q, byrow=T)
A34<-A[-1,]

## Update the cutpoints for latent variables Z corresponding to ordinal outcomes. 
# Note that one of cutpoints or the intercept should be fixed for identifiability. In the simulation study, the first finite element of the cutpoints is set to be zero.
for (i in (Q1+1):Q){
for (j in 2:(length(lvl.ys[[i-Q1]])-1)){
gam[j+1,i-Q1]<-runif(1,min=max(Z.all[ys[,i]==lvl.ys[[i-Q1]][j],i], gam[j,i-Q1]), max=min(Z.all[ys[,i]==lvl.ys[[i-Q1]][j+1],i], gam[j+2,i-Q1]))
}
}

## Update degrees of freedom nu.t1 of a multivarite t distribution for continuous outcomes.
den.nu.t1<-apply(dgamma(phis1.mat, shape=can.nu.t.mat/2, rate=can.nu.t.mat/2, log=T), 2, sum)
nu.t1<-min(can.nu.t[runif(1)<=cumsum(exp(den.nu.t1-max(den.nu.t1))/sum(exp(den.nu.t1-max(den.nu.t1))))])

## Update degrees of freedom nu.t3 of a multivarite t distribution for mediators.
den.nu.t3<-apply(dgamma(phis3.mat, shape=can.nu.t.mat/2, rate=can.nu.t.mat/2, log=T), 2, sum)
nu.t3<-min(can.nu.t[runif(1)<=cumsum(exp(den.nu.t3-max(den.nu.t3))/sum(exp(den.nu.t3-max(den.nu.t3))))])


## save MCMC samples
post.Ws[,,ni]<-Ws
post.A2[,,ni]<-A2
post.A34[,,ni]<-A34
post.Bs2[,,ni]<-Ws%*%A2
post.Bs3[,,ni]<-Ws%*%A34[1:K,]		
post.bt02[ni,]<-bt02s
post.bt03[ni,]<-bt03s
post.R2[ni,]<-Rs2[lower.tri(Rs2, diag=T)]
post.R3[ni,]<-Rs3[lower.tri(Rs3, diag=T)]
post.gam[,,ni]<-gam[-c(1,5),]
post.phis1[ni,]<-phis1
post.phis2[ni,]<-phis2
post.phis3[ni,]<-phis3
post.nu1[ni]<-nu.t1
post.nu3[ni]<-nu.t3
}

