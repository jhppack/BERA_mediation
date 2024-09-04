BERA.med<-function(ys, Xs, ms, indx.Ws, no.cont.y, nitr=50000, hyper.parm=list(E.A2=0, V.A2=100, E.A3=0, V.A3=100, E.A4=0, V.A4=100, E.W=0, V.W=10, df.C=0.01, alp.MH=4)){
#### This function fits the BERA-mediation model. The following function arguments should be specified:

### Arguments
## - ys: An outcome matrix. Continuous and ordinal outcomes should be entered in columns, with continuous outcomes first. The number of ordinal outcomes must be at least two.
## - Xs: A matrix of correlated predictors. The column corresponding to the intercept term should not be included.
## - ms: A mediator matrix.
## - indx.Ws: An indicator matrix specifying which elements of the component weight matrix W are non-zero.
## - no.cont.y: The number of continuous outcomes in `ys`.
## - nitr: The number of iterations for the proposed MCMC algorithm. The default value is 50,000.
## - hyper.parm: Values for the hyperparameters of the set of proposed prior distributions. "E" and "V" before the period (.) denote the mean and variance of a normal prior 
##	distribution, respectively. The letters following the period correspond to the matrices in the proposed model. The default values are set to make the prior 
##	distribution non-informative. These are the same values used for Simulation Setting 1, such as `list(E.A2=0, V.A2=100, E.A3=0, V.A3=100, E.A4=0, V.A4=100, E.W=0, V.W=10, df.C=0.01, alp.MH=4)`.

### Value
## The output of this function is returned as a list with components containing posterior samples for parameters of interest. The component names are as follows:
## - Ws: Non-zero component weights of W.
## - A2: Regression coefficients of XW for mediators.
## - A3: Regression coefficients of XW for outcomes.
## - A4: Regression coefficients of mediators for outcomes.
## - Bs2: The product of W and A2.
## - Bs3: The product of W and A3.
## - bt02: The intercept vector for mediators.
## - bt03: The intercept vector for outcomes.
## - R2: Residual covariance matrix for mediators.
## - R3: Residual covariance matrix for outcomes.
## - nu1: Degrees of freedom for a t-distribution for continuous outcomes.
## - nu2: Degrees of freedom for a t-distribution for mediators.
 


library(mvtnorm)
library(tmvtnorm)
library(MCMCpack)

## Denote variables used in the algorithm.
N<-nrow(ys)			# sample size
K<-ncol(indx.Ws)		# no. of latent variables
Q<-ncol(ys)			# no. of outcomes
Q1<-no.cont.y		# no. of continuous outcomes
Q2<-Q-Q1			# no. of ordinal outcomes
S<-ncol(ms)			# no. of mediators
P<-ncol(Xs)			# no. of predictors
E.A2<-hyper.parm[[1]]			# prior mean vector for regression coefficents A2
V.A2<-as.matrix(hyper.parm[[2]])	# prior variance matrix for regression coefficents A2
E.A3<-hyper.parm[[3]]			# prior mean vector for regression coefficents A3
V.A3<-as.matrix(hyper.parm[[4]])	# prior variance matrix for regression coefficents A3
E.A4<-hyper.parm[[5]]			# prior mean vector for regression coefficents A4
V.A4<-as.matrix(hyper.parm[[6]])	# prior variance matrix for regression coefficents A4
E.W<-hyper.parm[[7]]			# prior mean vector for component weights W
V.W<-as.matrix(hyper.parm[[8]])	# prior variance matrix for component weights W
df.C<-hyper.parm[[9]]			# degrees of freedom for a Wishart distribution for the inverse of residual covariance matrices
V.R<-1/0.01*diag(Q)			# scale matrix for a Wishart distribution for the inverse of residual covariance matrices	
alp.MH<-hyper.parm[[10]]		# width for a candidate value for a M-H algorithm			

if (length(E.A2)!=nrow(V.A2)) stop("The dimensions of the prior mean and covariance matrix for regression coefficients A2 do not match.")
if (length(E.A3)!=nrow(V.A3)) stop("The dimensions of the prior mean and covariance matrix for regression coefficients A3 do not match.")
if (length(E.A4)!=nrow(V.A4)) stop("The dimensions of the prior mean and covariance matrix for regression coefficients A4 do not match.")
if (length(E.W)!=nrow(V.W)) stop("The dimensions of the prior mean and covariance matrix for weights W do not match.")
if (ncol(Xs)!=nrow(indx.Ws)) stop("The number of preditors and the numbers of rows in indx.Ws do not match.")

## Fix a degree of freedom for the t-distributions underlying ordinal outcomesww
nu.t2<-7.3			

## Initial values
phis1<-rgamma(N, shape=nu.t2/2, rate=nu.t2/2)		# phi for continuous outcomes
phis2<-rgamma(N, shape=nu.t2/2, rate=nu.t2/2)		# phi for ordinal outcomes
phis3<-rgamma(N, shape=nu.t2/2, rate=nu.t2/2)		# phi for mediators
Rs2<-diag(S)*1							# residual covariance matrix for mediators						
Rs2.inv<-solve(Rs2)
Rs3<-diag(Q)*1							# residual covariance matrix for continuous and ordinal outcomes						
diag(Rs3)[(Q1+1):Q]<-1							# scaled for the ordinal variables
Rs3.inv<-solve(Rs3)

if (length(E.A2)==1){E.A2<-rep(E.A2, S); V.A2<-diag(S)*as.numeric(V.A2)}
if (length(E.A3)==1){E.A3<-rep(E.A3, K); V.A3<-diag(K)*as.numeric(V.A3)}
if (length(E.A4)==1){E.A4<-rep(E.A4, S); V.A4<-diag(S)*as.numeric(V.A4)}
if (length(E.W)==1){E.W<-rep(E.W, P); V.W<-diag(P)*as.numeric(V.W)}

A2<-t(rmvnorm(S, mean=E.A2, sigma=V.A2)*0.01)		# regression coefficents A2
bt02s<-rnorm(S, 0, 100)*0.01					# intercept vector for mediators
bt02s.mat<-matrix(bt02s, N, S, byrow=T)
# Define the mean vector and covariance matrix for the intercept term and A2
E.A2<-c(0, E.A2)
V.A2.tmp<-V.A2; V.A2<-diag(S+1)*100; V.A2[-1,-1]<-V.A2.tmp	

A3<-t(rmvnorm(Q, mean=E.A3, sigma=V.A3)*0.01)		# regression coefficents A3
A4<-t(rmvnorm(Q, mean=E.A4, sigma=V.A4)*0.01)		# regression coefficents A4
A34<-rbind(A3, A4)						# regression coefficents A3 and A4
bt03s<-rnorm(Q, 0, 100)*0.01					# intercept vector for outcomes
bt03s.mat<-matrix(bt03s, N, Q, byrow=T)
# Define the mean vector and covariance matrix for the intercept term, A3, and A4 
E.A<-c(0, E.A3, E.A4)
V.A.tmp<-matrix(0, K+S, K+S); V.A.tmp[1:K, 1:K]<-V.A3; V.A.tmp[(1:S)+K, (1:S)+K]<-V.A4;V.A<-diag(K+S+1)*100; V.A[-1,-1]<-V.A.tmp

Ws<-t(rmvnorm(K, mean=E.W, sigma=V.W))			# component weights W
Ws[indx.Ws==0]<-0

lvl.ys<-gam<-vector(length=Q2, mode='list')
gam<-matrix(NA, 20, Q2)		# cutpoints for ordinal outcomes. The maximum of 20 levels for each categorical variable is assumed. 
can.gam<-(0:20)/20						
for (j in 1:Q2){ 
lvl.ys[[j]]<-sort(unique(ys[,Q1+j]))			# unique values for ordinal outcomes
gam[1:(length(lvl.ys[[j]])+1),j]<-c(-Inf, (0:(length(lvl.ys[[j]])-2))/(length(lvl.ys[[j]])), Inf)
}

## Grid points for degrees of freedom for multivariate t distributions
can.nu.t<-seq(0.5, 30, by=0.5)
can.nu.t.mat<-matrix(can.nu.t, ns, length(can.nu.t), byrow=T)

## Define matrices and vectors for posterior samples for parameters of interest
post.Ws<-array(0, dim=c(P,K,nitr))
post.A2<-array(0, dim=c(K,S,nitr))
post.A34<-array(0, dim=c(K+S,Q,nitr))
post.Bs2<-array(0, dim=c(P,S,nitr))
post.Bs3<-array(0, dim=c(P,Q,nitr))
post.bt02<-matrix(0, nitr, S)
post.bt03<-matrix(0, nitr, Q)
post.R2<-matrix(0, nitr, S*(S+1)/2)
post.R3<-matrix(0, nitr, Q*(Q+1)/2)
post.nu1<-rep(0, nitr)
post.nu3<-rep(0, nitr)

### MCMC algorithm begins here.
for (ni in 1:nitr){

## Update the latent variables Zs for ordinal outcomes.
if(ni>1){tmp<-Zs}
Zs<-t(sapply(1:N, function(x){
mu.x<-as.numeric(bt03s+Xs[x,]%*%(Ws%*%A34[1:K,])+ms[x,]%*%A34[K+(1:S),])
Sig.x<-diag(c(rep(sqrt(sig2.1/phis1[x]), Q1), rep(sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3%*%diag(c(rep(sqrt(sig2.1/phis1[x]), Q1), rep(sqrt(sig2.2/phis2[x]), Q2)))
## conditional mean and covariance
mu.i<-as.numeric(mu.x[Q1+1:Q2]+Sig.x[Q1+(1:Q2),1:Q1]%*%solve(Sig.x[1:Q1,1:Q1])%*%(ys[x,1:Q1]-mu.x[1:Q1]))
Sig.i<-Sig.x[Q1+(1:Q2),Q1+(1:Q2)]-Sig.x[Q1+(1:Q2),1:Q1]%*%solve(Sig.x[1:Q1,1:Q1])%*%Sig.x[1:Q1, Q1+(1:Q2)]
rtmvnorm(1, lower=gam[ys[x,Q1+1:Q2]+c(0,nrow(gam))], upper=gam[ys[x,Q1+1:Q2]+c(0,nrow(gam))+1], 
mean=mu.i, sigma=Sig.i, algorithm="gibbs", burn.in.samples=1000)}))

## additional step to avoid NaN due to a bad combination of parameters
if(sum(is.na(Zs))>0){
Zs[apply(is.na(Zs),1,any),]<-tmp[apply(is.na(Zs),1,any),]
}

## Define all latent variables for both continuous and ordinal outcomes.
Zs.all<-cbind(ys[, 1:Q1], Zs)

## Update the component weight Wk.
for (k in 1:K){
if (sum(indx.Ws[,k]==1)>1){
post.var.Wk<-solve(
matrix(apply(sapply(1:N, function(x){outer(Xs[x, indx.Ws[,k]==1], Xs[x, indx.Ws[,k]==1])*(as.numeric(A2[k,]%*%diag(rep(sqrt(phis3[x]), S))%*%Rs2.inv%*%diag(rep(sqrt(phis3[x]), S))%*%A2[k,])
+as.numeric(A34[k,]%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3.inv%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%A34[k,]))
}), 1, sum), sum(indx.Ws[,k]), sum(indx.Ws[,k]))+solve(V.W[indx.Ws[,k]==1,indx.Ws[,k]==1]))

if (K>2){
post.mean.Wk<-post.var.Wk%*%(
apply(sapply(1:N, function(x){Xs[x, indx.Ws[,k]==1]*(as.numeric(as.numeric(A2[k,]%*%diag(rep(sqrt(phis1[x]), S))%*%Rs2.inv%*%diag(rep(sqrt(phis1[x]), S)))%*%
	as.numeric(ms[x,]-bt02s-Xs[x,]%*%(Ws[,-k]%*%A2[-k,])))
+as.numeric(as.numeric(A34[k,]%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3.inv%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2))))%*%
	as.numeric(Zs.all[x,]-bt03s-Xs[x,]%*%(Ws[,-k]%*%A34[1:K,][-k,])-ms[x,]%*%A34[K+(1:S),])))
}), 1, sum))
}else{
post.mean.Wk<-post.var.Wk%*%(
apply(sapply(1:N, function(x){Xs[x, indx.Ws[,k]==1]*(as.numeric(as.numeric(A2[k,]%*%diag(rep(sqrt(phis1[x]), S))%*%Rs2.inv%*%diag(rep(sqrt(phis1[x]), S)))%*%
	as.numeric(ms[x,]-bt02s-Xs[x,]%*%(Ws[,-k]%o%A2[-k,])))
+as.numeric(as.numeric(A34[k,]%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3.inv%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2))))%*%
	as.numeric(Zs.all[x,]-bt03s-Xs[x,]%*%(Ws[,-k]%o%A34[1:K,][-k,])-ms[x,]%*%A34[K+(1:S),])))
}), 1, sum))
}
Ws[indx.Ws[,k]==1,k]<-t(rmvnorm(1, mean=as.vector(post.mean.Wk), sigma=post.var.Wk))
# normalization to satisfy the constrant that t(Xs%*%Ws[,k])%*%(Xs%*%Ws[,k])=N.
Ws[,k]<-Ws[,k]*as.numeric(sqrt(N)/sqrt(t(Xs%*%Ws[,k])%*%(Xs%*%Ws[,k])))
}else{Ws[indx.Ws[,k]==1,k]<-1}			# if there is only one predictor for a latent variable.
}


## Update the regression coefficients A2 including the intercept for mediators ms.
Fs2<-cbind(1, Xs%*%Ws)
tmp.XsXs2<-matrix(apply(sapply(1:ns, function(x){t(kronecker(diag(S),matrix(Fs2[x,],1)))%*%diag(rep(sqrt(phis3[x]),S))%*%Rs2.inv%*%diag(rep(sqrt(phis3[x]),S))%*%(kronecker(diag(S),matrix(Fs2[x,],1)))}), 1, sum), dim(Fs2)[2]*S) 
tmp.XsZs2<-apply(sapply(1:ns, function(x){t(kronecker(diag(S),matrix(Fs2[x,],1)))%*%diag(rep(sqrt(phis3[x]),S))%*%Rs2.inv%*%diag(rep(sqrt(phis3[x]),S))%*%ms[x,]}),1,sum)
A2<-matrix(rmvnorm(1, mean=as.numeric(solve(tmp.XsXs2+solve(kronecker(diag(S),V.A2)))%*%(tmp.XsZs2+solve(kronecker(diag(S),V.A2))%*%rep(E.A2, S))), sigma=solve(tmp.XsXs2+solve(kronecker(diag(S),V.A2)))), K+1)

## Update the residual covariance matrix  for mediators ms.
tmp.ee2<-matrix(apply(sapply(1:ns, function(x){t((ms[x,]-Fs2[x,]%*%A2)%*%diag(rep(sqrt(phis3[x]), S)))%*%(ms[x,]-Fs2[x,]%*%A2)%*%diag(rep(sqrt(phis3[x]), S))}), 1, sum), S) 
Rs2.inv<-rWishart(1, df.C+ns, solve(solve(V.R[1:S, 1:S])+tmp.ee2))[,,1]
Rs2<-chol2inv(chol(Rs2.inv))

## Update the regression coefficients A including the intercept for latent outcomes Zs.all.
Fs3<-cbind(1, Xs%*%Ws, ms)
tmp.XsXs3<-matrix(apply(sapply(1:ns, function(x){t(kronecker(diag(Q),matrix(Fs3[x,],1)))%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3.inv%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%(kronecker(diag(Q),matrix(Fs3[x,],1)))}), 1, sum), dim(Fs3)[2]*Q) 
tmp.XsZs3<-apply(sapply(1:ns, function(x){t(kronecker(diag(Q),matrix(Fs3[x,],1)))%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3.inv%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))%*%Zs.all[x,]}),1,sum)
A<-matrix(rmvnorm(1, mean=as.numeric(solve(tmp.XsXs3+solve(kronecker(diag(Q),V.A)))%*%(tmp.XsZs3+solve(kronecker(diag(Q),V.A))%*%rep(E.A, Q))), sigma=solve(tmp.XsXs3+solve(kronecker(diag(Q),V.A)))), K+S+1)

## Update the residual covariance matrix for latent outcomes Zs.all.
tmp.ee3<-matrix(apply(sapply(1:ns, function(x){t((Zs.all[x,]-Fs3[x,]%*%A)%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2))))%*%(Zs.all[x,]-Fs3[x,]%*%A)%*%diag(c(rep(1/sqrt(sig2.1/phis1[x]), Q1), rep(1/sqrt(sig2.2/phis2[x]), Q2)))}), 1, sum), Q) 
Rs3.inv<-rWishart(1, df.C+ns, solve(solve(V.R)+tmp.ee3))[,,1]
Rs3<-chol2inv(chol(Rs3.inv))

## M-H step for subject-specific precision parameter phis1 for continuous outcomes
lower.bnd1<-apply(cbind(phis1-alp.MH, 0), 1, max)
upper.bnd1<-phis1+alp.MH
can.phis1<-runif(N, lower.bnd1, upper.bnd1) 

target.ratio1<-sapply(1:N, function(x){exp(dmvnorm(Zs.all[x,], mean=(Fs3%*%A)[x,], sigma=diag(c(rep(sqrt(sig2.1/can.phis1[x]),Q1), rep(sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3%*%diag(c(rep(sqrt(sig2.1/can.phis1[x]),Q1), rep(sqrt(sig2.2/phis2[x]), Q2))), log=T)+
dgamma(can.phis1[x], shape=nu.t1/2, rate=nu.t1/2, log=T)-
dmvnorm(Zs.all[x,], mean=(Fs3%*%A)[x,], sigma=diag(c(rep(sqrt(sig2.1/phis1[x]),Q1), rep(sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3%*%diag(c(rep(sqrt(sig2.1/phis1[x]),Q1), rep(sqrt(sig2.2/phis2[x]), Q2))), log=T)-
dgamma(phis1[x], shape=nu.t1/2, rate=nu.t1/2, log=T))})
phis1<-ifelse(runif(N)>target.ratio1, phis1, can.phis1)
phis1.mat<-matrix(phis1, N, length(can.nu.t))

## M-H step for subject-specific precision parameter phis2 for latent variables for ordinal outcomes
lower.bnd<-apply(cbind(phis2-alp.MH, 0), 1, max)
upper.bnd<-phis2+alp.MH
can.phis2<-runif(N, lower.bnd, upper.bnd) 

target.ratio<-sapply(1:N, function(x){exp(dmvnorm(Zs.all[x,], mean=(Fs3%*%A)[x,], sigma=diag(c(rep(sqrt(sig2.1/phis1[x]),Q1), rep(sqrt(sig2.2/can.phis2[x]), Q2)))%*%Rs3%*%diag(c(rep(sqrt(sig2.1/phis1[x]),Q1), rep(sqrt(sig2.2/can.phis2[x]), Q2))), log=T)+
dgamma(can.phis2[x], shape=nu.t2/2, rate=nu.t2/2, log=T)-
dmvnorm(Zs.all[x,], mean=(Fs3%*%A)[x,], sigma=diag(c(rep(sqrt(sig2.1/phis1[x]),Q1), rep(sqrt(sig2.2/phis2[x]), Q2)))%*%Rs3%*%diag(c(rep(sqrt(sig2.1/phis1[x]),Q1), rep(sqrt(sig2.2/phis2[x]), Q2))), log=T)-
dgamma(phis2[x], shape=nu.t2/2, rate=nu.t2/2, log=T))})
phis2<-ifelse(runif(N)>target.ratio, phis2, can.phis2)
phis2.mat<-matrix(phis2, N, length(can.nu.t))

## Update subject-specific precision parameter phis3 for mediators ms.
phis3<-rgamma(ns, shape=(nu.t3+S)/2, rate=nu.t3/2+1/2*sapply(1:ns, function(x){(ms[x,]-Fs2[x,]%*%A2)%*%Rs2.inv%*%t(ms[x,]-Fs2[x,]%*%A2)}))
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
gam[j+1,i-Q1]<-runif(1,min=max(Zs.all[ys[,i]==lvl.ys[[i-Q1]][j],i], gam[j,i-Q1]), max=min(Zs.all[ys[,i]==lvl.ys[[i-Q1]][j+1],i], gam[j+2,i-Q1]))
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
post.nu1[ni]<-nu.t1
post.nu3[ni]<-nu.t3
}
return(list(Ws=post.Ws, A2=post.A2, A3=post.A34[1:K,,], A4=post.A34[-(1:K),,], Bs2=post.Bs2, Bs3=post.Bs3, bt02=post.bt02, bt03=post.bt03,
	R2=post.R2, R3=post.R3, nu1=post.nu1, nu3=post.nu3))
}


BERA.med.summary<-function(BERA.med.out, indx.W, burnin=2000, thinning=10){
### This function organizes and summarizes the output of BERA.med. The following function arguments must be specified:

### Arguments
## BERA.med.out: The output from the BERA.med function.
## indx.W: An indicator matrix specifying which elements of the component weight matrix W are non-zero.
## burnin: The number of iterations to discard as part of the burn-in period. The default value is 2,000.
## thinning: The thinning interval for the MCMC samples. The default value is 10.

### Value
## The output of this function is returned as a list with components that contain the posterior mean and 95% credible intervals for the parameters of interest. The component names are as follows:
## - Ws: Non-zero component weights of W.
## - A2: Regression coefficients of XW for mediators.
## - A3: Regression coefficients of XW for outcomes.
## - A4: Regression coefficients of mediators for outcomes.
## - bt02: The intercept vector for mediators.
## - bt03: The intercept vector for outcomes.
## - R2: Residual covariance matrix for mediators.
## - R3: Residual covariance matrix for outcomes.
	
S<-dim(BERA.med.out$A2)[2]		# no. of mediators
Q<-dim(BERA.med.out$A3)[2]		# no. of outcomes
K<-ncol(indx.W)				# no. of latent variables

post.seq<-seq(burnin,dim(BERA.med.out$A2)[3], by=thinning)

summary.Ws<-cbind(apply(BERA.med.out$Ws[,,post.seq],c(1,2), mean)[which(indx.W==1)],
apply(BERA.med.out$Ws[,,post.seq],c(1,2), quantile, prob=c(0.025))[which(indx.W==1)],
apply(BERA.med.out$Ws[,,post.seq],c(1,2), quantile, prob=c(0.975))[which(indx.W==1)])
rownames(summary.Ws)<-paste0('W', matrix(1:nrow(indx.W), nrow(indx.W), ncol(indx.W)),
matrix(1:ncol(indx.W), nrow(indx.W), ncol(indx.W), byrow=T))[indx.W==1]
colnames(summary.Ws)<-c("post.mean", "95%CI.lower", "95%CI.upper")

summary.A2<-cbind(apply(BERA.med.out$A2[,,post.seq],c(1,2), mean)[1:(S*K)],
apply(BERA.med.out$A2[,,post.seq],c(1,2), quantile, prob=c(0.025))[1:(S*K)],
apply(BERA.med.out$A2[,,post.seq],c(1,2), quantile, prob=c(0.975))[1:(S*K)])
rownames(summary.A2)<-paste0('A2.', matrix(1:dim(BERA.med.out$A2)[1], dim(BERA.med.out$A2)[1], dim(BERA.med.out$A2)[2]),
matrix(1:dim(BERA.med.out$A2)[2], dim(BERA.med.out$A2)[1], dim(BERA.med.out$A2)[2], byrow=T))
colnames(summary.A2)<-c("post.mean", "95%CI.lower", "95%CI.upper")


summary.A3<-cbind(apply(BERA.med.out$A3[,,post.seq],c(1,2), mean)[1:(K*Q)],
apply(BERA.med.out$A3[,,post.seq],c(1,2), quantile, prob=c(0.025))[1:(K*Q)],
apply(BERA.med.out$A3[,,post.seq],c(1,2), quantile, prob=c(0.975))[1:(K*Q)])
rownames(summary.A3)<-paste0('A3.', matrix(1:dim(BERA.med.out$A3)[1], dim(BERA.med.out$A3)[1], dim(BERA.med.out$A3)[2]),
matrix(1:dim(BERA.med.out$A3)[2], dim(BERA.med.out$A3)[1], dim(BERA.med.out$A3)[2], byrow=T))
colnames(summary.A3)<-c("post.mean", "95%CI.lower", "95%CI.upper")

summary.A4<-cbind(apply(BERA.med.out$A4[,,post.seq],c(1,2), mean)[1:(S*Q)],
apply(BERA.med.out$A4[,,post.seq],c(1,2), quantile, prob=c(0.025))[1:(S*Q)],
apply(BERA.med.out$A4[,,post.seq],c(1,2), quantile, prob=c(0.975))[1:(S*Q)])
rownames(summary.A4)<-paste0('A4.', matrix(1:dim(BERA.med.out$A4)[1], dim(BERA.med.out$A4)[1], dim(BERA.med.out$A4)[2]),
matrix(1:dim(BERA.med.out$A4)[2], dim(BERA.med.out$A4)[1], dim(BERA.med.out$A4)[2], byrow=T))
colnames(summary.A4)<-c("post.mean", "95%CI.lower", "95%CI.upper")

summary.bt02<-cbind(apply(BERA.med.out$bt02[post.seq,],c(2), mean),
apply(BERA.med.out$bt02[post.seq,],c(2), quantile, prob=c(0.025)),
apply(BERA.med.out$bt02[post.seq,],c(2), quantile, prob=c(0.975)))
rownames(summary.bt02)<-paste0('bt02.', 1:dim(BERA.med.out$bt02)[2])
colnames(summary.bt02)<-c("post.mean", "95%CI.lower", "95%CI.upper")

summary.bt03<-cbind(apply(BERA.med.out$bt03[post.seq,],c(2), mean),
apply(BERA.med.out$bt03[post.seq,],c(2), quantile, prob=c(0.025)),
apply(BERA.med.out$bt03[post.seq,],c(2), quantile, prob=c(0.975)))
rownames(summary.bt03)<-paste0('bt03.', 1:dim(BERA.med.out$bt03)[2])
colnames(summary.bt03)<-c("post.mean", "95%CI.lower", "95%CI.upper")

summary.R2<-cbind(apply(BERA.med.out$R2[post.seq,],c(2), mean),
apply(BERA.med.out$R2[post.seq,],c(2), quantile, prob=c(0.025)),
apply(BERA.med.out$R2[post.seq,],c(2), quantile, prob=c(0.975)))
rownames(summary.R2)<-matrix(paste0('R2.', matrix(1:S, S, S), matrix(1:S, S, S, byrow=T)), ncol=S)[lower.tri(matrix(1, S, S), diag=T)]
colnames(summary.R2)<-c("post.mean", "95%CI.lower", "95%CI.upper")

summary.R3<-cbind(apply(BERA.med.out$R3[post.seq,],c(2), mean),
apply(BERA.med.out$R3[post.seq,],c(2), quantile, prob=c(0.025)),
apply(BERA.med.out$R3[post.seq,],c(2), quantile, prob=c(0.975)))
rownames(summary.R3)<-matrix(paste0('R3.', matrix(1:Q, Q, Q), matrix(1:Q, Q, Q, byrow=T)), ncol=2)[lower.tri(matrix(1, Q, Q), diag=T)]
colnames(summary.R3)<-c("post.mean", "95%CI.lower", "95%CI.upper")

return(list(Ws=summary.Ws, A2=summary.A2, A3=summary.A3, A4=summary.A4, bt02=summary.bt02, bt03=summary.bt03,
	R2=summary.R2, R3=summary.R3))
}
