rm(list=ls())
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
#¤¤¤ make sure the data is generated using the data function (AssortativeMatingSim put in a function)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
SimAssMating<-function(
n.pairs=100, n.matings=3, n.repeated.measures=3,mean.mal=2,mean.fem=3,
var.i.mal=2, var.i.fem=2,rho.i.m.f=0.5
, var.CommEnv.mal=1.1, var.CommEnv.fem=1.1,rho.CommEnv.m.f=0 
,var.res.mal=1, var.res.fem=1)
{
## assortative mating
#simulation
#remating is implemented by having the old partners removed
# to avoid the challnge of remating a pool of individuals so that 
# they mate assortative (according to some covariance)
#mean trait values
#variances and covariances
#correlation induces by the common environment of the pair 
cov.CommEnv.m.f=rho.CommEnv.m.f*sqrt(var.CommEnv.mal*var.CommEnv.fem) #comm env may be correlated
Sigma.CommEnv <- matrix(c(var.CommEnv.mal,cov.CommEnv.m.f,cov.CommEnv.m.f,var.CommEnv.fem), ncol=2)
#

#######################################################################
## 1. create data
#######################################################################
require(mvtnorm)
data=data.frame()
for (m in 1:n.matings) {
	if (m==1) { #first mating round
		i.mal=rnorm(n=n.pairs,mean=mean.mal,sd=sqrt(var.i.mal))
		i.fem=mean.fem+(rho.i.m.f*sqrt(var.i.fem)/sqrt(var.i.mal))*(i.mal-mean.mal)
		i.fem=rnorm(n=length(i.fem),mean=i.fem,sd=sqrt(var.i.fem*(1-rho.i.m.f^2)))
		#store data
		id.n.i.mal=cbind(m,seq(1:n.pairs),i.mal)
		id.n.i.fem=cbind(m,seq(1:n.pairs),i.fem)
	} #if m==1
	if (m==2) { #second mating round: all females die and are replaced by new females who mate according to assortative mating rule
		i.fem=mean.fem+(rho.i.m.f*sqrt(var.i.fem)/sqrt(var.i.mal))*(i.mal-mean.mal)
		i.fem=rnorm(n=length(i.fem),mean=i.fem,sd=sqrt(var.i.fem*(1-rho.i.m.f^2)))
		#store data
		id.n.i.mal=rbind(id.n.i.mal,cbind(m,seq(1:n.pairs),i.mal))
		id.n.i.fem=rbind(id.n.i.fem,cbind(m,n.pairs+seq(1:n.pairs),i.fem))
	} #if m==3
	if (m==3) { #third mating round: all males die and are replaced by new males who mate according to assortative mating rule
		i.mal=mean.mal+(rho.i.m.f*sqrt(var.i.mal)/sqrt(var.i.fem))*(i.fem-mean.fem)
		i.mal=rnorm(n=length(i.mal),mean=i.mal,sd=sqrt(var.i.mal*(1-rho.i.m.f^2)))
		#store data
		id.n.i.mal=rbind(id.n.i.mal,cbind(m,n.pairs+seq(1:n.pairs),i.mal))
		id.n.i.fem=rbind(id.n.i.fem,cbind(m,n.pairs+seq(1:n.pairs),i.fem))
	} #if m==3


	for (r in 1:n.repeated.measures) {
		# the common environment and residuals are added to the individual specific values
		tmp=rmvnorm(n=n.pairs, mean=c(0,0), sigma=Sigma.CommEnv)
		c.e.mal=tmp[,1]
		c.e.fem=tmp[,2]
		res.mal=rnorm(n=n.pairs,mean=0,sd=sqrt(var.res.mal))
		res.fem=rnorm(n=n.pairs,mean=0,sd=sqrt(var.res.fem))
		data=rbind(data,
		cbind(id.n.i.mal[id.n.i.mal[,1]==m,],id.n.i.fem[id.n.i.fem[,1]==m,2:3], # mating, ID and i for mal and fem
		id.n.i.mal[id.n.i.mal[,1]==m,3]+c.e.mal+res.mal,id.n.i.fem[id.n.i.fem[,1]==m,3]+c.e.fem+res.fem #data points for mal and fem
		))
	} #for r
} #for a
## some processing
i.data<-cbind(id.n.i.mal,id.n.i.fem)
names(data)[c(2,4,6,7)]=c("ID.mal","ID.fem","trait.mal","trait.fem")
# create the sets needed for analysis
N.obs=dim(data)[1]
M=data[,"ID.mal"]
F=data[,"ID.fem"]
Mm=Fs=seq(n.pairs)
Fm=Ms=seq((1+n.pairs),2*n.pairs)
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
#package all and return as list
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
return(list(data=data,n.pairs=n.pairs, n.matings=n.matings, n.repeated.measures=n.repeated.measures
,mean.mal=mean.mal,mean.fem=mean.fem,
var.i.mal=var.i.mal, var.i.fem=var.i.fem,rho.i.m.f=rho.i.m.f
, var.CommEnv.mal=var.CommEnv.mal, var.CommEnv.fem=var.CommEnv.fem,rho.CommEnv.m.f=rho.CommEnv.m.f 
,var.res.mal=var.res.mal, var.res.fem=var.res.fem
, N.obs=N.obs, M=M, F=F, Mm=Mm, Fs=Fs, Fm=Fm,Ms=Ms))
} #end of function


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
# function for plotting the output and calculating basic stats
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
####################################################################################################################
plot.sim.data<-function(sim.data){
data=sim.data$data
# some simple plotting etc
#plot(data$trait.mal, data$trait.fem)
#select only one observation per pair (traditional analysis)
pairs=as.matrix(unique(data[,c("ID.mal","ID.fem")]))
pair.dat=data.frame()
for (p in 1:dim(pairs)[1]) {
	index=which((data[,"ID.mal"]==pairs[p,1])&(data[,"ID.fem"]==pairs[p,2]))
	pair.dat=rbind(pair.dat,data[index[1],])
}
#traditional correlation for assortative mating
assortative.trad=cor(pair.dat$trait.fem,pair.dat$trait.mal)
layout(matrix(c(1,2),1,2))
#and a plot of this
plot(trait.fem~trait.mal,data=pair.dat, main=paste("Traditional, cor = ",as.character(format(assortative.trad,digits=2))))
#assortative mating with respect to the individual specific values
assortative.i=cor(pair.dat$i.fem,pair.dat$i.mal)
#and a plot of this
plot(i.fem~i.mal,data=pair.dat, main=paste("Based on i, cor = ",as.character(format(assortative.i,digits=2))))
} #end of plotting function


#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
#examples
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
#1
#no residual variance (all c.e.)
#no ass.mating on the basis of i.mal,i.fem
#repeatability of trait=0.5
sim.data=SimAssMating(var.i.mal=2, var.i.fem=2,rho.i.m.f=0
, var.CommEnv.mal=4, var.CommEnv.fem=4,rho.CommEnv.m.f=0 
,var.res.mal=0, var.res.fem=0)
#output mod1 ok
            TRUTH  mean   2.5% 97.5%
mean.fem        3 2.905  2.642  3.16
mean.mal        2 1.907  1.674  2.15
rho.res         0 0.042 -0.031  0.12
var.i.fem       2 2.238  1.675  2.99
var.i.mal       2 1.934  1.428  2.58
var.res.fem     4 4.165  3.763  4.61
var.res.mal     4 4.394  3.948  4.88

#2
#no residual variance (all c.e.)
#strong ass.mating on the basis of common environment
#repeatability of trait=0.5
sim.data=SimAssMating(var.i.mal=2, var.i.fem=2,rho.i.m.f=0
, var.CommEnv.mal=4, var.CommEnv.fem=4,rho.CommEnv.m.f=0.8 
,var.res.mal=0, var.res.fem=0)
#output mod1 ok
            TRUTH  mean   2.5% 97.5%
mean.fem      3.0  2.9 2.63  3.14
mean.mal      2.0  1.9 1.73  2.15
rho.res       0.8  0.8 0.78  0.83
var.i.fem     2.0  2.4 1.88  3.07
var.i.mal     2.0  1.6 1.21  2.07
var.res.fem   4.0  4.3 3.91  4.78
var.res.mal   4.0  4.1 3.68  4.52

#output mod2 ok

            TRUTH  mean   2.5% 97.5%
mean.fem      3.0 2.8708  2.60  3.15
mean.mal      2.0 1.7719  1.53  2.02
rho.i         0.0 0.0062 -0.15  0.17
rho.res       0.8 0.8290  0.80  0.85
var.i.fem     2.0 2.4053  1.86  3.11
var.i.mal     2.0 1.8902  1.43  2.44
var.res.fem   4.0 4.3971  3.98  4.87
var.res.mal   4.0 4.1586  3.75  4.60



#3
#no residual variance (all c.e.)
#strong ass.mating on the basis of i.mal, i.fem
#repeatability of trait=0.5
sim.data=SimAssMating(var.i.mal=2, var.i.fem=2,rho.i.m.f=0.8
, var.CommEnv.mal=4, var.CommEnv.fem=4,rho.CommEnv.m.f=0 
,var.res.mal=0, var.res.fem=0)
#output mod1: all the variances correct, but no covariance detected 

            TRUTH  mean   2.5% 97.5%
mean.fem        3 2.989 2.77  3.21
mean.mal        2 2.005 1.80  2.21
rho.res         0 0.094 0.02  0.17
var.i.fem       2 1.604 1.16  2.13
var.i.mal       2 1.443 1.01  2.01
var.res.fem     4 3.782 3.40  4.20
var.res.mal     4 3.812 3.44  4.24

# the above is totall off, since the traditional assortative mating implies r>=0.2
plot.sim.data(sim.data)
#model 2 does the job!

            TRUTH    mean   2.5% 97.5%
mean.fem      3.0  3.0480  2.766 3.326
mean.mal      2.0  1.8721  1.599 2.154
rho.i         0.8  0.7737  0.669 0.860
rho.res       0.0 -0.0049 -0.078 0.067
var.i.fem     2.0  2.1378  1.541 2.883
var.i.mal     2.0  2.2541  1.636 2.985
var.res.fem   4.0  4.1829  3.759 4.628
var.res.mal   4.0  3.9146  3.523 4.334

#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
#code for analysis of assortative mating
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
ass.mating.mod1<-function() {
#priors and constraints
	#means
	mean.mal~dunif(-10,10)
	mean.fem~dunif(-10,10)
	#SD
	sigma.i.mal~dunif(0,10)
	sigma.i.fem~dunif(0,10)
	sigma.res.mal~dunif(0,10)
	sigma.res.fem~dunif(0,10)
	#variances
	var.i.mal<-pow(sigma.i.mal,2)
	var.i.fem<-pow(sigma.i.fem,2)
	var.res.mal<-pow(sigma.res.mal,2)
	var.res.fem<-pow(sigma.res.fem,2)
	#correlation
	#rho.i~dunif(-1,1) #correlation individual-specific value mal and fem
	rho.res~dunif(-1,1) #correlation common environment mal and fem
	#precision
	tau.i.mal<-pow(sigma.i.mal,-2)
	tau.i.fem<-pow(sigma.i.fem,-2)
	#rescaled precision for correlated traits
	#tau2.i.mal<-tau.i.mal/(1.-pow(rho.i,2))
	#tau2.i.fem<-tau.i.fem/(1.-pow(rho.i,2))
	# Constructing the covariance residual matrix and the corresponding precision matrix.
    	cov[1,1] <- sigma.res.mal * sigma.res.mal
    	cov[1,2] <- sigma.res.mal * sigma.res.fem * rho.res
    	cov[2,1] <- sigma.res.fem * sigma.res.mal * rho.res
    	cov[2,2] <- sigma.res.fem * sigma.res.fem
    	prec[1:2,1:2] <- inverse(cov[,])
	#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
	# model
	# monogamous pairs (not implemented)
	# males mated with multiple partners: Mm is the set of ID for these males, Fs female IDs of fems mated once in their life
	for (m in 1:Nfs) {
		i.mal[Mm[m]]~dnorm(mean.mal,tau.i.mal)
		#females mated with these males
		i.fem[Fs[m]]~dnorm(mean.fem,tau.i.fem)
	}
	# males mated with multiple partners: Mm is the set of ID for these males, Fm female IDs of fems mated multiple times
	for (n in 1:Nfs) {
		#females mated with these males
		i.fem[Fm[n]]~dnorm(mean.fem,tau.i.fem)
	}
	# females mated with multiple partners: Fm is the set of ID for these females, Ms male IDs single mated males
	for (f in 1:Nms) {
		#males mated with these females
		i.mal[Ms[f]]~dnorm(mean.mal,tau.i.mal)
	}
	# data
	for (k in 1:N.obs) {
		mu[k,1]<-i.mal[M[k]]
		mu[k,2]<-i.fem[F[k]]
		y[k,1:2] ~ dmnorm(mu[k,1:2], prec[ , ])
	}
} #model end

require(R2jags)

# prepare the data for analysis
#bundle data
jags.data<-list(
y=as.matrix(sim.data$data[,c("trait.mal","trait.fem")])
, N.obs=sim.data$N.obs
, M=sim.data$M
, F=sim.data$F
, Ms=sim.data$Ms
, Mm=sim.data$Mm
, Fs=sim.data$Fs
, Fm=sim.data$Fm
, Nfs=length(sim.data$Fs)
, Nms=length(sim.data$Ms)
)

#initial values
inits<-function(){list(
sigma.i.mal=runif(3,0,10), sigma.i.fem=runif(3,0,10),
sigma.res.mal=runif(3,0,10), sigma.res.fem=runif(3,0,10)
)}
#mean.mal=runif(3,-10,10), mean.fem=runif(3,-1,1),
#,rho.i=runif(3,-1,1), rho.res=runif(3,-1,1)

#monitored parameters
params<-c("mean.mal","mean.fem","var.i.mal","var.i.fem","var.res.mal","var.res.fem","rho.res")

#params<-c("mean.mal","mean.fem","var.i.mal","var.i.fem","rho.i","var.res.mal","var.res.fem","rho.res")
#MCMC setting
ni<-10
nt<-2
nb<-2
nc<-3

#Call jags auto
as.out1<-jags(jags.data,,params,model.file=ass.mating.mod1,n.chains=nc,n.iter=ni)
#updating (burn in)
as.out1.upd<-update(as.out1, n.iter=100)
#some autojags rounds
as.out1.auto<- autojags(as.out1.upd)
TRUTH=c(sim.data$mean.fem,sim.data$mean.mal,sim.data$rho.CommEnv.m.f,sim.data$var.i.fem,sim.data$var.i.mal,
(sim.data$var.res.fem+sim.data$var.CommEnv.fem),(sim.data$var.res.mal+sim.data$var.CommEnv.mal))
print(cbind(TRUTH,as.out1.auto$BUGSoutput$summary[2:8,c(1,3,7)]),digits=2)



#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
#code for analysis of assortative mating: including the cov between i.mal and i.fem
#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
ass.mating.mod2<-function() {
#priors and constraints
	#means
	mean.mal~dunif(-10,10)
	mean.fem~dunif(-10,10)
	#SD
	sigma.i.mal~dunif(0,10)
	sigma.i.fem~dunif(0,10)
	sigma.res.mal~dunif(0,10)
	sigma.res.fem~dunif(0,10)
	#variances
	var.i.mal<-pow(sigma.i.mal,2)
	var.i.fem<-pow(sigma.i.fem,2)
	var.res.mal<-pow(sigma.res.mal,2)
	var.res.fem<-pow(sigma.res.fem,2)
	#correlation
	rho.i~dunif(-1,1) #correlation individual-specific value mal and fem
	rho.res~dunif(-1,1) #correlation common environment mal and fem
	#precision
	tau.i.mal<-pow(sigma.i.mal,-2)
	tau.i.fem<-pow(sigma.i.fem,-2)
	#rescaled precision for correlated traits
	tau2.i.mal<-tau.i.mal/(1.-pow(rho.i,2))
	tau2.i.fem<-tau.i.fem/(1.-pow(rho.i,2))
	# Constructing the covariance residual matrix and the corresponding precision matrix.
    	cov[1,1] <- sigma.res.mal * sigma.res.mal
    	cov[1,2] <- sigma.res.mal * sigma.res.fem * rho.res
    	cov[2,1] <- sigma.res.fem * sigma.res.mal * rho.res
    	cov[2,2] <- sigma.res.fem * sigma.res.fem
    	prec[1:2,1:2] <- inverse(cov[,])
	#¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤
	# model
	# monogamous pairs (not implemented)
	# males mated with multiple partners: Mm is the set of ID for these males, Fs female IDs of fems mated once in their life
	for (m in 1:Nfs) {
		i.mal[Mm[m]]~dnorm(mean.mal,tau.i.mal)
		#females mated with these males
		tmp.fem[Fs[m]]<-mean.fem+(rho.i*sigma.i.fem/sigma.i.mal)*(i.mal[Mm[m]]-mean.mal)
		i.fem[Fs[m]]~dnorm(tmp.fem[Fs[m]],tau2.i.fem)
	}
	# males mated with multiple partners: Mm is the set of ID for these males, Fm female IDs of fems mated multiple times
	for (n in 1:Nfs) {
		#females mated with these males
		tmp.fem[Fm[n]]<-mean.fem+(rho.i*sigma.i.fem/sigma.i.mal)*(i.mal[Mm[n]]-mean.mal)
		i.fem[Fm[n]]~dnorm(tmp.fem[Fm[n]],tau2.i.fem)
	}
	# females mated with multiple partners: Fm is the set of ID for these females, Ms male IDs single mated males
	for (f in 1:Nms) {
		#males mated with these females
		tmp.mal[Ms[f]]<-mean.mal+(rho.i*sigma.i.mal/sigma.i.fem)*(i.fem[Fm[f]]-mean.fem)
		i.mal[Ms[f]]~dnorm(tmp.mal[Ms[f]],tau2.i.mal)
	}
	# data
	for (k in 1:N.obs) {
		mu[k,1]<-i.mal[M[k]]
		mu[k,2]<-i.fem[F[k]]
		y[k,1:2] ~ dmnorm(mu[k,1:2], prec[ , ])
	}
} #model end

require(R2jags)

# prepare the data for analysis
#bundle data
jags.data<-list(
y=as.matrix(sim.data$data[,c("trait.mal","trait.fem")])
, N.obs=sim.data$N.obs
, M=sim.data$M
, F=sim.data$F
, Ms=sim.data$Ms
, Mm=sim.data$Mm
, Fs=sim.data$Fs
, Fm=sim.data$Fm
, Nfs=length(sim.data$Fs)
, Nms=length(sim.data$Ms)
)

#initial values

#monitored parameters
params<-c("mean.mal","mean.fem","var.i.mal","var.i.fem","rho.i","var.res.mal","var.res.fem","rho.res")
#MCMC setting
ni<-10
nt<-2
nb<-2
nc<-3

#Call jags auto
as.out2<-jags(jags.data,,params,model.file=ass.mating.mod2,n.chains=nc,n.iter=ni)
#updating (burn in)
as.out2.upd<-update(as.out2, n.iter=100)
#some autojags rounds
as.out2.auto<- autojags(as.out2.upd)
TRUTH=c(sim.data$mean.fem,sim.data$mean.mal,sim.data$rho.i,sim.data$rho.CommEnv.m.f,sim.data$var.i.fem,sim.data$var.i.mal,
(sim.data$var.res.fem+sim.data$var.CommEnv.fem),(sim.data$var.res.mal+sim.data$var.CommEnv.mal))
print(cbind(TRUTH,as.out2.auto$BUGSoutput$summary[2:9,c(1,3,7)]),digits=2)


