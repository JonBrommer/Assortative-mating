rm(list=ls())
setwd("C:\\Users\\joegbr\\Documents\\data\\Tutkimus\\projects\\NielsDingemanse\\AssortativeMating")
## assortative mating
#simulation
#remating is implemented by having the old partners removed
# to avoid the challnge of remating a pool of individuals so that 
# they mate assortative (according to some covariance)
n.pairs=1000
n.matings=3
n.repeated.measures=3
#mean trait values
mean.mal=10
mean.fem=10
#variances and covariances
var.i.mal=1
var.i.fem=1
rho.i.m.f=0.5 #correlation between individual specific values that mate
var.CommEnv.mal=1 #common environment
var.CommEnv.fem=1
rho.CommEnv.m.f=0.8 #correlation induces by the common environment of the pair 
cov.CommEnv.m.f=rho.CommEnv.m.f*sqrt(var.CommEnv.mal*var.CommEnv.fem) #comm env may be correlated
Sigma.CommEnv <- matrix(c(var.CommEnv.mal,cov.CommEnv.m.f,cov.CommEnv.m.f,var.CommEnv.fem), ncol=2)
var.res.mal=1
var.res.fem=1 #residuals are not correlated by definition
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
		i.fem=rnorm(n=length(i.fem),mean=i.fem,sd=sqrt(var.i.fem))
		#store data
		id.n.i.mal=cbind(m,seq(1:n.pairs),i.mal)
		id.n.i.fem=cbind(m,seq(1:n.pairs),i.fem)
	} #if m==1
	if (m==2) { #second mating round: all females die and are replaced by new females who mate according to assortative mating rule
		i.fem=mean.fem+(rho.i.m.f*sqrt(var.i.fem)/sqrt(var.i.mal))*(i.mal-mean.mal)
		i.fem=rnorm(n=length(i.fem),mean=i.fem,sd=sqrt(var.i.fem))
		#store data
		id.n.i.mal=rbind(id.n.i.mal,cbind(m,seq(1:n.pairs),i.mal))
		id.n.i.fem=rbind(id.n.i.fem,cbind(m,n.pairs+seq(1:n.pairs),i.fem))
	} #if m==3
	if (m==3) { #third mating round: all males die and are replaced by new males who mate according to assortative mating rule
		i.mal=mean.mal+(rho.i.m.f*sqrt(var.i.mal)/sqrt(var.i.fem))*(i.fem-mean.fem)
		i.mal=rnorm(n=length(i.mal),mean=i.mal,sd=sqrt(var.i.mal))
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
####################################################################################################################
## analysis
####################################################################################################################
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



