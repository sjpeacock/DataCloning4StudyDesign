###################################################################################
# This code accompanies the paper:
#
# Study design and parameter estimability for spatial and temporal ecological models
#
# Submitted May 11, 2016 to Methods in Ecology and Evolution by
# Stephanie J Peacock, Martin Krkosek, Mark Lewis, and Subhash Lele
# Questions should be directed to corresponding author SJ Peacock at
# stephanie.j.peacock at gmail.com
###################################################################################

#########################################################
#########################################################

rm(list=ls())

library(dclone)
library(gplots)
library(parallel)
library(boot)

#########################################################
# 1) Read in data and define parameters
########################################################
yr<-2003; r<-1
data<-read.delim("Data/Leps.txt", header=TRUE)
data<-subset(data, year==yr&rep==r&route!="Lower Knight")

#Calculate vector of distance for each data point
sum.data<-read.delim("Data/Summary.txt", header=TRUE)
sum.data<-subset(sum.data, year==yr&rep==r&route!="Lower Knight")

# #------------------------------------------------------
# # November 1, 2015: Try changing the location of some points to spread out data spatially
# # Might help estimability of f vs k?
data$distance<-sort(sum.data$distance)[data$site_num]

set.seed(289365) #get reproducible data for the moved sites
# #change site #13 (29.5 kms) to -30 kms
ind<-which(data$site_num==13)
data$C[ind]<-rpois(length(ind), lambda=0.007178194)
data$H[ind]<-rpois(length(ind), lambda=0.03365534)
data$M[ind]<-rpois(length(ind), lambda=0.02384308)
data$distance[ind]<-rep(-30, length(ind))
data$site_num[ind]<-0

# #changes site #5 to -40 kms
ind<-which(data$site_num==5)
data$C[ind]<-rpois(length(ind), lambda=0.003922110)
data$H[ind]<-rpois(length(ind), lambda=0.02263692)
data$M[ind]<-rpois(length(ind), lambda=0.02233066)
data$distance[ind]<-rep(-40, length(ind))
data$site_num[ind]<--1

data$site_num2<-rep(0, dim(data)[1])
I<-1
for(i in 1:length(sort(sum.data$distance))){
	data$site_num2[data$site_num==sort(unique(data$site_num))[i]]<-I
	I<-I+1
}

data$site_num<-data$site_num2

dist<-sort(unique(data$distance))

prior.mean<-rbind(c(-7, 2, 1, 0, 3.5, 3), c(-5, 0, 5, 4, 5, 1.2), c(-10, 3.5, 2.5, -3, 2, 5), c(-2, -1, -1, 2, 0, 0))
prior.sd<-rbind(rep(0.5, 6), c(0.3, 0.7, 1.5, 0.3, 0.3, 0.5), c(1, 0.3, 0.3, 1.2, 1, 1), rep(0.5, 6))

# par(mfrow=c(2,3), mar=c(3,3,2,1), oma=c(2,2,1,0))
# for(i in 1:6){
	# xmin<-min(prior.means[,i]-1.96)
	# xmax<-max(prior.means[,i]+1.96)
	# x<-seq(xmin, xmax, length.out=1000)	
	# plot(x, dnorm(x, mean=prior.means[1,i], sd=prior.sd[1,i]), "l", bty="l", main=pnames[i], yaxt="n")
	# for(j in 2:4) lines(x, dnorm(x, mean=prior.means[j,i], sd=prior.sd[j,i]), col=j)	
# }

set.seed(3495)		
init.fun<-function(){
	return(list(
		kappa.log=rnorm(1, prior.mean[1,1], 1), 
		sc.logit=rnorm(1, prior.mean[1,2], 1), 
		sh.logit=rnorm(1, prior.mean[1,2], 1), 
		Lc.log=rnorm(1, prior.mean[1,3], 1), 
		alpha.log=rnorm(1, prior.mean[1,4], 1), 
		D.log=rnorm(1, prior.mean[1,5], 1),
		Lh_R.log=rnorm(1, prior.mean[1,6], 1),
		Lm_R.log=rnorm(1, prior.mean[1,6], 1)))
		}

init.vals<-list(
	list(init.fun(), init.fun(), init.fun()),
	list(init.fun(), init.fun(), init.fun()),
	list(init.fun(), init.fun(), init.fun()),
	list(init.fun(), init.fun(), init.fun()))

#########################################################
# 5) Data cloning
#########################################################
source("model_2003.R")
setwd("allPriors")

# Need to export data, prior.mean, prior.sd, init.vals, model
fitPrior<-function(x){
	
	library(dclone)
	
	dat<-list(y=0, gamma=1.56, u_n=4/5, u_c=1/5, x=dist, L=c(data$C, data$H, data$M), site=rep(data$site_num, 3), stage=c(rep(1, dim(data)[1]), rep(2, dim(data)[1]), rep(3, dim(data)[1])), prior=rbind(prior.mean[x,], 1/(prior.sd[x,])^2))

	cl<-makeCluster(3, type="SOCK")
	K<-c(10,15,20)	
	t.start<-proc.time()

	dc.mod<-dc.parfit(cl,		#cluster
					data=dat, 
					params=c("kappa.log", "sc.logit", "sh.logit", "Lc.log", "alpha.log", "D.log", "Lh_R.log", "Lm_R.log"),
					inits=init.vals[[x]],
					model=model, n.adapt=5000, n.update=40000,n.iter=20000, thin=20,
					n.chains=3, unchanged=c("y", "gamma", "u_n", "u_c", "prior"),
					n.clones=K)
	all.time<-(proc.time()-t.start)[3]/60
	
	stopCluster(cl)

	return(list(dc.mod, all.time))
	
	} #end fitPrior function

#---------------------------------------------------
t0<-proc.time()
cl0<-makeCluster(4)
clusterExport(cl0, varlist=list("data", "dist", "init.vals", "model", "prior.mean", "prior.sd"))
X<-clusterApply(cl0, x=c(1:4), fun=fitPrior)
stopCluster(cl0)
cat("Process time (minutes) = ", (proc.time()-t0)[3]/60, "\n")
	
save.image("20160113_allPriors_moved.RData")