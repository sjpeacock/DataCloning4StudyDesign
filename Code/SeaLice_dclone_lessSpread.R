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

# rm(list=ls())

library(dclone)
library(gplots)
library(parallel)
library(boot)

#########################################################
# 1) Read in data and define parameters
########################################################
data<-read.delim("Data/Leps.txt", header=TRUE)

#Calculate vector of distance for each data point
sum.data<-read.delim("Data/Summary.txt", header=TRUE)

# #------------------------------------------------------
# Try moving some of the points closer together to see if estimability gets worse
data$distance<-sort(sum.data$distance)[data$site_num]

set.seed(289365) #get reproducible data for the moved sites
# change site #1 (-20 kms) to 0 kms
ind<-which(data$site_num==1)
data$C[ind]<-rpois(length(ind), lambda=0.055059605)
data$H[ind]<-rpois(length(ind), lambda=0.25548978)
data$M[ind]<-rpois(length(ind), lambda=0.06294714)
data$distance[ind]<-rep(0, length(ind))
data$site_num[ind]<-6.5
# changes site #2 (-15.5 km) to 5 kms
ind<-which(data$site_num==2)
data$C[ind]<-rpois(length(ind), lambda=0.057602022)
data$H[ind]<-rpois(length(ind), lambda=0.33943431)
data$M[ind]<-rpois(length(ind), lambda=0.08682123)
data$distance[ind]<-rep(5, length(ind))
data$site_num[ind]<-8.5

# changes site #3 (11.5 km) to 16 kms
ind<-which(data$site_num==3)
data$C[ind]<-rpois(length(ind), lambda=0.043019509)
data$H[ind]<-rpois(length(ind), lambda=0.45076160)
data$M[ind]<-rpois(length(ind), lambda=0.18972329)
data$distance[ind]<-rep(16, length(ind))
data$site_num[ind]<-10.5

data$site_num2<-rep(0, dim(data)[1])
I<-1
for(i in 1:length(sort(sum.data$distance))){
	data$site_num2[data$site_num==sort(unique(data$site_num))[i]]<-I
	I<-I+1
}

data$site_num<-data$site_num2
x<-sort(unique(data$distance))

# #------------------------------------------------------


# Physical parameters
y<-0				# 2003 Farm location
gamma<-1.56			# Advection parameter
u_n<-4/5			# Nauplii mortality
u_c<-1/5			# Cope mortality			

		
prior.mean<-c(-7, 2, 1, 0, 3.5, 3)
prior.sd<-0.5

set.seed(3495)		
init.fun<-function(){
	return(list(
		kappa.log=rnorm(1, prior.mean[1], 1), 
		sc.logit=rnorm(1, prior.mean[2], 1), 
		sh.logit=rnorm(1, prior.mean[2], 1), 
		Lc.log=rnorm(1, prior.mean[3], 1), 
		alpha.log=rnorm(1, prior.mean[4], 1), 
		D.log=rnorm(1, prior.mean[5], 1),
		Lh_R.log=rnorm(1, prior.mean[6], 1),
		Lm_R.log=rnorm(1, prior.mean[6], 1)))
		}
init.vals<-list(init.fun(), init.fun(), init.fun())

#########################################################
# 5) Data cloning
#########################################################
source("Code/model.R")

dat<-list(y=y, gamma=1.56, u_n=4/5, u_c=1/5, x=x, L=c(data$C, data$H, data$M), site=rep(data$site_num, 3), stage=c(rep(1, dim(data)[1]), rep(2, dim(data)[1]), rep(3, dim(data)[1])), prior=rbind(prior.mean, 1/(prior.sd)^2))


K<-c(1:25)
cl<-makeCluster(4, type="SOCK")
t.start<-proc.time()

dc.mod<-dc.parfit(cl,		#cluster
				data=dat, 
				params=c("kappa.log", "sc.logit", "sh.logit", "Lc.log", "alpha.log", "D.log", "Lh_R.log", "Lm_R.log"),
				#params=c("k", "sc", "sh", "Lc", "f", "D", "Lh_R", "Lm_R"),
				model=model,n.adapt=5000, n.update=40000,n.iter=20000,  thin=20,
				n.chains=3, unchanged=c("y", "gamma", "u_n", "u_c", "prior"),
				n.clones=K)
all.time<-(proc.time()-t.start)[3]/60
cat("Process time (minutes) = ", all.time)

stopCluster(cl)
