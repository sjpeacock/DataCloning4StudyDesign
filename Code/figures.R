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

# setwd("DataCloning4StudyDesign")

library(dclone)
library(gplots)
library(gtools)

#Color scheme for all figures
colR<-c("#d7191c", 1, "#2c7bb6")

#########################################################
# Data cloning analysis
######################################################### 

# List to store three different data sets and 3 different priors
allMod<-list(); length(allMod)<-3
allDat<-list();length(allX)<-3
allX<-list();length(allX)<-3

# 1) Less-spread data
load('Code/Workspaces/20160125_lessSpreadData_prior1.RData')
allMod[[1]]<-dc.mod
allDat[[1]]<-data
allX[[1]]<-x

# 2) Original data
load('Code/Workspaces/20160125_originalData_prior1.RData')
allMod[[2]]<-dc.mod
allDat[[2]]<-data
allX[[2]]<-x

# 3) More-spread data
load('Code/Workspaces/20160125_moreSpreadData_prior1.RData')
allMod[[3]]<-dc.mod
allDat[[3]]<-data
allX[[3]]<-x

#---------------------------------------------------------
# Fig. 2: Data cloning results: parameters of interest
#---------------------------------------------------------

DC<-list(); length(DC)<-3
for(i in 1:3){ DC[[i]]<-dctable(allMod[[i]])}

# Plot of varaiance over number of clones
par.names<-c(expression(italic(D)), expression(lambda[c]), expression(italic(L[h])), expression(italic(L[m])), expression(alpha), expression(kappa), expression(italic(s[c])), expression(italic(s[h]))) 

quartz(width=3.2, height=4, pointsize=9)
par(mfrow=c(2,1), mar=c(1,1,2,1), oma=c(2,3,0,0));
k<-0
for(i in 1:2){
	k<-k+1
	plot(1:10, 1/c(1:10), "l", lty=2, bty="l", xlab="", ylab="", ylim=c(0,1), xaxt="n", las=1)
	if(i==1) mtext(side=3, adj=0, expression(paste("a) Ambient source strength (", kappa*beta*v^-1, ")")), line=0.5)
	if(i==2) mtext(side=3, adj=0, expression(paste("b) Farm source strength (", alpha*beta*v^-1, ")")), line=0.5)
	for(j in 1:3){
		lines(DC[[j]][[c(6,5)[i]]]$n.clones, DC[[j]][[c(6,5)[i]]]$sd^2/(DC[[j]][[c(6,5)[i]]]$sd[1]^2), col=colR[j])
		points(DC[[j]][[c(6,5)[i]]]$n.clones, DC[[j]][[c(6,5)[i]]]$sd^2/(DC[[j]][[c(6,5)[i]]]$sd[1]^2), pch=c(21:24)[j], bg=c("white",colR[j])[1*(round(DC[[j]][[c(6,5)[i]]]$r.hat,1)==1)+1], col=colR[j], cex=1.2)
	}
	
	if(i==2) axis(side=1) else axis(side=1, labels=FALSE)

	if(i==2) legend(6, 1.1, pch=c(21:24), lwd=1, lty=1, col=colR, c("Less-spread", "Original", "More-spread"), bty="n", pt.bg=colR, xpd=NA)
	
	}
mtext(side=1, outer=TRUE, "Number of clones (K)", line=1)
mtext(side=2, outer=TRUE, "Scaled variance", line=2)

#---------------------------------------------------------
# Fig. 3: Data cloning results: other parameters
#---------------------------------------------------------
quartz(width=6, height=3.5, pointsize=9)
par(mfrow=c(2,3), mar=c(1,1,2,0), oma=c(4,4,0,1));
j<-1
for(i in 1:6){
	I<-c(1,2,7:8,3:4)[i]
	plot(1:10, 1/c(1:10), "l", lty=2, bty="l", xlab="", ylab="", ylim=c(0,1), xaxt="n", yaxt="n")
	if(i==1) mtext(side=3, adj=0, expression(paste("a) ", italic(D))), line=0.5)
	if(i==2) mtext(side=3, adj=0, expression(paste("b) ", lambda[c])), line=0.5)
	if(i==3) mtext(side=3, adj=0, expression(paste("c) ", italic(s[c]))), line=0.5)
	if(i==4) mtext(side=3, adj=0, expression(paste("d) ", italic(s[h]))), line=0.5)
	if(i==5) mtext(side=3, adj=0, expression(paste("e) ", italic(L[h]))), line=0.5)
	if(i==6) mtext(side=3, adj=0, expression(paste("f) ", italic(L[m]))), line=0.5)
	
	for(j in 1:3){
		lines(DC[[j]][[I]]$n.clones, DC[[j]][[I]]$sd^2/(DC[[j]][[I]]$sd[1]^2), col=colR[j])
		points(DC[[j]][[I]]$n.clones, DC[[j]][[I]]$sd^2/(DC[[j]][[I]]$sd[1]^2), pch=c(21:24)[j], bg=c("white",colR[j])[1*(round(DC[[j]][[I]]$r.hat,1)==1)+1], col=colR[j], cex=1.2)
	}
	
	if(is.element(i, c(4:6))) axis(side=1, cex.axis=1.4) else axis(side=1, labels=FALSE)
	if(is.element(i, c(1,4))) axis(side=2, cex.axis=1.4, las=1) else axis(side=2, labels=FALSE)

	}
mtext(side=1, outer=TRUE, "Number of clones (K)", line=2)
mtext(side=2, outer=TRUE, "Scaled variance", line=2.5)


######################################################### 
# Fig 4: Bivariate plots
######################################################### 
par.names<-c(expression(italic(D)), expression(lambda[c]), expression(italic(L[h])), expression(italic(L[m])), expression(alpha), expression(kappa), expression(italic(s[c])), expression(italic(s[h]))) 

n.points<-100
n.seq<-round(seq(1,1000,length.out=n.points))

ranges<-matrix(nrow=8, ncol=2)
for(i in 1:8){
	ranges[i,]<-range(c(allMod[[1]][[1]][n.seq,i], allMod[[1]][[2]][n.seq,i], allMod[[1]][[3]][n.seq,i], allMod[[2]][[1]][n.seq,i], allMod[[2]][[2]][n.seq,i], allMod[[2]][[3]][n.seq,i], allMod[[3]][[1]][n.seq,i], allMod[[3]][[2]][n.seq,i], allMod[[3]][[3]][n.seq,i]))
	}

# All parameters
quartz(width=6.3, height=6)
par(mfrow=c(7,7), mar=c(0,0,0,0), mgp=c(2,1,0), oma=c(1,4,3,1))
for(j in 1:7){
	for(i in 2:8){	
		if(i>j){
			plot(1,1,"n",  xlab="", ylab="", xaxt="n", yaxt="n", xlim=ranges[i,], ylim=ranges[j,])
			for(k in 1:3){
				for(m in 1:3){ 
					points(allMod[[k]][[m]][n.seq,i], allMod[[k]][[m]][n.seq,j], col=colR[k], pch=c(1,0,5)[k])
					}}
			if(i==j+1){
				mtext(side=2, par.names[j], line=1, las=1)
				}
			if(j==1) mtext(side=3, line=0.5, par.names[i])
			
		} else plot(1,1,"n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
	}
}

legend(-27.2, 27.2, pch=c(1,0,5), col=colR, c("Less spread", "Original data", "More spread"), xpd=NA, cex=1.2)



######################################################### 
# Fig 5: Model fits to data
######################################################### 
source("Code/sim_model.R")

boot<-function(x, n=1000){
	xi<-matrix(sample(x, n*length(x), replace=TRUE), nrow=n, ncol=length(x))
	mean.x<-apply(xi, 1, mean, na.rm=TRUE)
	ret<-c(mean(x, na.rm=TRUE), quantile(mean.x, c(0.025, 0.975), na.rm=TRUE))
	return(ret)
}

Z<-list();length(Z)<-9; dim(Z)<-c(3,3)
for(i in 1:3){for(j in 1:3) Z[[j,i]]<-matrix(nrow=length(x), ncol=3)}
for(j in 1:3){# for each dataset
	for(i in 1:length(x)){
		Z[[j,1]][i,]<-boot(allDat[[j]]$C[which(allDat[[j]]$site_num==i)])
		Z[[j,2]][i,]<-boot(allDat[[j]]$H[which(allDat[[j]]$site_num==i)])
		Z[[j,3]][i,]<-boot(allDat[[j]]$M[which(allDat[[j]]$site_num==i)])
		}
	}	
	
allParams<-list(); length(allParams)<-3
for(i in 1:3){ allParams[[i]]<-summary(allMod[[i]])[[1]][,1]}

allPred<-list(); length(allPred)<-3
for(i in 1:3){ allPred[[i]]<-sim(allParams[[i]], x=c(-60:60)) }

col2<-rep(colR[2], 16); col2[c(1,2,3)]<-colR[1]; col2[c(5,13)]<-colR[3]

quartz(width=6.3, height=4, pointsize=11)
par(mfrow=c(3,3), mar=c(1,1,1,1), oma=c(3,4,2,1))
for(j in 1:3){
	for(i in 1:3){			
		plotCI(allX[[i]], Z[[i,j]][,1], "n", li= Z[[i,j]][,2], ui= Z[[i,j]][,3], xlab="", ylab="", bty="l", gap=0.3, xlim=c(-40,40), ylim=c(0, c(0.2, 0.7, 0.8)[j]), yaxt="n", xaxt="n", col=NA)
		
		lines(c(-60:60), allPred[[i]][,j], col=colR[i])
		abline(v=0, lty=2)
		
		plotCI(allX[[i]], Z[[i,j]][,1], li= Z[[i,j]][,2], ui= Z[[i,j]][,3], gap=0.3, pch=21, pt.bg="white", col=colR[i], add=TRUE)
		if(j==1) mtext(side=3, c("Less-spread", "Original", "More-spread")[i], line=1)
		if(j==3) axis(side=1) else axis(side=1, labels=FALSE)
		if(i==1) axis(side=2, las=1) else axis(side=2, labels=FALSE)
		if(i==1) mtext(side=2, c("C(x)", "H(x)", "M(x)")[j], line=3)
		
		if(i==1) points(allX[[i]][c(4,7,10)], Z[[i,j]][c(4,7,10),1], col=colR[1], pch=19)
		if(i==2) points(allX[[i]][c(1,2,3)], Z[[i,j]][c(1,2,3),1], col=1, pch=21, bg=colR[1])
		if(i==2) points(allX[[i]][c(5,13)], Z[[i,j]][c(5,13),1], col=1, pch=21, bg=colR[3])
		if(i==3) points(allX[[i]][c(1,2)], Z[[i,j]][c(1,2),1], col=colR[3], pch=19)
		
		}
}
mtext(side=1, outer=TRUE, "Distance along migration (km)", line=1.5)

######################################################### 
# Fig 6: Parameter estimates
######################################################### 

allParams<-list(); length(allParams)<-3*3; dim(allParams)<-c(3,3) # dataset, metric (mean, lci, uci)
for(i in 1:3){
	if(length(allMod[[i]])>0){
			allParams[[i,1]]<-summary(allMod[[i]])[[1]][,1]
			allParams[[i,2]]<-summary(allMod[[i]])[[2]][,1]
			allParams[[i,3]]<-summary(allMod[[i]])[[2]][,5]
		}
	}

rangeP<-matrix(0, nrow=8, ncol=2)
for(k in 1:8){
	rangeP[k,1]<-min(allParams[[1,2]][k])
	for(i in 2:3){
		rangeP[k,1]<-min(c(rangeP[k,1], allParams[[i,2]][k]))
		}
	rangeP[k,2]<-max(allParams[[1,3]][k])
	for(i in 2:3){
		rangeP[k,2]<-max(c(rangeP[k,2], allParams[[i,3]][k]))
		}
	}

quartz(width=5, height=2.5, pointsize=11)
par(mar=c(3,1,1,1), oma=c(0,0,0,0))
for(k in 1:8){
	plotCI(k, allParams[[1,1]][k], col=NA, xlim=c(1,8.8), ylim=rangeP[k,], yaxt="n", xaxt="n", bty="n")
	plotCI(k+0.2, allParams[[1,1]][k], li=allParams[[1,2]][k], ui= allParams[[1,3]][k], col=colR[1], gap=0.3, add=TRUE, pch=21, pt.bg="white", sfrac=0.008)
		plotCI(k+0.375, allParams[[2,1]][k], li=allParams[[2,2]][k],ui= allParams[[2,3]][k], col=colR[2], gap=0.3, add=TRUE, pch=22, pt.bg="white", sfrac=0.008)
		plotCI(k+0.55, allParams[[3,1]][k], li=allParams[[3,2]][k],ui= allParams[[3,3]][k], col=colR[3], gap=0.3, add=TRUE, pch=23, pt.bg="white", xpd=NA, sfrac=0.008)
	if(k<8) par(new=TRUE)
	}
axis(side=1, at=c(1:8+0.3), labels=par.names, tick=FALSE, line=-0.5)	
axis(side=1, at=1:9-0.15, labels=FALSE)
#legend(8.8, 5, col=colR, lwd=1, pch=21:23, pt.bg=colR, c("Less spread", "True", "More spread"), xpd=NA, bty="n", pt.cex=0.8, cex=0.8)
mtext(side=2, "Parameter estimate")


######################################################### 
# Fig S1: Prior comparison
######################################################### 
priorInd<-c(5,3,6,6,4,1,2,2)

rangeP2<-matrix(0, nrow=8, ncol=2)
for(k in 1:8){
	rangeP2[k,1]<-min(allParams[[1,1,2]][k])
	for(i in 2:3){for(j in 2:3){ rangeP2[k,1]<-min(c(rangeP2[k,1], allParams[[i,j,2]][k]))}}
	for(i in 1:2){rangeP2[k,1]<-min(c(rangeP2[k,1], allPrior[[i]][1,priorInd[k]]))}
		
	rangeP2[k,2]<-max(allParams[[1,1,3]][k])
	for(i in 2:3){for(j in 2:3){ rangeP2[k,2]<-max(c(rangeP2[k,2], allParams[[i,j,3]][k]))}}
	for(i in 1:2){rangeP2[k,2]<-max(c(rangeP2[k,2], allPrior[[i]][1,priorInd[k]]))}
	
}

quartz(width=6.3, height=4, pointsize=11)
par(mfrow=c(2,4), mar=c(3,3,2,0), oma=c(2,2,2,1))
for(k in 1:8){
	P<-seq(rangeP2[k,1], rangeP2[k,2], length.out=100)
	plot(P, dnorm(P, mean=allPrior[[1]][1,priorInd[k]], sd=allPrior[[1]][2,priorInd[k]]), "l", bty="l", las=1, xlab="", ylab="", col=grey(0.8), lwd=2)
	lines(P, dnorm(P, mean=allPrior[[2]][1,priorInd[k]], sd=allPrior[[2]][2,priorInd[k]]), lty=2, col=grey(0.8), lwd=2)
	# lines(P, dnorm(P, mean=allPrior[[3]][1,priorInd[k]], sd=allPrior[[3]][2,priorInd[k]]), lty=3, col=grey(0.8), lwd=2)
	
	for(i in 1:3){for(j in 1:3){abline(v=allParams[[i,j,1]][k], col=colR[i], lty=j, lwd=1.2)}}
	
	if(k==1) mtext(side=3, adj=0, expression(paste("a) ", italic(D))), line=0.5)
	if(k==2) mtext(side=3, adj=0, expression(paste("b) ", lambda[c])), line=0.5)
	if(k==3) mtext(side=3, adj=0, expression(paste("c) ", italic(L[h]))), line=0.5)
	if(k==4) mtext(side=3, adj=0, expression(paste("d) ", italic(L[m]))), line=0.5)
	if(k==5) mtext(side=3, adj=0, expression(paste("e) ", kappa*beta*v^-1)), line=0.5)
	if(k==6) mtext(side=3, adj=0, expression(paste("f) ", alpha*beta*v^-1)), line=0.5)
	if(k==7) mtext(side=3, adj=0, expression(paste("g) ", italic(s[c]))), line=0.5)
	if(k==8) mtext(side=3, adj=0, expression(paste("h) ", italic(s[h]))), line=0.5)
	
	
	}

