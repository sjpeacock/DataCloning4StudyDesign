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

#################################################################################################
# This file plots the supplemental figures looking at the parameter estimates under
# four different prior assumptions.
#################################################################################################
library(dclone)
library(gplots)

#setwd("~/Google Drive/Data cloning/Sea lice model/MEE paper/Code/DataCloning4StudyDesign/Supplement")

# Select which dataset to plot (1 = less spread, 2 = Original, 3 = more spread)
dataset<-3

if(dataset==1){
	load('20160113_lessSpreadData_allPriors.RData')
	titles<-"Less-spread simulated data"
	filenames<-"lessSpread_allPriors.pdf"
	}
if(dataset==2){
	load('20160112_originalData_allPriors.RData')
	titles<-"Original data"
	filenames<-"original_allPriors.pdf"
	}
if(dataset==3){
	load('20160113_moreSpreadData_allPriors.RData')
	titles<-"More-spread simulated data"
	filenames<-"moreSpread_allPriors.pdf"
	}
pnames<-c(expression(paste("log(", italic(D), ")")), expression(paste("log(", lambda[c], ")")), expression(paste("log(", italic(L[h]), ")")), expression(paste("log(", italic(L[m]), ")")), expression(paste("log(", alpha, ")")), expression(paste("log(", kappa, ")")), expression(paste("logit(", italic(s[c]), ")")), expression(paste("logit(", italic(s[h]), ")")))

DC<-list(); length(DC)<-4
for(p in 1:4) DC[[p]]<-dctable(X[[p]][[1]])

pdf(file=filenames, width=7, height=7)
par(mfrow=c(4,2), mar=c(3,3,2,1), oma=c(2,2,2,0))
for(i in c(5,6,1:4,7:8)){
	plot(DC[[1]][[i]]$n.clones, DC[[1]][[i]]$mean, "n", ylim=range(DC[[1]][[i]][,c(4,8)], DC[[2]][[i]][,c(4,8)], DC[[3]][[i]][,c(4,8)], DC[[4]][[i]][,c(4,8)]), xaxt="n", las=1, bty="l", xlim=c(8,22))
	for(p in 1:4) lines(DC[[p]][[i]]$n.clones, DC[[p]][[i]]$mean, col=p)
	axis(side=1, at=c(10,15,20))
	for(p in 1:4) plotCI(DC[[p]][[i]]$n.clones, DC[[p]][[i]]$mean, li=DC[[p]][[i]]$'2.5%', ui=DC[[p]][[i]]$'97.5%', gap=0.3, pch=c(19,21)[as.numeric(DC[[p]][[i]]$r.hat>1.1)+1], pt.bg="white", add=TRUE, col=p)
	#mtext(side=3, names(DC[[1]])[i], line=0.5)
	mtext(side=3, pnames[i], line=0.5)
}
mtext(side=1, outer=TRUE, "Number of clones (K)")
mtext(side=2, outer=TRUE, "Posterior mean and 95% CI")
mtext(side=3, outer=TRUE, titles, line=0.5, font=2)
dev.off()
