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

sim<-function(params,x){

# Parameters
kappa<-exp(as.numeric(params['kappa.log']))
sc<- exp(as.numeric(params['sc.logit']))/(1+exp(as.numeric(params['sc.logit'])))
sh<- exp(as.numeric(params['sh.logit']))/(1+exp(as.numeric(params['sh.logit'])))
Lc<-exp(as.numeric(params['Lc.log']))
alpha <- exp(as.numeric(params['alpha.log']))
D<-exp(as.numeric(params['D.log']))
Lh_R<-exp(as.numeric(params['Lh_R.log']))
Lm_R<-exp(as.numeric(params['Lm_R.log']))

	
# Define parameters for larval distribution
a1<-(gamma+(gamma^2+4*u_n*D)^0.5)/(2*D); a2<-abs((gamma-(gamma^2+4*u_n*D)^0.5)/(2*D));
b1<-(gamma+(gamma^2+4*u_c*D)^0.5)/(2*D); b2<-abs((gamma-(gamma^2+4*u_c*D)^0.5)/(2*D));

# Larval constant (integral from -infty to infty)
c1<-1/(a1*(a1+b2))+1/(b1*(a1-b1))-1/(a1*(a1-b1))+1/(b1*(a2+b1))+1/(b2*(a1+b2))+1/(b2*(a2-b2))-1/(a2*(a2-b2))+1/(a2*(a2+b1))

#----------------------------------------------------------------------------------------
## Farm Footprints
# Create indices of which solution to use (left hand, right hand, or middle)
# The solution is in three parts, and we calculate all three and then the index tells us which of the parts we keep.
index<-matrix();length(index)<-3*length(x)*3;dim(index)<-c(3,length(x), 3)
for(i in 1:length(x)){ # for each site		
		# Copepidite
		index[1,i,1]<-(x[i]<=y)
		index[1,i,2]<-((x[i]-Lc) <= y) && (x[i]>y)
		index[1,i,3]<-((x[i]-Lc) > y)
		# Chalimus
		index[2,i,1]<-((x[i]-Lc) <= y)
		index[2,i,2]<-((x[i]-Lc*(1+Lh_R)) <= y && (x[i]-Lc) > y)
		index[2,i,3]<-((x[i]-Lc*(1+Lh_R)) > y)
		# Motile
		index[3,i,1]<-((x[i]-Lc*(1+Lh_R)) <= y)
		index[3,i,2]<-((x[i]-(Lc*(1+Lh_R+Lm_R))) <= y && (x[i]-Lc*(1+Lh_R)) > y)
		index[3,i,3]<-((x[i]-(Lc*(1+Lh_R+Lm_R))) > y)
		}# end i sites
		
#----------------------------------------------------------------------------------------
# Solve all three solutions for each x
all.pred<-matrix();length(all.pred)<-3*length(x)*3;dim(all.pred)<-c(3, length(x), 3)

for(i in 1:length(x)){ #for each site	
		#--------------------------------------
		#Copepodite
		#--------------------------------------
		# Left-hand solution
		all.pred[1,i,1]<-exp(b1*(x[i]-y))*(1-exp(-b1*Lc))*(a2+a1)/(b1*(a1-b1)*(a2+b1)) -exp(a1*(x[i]-y))*(1-exp(-a1*Lc))*((b1+b2)/(a1*(a1+b2)*(a1-b1)))

		# Middle solution
		all.pred[1,i,2]<--(1-exp(a1*(x[i]-Lc-y)))*(((b1+b2))/(a1*(a1+b2)*(a1-b1)))+(1-exp(b1*(x[i]-Lc-y)))*((a2+a1)/(b1*(a1-b1)*(a2+b1)))+(1-exp(-b2*(x[i]-y)))*((a2+a1)/(b2*(a1+b2)*(a2-b2)))-(1-exp(-a2*(x[i]-y)))* ((b2+b1)/(a2*(a2+b1)*(a2-b2)))
		
		# Right hand solution
	all.pred[1,i,3]<-exp(-a2*(x[i]-y))*(1-exp(a2*Lc))*((b2+b1)/(a2*(a2-b2)*(a2+b1)))-exp(-b2*(x[i]-y))*(1-exp(b2*Lc))*((a2+a1)/(b2*(a1+b2)*(a2-b2)))
		
		#--------------------------------------
		#Chalimus
		#--------------------------------------
		# Left-hand solution
		all.pred[2,i,1]<-exp(b1*(x[i]-y))*(exp(-b1*Lc)-exp(-b1*Lc*(1+Lh_R)))*(a2+a1)/(b1*(a1-b1)*(a2+b1)) -exp(a1*(x[i]-y))*(exp(-a1*Lc)-exp(-a1*Lc*(1+Lh_R)))*((b1+b2)/(a1*(a1+b2)*(a1-b1)))

		# Middle solution
		all.pred[2,i,2]<--(1-exp(a1*(x[i]-Lc*(1+Lh_R)-y)))*(((b1+b2))/(a1*(a1+b2)*(a1-b1)))+(1-exp(b1*(x[i]-Lc*(1+Lh_R)-y)))*((a2+a1)/(b1*(a1-b1)*(a2+b1)))+(1-exp(-b2*(x[i]-Lc-y)))*((a2+a1)/(b2*(a1+b2)*(a2-b2)))-(1-exp(-a2*(x[i]-Lc-y)))* ((b2+b1)/(a2*(a2+b1)*(a2-b2)))
		
		# Right hand solution
	all.pred[2,i,3]<-exp(-a2*(x[i]-y))*(exp(a2*Lc)-exp(a2*Lc*(1+Lh_R)))*((b2+b1)/(a2*(a2-b2)*(a2+b1)))-exp(-b2*(x[i]-y))*(exp(b2*Lc)-exp(b2*Lc*(1+Lh_R)))*((a2+a1)/(b2*(a1+b2)*(a2-b2)))
		
		#--------------------------------------
		#Motile
		#--------------------------------------
		# Left-hand solution
		all.pred[3,i,1]<-exp(b1*(x[i]-y))*(exp(-b1*Lc*(1+Lh_R))-exp(-b1*Lc*(1+Lh_R+Lm_R)))*(a2+a1)/(b1*(a1-b1)*(a2+b1))-exp(a1*(x[i]-y))*(exp(-a1*Lc*(1+Lh_R))-exp(-a1*Lc*(1+Lh_R+Lm_R)))*((b1+b2)/(a1*(a1+b2)*(a1-b1)))

		# Middle solution
		all.pred[3,i,2]<--(1-exp(a1*(x[i]-Lc*(1+Lh_R+Lm_R)-y)))*(((b1+b2))/(a1*(a1+b2)*(a1-b1)))+(1-exp(b1*(x[i]-Lc*(1+Lh_R+Lm_R)-y)))*((a2+a1)/(b1*(a1-b1)*(a2+b1)))+(1-exp(-b2*(x[i]-Lc*(1+Lh_R)-y)))*((a2+a1)/(b2*(a1+b2)*(a2-b2)))-(1-exp(-a2*(x[i]-Lc*(1+Lh_R)-y)))*((b2+b1)/(a2*(a2+b1)*(a2-b2)))
		
		# Right hand solution
	all.pred[3,i,3]<-exp(-a2*(x[i]-y))*(exp(a2*Lc*(1+Lh_R))-exp(a2*Lc*(1+Lh_R+Lm_R)))*((b2+b1)/(a2*(a2-b2)*(a2+b1)))-exp(-b2*(x[i]-y))*(exp(b2*Lc*(1+Lh_R))-exp(b2*Lc*(1+Lh_R+Lm_R)))*((a2+a1)/(b2*(a1+b2)*(a2-b2)))

} #end site i

#----------------------------------------------------------------------------------------
# Multiple by the index and sum all three solutions together
# should only be one solution that is non-zero for each location
int.solns<-matrix();length(int.solns)<-3*length(x); dim(int.solns)<-c(3,length(x))
for(m in 1:3){ # for each stage
	for(i in 1:length(x)){ #for each site
		int.solns[m,i]<-all.pred[m,i,1]*index[m,i,1] + all.pred[m,i,2]*index[m,i,2] + all.pred[m,i,3]*index[m,i,3]
			}}
		
#----------------------------------------------------------------------------------------
# Add solutions from all three farms, scaled by their input
pred<-matrix(nrow=length(x), ncol=3)
for(i in 1:length(x)){ # for each site
	pred[i,1]<-kappa*Lc + alpha/c1*int.solns[1,i]
	pred[i,2]<-sc*(kappa*Lc*Lh_R + alpha/c1*int.solns[2,i])
	pred[i,3]<-sh*sc*(kappa*Lc*Lm_R + alpha/c1*int.solns[3,i])
	}

#-----------------------------------------------------------------------------------------------
return(pred)
#-----------------------------------------------------------------------------------------------	
	}