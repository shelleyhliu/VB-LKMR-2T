##############
# Relative importance of each metal at two time windows
##############

source("VB_Postmean2T.R")

doginv = FALSE
par(mfrow=c(1,2))

ylim = c(-2, 2)
mat.time = matrix(NA, nrow=numtime, ncol=grpsize)
mat.sd = matrix(NA, nrow=numtime, ncol=grpsize)
qs = c(0.25, 0.75, 0.5)

true.mat	= matrix(rep(0,numtime*grpsize), nrow=numtime, ncol=grpsize)

##############
# Time 1
##############

Z.Time = Z.1
true.mat[1,] = rep(0, grpsize)
for (j in 1:grpsize) {
	cross.sec = rbind(apply(Z.Time, 2, median), apply(Z.Time, 2, median))
	cross.sec[,j] = c(quantile(Z.Time[,j], qs[2]), quantile(Z.Time[,j], qs[1]))
	
	hgrid.cross.sec <- newh.postmean.gfl.1(Time = 1:5, g = 1, numtime = numtime, Znew=cross.sec)
	mat.time[1,j] = hgrid.cross.sec$postmean[1] - hgrid.cross.sec$postmean[2]
	
	mat.sd[1,j] = sqrt(hgrid.cross.sec$postvar[2,2] + hgrid.cross.sec$postvar[1,1] - 2*hgrid.cross.sec$postvar[1,2])

	if (j==1 && Scenario2 == TRUE) {
	  true.mat[1,j] = 0.5*(cross.sec[1,j]^2 + cross.sec[1,j+1]^2 - cross.sec[1, j+3]^2 - cross.sec[1,j+4]^2 + 0.5*cross.sec[1,j+3]*cross.sec[1,j+4] + cross.sec[1,j+3] + cross.sec[1,j+4]) - 0.5*(cross.sec[2,j]^2 + cross.sec[2,j+1]^2 - cross.sec[2, j+3]^2 - cross.sec[2,j+4]^2 + 0.5*cross.sec[2,j+3]*cross.sec[2,j+4] + cross.sec[2,j+3] + cross.sec[2,j+4])
	  
	}
	
	if (j==2 && Scenario2 == TRUE) {
	  j = j - 1
	  true.mat[1,j+1] = 0.5*(cross.sec[1,j]^2 + cross.sec[1,j+1]^2 - cross.sec[1, j+3]^2 - cross.sec[1,j+4]^2 + 0.5*cross.sec[1,j+3]*cross.sec[1,j+4] + cross.sec[1,j+3] + cross.sec[1,j+4]) - 0.5*(cross.sec[2,j]^2 + cross.sec[2,j+1]^2 - cross.sec[2, j+3]^2 - cross.sec[2,j+4]^2 + 0.5*cross.sec[2,j+3]*cross.sec[2,j+4] + cross.sec[2,j+3] + cross.sec[2,j+4])
	}
	
	if (j==3 && Scenario2 == TRUE) {
	  j = j - 2
	  true.mat[1,j+2] = 0.5*(cross.sec[1,j]^2 + cross.sec[1,j+1]^2 - cross.sec[1, j+3]^2 - cross.sec[1,j+4]^2 + 0.5*cross.sec[1,j+3]*cross.sec[1,j+4] + cross.sec[1,j+3] + cross.sec[1,j+4]) - 0.5*(cross.sec[2,j]^2 + cross.sec[2,j+1]^2 - cross.sec[2, j+3]^2 - cross.sec[2,j+4]^2 + 0.5*cross.sec[2,j+3]*cross.sec[2,j+4] + cross.sec[2,j+3] + cross.sec[2,j+4])
	}
	
	if (j==4 && Scenario2 == TRUE) {
	  j = j - 3
	  true.mat[1,j+3] = 0.5*(cross.sec[1,j]^2 + cross.sec[1,j+1]^2 - cross.sec[1, j+3]^2 - cross.sec[1,j+4]^2 + 0.5*cross.sec[1,j+3]*cross.sec[1,j+4] + cross.sec[1,j+3] + cross.sec[1,j+4]) - 0.5*(cross.sec[2,j]^2 + cross.sec[2,j+1]^2 - cross.sec[2, j+3]^2 - cross.sec[2,j+4]^2 + 0.5*cross.sec[2,j+3]*cross.sec[2,j+4] + cross.sec[2,j+3] + cross.sec[2,j+4])
	}
	
	if (j==5 && Scenario2 == TRUE) {
	  j = j - 4
	  true.mat[1,j+4] = 0.5*(cross.sec[1,j]^2 + cross.sec[1,j+1]^2 - cross.sec[1, j+3]^2 - cross.sec[1,j+4]^2 + 0.5*cross.sec[1,j+3]*cross.sec[1,j+4] + cross.sec[1,j+3] + cross.sec[1,j+4]) - 0.5*(cross.sec[2,j]^2 + cross.sec[2,j+1]^2 - cross.sec[2, j+3]^2 - cross.sec[2,j+4]^2 + 0.5*cross.sec[2,j+3]*cross.sec[2,j+4] + cross.sec[2,j+3] + cross.sec[2,j+4])
	}
	}
	
#Plot
plot(1:5, mat.time[1,], xaxt="n",
    ylim=ylim,
    pch=15, xlab="Metal", ylab="Main effect of each metal",
    main="Time 1"
)

arrows(1:5, mat.time[1,] -1.96*mat.sd[1,], 1:5, mat.time[1,] +1.96*mat.sd[1,], length=0.05, angle=90, code=3)
axis(1, at=1:5, labels=c("Z1", "Z2", "Z3", "Z4", "Z5"))
abline(h=0)
points(1:grpsize, true.mat[1,], col="red", pch=19)



##############
# Time 2
##############

Z.Time = Z.2
true.mat[2,] = rep(0, grpsize)
for (j in 1:grpsize) {
	cross.sec = rbind(apply(Z.Time, 2, median), apply(Z.Time, 2, median))
	cross.sec[,j] = c(quantile(Z.Time[,j], qs[2]), quantile(Z.Time[,j], qs[1]))
	
	hgrid.cross.sec <- newh.postmean.gfl.2(Time = 6:10, g = 2, numtime = numtime, Znew=cross.sec)
	
	mat.time[2,j] = hgrid.cross.sec$postmean[1] - hgrid.cross.sec$postmean[2]
	mat.sd[2,j] = sqrt(hgrid.cross.sec$postvar[2,2] + hgrid.cross.sec$postvar[1,1] - 2*hgrid.cross.sec$postvar[1,2])
	if (Scenario1){
	  if (j==1 ) {
	    true.mat[2,j] = 0.5*(cross.sec[1,j]^2 - cross.sec[1,(j+1)]^2 + 0.5*cross.sec[1,j]*cross.sec[1,(j+1)] + cross.sec[1,j] + cross.sec[1,(j+1)] - (cross.sec[2,j]^2 - cross.sec[2,(j+1)]^2 + 0.5*cross.sec[2,j]*cross.sec[2,(j+1)] + cross.sec[2,j] + cross.sec[2,(j+1)]) )
	  }
	  if (j==2) {
	    true.mat[2,j] = 0.5*(cross.sec[1,j]^2 - cross.sec[1,(j-1)]^2 + 0.5*cross.sec[1,j]*cross.sec[1,(j-1)] + cross.sec[1,j] + cross.sec[1,(j-1)] - (cross.sec[2,j]^2 - cross.sec[2,(j-1)]^2 + 0.5*cross.sec[2,j]*cross.sec[2,(j-1)] + cross.sec[2,j] + cross.sec[2,(j-1)]) )
	  }
	}
	if (j==1 && Scenario2 == TRUE) {

	  true.mat[2,j] = 0.5*(cross.sec[1,j]^2 - cross.sec[1,j+1]^2 - cross.sec[1,j]*cross.sec[1,j+1] + cross.sec[1,j+2]^2 + cross.sec[1,j+3] - cross.sec[1,j+4]) - 0.5*(cross.sec[2,j]^2 - cross.sec[2,j+1]^2 - cross.sec[2,j]*cross.sec[2,j+1] + cross.sec[2,j+2]^2 + cross.sec[2,j+3] - cross.sec[2,j+4])
	}
	
	if (j==2 && Scenario2 == TRUE) {
	  j = j - 1
	  true.mat[2,j+1] = 0.5*(cross.sec[1,j]^2 - cross.sec[1,j+1]^2 - cross.sec[1,j]*cross.sec[1,j+1] + cross.sec[1,j+2]^2 + cross.sec[1,j+3] - cross.sec[1,j+4]) - 0.5*(cross.sec[2,j]^2 - cross.sec[2,j+1]^2 - cross.sec[2,j]*cross.sec[2,j+1] + cross.sec[2,j+2]^2 + cross.sec[2,j+3] - cross.sec[2,j+4])
	}
	
	if (j==3 && Scenario2 == TRUE) {
	  j = j - 2
	  true.mat[2,j+2] = 0.5*(cross.sec[1,j]^2 - cross.sec[1,j+1]^2 - cross.sec[1,j]*cross.sec[1,j+1] + cross.sec[1,j+2]^2 + cross.sec[1,j+3] - cross.sec[1,j+4]) - 0.5*(cross.sec[2,j]^2 - cross.sec[2,j+1]^2 - cross.sec[2,j]*cross.sec[2,j+1] + cross.sec[2,j+2]^2 + cross.sec[2,j+3] - cross.sec[2,j+4])
	}
	
	if (j==4 && Scenario2 == TRUE) {
	  j = j - 3
	  true.mat[2,j+3] = 0.5*(cross.sec[1,j]^2 - cross.sec[1,j+1]^2 - cross.sec[1,j]*cross.sec[1,j+1] + cross.sec[1,j+2]^2 + cross.sec[1,j+3] - cross.sec[1,j+4]) - 0.5*(cross.sec[2,j]^2 - cross.sec[2,j+1]^2 - cross.sec[2,j]*cross.sec[2,j+1] + cross.sec[2,j+2]^2 + cross.sec[2,j+3] - cross.sec[2,j+4])
	}
	
	if (j==5 && Scenario2 == TRUE) {
	  j = j - 4
	  true.mat[2,j+4] = 0.5*(cross.sec[1,j]^2 - cross.sec[1,j+1]^2 - cross.sec[1,j]*cross.sec[1,j+1] + cross.sec[1,j+2]^2 + cross.sec[1,j+3] - cross.sec[1,j+4]) - 0.5*(cross.sec[2,j]^2 - cross.sec[2,j+1]^2 - cross.sec[2,j]*cross.sec[2,j+1] + cross.sec[2,j+2]^2 + cross.sec[2,j+3] - cross.sec[2,j+4])
	}
}
	
#Plot
plot(1:5, mat.time[2,], xaxt="n",
    ylim=ylim,
    pch=15, xlab="Metal", ylab="Main effect of each metal",
    main="Time 2"
)

arrows(1:5, mat.time[2,] -1.96*mat.sd[2,], 1:5, mat.time[2,] +1.96*mat.sd[2,], length=0.05, angle=90, code=3)
axis(1, at=1:5, labels=c("Z1", "Z2", "Z3", "Z4", "Z5"))
abline(h=0)

points(1:grpsize, true.mat[2,], col="red", pch=19)




