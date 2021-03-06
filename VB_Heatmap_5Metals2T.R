######################
# VB Heatmap
######################

source("VB_Postmean2T.R")

doginv	=FALSE
ngrid 	= 30
j 		= 2 
min.plot.dist = 0.5 ## only color level plot for points within this distance from an observed data point
ylim = c(-5, 5)

##################################
# Time 1: 
##################################

Z.Time = Z.1

qs <- c(0.1, 0.5, 0.9, 0.05, 0.95)
z1.1 <- seq(quantile(Z.Time[,1], qs[4]), quantile(Z.Time[,1], qs[5]), length=ngrid)
z1.2 <- seq(quantile(Z.Time[,2], qs[4]), quantile(Z.Time[,2], qs[5]), length=ngrid)

cross.sec.1 = cbind(expand.grid(z1=z1.1, z2 = z1.2), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]))

grid = expand.grid(z1=z1.1, z2 = z1.2)
mindists = rep(NA, nrow(grid))
for(i in seq_along(mindists)) {
		pt <- as.numeric(grid[i,])
		dists <- as.matrix(dist(rbind(pt, Z.Time[,1:2])))["pt",-1]
		mindists[i] <- min(dists)
}
rm(grid, pt, dists)

hgrid.T.1 <- newh.postmean.gfl.1(Time = 1:5, g = 1, numtime = 2, Znew=cross.sec.1, sel=sel)$postmean
hgrid.T.1[mindists > min.plot.dist] <- NA
hgrid.T.1 = matrix(hgrid.T.1, nrow=ngrid)

if (Scenario1){
  true.1 = matrix(rep(0, ngrid*ngrid), nrow=ngrid)
}
if (Scenario2) {
  true.1 = matrix(0.5*(cross.sec.1[,1]^2 + cross.sec.1[,2]^2 - cross.sec.1[,4]^2 - cross.sec.1[,5]^2 + 0.5*cross.sec.1[,4]*cross.sec.1[,5] + cross.sec.1[,4] + cross.sec.1[,5]), nrow=ngrid)
}

##################################
# Time 2: 
##################################

Z.Time = Z.2

z2.1 <- seq(quantile(Z.Time[,1], qs[4]), quantile(Z.Time[,1], qs[5]), length=ngrid)
z2.2 <- seq(quantile(Z.Time[,2], qs[4]), quantile(Z.Time[,2], qs[5]), length=ngrid)

cross.sec.2 = cbind(expand.grid(z1=z1.1, z2 = z1.2), z3=median(Z.Time[,3]), z4=median(Z.Time[,4]), z5=median(Z.Time[,5]))

grid = expand.grid(z1=z2.1, z2 = z2.2)
mindists = rep(NA, nrow(grid))
for(i in seq_along(mindists)) {
		pt <- as.numeric(grid[i,])
		dists <- as.matrix(dist(rbind(pt, Z.Time[,1:2])))["pt",-1]
		mindists[i] <- min(dists)
}

hgrid.T.2 <- newh.postmean.gfl.2(Time = 6:10, g = 2, numtime = 2, Znew=cross.sec.2, sel=sel)$postmean
hgrid.T.2[mindists > min.plot.dist] <- NA
hgrid.T.2 = matrix(hgrid.T.2, nrow=ngrid)

if (Scenario1){
  true.2 = matrix(0.5*(cross.sec.1[,1]^2 - cross.sec.1[,2]^2 + 1/2*cross.sec.1[,1]*cross.sec.1[,2] + cross.sec.1[,1] + cross.sec.1[,2]), nrow=ngrid)
}
if (Scenario2) {
  true.2 = matrix(0.5*(cross.sec.1[,1]^2 - cross.sec.1[,2]^2 - cross.sec.1[,1]*cross.sec.1[,2] + cross.sec.1[,3]^2 + cross.sec.1[,4] - cross.sec.1[,5]), nrow=ngrid)
}

###############################

# Multiple images with a common legend

set.panel()
zlim = c(-2.5,2.5)

par(oma=c( 0,0,0,4)) # margin of 4 spaces width at right hand side
set.panel( 2,2) 

# now draw all your plots using usual image command

image(z1.1, z1.2, hgrid.T.1, xlab="Z1", ylab="Z2", col=tim.colors(), main="Time 1: VB-LKMR Estimated", zlim=zlim)
points(Z[,1:5])
image(z2.1, z2.2, hgrid.T.2, xlab="Z1", ylab="Z2", col=tim.colors(), main="Time 2: VB-LKMR Estimated", zlim=zlim)
points(Z[,6:10])
image(z1.1, z1.2, true.1, xlab="Z1", ylab="Z2", col=tim.colors(), main="Time 1: True", zlim=zlim)
points(Z[,1:5])
image(z2.1, z2.2, true.2, xlab="Z1", ylab="Z2", col=tim.colors(), main="Time 2: True", zlim=zlim)
points(Z[,6:10])


par(oma=c( 0,0,0,1)) # reset margin to be much smaller.
image.plot( legend.only=TRUE, zlim=zlim) 

# image.plot tricked into  plotting in margin of old setting 

set.panel() # reset plotting device
#dev.off()

###############################


