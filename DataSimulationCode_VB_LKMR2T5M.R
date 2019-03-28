#############################
# Shelley H. Liu 
# Data Simulation File for VB-LKMR
# Simulating gradual quadratic effects model
#############################

set.seed(1)
library(lars); library(mvtnorm); library(SuppDists); library(MCMCpack); library(grplasso); library(magic); library(kernlab); library(MASS); library(fields)

N 				= 200 # Number of subjects
n         = N
T 				= 2 # Number of time points
P 				= N*T
M				  = 5 # Number of metals
grpsize   = M
numtime   = T

TauG = N #Group size for tau. For the kernel case it should be N.

poly = polydot(degree = 2, offset = 1)

offDiagonal = function(x) {
		diagBelow = diag(x)
		i = 1
		while (i <= TauG) { 
			diagBelow=rbind(rep(0,length(x)	+i),cbind(diagBelow,rep(0,length(x) + i - 1)))
			i = i + 1
		}
		mat <- diagBelow + t(diagBelow) - diag(diag(diagBelow))
		return(mat)
}

#############################

cofactor		= 1e-7

ar 			= 0.8

time_mat <- matrix( ar, 2,2 ) 
diag(time_mat) <- 1

metal_mat <- matrix( c(1, 0.19, 0.32, 0.15, 0.4, 0.19, 1, 0.11, 0.25, 0.4, 0.32, 0.11, 1, 0.17, 0.24, 0.15, 0.25, 0.17, 1, 0.35, 0.4, 0.4, 0.24, 0.35, 1), 5,5, byrow=TRUE) 

cmat = kronecker(time_mat, metal_mat)
obs <- matrix( mvrnorm(N, mu=rep(0, 5*2), Sigma=cmat), nrow = N)

Z = obs
Z 		= apply(Z, 2, scale)

Z.1 = Z[,1:5]
Z.2 = Z[,6:10]

X				= cbind(scale(matrix(rnorm(N, mean = 10, sd = 1))), scale(matrix(sample(1:2, N, replace = TRUE))))
beta			= c(1,1)

res.sd = 1 
if (Scenario1){
  Y = matrix(rnorm(n=N, mean=0, sd=res.sd)) + X %*% beta + 0.5*(Z.2[,1]^2 - Z.2[,2]^2 + 1/2*Z.2[,1]*Z.2[,2] + Z.2[,1] + Z.2[,2]) 
}
if (Scenario2){
  Y = matrix(rnorm(n=N, mean=0, sd=res.sd)) + X %*% beta + 0.5*(Z.1[,1]^2 + Z.1[,2]^2 - Z.1[,4]^2 - Z.1[,5]^2 + 0.5*Z.1[,4]*Z.1[,5] + Z.1[,4] + Z.1[,5]) + 0.5*(Z.2[,1]^2 - Z.2[,2]^2 - Z.2[,1]*Z[,2] + Z.2[,3]^2 + Z.2[,4] - Z.2[,5])
}
# Make W matrix

W 				= diag(N)
counter 		= 1
while(counter < T) {
	W 		= cbind(W, diag(N))
	counter = counter+1
}
