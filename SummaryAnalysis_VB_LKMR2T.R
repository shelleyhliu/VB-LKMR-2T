#############################
# Shelley H. Liu
# Analysis of output
# Updated 3-25-2019
#############################

Summary.MFVB 	= matrix(NA, nrow=T, ncol=4)

# Time point 1
ind 			= 1:N
h_hat.1			= E.q.h[ind]
if (Scenario1){
  h_true.1		= rep(0, N)
}
if (Scenario2){
  h_true.1 = 0.5*(Z.1[,1]^2 + Z.1[,2]^2 - Z.1[,4]^2 - Z.1[,5]^2 + 0.5*Z.1[,4]*Z.1[,5] + Z.1[,4] + Z.1[,5])
}

model.1 = lm(h_hat.1 ~ h_true.1)

if (Scenario1){
  Summary.MFVB[1,] = c(summary(model.1)$coef[1,1], 0, 0, sqrt( sum( (h_hat.1 - h_true.1)^2) / N))
}
if (Scenario2){
  Summary.MFVB[1,] = c(summary(model.1)$coef[1,1], summary(model.1)$coef[2,1], summary(model.1)$r.sq, sqrt( sum( (h_hat.1 - h_true.1)^2) / N))
}

# Time point 2
ind 			= (N+1):(2*N)
h_hat.2 		= E.q.h[ind]
if (Scenario1){
  h_true.2 		= 0.5*(Z.2[,1]^2 - Z.2[,2]^2 + 1/2*Z.2[,1]*Z.2[,2] + Z.2[,1] + Z.2[,2])
}
if (Scenario2){
  h_true.2 = 0.5*(Z.2[,1]^2 - Z.2[,2]^2 - Z.2[,1]*Z.2[,2] + Z.2[,3]^2 + Z.2[,4] - Z.2[,5])
}
model.2 = lm(h_hat.2 ~ h_true.2)

Summary.MFVB[2,] = c(summary(model.2)$coef[1,1], summary(model.2)$coef[2,1], summary(model.2)$r.sq, sqrt( sum( (h_hat.2 - h_true.2)^2) / N))

print("VB: Intercept, Slope, R^2, RMSE")
print( round(Summary.MFVB, 2) )


