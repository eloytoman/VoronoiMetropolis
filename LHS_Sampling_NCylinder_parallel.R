library(lhs)
source("~/ProjectVoronoiMetropolis/VoronoiMetropolis/N-cylinder/Parallel_N_cylinder.R")

A <- randomLHS(100,3)

#First column: contractile (gamma)
#Second column: tension (lambda)
#Third column: contractile parameter (s0)
#Fourth column: bending parameter alpha

B <- matrix(nrow = nrow(A), ncol = ncol(A))
B[,1] <- qunif(A[,1], min = 0, max = 0.1)
B[,2] <- qunif(A[,2], min = -0.8, 0.3)
B[,3] <- qlnorm(A[,3], meanlog = 0, sdlog = 1.9)


cl <- makeCluster(3)
registerDoParallel(cl)

results <- foreach(i=1:3,.combine = rbind, .packages = "deldir") %dopar% {
  do.call(metropolisad, list(seed = 100, steps = 1, L=10, Ratio = 10,
                             gamma_ad = B[i,1], lambda_ad = B[i,2], s0 = B[i,3] ))
}

stopCluster(cl)

stopImplicitCluster()