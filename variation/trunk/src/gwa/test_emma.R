#!/usr/bin/Rscript

#2008-10-08 code to test emma, my beta estimate and variance component. (check emma/R/emma.R)

load('./emma/data/emmadat.rda')
str(emmadat)
source('./emma/R/emma.R')
objects()

#debug(emma.ML.LRT)

#rs <- emma.ML.LRT(emmadat$ys[1,], emmadat$xs,emmadat$K)
rs <- emma.REML.t(emmadat$ys[1,], emmadat$xs,emmadat$K)
#traceback()

#warnings()

#str(rs)
cat("rs$ML1s: ", rs$ML1s, "\n")
cat("rs$ML0s: ", rs$ML0s, "\n")


output_variance = function(marker_index)
{
cat("pvalue=", rs$ps[marker_index], " stats=", rs$stats[marker_index],  " vgs=", rs$vgs[marker_index], " ves=", rs$ves[marker_index], " beta0_est=", rs$beta0_est[marker_index], " beta1_est=", rs$beta1_est[marker_index], " beta0_est1=", rs$beta0_est1[marker_index], " beta1_est1=", rs$beta1_est1[marker_index], "\n")
no_of_individuals = length(emmadat$ys[1,])

phenotype_var = var(emmadat$ys[1,])
cat("variance of phenotype is:", phenotype_var, "\n")

X0 <- matrix(1,no_of_individuals,1)
beta = matrix(c( rs$beta0_est[marker_index], rs$beta1_est[marker_index]), 2,1)
x_beta_est = cbind(X0, emmadat$xs[marker_index,]) %*% beta
#cat("x_beta_est: ", x_beta_est, "\n")
x_beta_var = var(x_beta_est)
cat("variance of X*beta_est at min pvalue:", x_beta_var, " percentage=", x_beta_var/phenotype_var, "=", rs$genotype_var_perc[marker_index], "\n")

x_beta_est_alt = emmadat$xs[marker_index,]*rs$beta1_est[marker_index]
cat("variance of X*beta_est (no X0) at min pvalue:", var(x_beta_est_alt), "\n")



mu = emmadat$K %*% (emmadat$K + diag(rs$ves[marker_index], no_of_individuals)) %*% (emmadat$ys[1,] - x_beta_est)
#cat("mu: ", mu, "\n")
cat("variance of mu: ", var(mu), "\n")
cat("covariance of (x_beta_est, mu): ", cov(x_beta_est, mu), "\n")

residual = emmadat$ys[1,] - x_beta_est - mu
#cat("residual: ", residual, "\n")
cat("variance of residual: ", var(residual), "\n")
cat("covariance of (x_beta_est, residual): ", cov(x_beta_est,residual), "\n")
cat("covariance of (mu, residual): ", cov(mu, residual), "\n")

cat("variance of x_beta_est + mu: ", var(x_beta_est + mu), "\n")
cat("variance of x_beta_est + mu+ residual: ", var(x_beta_est + mu+ residual), "\n")
cat("\n")
}

i = which.max(rs$vgs)
vgs_max_index = i
cat("max vgs at position:", i, "\n")
output_variance(i)


i = which.max(rs$ves)
ves_max_index = i
cat("max ves at position:", i, "\n")
output_variance(i)

i = which.min(rs$ps)
cat("min pvalue at position:", i, "\n")
output_variance(i)

