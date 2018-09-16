
# tests for OTM functions

rm(list=ls())

library(Rgeogram)

set.seed(12345)

# 2D

chi1 <- c(0.2, 0.7, .4, .5)
chi2 <- c(0.1, 0.3, .9, .4)

chi_kj = cbind(chi1,chi2)/2

res <- otm2D(chi_kj)
otm2D(chi_kj,cbind(rep(0.25,4)))

vtilde <- res$weights

v <- - vtilde + (chi_kj[,1]^2 + chi_kj[,2]^2)
delta_j <- v[4] - v[1:3] 

n_samp <- 1000

chi_kj = matrix(runif(n_samp*2),ncol=2)/2

res <- otm2D(chi_kj)

#
# 3D

# need dim + 1 sample points

n_samp <- 1000

chi_kj = matrix(runif(n_samp*3),ncol=3)/2

res <- otm3D(chi_kj)
res$elapsed_time

#
