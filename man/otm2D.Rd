\name{otm2D}

\alias{otm2D}

\title{2D Semi-Discrete Optimal Transport}

\description{
    Solve a two-dimensional semi-discrete optimal transport problem.
}

\usage{
otm2D(chi_mat, weights_in = NULL)
}

\arguments{
    \item{chi_mat}{
    an nx2 matrix of input values.}
    \item{weights_in}{
    an n-length vector of weights (e.g., the market shares).}
}

\details{
    See the references below.
}

\value{
The function returns a list of objects:
    \item{weights_out}{
    weights of the power diagram.}
    \item{elapsed_time}{
    running time of algorithm.}
}

\author{Keith O'Hara}

\references{
    L\'{e}vy, B. (2015): ``A numerical algorithm for L2 semi-discrete optimal transport in 3D,'' ESAIM: Mathematical Modelling and Numerical Analysis, 49(6), 1693--1715.
}

\examples{
\dontrun{
chi1 <- c(0.2, 0.7, .4, .5)
chi2 <- c(0.1, 0.3, .9, .4)

chi_kj = cbind(chi1,chi2)/2

res <- otm2D(chi_kj)
otm2D(chi_kj,cbind(rep(0.25,4)))

vtilde <- res$weights

v <- - vtilde + (chi_kj[,1]^2 + chi_kj[,2]^2)
delta_j <- v[4] - v[1:3] 
}
}

\keyword{otm2D}
\concept{otm2D}