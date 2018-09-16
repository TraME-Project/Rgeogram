\name{otm3D}

\alias{otm3D}

\title{3D Semi-Discrete Optimal Transport}

\description{
    Solve a three-dimensional semi-discrete optimal transport problem.
}

\usage{
otm3D(chi_mat, weights_in = NULL)
}

\arguments{
    \item{chi_mat}{
    an nx3 matrix of input values.}
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
n_samp <- 1000

chi_kj = matrix(runif(n_samp*3),ncol=3)/2

res <- otm3D(chi_kj)
}
}

\keyword{otm3D}
\concept{otm3D}