\name{bellman_ford}

\alias{bellman_ford}

\title{Bellman-Ford Algorithm}

\description{
    Solve the shortest path problem using the Bellman-Ford algorithm.
}

\usage{
bellman_ford(n_nodes,source_ind,arcs_matrix)
}

\arguments{
    \item{n_nodes}{
    The number of nodes in the network}
	\item{source_ind}{
	The source node index.}
	\item{arcs_matrix}{
	A K x 3 matrix of arcs for the network.}
}

\details{
    See \url{https://en.wikipedia.org/wiki/Bellman\%E2\%80\%93Ford_algorithm}.
}

\value{
The function returns a list of objects:

    \item{min_dist}{
	A numeric vector containing the minimum distance from the source node to other nodes in the network.}
	\item{path_list}{
	The minimum distance path for each node.}
    \item{elapsed_time}{
	Runtime of the algorithm.}
}

\author{Keith O'Hara}

\seealso{
  \code{\link{get_shortest_path}}.
}

\examples{
\dontrun{
arcs <- rbind(c(1,2,7),
              c(1,3,9),
              c(1,6,14),
              c(2,3,10),
              c(2,4,15),
              c(3,4,11),
              c(3,6,2),
              c(4,5,6),
              c(6,5,9))

nbNodes = max(max(arcs[,1]),max(arcs[,2]))

# set source node

s_n <- 1

# run solver

sol <- bellman_ford(nbNodes,s_n,arcs)
}
}

\keyword{bellman_ford}
\concept{bellman_ford}