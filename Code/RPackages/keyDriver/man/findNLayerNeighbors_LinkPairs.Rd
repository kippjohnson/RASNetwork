\name{findNLayerNeighborsLinkPairs}
\alias{findNLayerNeighborsLinkPairs}
\title{
Find new neighbors N degrees of separation from original linked list.
}
\description{
Finds new neighbors separated by N degrees.
}
\usage{
findNLayerNeighborsLinkPairs( linkpairs , subnetNodes , nlayers = 1 , directed = FALSE )
}
\arguments{
  \item{linkpairs}{
List of linked paired entities.
}
  \item{subnetNodes}{
Set of nodes of interest.  Could be a signature.
}
  \item{nlayers}{
Degrees of separation to be investigated.
}
  \item{directed}{
Is the graph directed or undirected.
}
}
\details{
Finds new list of paired entities of neighbors separated by N degrees of separation. 
}
\value{
%% To be added
}
\author{
Bin Zhang
}
\examples{
## To be added
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ networks }
