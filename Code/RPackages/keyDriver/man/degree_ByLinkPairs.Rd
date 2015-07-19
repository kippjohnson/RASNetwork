\name{degreeByLinkPairs}
\alias{degreeByLinkPairs}
\title{
Find links of pairs by degree of relationship
}
\description{
Search through a set of linked pairs for additional paired relationships.  Add
until you cannot find any additional pairs.
}
\usage{
degreeByLinkPairs(linkpairs , directed = FALSE , cleangarbage = FALSE )
}
\arguments{
  \item{linkpairs}{
    Vector of linked pairs.
  }
  \item{directed}{
    Is this a directed or undirected graph.
  }
  \item{cleangarbage}{
    Collect garbage each round.
  }
}
\details{
Loop through a linked list of paired entities to find additional relationships
not listed.
}
\value{
Vector of linked entities.
}
\author{
Bin Zhang
}
\examples{
## Needs to be added.
}
\keyword{ networks }
