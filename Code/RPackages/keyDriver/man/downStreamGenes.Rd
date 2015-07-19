\name{downStreamGenes}
\alias{downStreamGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{something
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
downStreamGenes(netpairs, seednodes, N = 100, directed = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{netpairs}{
%%     ~~Describe \code{netpairs} here~~
}
  \item{seednodes}{
%%     ~~Describe \code{seednodes} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{directed}{
%%     ~~Describe \code{directed} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(netpairs, seednodes, N= 100, directed=T)
{
   prenodes = seednodes
   cnt = N
   while(T) {
      retlinks = findNLayerNeighbors_LinkPairs(linkpairs=netpairs, subnetNodes=prenodes, 
                      nlayers=1, directed=directed)

      if(is.null(retlinks)){return (NULL); }

      curnodes = union(retlinks[,1],retlinks[,2]) 
      pcdiff   = setdiff(curnodes, prenodes)
      prenodes = curnodes

      if(length(pcdiff)==0){break}

      cnt= cnt-1
      if (cnt==0) {break}
   }

   if(is.null(retlinks)){return (NULL)}
   else {
      return(curnodes)
   }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
