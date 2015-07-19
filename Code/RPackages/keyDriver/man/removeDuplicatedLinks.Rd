\name{removeDuplicatedLinks}
\alias{removeDuplicatedLinks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Something
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
removeDuplicatedLinks(linkpairs, directed = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{linkpairs}{
%%     ~~Describe \code{linkpairs} here~~
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
function(linkpairs, directed=F){

    if(dim(linkpairs)[1]==1){
       return (linkpairs)
    }

    links = paste(linkpairs[,1], linkpairs[,2],sep="\t")

    # 1. remove duplications 
    #
    cleanedlinkMatrix = union(links, NULL)
    length(cleanedlinkMatrix)

    #ofname ="tmp.txt"
    #write.table( cleanedlinkMatrix, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)

    # 2. remove inversed duplications
    #
    #linkMatrix <- read.delim(ofname, sep="\t", header=F) # not good for batch operation, ie, many of the same jobs running
    linkMatrix  <- getAllParts(cleanedlinkMatrix, "\t")
    dim(linkMatrix)
    linkMatrix  = as.matrix(linkMatrix)

    if(directed){
       return(linkMatrix)
    }

    if (dim(linkMatrix)[1]==1) {return(linkMatrix);}

    #  first, remove self-interactions are also removed
    #
    selSelfLinks      = linkMatrix[,1]==linkMatrix[,2]
    linkMatrix        = linkMatrix[!selSelfLinks,]
    cleanedlinkMatrix = cleanedlinkMatrix[!selSelfLinks]

    reversedLinks = paste(linkMatrix[,2], linkMatrix[,1],sep="\t")

    no.links = length(reversedLinks)
    reversedLinksWithIndex = cbind(reversedLinks, c(1:no.links) )

    merged = merge(cleanedlinkMatrix, reversedLinksWithIndex, by.x=1,by.y=1,all=F)
    dim(merged)

    if (dim(merged)[1]>0) {
        merged      = as.matrix(merged)
        removedCols = as.integer(merged[,2])

        # cosntruct non-duplicated interactions
        #
        dupLinks    = cleanedlinkMatrix[removedCols]
        dupLinksRev = reversedLinksWithIndex[removedCols]
        uniques     = NULL
        for(i in c(1:length(dupLinks)) )
        {
           found    = is.element(dupLinks[i], uniques)
           foundRev = is.element(dupLinksRev[i], uniques)
           combined = found | foundRev
           if (!combined){
               uniques=c(uniques,dupLinks[i])
           }
        }
        length(uniques)
        xlinkMatrix  = c(cleanedlinkMatrix[-removedCols], uniques)
        #write.table( cleanedlinkMatrix[-removedCols], ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
        #write.table( uniques, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE,append=T)
    }else{
        #write.table( cleanedlinkMatrix, ofname, sep="\t", quote=FALSE, col.names=F, row.names=FALSE)
        xlinkMatrix  = cleanedlinkMatrix
    }

    #linkMatrix <- read.delim(ofname, sep="\t", header=F)
    #dim(linkMatrix)
    #linkMatrix  = as.matrix(linkMatrix)

    linkMatrix = getAllParts(xlinkMatrix, "\t")

    return(linkMatrix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
