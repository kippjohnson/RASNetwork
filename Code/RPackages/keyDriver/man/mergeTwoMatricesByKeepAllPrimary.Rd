\name{mergeTwoMatricesByKeepAllPrimary}
\alias{mergeTwoMatricesByKeepAllPrimary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Something
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mergeTwoMatricesByKeepAllPrimary(primaryMatrix, minorMatrix, missinglabel = "", keepAllPrimary = T, keepPrimaryOrder = T, keepAll = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{primaryMatrix}{
%%     ~~Describe \code{primaryMatrix} here~~
}
  \item{minorMatrix}{
%%     ~~Describe \code{minorMatrix} here~~
}
  \item{missinglabel}{
%%     ~~Describe \code{missinglabel} here~~
}
  \item{keepAllPrimary}{
%%     ~~Describe \code{keepAllPrimary} here~~
}
  \item{keepPrimaryOrder}{
%%     ~~Describe \code{keepPrimaryOrder} here~~
}
  \item{keepAll}{
%%     ~~Describe \code{keepAll} here~~
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
function(primaryMatrix, minorMatrix, missinglabel="", keepAllPrimary=T, keepPrimaryOrder=T, keepAll=F)
{
  no.promarycols <- dim(primaryMatrix)[2]
  no.mustbegenes <- dim(primaryMatrix)[1]

  # we add in one more column to indicate which genes are mustbeincluded after being merged with mcg
  keyword="mustbeused"
  mustbeGenesMatrix = cbind(primaryMatrix, c(1:no.mustbegenes), rep(keyword, no.mustbegenes) )

  if (is.null(colnames(primaryMatrix)) ){
    colnames(mustbeGenesMatrix) <- c( c(1:no.promarycols), "primorder", keyword)
  }else{
    colnames(mustbeGenesMatrix) <- c( colnames(primaryMatrix), "primorder", keyword)
  }
  dim(mustbeGenesMatrix)

  if(is.null(keepAllPrimary) ){ #normal merge: to have the common elements
    myMatrix = merge(mustbeGenesMatrix, minorMatrix, by.x=1, by.y=1,all.x=F,sort=F,all=F)
  }else{
    myMatrix = merge(mustbeGenesMatrix, minorMatrix, by.x=1, by.y=1,all.x=T,sort=F,all=T)
  }

  dim(myMatrix)
  nocols.mymatrix <- dim(myMatrix)[2]

  #the mustbeused genes which are not included in minor have NAs in the column $mustbeused
  #so we can use this information to figure out which mustbeused genes missing in minorMatrix
  myMatrix[,nocols.mymatrix] = ifelse( is.na(myMatrix[,nocols.mymatrix]), missinglabel, as.character(myMatrix[,nocols.mymatrix]) )

  orders = order( as.numeric(as.matrix(myMatrix[, no.promarycols+1])))
  if (keepPrimaryOrder)
      myMatrix = myMatrix[orders,]

  if (is.null(keepAllPrimary) ){
     selected = rep(T, dim(myMatrix)[1])
  }else{
     if (keepAllPrimary)
       selected = !(is.na(myMatrix[, no.promarycols+2]))
     else #return the row-elements in minor which are missed in primary
       selected = is.na(myMatrix[, no.promarycols+2])
  }

  sum(selected)

  #keep the primary matrix and remove the mustbeused column
  myMatrix[selected, -c(no.promarycols+1, no.promarycols+2)]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
