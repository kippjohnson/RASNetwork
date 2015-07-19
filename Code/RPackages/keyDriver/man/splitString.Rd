\name{splitString}
\alias{splitString}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Something
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
splitString(mystring, separator = "; ")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mystring}{
%%     ~~Describe \code{mystring} here~~
}
  \item{separator}{
%%     ~~Describe \code{separator} here~~
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
function(mystring, separator="; "){
  splitted = NULL
  for (each in mystring){
     if (is.na(each) | is.null(each)){
        next
     }
     a=unlist( strsplit(each, separator) )
     splitted =c(splitted, a)
  }
  #a=unlist( strsplit(mystring, separator) )
  return(splitted )
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
