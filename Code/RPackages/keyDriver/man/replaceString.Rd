\name{replaceString}
\alias{replaceString}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Something
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
replaceString(fullfnames, oldstr, newstr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fullfnames}{
%%     ~~Describe \code{fullfnames} here~~
}
  \item{oldstr}{
%%     ~~Describe \code{oldstr} here~~
}
  \item{newstr}{
%%     ~~Describe \code{newstr} here~~
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
function(fullfnames, oldstr, newstr){

  no.files = length(fullfnames)
  res = NULL
  for(each in fullfnames) {
    #print(paste(i, "/", no.files, ":", each) )

    each2 = paste(each, oldstr, sep="")
    splitted = splitString(each2, oldstr)

    neweach = concatenate(splitted, newstr)

    if(F) {
    neweach  =""
    for (is in splitted) {
       neweach= paste(neweach, newstr, sep=is)
    }
    }

    #oldeach  = paste(pathnet, each,   sep="")
    #neweach  = paste(pathnet, newstr, splitted[2], sep="")
    #a=file.rename(from=oldeach, to=neweach)
    #print(a)

    res = c(res, neweach)
  }
  return (res)  
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
