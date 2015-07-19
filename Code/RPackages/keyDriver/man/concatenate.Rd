\name{concatenate}
\alias{concatenate}
\title{Concatenate a vector of strings, separated by a separator string}
\description{
This function takes a vector of strings as input and a string separator.  It
then pastes the vector of strings together using the separator as a filler
between vector entries.}
\usage{concatenate( myvect , mysep = "" )}
\arguments{
  \item{myvect}{vector of strings to be pasted together}
  \item{mysep}{string used as filler separating vector entries}
}
\details{
This function takes a vector of strings as input and a string separator.  It
then pastes the vector of strings together using the separator as a filler
between vector entries.}
\value{
Returns a single string consisting of each entry of myvect separated by mysep.}
\author{
Bin Zhang
}
\examples{
a <- c( "The" , "quick" , "brown" , "fox" , "jumped" , "over" , "the" , "lazy" , "dog" )
concatenate( a , " " )
}
\keyword{paste}
