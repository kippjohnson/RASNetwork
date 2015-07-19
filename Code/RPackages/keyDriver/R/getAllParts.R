getAllParts <- function( fullfnames , sep = "-" , retLen = FALSE )
{
  unSplit <- function( x ) { unlist( strsplit( x , sep ) ) }
  ReturnMatrix <- function( x , returnMatrixEnv , unSplit )
  {
	  usplit <- unSplit( x )
	  returnMatrixEnv$returnMatrixList[[returnMatrixEnv$i]] <- usplit
	  returnMatrixEnv$i <- returnMatrixEnv$i + 1
	  return( length( usplit ) )
  }
  makeMatrix <- function( x , returnMatrixEnv , rowLengths )
  {
	  returnMatrixEnv$returnMatrix[returnMatrixEnv$i,1:rowLengths[returnMatrixEnv$i]] <- x
	  returnMatrixEnv$i <- returnMatrixEnv$i + 1
	  return( rowLengths[returnMatrixEnv$i-1] )
  }
  returnLengths <- function( x , unSplit ) { return( length( unSplit( x ) ) ) }
  if ( retLen )
  {
	return( unlist( lapply( fullfnames , returnLengths , unSplit ) ) )
  }
  else
  {
    returnMatrixEnv <- new.env()
	returnMatrixEnv$returnMatrixList <- vector( "list" , length = length( fullfnames ) )
	returnMatrixEnv$i <- 1
    rowLengths <- unlist( lapply( fullfnames , ReturnMatrix , returnMatrixEnv , unSplit ) )
	returnMatrixEnv$returnMatrix <- matrix( "" , nrow = length( rowLengths ) , ncol = rowLengths[which.max( rowLengths )] )
	returnMatrixEnv$i <- 1
	rowLengths <- unlist( lapply( returnMatrixEnv$returnMatrixList , makeMatrix , returnMatrixEnv , rowLengths ) )
	return( returnMatrixEnv$returnMatrix )
  }
}

