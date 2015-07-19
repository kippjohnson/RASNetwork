removeDuplicatedLinks <- function( linkpairs , directed = FALSE )
{
    if ( isTRUE( all.equal( dim( linkpairs )[1] , 1 ) ) )
	{
       return( linkpairs )
    }
    links <- paste( linkpairs[,1] , linkpairs[,2] , sep = "\t" )
    # 1. remove duplications 
    #
    cleanedlinkMatrix <- union( links , NULL )
    # 2. remove inversed duplications
    #
    linkMatrix  <- as.matrix( getAllParts( cleanedlinkMatrix , "\t" ) )
    if ( directed || isTRUE( all.equal( dim( linkMatrix )[1] , 1 ) ) )
	{
		return( linkMatrix )
	}
    #  first, remove self-self interactions
    #
    selfSelfLinks <- linkMatrix[,1] == linkMatrix[,2]
    linkMatrix <- linkMatrix[!selfSelfLinks,]
    cleanedlinkMatrix <- cleanedlinkMatrix[!selfSelfLinks]
	# Now, create reverse links
    reversedLinks <- cbind( paste( linkMatrix[,2] , linkMatrix[,1] , sep = "\t" ) , c( 1:length( cleanedlinkMatrix ) ) )
	removedCols <- as.integer( merge( cleanedlinkMatrix , reversedLinks , by.x = 1 , by.y = 1 , all = FALSE )[,2] )
	if ( length( removedCols ) > 0 )
	{
        # construct non-duplicated interactions
        #
        dupLinks <- cleanedlinkMatrix[removedCols]
        dupLinksRev <- reversedLinks[removedCols]
        uniques <- NULL
        for ( i in c( 1:length( dupLinks ) ) )
        {
           if ( !( is.element( dupLinks[i] , uniques ) | is.element( dupLinks[i] , uniques ) ) )
		   {
               uniques <- c( uniques , dupLinks[i] )
           }
        }
        xlinkMatrix <- c( cleanedlinkMatrix[-removedCols] , uniques )
    }
	else
	{
        xlinkMatrix <- cleanedlinkMatrix
    }
    return( getAllParts( xlinkMatrix , "\t" ) )
}

