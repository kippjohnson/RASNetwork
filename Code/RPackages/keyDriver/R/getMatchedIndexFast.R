getMatchedIndexFast <- function( cvector , subvect )
{
	orgIndex <- cbind( cvector , c( 1:length( cvector ) ) )
	subIndex <- cbind( subvect , c( 1:length( subvect ) ) )
	merged <- as.data.frame( merge( subIndex , orgIndex , by.x = 1 , by.y = 1 , all.x = TRUE ) )
	if ( dim( merged )[1] > 1 )
	{
		od <- order( merged[,2] )  # restore the original order of subvect
		merged <- merged[od,]
	}
	return( merged[,3] )
}
