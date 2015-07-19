downStreamGenes <- function( netpairs , seednodes , N = 100 , directed = TRUE )
{
   prenodes <- seednodes
   cnt <- N
   pcdiff <- 1
   while( length( pcdiff ) > 0 && cnt > 0 )
   {
      retlinks <- findNLayerNeighborsLinkPairs( linkpairs = netpairs , subnetNodes = prenodes ,
               nlayers = 1 , directed = directed )
      if( is.null( retlinks ) )
	  {
		  return( NULL )
	  }
      curnodes <- union( retlinks[,1] , retlinks[,2] ) 
      pcdiff <- setdiff( curnodes , prenodes )
      prenodes <- curnodes
      cnt <- cnt - 1
   }

   if ( is.null( retlinks ) )
   {
	   return( NULL )
   }
   else
   {
      return( curnodes )
   }
}

