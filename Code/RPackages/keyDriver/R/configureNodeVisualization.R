configureNodeVisualization <- function( allnodes , signature , kdaMatrix , bNodeSz = 40 , bFontSz = 12 )
{
	
	# SIG--signature; NSIG--not signature; GKD--Global KeyDriver; LKD--Local KeyDriver; NKD--Not KeyDriver
	#
	xcategories <- c( "SIG_GKD" , "SIG_LKD" , "SIG_NKD" , "NSIG_GKD" , "NSIG_LKD" , "NSIG_NKD" )
	xcolors <- c( "red" , "blue" , "lightgreen" , "red" , "blue" , "grey" )
	names( xcolors ) <- xcategories
	xshapes <- c( "square" , "square" , "circle" , "circle" , "circle", "circle" )
	names( xshapes ) <- xcategories
	xsizes <- c( 3 * bNodeSz , 2 * bNodeSz , bNodeSz , 3 * bNodeSz , 2 * bNodeSz , bNodeSz )
	names(xsizes) <- xcategories
	xfontsz <- c( 3 * bFontSz , 2 * bFontSz , bFontSz , 3 * bFontSz , 2 * bFontSz , bFontSz )
	names( xfontsz ) <- xcategories
	
	no.nodes <- length( allnodes )
	
	# legend table 
	legendtb <- cbind( xcategories , xshapes , xcolors , xcolors , xsizes , xfontsz )
	colnames( legendtb ) <- c( "label" , "shape" , "color" , "border" , "node_size" , "font_size" )
	
	sigInNet <- intersect( allnodes , signature )
	sigStatus <- rep( "NSIG" , no.nodes )
	names( sigStatus ) <- allnodes
	sigStatus[sigInNet] <- "SIG"
	kdrStatus <- rep( "NKD" ,  no.nodes )
	names( kdrStatus ) <- allnodes

	nf.cols <- dim( kdaMatrix )[2]
	nf.rows <- dim( kdaMatrix )[1]
	keydrvNames <- NULL
	if ( nf.rows > 0 )
	{
		keydrv <- as.integer( kdaMatrix[,nf.cols] )
		# global driver
		keysel <- which( keydrv == 1 )
		keydrvNames <- kdaMatrix[keysel,1]
		kdrStatus[keydrvNames] <- "GKD"

		# local driver
		if ( any( keydrv == 0 ) )
		{
			keysel <- which( keydrv == 0 )
			keydrvNames <- kdaMatrix[keysel,1]
			kdrStatus[keydrvNames] <- "LKD"
		}

		# combined signature-keydriver status
		#
		sigkdrStatus <- paste( sigStatus , kdrStatus , sep = "_" )
		hnList <- tapply( allnodes , sigkdrStatus , list ) # make a list for each category
		sigkdrNames <- names( hnList )

		isNonSig <- intersect( xcategories[4:6] , sigkdrNames ) # if all nodes are signatures, we use only circle for display
		if ( isTRUE( all.equal( length( isNonSig ) , 0 ) ) )
		{
			xshapes <- c( "circle" , "circle" , "circle" , "circle" , "circle" , "circle" )
			names( xshapes ) <- xcategories
		}

		# set up actual visualization properties
		yHighColor <- xcolors[sigkdrNames]
		yHighShape <- xshapes[sigkdrNames]
		yHighSize <- xsizes[sigkdrNames]
		yHighFontSZ <- xfontsz[sigkdrNames]
	}
	else
	{
		hnList <- list( sigInNet ) # highlight only signature
		yHighColor <- c( "brown" )
		yHighShape <- c( "circle" )
		yHighSize <- c( "1" )
		yHighFontSZ <- c( "1" )
	}
	return( list( hnList , cbind( yHighColor , yHighShape , yHighSize , yHighFontSZ ) , legendtb ) )
}
