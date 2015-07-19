keyDriverAnalysis <- function( inputnetwork , signature , directed = TRUE , nlayerExpansion = 1 ,
		nlayerSearch = 6 , enrichedNodesPercentCut = -1 , boostHubs = TRUE , dynamicSearch = TRUE ,
		FETpValueCut = 0.05 , useCorrectedpValue = TRUE , outputfile = NULL )
{
	if ( !is.null( outputfile ) )
	{
		#onetFname <- paste( outputfile , outputDir , key2 , ".pair" , sep = '' )
		snpFname <- paste( outputfile , ".snp" ,  sep = '' )
		kdFname <- paste( outputfile , "_keydriver.xls" , sep = '' )
	}

	# overlap between network & signature
	wholenodes <- union( inputnetwork[,1] , inputnetwork[,2] )
	no.wholenodes <- length( wholenodes )
	wholeOvlp <- intersect( wholenodes , signature )
	no.wholeOvlp <- length( wholeOvlp )

	if ( nlayerExpansion >= 1 )
	{
		# expand network by n-layer nearest neighbors
		expandNet <- findNLayerNeighborsLinkPairs( linkpairs = inputnetwork ,
				subnetNodes = signature , nlayers = nlayerExpansion , directed = directed )
	}
	else if ( isTRUE( all.equal( nlayerExpansion , 0 ) ) )
	{
		# no expansion
		expandNet <- getSubnetworkLinkPairs( linkpairs = inputnetwork ,
				                             subnetNodes = signature )
	}
	else
	{
		expandNet <- inputnetwork
	}

	print( paste( "dim(expandNet): " , dim( expandNet ) ) )

	allnodes <- sort( union( expandNet[,1] , expandNet[,2] ) )
	no.nodes <- length( allnodes )
	#write.table(expandNet, onetFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

	# convert IDs into indices
	netIdxSrc <- getMatchedIndexFast( allnodes , expandNet[,1] )
	netIdxDst <- getMatchedIndexFast( allnodes , expandNet[,2] )
	signatIdx <- getMatchedIndexFast( allnodes , intersect( allnodes , signature ) )
	expandNetIdx <- cbind( netIdxSrc , netIdxDst )

	################################################################################################
	# 4. keydriver for a given network
	#
	#linkpairs=expandNet;signature=genes;directed=directed; nlayers=6; min_downstreamnodes=min_ds_cut; FETpValueCut=0.05; boostHubs=T; dynamicSearch=dynamicSearch
	
	if ( directed )
	{
		ret <- keydriverInSubnetwork( linkpairs = expandNetIdx , signature=signatIdx ,
				background = c( no.wholenodes , no.wholeOvlp ) , directed = directed ,
				nlayers = nlayerSearch , enrichedNodesPercentCut = enrichedNodesPercentCut ,
				FETpValueCut = FETpValueCut , boostHubs = boostHubs ,
				dynamicSearch = dynamicSearch , bonferroniCorrection = useCorrectedpValue )
	}
	else
	{
		ret <- keydriverInSubnetwork( linkpairs = expandNetIdx , signature = signatIdx ,
				background = c( no.wholenodes , no.wholeOvlp ) , directed = directed ,
				nlayers = nlayerSearch , enrichedNodesPercentCut = enrichedNodesPercentCut ,
				FETpValueCut = FETpValueCut , boostHubs = boostHubs ,
				dynamicSearch = dynamicSearch , bonferroniCorrection = useCorrectedpValue )
	}

	if ( is.null( ret ) )
	{
		return( NULL )
	}

	# retrieve results
	#
	fkd <- ret[[1]]
	parameters <- ret[[2]]

	print(fkd)
	fkd[,1] <- allnodes[as.integer( fkd[,1] )]   

	if ( !is.null( outputfile ) )
	{
		write.table( fkd , kdFname , sep="\t" , quote = FALSE , col.names = TRUE ,
				     row.names = FALSE )

		################################################################################################
		#  output networks & key drivers for visualization
		#
		#     Cytoscape output: 1) network file - *_cys.txt 2) node property file: *_cys-nodes.txt
		#

		nodeprop <- configureNodeVisualization( allnodes = allnodes , signature = genes ,
				                                kdaMatrix = fkd )

		hnList <- nodeprop[[1]] # node subcategpries
		listprop <- nodeprop[[2]] # visual properties for each subcategory
		legend <- nodeprop[[3]] # legend table for visual properties

		resf <- makeSNP( netpairsWtype = expandNet , edgecolorlevels = c( "grey" ) ,
				  highlightNodes = hnList , normColor = "grey" , highColor = listprop[,1] ,
				  normShape = "circle" , highShape = listprop[,2] , normNodeSize = "40" ,
				  highNodeSize = listprop[,3] , normFontSize = "12" ,
				  highFontSize = listprop[,4] , legendtable = legend , snafile = snpFname )

		result <- list( expandNet , fkd , ret[[2]] , getFileFullNameNopath( resf ) )
#		result <- list( expandNet , fkd , ret[[2]] , getFileName( snpFname ) )
		names( result ) <- c( "subnetwork" , "keydrivers" , "parameters" , "files" )
	}
	else
	{
		result <- list( expandNet , fkd , ret[[2]] )
		names( result ) <- c( "subnetwork" , "keydrivers" , "parameters" )
	}

	return( result )
}
