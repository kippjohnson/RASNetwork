keydriverInSubnetwork <- function( linkpairs , signature , background = NULL ,
		directed = TRUE , nlayers = 6 , enrichedNodesPercentCut = -1 , FETpValueCut = 0.05 ,
		boostHubs = TRUE , dynamicSearch = TRUE , bonferroniCorrection = TRUE )
{
	allnodes <- union( linkpairs[,1] , linkpairs[,2] )
	no.subnetsize <- length( allnodes )

	# whole network nodes as the signature
	networkAsSignature <- isTRUE( all.equal( length( setdiff( allnodes , signature ) ) , 0 ) )
	
	overlapped <- intersect( allnodes , signature )
	no.overlapped <- length( overlapped ) # within the subnetwork

	if ( is.null( background ) )
	{
		background2 <- c( no.subnetsize , no.overlapped ) 
	}
	else
	{
		background2 <- background
	}

	keydrivers <- NULL
	kdMatrix <- NULL
	kdIndex <- NULL # indices of keydrivers in dsnodesList

	dsnodesList <- as.list( rep( 0 , no.subnetsize ) )
	no.dsnodes <- rep( 0 , no.subnetsize )
	cnt <- 1

	intv <- as.integer( no.subnetsize / 10 )
	print( "find downstream genes" )

	# set up searching range
	if ( dynamicSearch )
	{ # dynamic search for optimal layer
		layers4Search <- c( 1:nlayers )
	} else{  # fixed layer for optimal layer
		layers4Search <- c( nlayers )
	}
	# if the network itself is the signature, no need for dynamic search
	if ( networkAsSignature )
	{
		layers4Search <- c( nlayers )
	}

	for ( i in c( 1:no.subnetsize ) )
	{
		if ( isTRUE( all.equal( i%%intv , 0 ) ) )
		{
			print( paste( i , "/" , no.subnetsize ) )
		}

		# initialization
		minpv <- 1
		minNoHits <- 0
		minNoIdn <- 0
		minLayer <- 0
		minDn <- 0
		minFc <- 0
		minpvW <- 1
		minFcW <- 0
		for ( y in layers4Search )
		{
			#netpairs=linkpairs; seednodes=allnodes[i]; N=nlayers; directed=directed
			idn <- downStreamGenes( netpairs = linkpairs , seednodes = allnodes[i] ,
					N = y , directed = directed )
			idn <- setdiff(idn, allnodes[i])
			no.idn <- length(idn)

			if ( !networkAsSignature )
			{# do enrichment test for only expanded subnetwork
				hits <- intersect( idn , overlapped )
				no.hits <- length( hits )

				if ( isTRUE( all.equal( no.hits , 0 ) ) )
				{
					next
				}

				foldchg <- ( no.hits / no.idn ) / ( no.overlapped / no.subnetsize )
				pv <- phyper( no.hits - 1 , no.idn , no.subnetsize - no.idn ,
						no.overlapped , lower.tail = FALSE )

				foldchgW <- ( no.hits / no.idn ) / ( background2[2] / background2[1] )
				pvW <- phyper( no.hits - 1 , no.idn , background2[1] - no.idn ,
						background2[2] , lower.tail = FALSE )

				if ( pv < minpv )
				{
					minpv <- pv
					minNoHits <- no.hits
					minNoIdn <- no.idn
					minLayer <- y
					minFc <- foldchg
					minpvW <- pvW
					minFcW <- foldchgW
				}
			}
			else
			{ # for non-expanded subnetwork
				no.hits <- no.idn
				minpv <- 0
				minNoHits <- no.idn
				minNoIdn <- no.idn
				minLayer <- y
				minFc <- 1
			}
		} #y
		
		# record the down stream genes for the biggest layer
		if ( no.idn > 0 )
		{
			dsnodesList[[i]] <- idn
			no.dsnodes[i] <- no.idn
		}

		correctMinPv <- minpv * no.subnetsize
		correctMinPv <- ifelse( correctMinPv > 1 , 1 , correctMinPv )
		
		res <- c( minNoHits , minNoIdn , no.overlapped , no.subnetsize , background2[2] ,
				background2[1] , length( signature ) , minLayer , minFcW , minpvW , minFc ,
				minpv , correctMinPv )
		kdMatrix <- rbind( kdMatrix , res )
		#print(res)
	}
	
	mymincut <- enrichedNodesPercentCut * no.overlapped
	if ( enrichedNodesPercentCut <= 0 )
	{
		mymincut = mean(no.dsnodes) + sd(no.dsnodes)
	}
	cutmatrix <- c( mean( no.dsnodes ) , sd( no.dsnodes ) , mymincut )

	# pick up key drivers by pvalue and no. of downstream genes
	ncols <- dim( kdMatrix )[2]
	
	if ( bonferroniCorrection )
	{ # use corrected pvalue
		kdSel <- ( kdMatrix[,ncols] < FETpValueCut ) & ( kdMatrix[,2] >= mymincut )
	}
	else
	{
		kdSel <- ( kdMatrix[,ncols-1] < FETpValueCut) & (kdMatrix[,2] >= mymincut )
	}

	if ( sum( kdSel ) > 0 )
	{
		keydrivers <- allnodes[kdSel]
		kdIndex <- c( 1:no.subnetsize )[kdSel]
		n.drivers <- length( keydrivers )

		#******************* local driver or not **************************************
		#
		# check whether a driver is in the downstream of other drivers
		keydrv <- rep( 0 , no.subnetsize )
		#if (!networkAsSignature) {
		for ( i in c( 1:n.drivers ) )
		{
			# Note that kdIndex[i] is the index of ith keydriver in kdMatrix  
			# restrict to only candidate drivers 
			iselA <- ( kdMatrix[,2] > kdMatrix[kdIndex[i],2] ) & kdSel
			isel <- c( 1:no.subnetsize )[iselA]
			if ( sum(isel) > 0 )
			{
				if ( directed )
				{
					ilocal <- setInSets( setC = allnodes[kdIndex[i]] ,
							setlist = dsnodesList[isel] )
				}
				else
				{
					ilocal <- setInSets( setC = dsnodesList[[kdIndex[i]]] ,
							setlist = dsnodesList[isel] )
				}
				keydrv[kdIndex[i]] <- !ilocal + 0
			}
			else
			{
				keydrv[kdIndex[i]] <- TRUE
			}
		}
		#}
	}
	
	# promote genes with many direct links to be key drivers
	#
	#              inlinks outlinks totallinks
	#0610031J06Rik       2        0          2
	#1110001J03Rik       0        1          1
	#
	if ( boostHubs )
	{
		
		if ( !networkAsSignature )
		{
			# for expanded network, restrict the boosted nodes to the key driver candidates
			kdSelB <- rep( FALSE , no.subnetsize )
			kdSelB[kdIndex] <- TRUE
			psel <- kdMatrix[,ncols-3] * no.subnetsize < 0.05
			kdPsd <-1
			kdPmean <- 1
			kdpvalues <- -log10( kdMatrix[,ncols-3] )
			kdpvalues <- ifelse( is.na( kdpvalues ) , 0 , kdpvalues )
			#histogram(kdpvalues)
			if ( sum( psel ) > 0 )
			{
				kdPmean <- mean( kdpvalues[psel] )
				kdPsd <- sd( kdpvalues[psel] )
				#kdPmean= median(kdpvalues[psel]); kdPsd= mad(kdpvalues[psel])
				print( as.numeric( signif( kdpvalues[psel] , 2 ) ) )
				directSel <- ( kdpvalues > ( kdPmean + kdPsd ) )
				directSel <- ifelse( is.na( directSel ) , FALSE , directSel )
				#print(directSel)
				if ( sum( directSel ) > 0 )
				{
					kdSel <- kdSel | directSel
					dIndex <- c( 1:no.subnetsize )[kdSel]
					keydrv <- rep( FALSE , no.subnetsize )
					keydrv[dIndex] <- TRUE
				}
			}
			cutmatrix <- rbind( c( mean( no.dsnodes ) , sd( no.dsnodes ) , mymincut ,
							kdPmean , kdPsd , kdPmean + kdPsd ) )

			colnames( cutmatrix ) <- c( "mean_downstream" , "sd_downstream" ,
					"enrichedNodes_cut" , "mean_logP" , "sd_logP" , "cut_logP" )
		}
		else
		{
			# for non-expanded network, consider all the nodes in the subnetwork
			kdSelB <- rep( TRUE , no.subnetsize )
			# align the degree with allnodes
			mydegree <- degreeByLinkPairs( linkpairs = linkpairs , directed = directed ,
					cleangarbage = FALSE )
			mIdx <- getMatchedIndexFast( rownames( mydegree ) , allnodes )
			mydegree <- mydegree[mIdx,]
			if ( directed )
			{
				directSel <- mydegree[,2] > mean( mydegree[,2] ) + 2 * sd( mydegree[,2] )
				cutmatrix <- rbind( c( mean( no.dsnodes ) , sd( no.dsnodes ) , mymincut ,
								mean( mydegree[,2] ) , sd( mydegree[,2] ) ,
								mean( mydegree[,2] ) + 2 * sd( mydegree[,2] ) ) )
			}
			else
			{
				directSel <- mydegree[,3] > mean( mydegree[,3] ) + 2 * sd( mydegree[,3] )
				cutmatrix <- rbind( c( mean( no.dsnodes ) , sd( no.dsnodes ) , mymincut ,
								mean( mydegree[,3] ) , sd( mydegree[,3] ) ,
								mean( mydegree[,3] ) + 2 * sd( mydegree[,3] ) ) )
			}
			directSel <- directSel & kdSelB

			directeHub <- rownames( mydegree )[directSel]
			isDirectHub <- setElementInSet( allnodes , directeHub )

			keydrv[isDirectHub] <- TRUE
			kdSel <- kdSel | isDirectHub
			colnames( cutmatrix ) <- c( "mean_downstream" , "sd_downstream" ,
					"cut_downstream" , "mean_degree" , "sd_degree" , "cut_degree" )
		}
	}
	else
	{
		cutmatrix <- rbind( c( mean( no.dsnodes ) , sd( no.dsnodes ) , mymincut , "F" ) )
		colnames( cutmatrix ) <- c( "mean_downstream" , "sd_downstream" , "cut_downstream" ,
				"boost_directhubs" )
	}
	
	if ( isTRUE( all.equal( sum( kdSel ) , 0 ) ) )
	{
		return( NULL )
	}
	
	##
	# in this case, signature is the network nodes themselves, so pvalue will be 0 for all nodes
	# so the driver will be the ones with most downstream genes
	#
	isSignature <- rep( 0 , no.subnetsize )
	names( isSignature ) <- allnodes
	isSignature[overlapped] <- 1

	fkd <- cbind( allnodes , isSignature , kdMatrix , keydrv + 0 )[kdSel,]

	if ( sum( kdSel ) > 1 )
	{
		nf.cols <- dim( fkd )[2]
		if ( networkAsSignature )
		{
			mo <- order( -as.integer( fkd[,3] ) )
		}
		else
		{
			mo <- order( as.numeric( fkd[,nf.cols - 1] ) )
		}

		fkd <- fkd[mo,]
		# put key driver on the top
		mo <- order( -as.integer( fkd[,nf.cols] ) )
		fkd <- fkd[mo,]
	}
	else
	{
		fkd <- rbind( fkd )
	}
	
	colnames( fkd ) <- c( "keydrivers" , "isSignature" , "hits" , "downstream" ,
			"signature_in_subnetwork" , "subnetwork_size" , "signature_in_network" ,
			"network_size" , "signature" , "optimal_layer" , "fold_change_whole" ,
			"pvalue_whole" , "fold_change_subnet" , "pvalue_subnet" ,
			"pvalue_corrected_subnet" , "keydriver" )

	print( fkd )

	ret <- as.list( c( 1:2 ) )
	ret[[1]] <- fkd
	ret[[2]] <- cutmatrix

	return( ret )
}
