###################################################################################################
#
#  This file is the main program for Key Driver Analysis (KDA)
#
#  Function: to identify the key regulators for a list of genes with regard to a given network.
#
#  Author: Bin Zhang, PhD, Sage Bionetworks
#  Contact: bin.zhang@sagebase.org
#
#  Release Date: Feb 5, 2010
#  
#  Modified by ______________  Date ____________
#
#
#  Reference: 
#    1. Jun Zhu, Bin Zhang, Erin Smith, Becky Drees, Rachel Brem, Roger Bumgarner, Eric E. Schadt. 
#       (2008) Integrating Large-Scale Functional Genomic Data to Dissect the Complexity of Yeast 
#       Regulatory Networks. Nature Genetics 40: 854-861
#    2. Bin Zhang, Linh M. Tran, Jonathan Derry, Jun Zhu, Stephen Friend. Breast Cancer Transcriptional
#       Networks: Pathways, Regulators, and Prediction. Manuscript
#

##----------------------------------------------------------------------------------------------------
# Input: 
#   1. gene list: finputlist ="input_gene_lists.xls"
#   2. causal network: fcausalnet="input_breastcancer_BN.xls"; 
#   3. directed or undirected network: directed = T
#   4. expand subnetwork based on L-layer neighbors: layer=0
#   5. gene annotation file(NULL if not available): fgeneinfo
#
# Output:
#    1. keydrivers for each subnetwork: "*_keydriver.xls"
#    2. Cytoscape network: *_cys.txt
#    3. Cytoscape node properties: *_cys-nodes.txt
#    4. combined keydrivers for all runs: "_KDx_combined.xls"
#    5. parameters used for key driver analysis: "_KDx_parameters.xls"
#
#
# ---------------------------------- Parameters to be changed -----------------------------------------

# example 2: key drivers for yeast eQTL hotspots
#
data( yeastcausalnet )
data( yeastinputlist )
data( yeastgeneinfo )

directed <- TRUE
nlayerExpansion <- 1
nlayerSearch <- 4
percentCut <- 0.10
dynamicSearch <- TRUE
useCorrectedpValue <- TRUE
boostHubs <- TRUE
pValueCut <- 0.05

# 2. specify the directory for holding analysis results
#
if ( directed )
{
	outputDir <- "KeyDrivers/"
} else
{
	outputDir <- "KeyDrivers-undirected/"
}

if ( is.na( ( finfo <- file.info( outputDir ) )["isdir"] ) )
{
	dir.create( outputDir )
} else if ( !finfo["isdir"] )
{
	error( "Output directory cannot be created as file exists with that name" )
}

#
# -----------------------------End of Parameters to be changed --------------------------------------

library( class )
library( cluster )
library( rpart )
library( sma ) # this is needed for plot.mat below
library( lattice ) # require is design for use inside functions 

if ( ( Sys.info() )["sysname"] == "Windows" )
{
	memory.size( TRUE )   # check the maximum memory that can be allocated
	memory.limit( size = 3800 )   # increase the available memory
}

################################################################################################
#    1. read in network

cnet <- as.matrix( fcausalnet )
dim( cnet )

totalnodes <- union( cnet[,1] , cnet[,2] )

fname <- "yeast"
fname <- paste( fname , "_L" , nlayerExpansion , sep = "" )

################################################################################################
# 2. read in gene lists

listMatrix <- finputlist
dim( listMatrix )
listMatrix <- as.matrix( listMatrix )
listMatrix[1:2,]
ncols <- dim( listMatrix )[2]

modules <- names( table( listMatrix[,ncols] ) )

# 1. X_KDx_combined.xls: a combination of all the tables "X_Y_keydriver.xls"
# 2. X_KDx_parameters.xls: the parameters used for all the keydriver analyses
# 3. X_KDy_cys.xls: 

xkdFall <- paste( outputDir , fname , "_KDx_combined.xls" , sep = '' )
xkdFpara <- paste( outputDir , fname , "_KDx_parameters.xls" , sep = '' )
ykdFres <- paste( outputDir , fname , "_KDy_cys.xls" ,  sep = '' )
xkdrMatrix <- NULL
paraMatrix <- NULL

################################################################################################
# 3. process each gene list
#
resfiles <- NULL
for ( em in modules )
{

	print( paste( "*****************" , em , "********************" ) )

	esel <- listMatrix[, ncols] == em

	# remove abnormal gene names
	#
	genes <- union( listMatrix[esel,1] , NULL )
	genes <- genes[genes!=""]
	genes <- genes[!is.na( genes )]
	no.genes <- length( genes )

	em2 <- replaceString( em , ":" , "" )
	em2 <- replaceString( em2 , " " , "-" )

	key2 <- paste( outputDir , fname , "_KD_" , em2 , sep = "" )

	ret <- keyDriverAnalysis( inputnetwork = cnet , signature = genes , directed = directed ,
			nlayerExpansion = nlayerExpansion , nlayerSearch = nlayerSearch ,
			enrichedNodesPercentCut = percentCut , boostHubs = boostHubs ,
			dynamicSearch = dynamicSearch , FETpValueCut = pValueCut ,
			useCorrectedpValue = useCorrectedpValue , outputfile = key2 )
	if ( is.null( ret ) )
	{
		next
	}

	fkd <- ret[[2]]
	parameters <- ret[[3]]

	fkd2 <- cbind( rep( em , dim( fkd )[1] ) , fkd )
	xkdrMatrix <- rbind( xkdrMatrix , fkd2 )
	paraMatrix <- rbind( paraMatrix , c( em2 , parameters ) )

	resfiles <- rbind( resfiles , ret[[4]] )

} #for (em in modules) {

# save all key drivers
colnames( xkdrMatrix ) <- c( "module" , colnames( xkdrMatrix )[-1] ) 
write.table( xkdrMatrix , xkdFall , sep = "\t" , quote = FALSE , col.names = TRUE ,
		row.names = FALSE )

# save all data files
write.table( resfiles , ykdFres , sep = "\t" , quote = FALSE , col.names = FALSE ,
		row.names = FALSE )

# save parameters used
#
colnames( paraMatrix ) <- c( "subnet" , colnames( parameters ) )
write.table( paraMatrix , xkdFpara , sep = "\t" , quote = FALSE , col.names = TRUE ,
		row.names = FALSE )

## ------------------------------------- END ------------------------------------------------
