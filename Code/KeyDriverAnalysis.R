b### Key Driver Analysis Script
### Kipp Johnson
### kipp.johnson@icahn.mssm.edu

### Adapted directly from help file of KeyDriver package

###
###
###
###
###

library(keyDriver);
rm(list=ls())

setwd("~/Projects/GelbRotation/NetworkData/KeyDriverAnalysis/Directed/")
################################################################################################
#### Set a global variable
var_nLayerExpansion <- 5 # Number of layers to expand network past signature (RASopathy) genes
####

# 1. Read in network and signature genes

# networkfilepath: directed network ( source {col1} --> target {col2 or col3} )
# listMatrix: gene list annotated (column 1) with module membership (column 2)

################################
#####                      #####
##### Network in 3 columns #####
#####                      #####
################################
### PF_AD_links_for_cytoscape
# networkfilepath <- '/Users/kwj/Projects/GelbRotation/NetworkData/Networks/annotated/PF_AD_links_for_cytoscape.anno.txt'
# listMatrix <- read.table('/Users/kwj/Projects/GelbRotation/NetworkData/Networks/nodes/PF_AD_links_for_cytoscape.nodes.txt',header=TRUE)
# colnum <- 3

### delimited_blood_network
# networkfilepath <- '/Users/kwj/Projects/GelbRotation/NetworkData/Networks/annotated/delimited_blood_network.anno.txt'
# listMatrix <- read.table('/Users/kwj/Projects/GelbRotation/NetworkData/Networks/nodes/delimited_blood_network.nodes.txt',header=TRUE)
# colnum <- 3

################################
#####                      #####
##### Network in 2 columns #####
#####                      #####
################################
### PF_normal_links_space_delimited
# networkfilepath <- '/Users/kwj/Projects/GelbRotation/NetworkData/Networks/annotated/PF_normal_links_space_delimited.anno.txt'
# listMatrix <- read.table('/Users/kwj/Projects/GelbRotation/NetworkData/Networks/nodes/PF_normal_links_space_delimited.nodes.txt',header=TRUE)
# colnum <- 2

### ileum_links_space_delimited
# networkfilepath <- '/Users/kwj/Projects/GelbRotation/NetworkData/Networks/annotated/ileum_links_space_delimited.anno.txt'
# listMatrix <- read.table('/Users/kwj/Projects/GelbRotation/NetworkData/Networks/nodes/ileum_links_space_delimited.nodes.txt',header=TRUE)
# colnum <- 2

### omental_links_space_delimited
# networkfilepath <- '/Users/kwj/Projects/GelbRotation/NetworkData/Networks/annotated/omental_links_space_delimited.anno.txt'
# listMatrix <- read.table('/Users/kwj/Projects/GelbRotation/NetworkData/Networks/nodes/omental_links_space_delimited.nodes.txt',header=TRUE)
# colnum <- 2

### delimited_risk_network
# networkfilepath <- '/Users/kwj/Projects/GelbRotation/NetworkData/Networks/annotated/delimited_risk_network.anno.txt'
# listMatrix <- read.table('/Users/kwj/Projects/GelbRotation/NetworkData/Networks/nodes/delimited_risk_network.nodes.txt',header=TRUE)
# colnum <- 2

### delimited_pan_intestine_network
networkfilepath <- '/Users/kwj/Projects/GelbRotation/NetworkData/Networks/annotated/delimited_pan_intestine_network.anno.txt'
listMatrix <- read.table('/Users/kwj/Projects/GelbRotation/NetworkData/Networks/nodes/delimited_pan_intestine_network.nodes.txt',header=TRUE)
colnum <- 2

net <- read.table(networkfilepath,header=TRUE)

if(colnum==2){
    # If network in 2 columns:
    cnet <- as.matrix(net[,1:2])}
  else{
    # If network in 3 columns:
    cnet <- as.matrix(net[,c(1,3)])
}

# Not necessary; done in KDA example
totalnodes <- union( cnet[,1] , cnet[,2] )

################################################################################################
# 2. Work with gene module lists

listMatrix <- listMatrix[listMatrix$Status=='RASopathy_gene',] # only use this one "module"

#dim( listMatrix )
listMatrix <- as.matrix( listMatrix ) # KDA requires matrix input
listMatrix[1:2,] # look at first two rows
ncols <- dim(listMatrix )[2]

modules <- names( table( listMatrix[,ncols] ) ) # modules to loop over

# 1. X_KDx_combined.xls: a combination of all the tables "X_Y_keydriver.xls"
# 2. X_KDx_parameters.xls: the parameters used for all the keydriver analyses
# 3. X_KDy_cys.xls:

################################################################################################
# 3. File output variables

## Function to get file name from some file path
## Probably not tolerant of edge cases, verify that it works for you
filename <- function(filepath){
  file_split <- strsplit(networkfilepath, split="/") # first split at /
  file_split <- strsplit(file_split[[1]][length(file_split[[1]])], "\\.") # then split at .
  return(file_split[[1]][1]) # return first item that was split at .
}

file_split <- filename(networkfilepath) # get file name

# Create a name to use as an output directory name for the analysis
outputDir <- paste("KD_",file_split,"_layer_",var_nLayerExpansion,"/", sep="")

# Create output directory of name outputdir
if ( is.na( ( finfo <- file.info( outputDir ) )["isdir"] ) )
{
  dir.create( outputDir )
} else if ( !finfo["isdir"] ){
  error( "Output directory already exists" )
}

fname <- file_split # copy variable for compatability with example
xkdFall    = paste(outputDir, fname, "_lay", var_nLayerExpansion, "_", "_KDx_combined.xls",  sep='')
xkdFpara   = paste(outputDir, fname, "_lay", var_nLayerExpansion, "_", "_KDx_parameters.xls",  sep='')
ykdFres    = paste(outputDir, fname, "_lay", var_nLayerExpansion, "_", "_KDy_cys.xls",  sep='')
xkdrMatrix = NULL; paraMatrix=NULL

################################################################################################
# 3. Process each gene list
#
resfiles = NULL
for (em in modules) { #looping over modules of interest

  print (paste("*****************", em, "********************")) # print module name

  esel <- listMatrix[, ncols] == em

  # remove abnormal gene names
  genes <- union(listMatrix[esel, 1], NULL)
  genes <- genes[genes!=""]
  genes <- genes[!is.na(genes)]
  no.genes <- length(genes)

  key2 <- paste(outputDir, fname, "_lay", var_nLayerExpansion, "_", "_KD_", em, sep="")

  # Using entirely default values, except for nLayerExpansion, which we set to **n**
  # Have also turned off useCorrectedpValue (It's a Bonferroni correction)
  ret <- keyDriverAnalysis(inputnetwork=cnet,
                           signature=genes,
                           directed=TRUE,
                           nlayerExpansion=var_nLayerExpansion,
                           nlayerSearch=6,
                           enrichedNodesPercentCut=0.1,
                           boostHubs=TRUE,
                           dynamicSearch=TRUE,
                           FETpValueCut=0.05,
                           useCorrectedpValue=FALSE,
                           outputfile=key2)

  # if ret is null, you have no key drivers (apparently)
  # this is poorly documented
  if(is.null(ret)) {print('ret is null')}

  fkd = ret[[2]]
  parameters = ret[[3]]

  fkd2 = cbind(rep(em, dim(fkd)[1]), fkd)
  xkdrMatrix = rbind(xkdrMatrix, fkd2)
  paraMatrix = rbind(paraMatrix, c(em, parameters) )

  resfiles = rbind(resfiles, ret[[4]])
}

for (em in modules) {

# save all key drivers
colnames(xkdrMatrix) = c("module", colnames(xkdrMatrix)[-1] )
write.table(xkdrMatrix, xkdFall, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

# save all data files
write.table(resfiles, ykdFres, sep="\t",quote=FALSE, col.names=F, row.names=FALSE)

# save parameters used
#
colnames(paraMatrix) <- c("subnet", colnames(parameters))
write.table(paraMatrix, xkdFpara, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)
}

