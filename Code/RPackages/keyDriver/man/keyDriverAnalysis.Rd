\name{keyDriverAnalysis}
\alias{keyDriverAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Top Function of Key Driver Analysis
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Identify key regulatory components, or key drivers, of sub-networks with respect to varying biological contexts.
}
\usage{
keyDriverAnalysis = function(inputnetwork, signature, directed=T, nlayerExpansion=1,
		nlayerSearch=6, enrichedNodesPercentCut=-1, boostHubs=T, dynamicSearch=T, 
		FETpValueCut=0.05, useCorrectedpValue=T, outputfile=NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputnetwork}{
  an input network in form of Nx2 matrix
%%     ~~Describe \code{inputnetwork} here~~
}
  \item{signature}{
  an input signature for which key drivers need be identified
%%     ~~Describe \code{signature} here~~
}
  \item{directed}{
  input network is directed or un-directed. TRUE for directed network and FALSE for un-directed one.
%%     ~~Describe \code{} here~~
}
  \item{nlayerExpansion=}{
  how many layers you want to expand to find a subnetwork centering around the signature.
%%     ~~Describe \code{} here~~
}
  \item{nlayerSearch}{
  the maximum layer to be searched for the enrichment test.
%%     ~~Describe \code{} here~~
}
  \item{enrichedNodesPercentCut}{
  minimum percent of number of neighbors a key driver should have.
%%     ~~Describe \code{} here~~
}
  \item{boostHubs}{
  should we boost the nodes with many directed links to be key drivers?
%%     ~~Describe \code{} here~~
}
  \item{dynamicSearch}{
  should we search different layers for the best enrichment?
%%     ~~Describe \code{} here~~
}
  \item{FETpValueCut}{
  Fisher exact test p-value cutoff
%%     ~~Describe \code{} here~~
}
  \item{useCorrectedpValue}{
  should we use Bonferroni corrected Pvalue for KDA?
%%     ~~Describe \code{} here~~
}
  \item{outputfile}{
  filename (no extension) for outputing results into files. By default (NULL), the output is not written to files.
%%     ~~Describe \code{} here~~
}
}
\details{

%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of four components: subnetwork, keydrivers, parameters, files:
\item{subnetwork}{a 2 column matrix representing the expanded subnetwork for identifying key drivers.}
\item{keydrivers}{a matrix of inferred key drivers with various properties.}
\item{parameters}{the paramters used in the analysis}
\item{files}{the files holding the results of the expanded subnetwork and keydrivers and they can be taken as input of a Cytoscape plugin for visualization.}
}
\references{
Bin Zhang, Justin Guinney, Dave A Henderson and Jun Zhu. Key Driver Analysis in Gene Regulatory Networks. Manuscript.
Zhu, J., B. Zhang, et al. (2008). "Integrating large-scale functional genomic data to dissect the complexity of yeast regulatory networks." Nat Genet 40(7): 854-861
%% ~put references to the literature/web site here ~
}
\author{
Bin Zhang, bin.zhang@sagebase.org
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
keydriver_in_subnetwork()
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library( class )
library( cluster )
library( rpart )
library( sma ) # this is needed for plot.mat below
library( lattice ) # require is design for use inside functions 

data( yeastcausalnet )
data( yeastinputlist )

outputDir <- "KeyDrivers/"
if ( is.na( ( finfo <- file.info( outputDir ) )["isdir"] ) )
{
	dir.create( outputDir )
} else if ( !finfo["isdir"] ){
	error( "Output directory cannot be created as file exists with that name" )
}

################################################################################################
#    1. read in network

cnet <- as.matrix( fcausalnet )
totalnodes <- union( cnet[,1] , cnet[,2] )

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

fname <- "kda_yeast"
xkdFall    = paste(outputDir, fname, "_KDx_combined.xls",  sep='')
xkdFpara   = paste(outputDir, fname, "_KDx_parameters.xls",  sep='')
ykdFres    = paste(outputDir, fname, "_KDy_cys.xls",  sep='')
xkdrMatrix = NULL; paraMatrix=NULL

################################################################################################
# 3. process each gene list
#
resfiles = NULL
for (em in modules) {
	
	print (paste("*****************", em, "********************"))
	
	esel = listMatrix[, ncols] == em
	
	# remove abnormal gene names
	#
	genes = union(listMatrix[esel, 1], NULL)
	genes <- genes[genes!=""]
	genes <- genes[!is.na(genes)]
	no.genes <- length(genes)
	
	em2 = replaceString(em, ":","")
	em2 = replaceString(em2, " ","-")
	
	key2 = paste(outputDir, fname, "_KD_", em2, sep="")
	
	ret=keyDriverAnalysis(inputnetwork=cnet, signature=genes, directed=directed, nlayerExpansion=nlayerExpansion,
			nlayerSearch=nlayerSearch, enrichedNodesPercentCut=percentCut, boostHubs=boostHubs, dynamicSearch=dynamicSearch, 
			FETpValueCut=pValueCut, useCorrectedpValue=useCorrectedpValue, outputfile=key2)
	if(is.null(ret)) {next}
	
	fkd = ret[[2]]
	parameters = ret[[3]]
	
	fkd2 = cbind(rep(em, dim(fkd)[1]), fkd)
	xkdrMatrix = rbind(xkdrMatrix, fkd2)
	paraMatrix = rbind(paraMatrix, c(em2, parameters) )
	
	resfiles = rbind(resfiles, ret[[4]])
	
} #for (em in modules) {

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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{key driver analysis}
\keyword{regulatory network}
\keyword{gene network}
\keyword{gene regulators}
