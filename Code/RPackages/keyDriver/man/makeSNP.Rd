\name{makeSNP}
\alias{makeSNP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Something
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeSNP(netpairsWtype, highlightNodes = NULL, edgecolorlevels, normColor = "grey", highColor = "red", normShape = "50", highShape = "50", directed = T, legendtable = NA, snafile = "tmp.sna")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{netpairsWtype}{
%%     ~~Describe \code{netpairsWtype} here~~
}
  \item{highlightNodes}{
%%     ~~Describe \code{highlightNodes} here~~
}
  \item{edgecolorlevels}{
%%     ~~Describe \code{edgecolorlevels} here~~
}
  \item{normColor}{
%%     ~~Describe \code{normColor} here~~
}
  \item{highColor}{
%%     ~~Describe \code{highColor} here~~
}
  \item{normShape}{
%%     ~~Describe \code{normShape} here~~
}
  \item{highShape}{
%%     ~~Describe \code{highShape} here~~
}
  \item{directed}{
%%     ~~Describe \code{directed} here~~
}
  \item{legendtable}{
%%     ~~Describe \code{legendtable} here~~
}
  \item{snafile}{
%%     ~~Describe \code{snafile} here~~
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
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(netpairsWtype, highlightNodes=NULL,
           edgecolorlevels,
           normColor="grey", highColor="red", 
           normShape="50",   highShape="50", directed= T, 
           legendtable=NA, snafile="tmp.sna")
{

    xfname = getFileName(snafile)
    fcys  = paste(xfname, "_cys.txt", sep="")
    fcysn = paste(xfname, "_cys-nodes.txt", sep="")
    write.table(netpairsWtype, fcys, sep="\t", quote=FALSE, col.names=T, row.names=FALSE)

    pcols = dim(netpairsWtype)[2]

    # consider both columns
    uniquenames    = union(as.character(netpairsWtype[,1]), as.character(netpairsWtype[,2]) )
    uniquenames    = sort(uniquenames)
    no.uniquenames = length(uniquenames)

    # initialize matrix with the size equal to no.uniquenames
    #
    name2idxMatrix = cbind(uniquenames, c(1:no.uniquenames) )
    

    # 0. make link index: A1 A2 T & A I
    # A1 A2 T & A I ==> A1 A2 T I1
    #
    leftIdx = merge(netpairsWtype, name2idxMatrix, by.x=1, by.y=1, all=F) 
    
    #  A1 A2 T I1 & A I ==> A2 A1 T I1 I2: I1 and I2 are the indices of A1 and A2 respectively
    #
    allIdx  = merge(leftIdx, name2idxMatrix, by.x=2, by.y=1, all=F)

    no.pairs = dim(allIdx)[1]

    if(pcols==2){
      no.elevels = length(edgecolorlevels)
      greyIdx    = c(1:no.elevels) [edgecolorlevels=="grey"]
      linksIdxMatrix = cbind( allIdx[,c(3,4)], rep(greyIdx,no.pairs) )
    } else{
      linksIdxMatrix = allIdx[,c(4,5,3)]
    }


    # 1. head string
    header = paste("*vertices ", as.character(no.uniquenames), sep="")

    # 2. vertices matrix 
    if ( is.null(highlightNodes) ){
        verticesMatrix = cbind(uniquenames, 
                               rep(normColor,no.uniquenames), 
                               rep(normShape,no.uniquenames) )
    } else{
      verticesMatrix = matrix("", no.uniquenames, 3)
      verticesMatrix[,1] = uniquenames
      verticesMatrix[,2] = rep(normColor,no.uniquenames)
      verticesMatrix[,3] = rep(normShape,no.uniquenames)
      xcolor= rep(normColor,no.uniquenames); names(xcolor) <- uniquenames
      xshape= rep(normShape,no.uniquenames); names(xshape) <- uniquenames

      # set color and shape for highlighted nodes
      #
      if ( !is.list(highlightNodes) ){
          highlightNodes2 = intersect(highlightNodes, uniquenames)
          xcolor[ highlightNodes2] = highColor
          xshape[ highlightNodes2] = highShape
      }else{
         no.highnodeSets = length(highlightNodes)
         for(il in c(1:no.highnodeSets) ) {
             highlightNodes2 = intersect(highlightNodes[[il]], uniquenames)
             if(length(highlightNodes2)==0){next};
             xcolor[highlightNodes2] = highColor[il]
             xshape[highlightNodes2] = highShape[il]
         }
      }
      verticesMatrix[,2] = as.character(xcolor)
      verticesMatrix[,3] = as.character(xshape)
    }

    #verticesMatrix <- as.matrix(verticesMatrix)
    colnames(verticesMatrix) <- c("nodename","color","shape")

    write.table(verticesMatrix, fcysn, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)


    #**************************************************************************
    #
    # 3. output indexed netpairs
    #
    # Legend
    if ( !is.na(legendtable) ) {
      mhead = paste("Legend", dim(legendtable)[1], sep=" ")
      write.table(as.matrix(mhead), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
      write.table(legendtable, snafile, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=T)
      # vertex
      write.table(as.matrix(header), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    } else{
      # vertex
      write.table(as.matrix(header), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=F)
    }

    write.table(verticesMatrix, snafile, sep="\t",quote=FALSE, col.names=T, row.names=FALSE, append=T)

    #edge color
    write.table(t(as.matrix(c("edge_color_levels", edgecolorlevels) )), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)

    # link pairs based index   
    #
    write.table(rbind(c("src","dst", "type")), snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
    write.table(linksIdxMatrix, snafile, sep="\t",quote=FALSE, col.names=F, row.names=FALSE, append=T)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
