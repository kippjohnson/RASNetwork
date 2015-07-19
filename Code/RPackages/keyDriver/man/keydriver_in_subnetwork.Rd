\name{keydriverInSubnetwork}
\alias{keydriverInSubnetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Something
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
keydriverInSubnetwork(linkpairs, signature, directed = T, nlayers = 6, min_downstreamnodes = -1, FETpValueCut = 0.05, boostHubs = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{linkpairs}{
%%     ~~Describe \code{linkpairs} here~~
}
  \item{signature}{
%%     ~~Describe \code{signature} here~~
}
  \item{directed}{
%%     ~~Describe \code{directed} here~~
}
  \item{nlayers}{
%%     ~~Describe \code{nlayers} here~~
}
  \item{min_downstreamnodes}{
%%     ~~Describe \code{min_downstreamnodes} here~~
}
  \item{FETpValueCut}{
%%     ~~Describe \code{FETpValueCut} here~~
}
  \item{boostHubs}{
%%     ~~Describe \code{boostHubs} here~~
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
function(linkpairs, signature, directed=T, nlayers=6, 
                                   min_downstreamnodes=-1, FETpValueCut=0.05, boostHubs=T) {

allnodes = union(linkpairs[,1], linkpairs[,2])
no.subnetsize = length(allnodes)

# whole network nodes as the signature
networkAsSignature = length(setdiff(allnodes, signature)) ==0

overlapped    = intersect(allnodes, signature)
no.overlapped = length(overlapped) # within the subnetwork

keydrivers= NULL
kdMatrix  = NULL
kdIndex   = NULL # indices of keydrivers in dsnodesList

dsnodesList = as.list(rep(0,no.subnetsize)); no.dsnodes = rep(0, no.subnetsize)
cnt = 1

intv = as.integer(no.subnetsize/10)
print("find downstream genes")

for(i in c(1:no.subnetsize) ) {
   if(i\%\%intv==0){ print(paste(i, "/", no.subnetsize)) }

   #netpairs=linkpairs; seednodes=allnodes[i]; N=nlayers; directed=directed
   idn = downStreamGenes(netpairs=linkpairs, seednodes=allnodes[i], N=nlayers, directed=directed)
   idn = setdiff(idn, allnodes[i])
   no.idn=length(idn); no.dsnodes[i] = no.idn
   if(no.idn>0) {
     dsnodesList[[i]] = idn # record the down stream genes
   }
  }

mymincut = min_downstreamnodes
if (min_downstreamnodes<=0){
   mymincut = mean(no.dsnodes) + sd(no.dsnodes)
  }
cutmatrix = c( mean(no.dsnodes), sd(no.dsnodes), mymincut)

print("Fisher exact test ...\n")

for(i in c(1:no.subnetsize) ) {
   
   no.idn=no.dsnodes[i];
   if (no.idn <mymincut) {next} #use the cut only for the unexpanded network

   idn = dsnodesList[[i]]
   
   pv = 0; no.hits = no.idn
   if(!networkAsSignature){
      hits    = intersect(idn, overlapped)
      no.hits = length(hits)
      if(no.hits==0){next}
      pv = phyper(no.hits - 1, no.idn, no.subnetsize-no.idn, no.overlapped, lower.tail=F)
   }

   res= c(no.hits, no.idn, no.overlapped, length(signature), no.subnetsize, pv)

   if(pv>FETpValueCut){next}

   kdMatrix  = rbind(kdMatrix, res)
   keydrivers= c(keydrivers, allnodes[i])
   kdIndex = c(kdIndex, i)

   cnt       = cnt + 1

   print(res)
  }
n.drivers = length(keydrivers)

if(is.null(kdMatrix)){return (NULL)}

#******************* local driver or not **************************************
#
# check whether a driver is in the downstream of other drivers
keydrv = rep(0, n.drivers)
#if (!networkAsSignature) {
for ( i in c(1:n.drivers) ) {

    iselA  = kdMatrix[,2] > kdMatrix[i,2]
    isel   = kdIndex[iselA] # restrict to only candidate drivers 

    if ( sum(isel)>0) {
       if(directed) {
          ilocal= setInSets(setC=keydrivers[i], setlist=dsnodesList[isel])
       } else{
          ilocal= setInSets(setC=dsnodesList[[ kdIndex[i] ]], setlist=dsnodesList[isel])
       }

       keydrv[i] = !ilocal + 0
    } else{
       keydrv[i] = T
    }
  }
#}

# promote genes with many direct links to be key drivers
#
#              inlinks outlinks totallinks
#0610031J06Rik       2        0          2
#1110001J03Rik       0        1          1
#
if(boostHubs) {
  
  if(!networkAsSignature){
     # for expanded network, restrict the boosted nodes to the key driver candidates
     kdSel = rep(F, no.subnetsize); kdSel[kdIndex]=T;
  } else{
     # for non-expanded network, consider all the nodes in the subnetwork
     kdSel = rep(T, no.subnetsize);
  }

  mydegree  = degree_ByLinkPairs(linkpairs=linkpairs, directed=directed, cleangarbage=F)
  if(directed) {
     directSel = mydegree[,2]> mean(mydegree[,2]) + 2*sd(mydegree[,2])
     cutmatrix = rbind( c(mean(no.dsnodes), sd(no.dsnodes), mymincut,
                       mean(mydegree[,2]),sd(mydegree[,2]), mean(mydegree[,2]) + 2*sd(mydegree[,2]) ))
  }else{
     directSel = mydegree[,3]> mean(mydegree[,3]) + 2*sd(mydegree[,3])
     cutmatrix = rbind( c(mean(no.dsnodes), sd(no.dsnodes), mymincut,
                       mean(mydegree[,3]),sd(mydegree[,3]), mean(mydegree[,3]) + 2*sd(mydegree[,3])))
  }
  directSel = directSel & kdSel

  directeHub= rownames(mydegree)[directSel]
  isDirectHub = setElementInSet(keydrivers, directeHub)
  keydrv[isDirectHub] = T

  colnames(cutmatrix) <- c("mean_downstream", "sd_downstream", "cut_downstream",
                           "mean_degree", "sd_degree", "cut_degree")

  } else{
  cutmatrix = rbind( c(mean(no.dsnodes), sd(no.dsnodes), mymincut, "F"))
  colnames(cutmatrix) <- c("mean_downstream", "sd_downstream", "cut_downstream", "boost_directhubs")
  }


# in this case, signature is the network nodes themselves, so pvalue will be 0 for all nodes
# so the driver will be the ones with most downsttream genes
#
if(networkAsSignature){
  mo = order(-kdMatrix[,1])
  } else{
  mo = order(kdMatrix[,6])
  }

fkd = cbind(keydrivers, kdMatrix, keydrv); nf.cols = dim(fkd)[2]

if(dim(fkd)[1] >1 ) {
 fkd = fkd[mo, ]
 # put key driver on the top
 mo  = order( -as.integer(fkd[,nf.cols]) )
 fkd = fkd[mo, ]
  }

colnames(fkd) <- c("keydrivers", "hits", "downstream", "signature_in_network", "signature", "network_size", "pvalue", "keydriver")

ret = as.list(c(1:2))
ret[[1]] = fkd
ret[[2]] = cutmatrix

return (ret)

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
