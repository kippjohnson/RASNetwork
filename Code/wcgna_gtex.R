
### GTEx analysis
rm(list=ls());
options(stringsAsFactors = FALSE);

library(data.table);
library(WGCNA);

## Read in and transform the data to a matrix
gte <- fread("~/Desktop/GTEx_Analysis_V4_eQTLInputFiles_geneLevelNormalizedExpression/Heart_Left_Ventricle.expr.txt")
gtemat <- as.matrix(gte[,-1, with=FALSE])
#heatmap(gtemat)

## WGCNA requires a data frame instead of a matrix
gtemat <- as.data.frame(gtemat)
gtemat <- t(gtemat)
#row.names(gtemat) <- gte$Id

## Preview the data
gtemat[1:5,1:5]
dim(gtemat)

## Verification that the sample data is good
gsg = goodSamplesGenes(gtemat, verbose = 3);
gsg$allOK

## Check data for outliers
sampleTree = hclust(dist(gtemat), method = "average");

#par(cex = 0.6);
#par(mar = c(0.3,0.3,0.3,0.3))
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)


### Automatic construction of the network
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(gtemat, powerVector = powers, verbose = 5)

### Plotting SFT
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90)

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

### Construction of the network, with scale parameter=6
net = blockwiseModules(gtemat, power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "gtex_tom",
                       maxBlockSize = 30000,
                       verbose = 3)

names(net)
mergedColors = labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

save(net, file='~/Projects/GelbRotation/GTEx-auto-network.Rdata')

propVarExplained(gtemat, net$colors, net$MEs)

ADJ1=abs(cor(gtemat,use="p"))^9
dissTOM=TOMdist(ADJ1)
