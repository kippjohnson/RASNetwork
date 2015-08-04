## Kipp Johnson
## kipp.johnson@icahn.mssm.edu

## Code to perform meta-analysis on p-values from
## Ben's cgea pipeline

rm(list=ls())

##########################################
## rJava options to facilitate xlsx input
##########################################
options(java.parameters = "-Xmx6000m");
jgc <- function(){ .jcall("java/lang/System", method = "gc") } ;
library(rJava);

########################################
## Load preliminary required packages
########################################
require(xlsx);
require(data.table);

########################################
## Read in required excel files from directory
########################################

combineXlsx <- function(filepath, filepattern="^[A-Za-z].*\\.xlsx"){
  oldwd <- getwd()
  setwd(filepath)

xlsxFiles <- dir(path=filepath,
                 pattern=filepattern) # Get all xlsx files in the directory

# print(xlsxFiles)
# }

outputdf <- data.frame(name=character(), # Initialize empty dataframe
                       compound=character(),
                       score=numeric(),
                       normScore=numeric(),
                       pvalue=numeric(),
                       adjP=numeric())

  for(file_i in 1:length(xlsxFiles)){
    name <- strsplit(xlsxFiles[file_i], split=".", fixed=TRUE)[[1]][1]
    sheet <- read.xlsx(xlsxFiles[file_i],sheetIndex=1, colIndex = 1:5)
    page <- cbind(name, sheet)
    outputdf <- rbind(outputdf, page)
    jgc() # garbage collection
  }
  return(outputdf)
  setwd(oldwd)
}

combinedSheets <- combineXlsx('/Users/kwj/Projects/GelbRotation/DrugScreen/cgea/xlsxFiles')

combinedSheets <- as.data.table(combinedSheets)
setkey(combinedSheets, name, score)

par(mfrow=c(2,4))
namelevels <- levels(combinedSheets$name)
for(i in 1:8){
  scores <- combinedSheets[name==namelevels[i]]$score
  hist(scores, main=paste(eval(namelevels[i])), cex=5)
}

namelevels <- levels(combinedSheets$name)
for(i in 1:8){
  scores <- combinedSheets[name==namelevels[i]]$pvalue
  hist(scores, main=paste(eval(namelevels[i])), cex=5)
}


combinedSheets <- combinedSheets[order(name, -score)]

topDrugs0 <- combinedSheets[, (.SD[score>0.30]), by=.(name)]
topDrugs <- combinedSheets[, (.SD[score>0.25]), by=.(name)]
topDrugs2 <- combinedSheets[, (.SD[score>0.20]), by=.(name)]


t0 <- table(topDrugs0$name)
t1 <- table(topDrugs$name)
t2 <- table(topDrugs2$name)

t012 <- cbind(t0,t1,t2)

kable(t012,
      col.names=c(
                  "Drugs (textgreater 0.30)",
                  "Drugs (textgreater 0.25)",
                  "Drugs (textgreater 0.20)"),
      format="latex")

t0 <- sort(table(as.character(as.factor(topDrugs0$compound))))
t1 <- sort(table(as.character(as.factor(topDrugs$compound))))
t2 <- sort(table(as.character(as.factor(topDrugs2$compound))))

x <- names(which(t0>=2))
for(i in 1:length(x)){ cat( c(x[i], "" ))}

x1 <- names(which(t1>=3))
for(i in 1:length(x1)){ cat( c(x1[i], ", " , sep=""))}

x2 <- names(which(t2>=5))
for(i in 1:length(x2)){ cat( c(x2[i], ", " , sep=""))}

x12 <- intersect(x1,x2)
x123 <- intersect(x12, x)
x123;

i12 <- union(x1,x2)
i123 <- union(i12, x)
i123;

match(tolower(i123), screenAll$Name)
length(which(!is.na(match(tolower(screenAll$Name), tolower(combinedSheets$compound)))))


# combinedSheetsDeg2 <- combineXlsx('/Users/kwj/Projects/GelbRotation/DrugScreen/cgea/xlsxFiles',
#                                   filepattern="^[a-z]*deg2.*\\.xlsx")
#
# combinedSheetsDeg3 <- combineXlsx('/Users/kwj/Projects/GelbRotation/DrugScreen/cgea/xlsxFiles',
#                                   filepattern="^[a-z]*deg3.*\\.xlsx")

######################################################
## Do p-value meta-analysis on combined excel sheets
######################################################

pvalMetaAnalysis <- function(combinedXlsxFile){
  sheets <- data.table(combinedXlsxFile)
  setkey(sheets, name, pvalue, compound)

  fisher.meta <- function(p) { ## Meta-analysis pvalue
    Xsq <- -2*sum(log(p))
    p.val <- pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
    return(c(meta.pvalue = p.val))
  }

  sheets[, fisher.meta.pval:=fisher.meta(adjP), by=.(compound)]
  sheets[,fisher.meta.pval.adj:=p.adjust(fisher.meta.pval, method="fdr"), by=.(name)]
  return(sheets)
}

sheets <- pvalMetaAnalysis(combinedSheets)
# sheetsDeg2 <- pvalMetaAnalysis(combinedSheetsDeg2)
# sheetsDeg3 <- pvalMetaAnalysis(combinedSheetsDeg3)

sheets[name==.(delimited_pan_intestine_deg2)]
length(which(sheets[name==.(delimited_pan_intestine_deg2), .(pvalue)]<0.000001))

######################################################
## Fixed effects model for score
######################################################

# http://www.meta-analysis.com/downloads/Meta-analysis%20fixed%20effect%20vs%20random%20effects.pdf

scoreFixedEffects <- function(combinedXlsxFile){
  sheets <- data.table(combinedXlsxFile)
  sheets[, inv.score.var:=1/var(score), by=.(name)] #assign weights by study
  sheets[, weighted.mean.score:=sum(score*inv.score.var)/sum(inv.score.var), by=.(compound)] #get weighted mean score
  sheets[, combined.effect.var:=1/sum(inv.score.var), by=.(compound)]
  sheets$combined.effect.se <- sqrt(sheets$combined.effect.var)
  sheets$z <- sheets$weighted.mean.score/sheets$combined.effect.se
  sheets$fixed_p <- 2*(1-pnorm(abs(sheets$z)))

  return(sheets)
  }

sheets <- scoreFixedEffects(combinedSheets)
hist(sheets[name==sheets$name[1],.(compound, fixed_p)]$fixed_p)

######################################################
## Random effects model for score
######################################################

# http://www.meta-analysis.com/downloads/Meta-analysis%20fixed%20effect%20vs%20random%20effects.pdf
#
# scoreRandomEffects <- function(combinedXlsxFile){
#   sheets <- data.table(combinedXlsxFile)
#   sheets <- data.table(combinedSheets)
#   sheets[, inv.score.var:=1/var(score), by=.(name)] #assign weights by study
#   sheets[, meanScore:=mean(score), by=.(compound)] # compute the mean score per compound
#   sheets[, Q:=sum(inv.score.var * (score - meanScore)^2), by=.(compound)] # compute Q statistic
#   sheets$df <- length(unique(sheets$name))-1
#   sheets[, C:=sum(inv.score.var)-sum(inv.score.var^2/inv.score.var), by=.(compound)]
#   sheets$Tsquare <- ifelse(sheets$Q>sheets$df, sheets$Q-sheets$df  ,0)
#   return(sheets)
# }

# sheets <- scoreFixedEffects(combinedSheets)

library(lme4);
score.lmer1 <- lmer(score ~ 1 + (1|name) + compound, data=combinedSheets)
out1 <- summary(score.lmer1)
p.val <- 1-pt(out1$coefficients[,3], df=Inf)
out1_coefficients <- cbind(out1$coefficients, p.val)
hist(out1_coefficients[,4])

length(which(out1_coefficients[,4]<0.01))
length(which(out1_coefficients[,4]<0.000001))


score.glm <- glm(score ~ name + compound, data=combinedSheets)
out2 <- summary(score.glm)
hist(out2$coefficients[,4])

######################################################
## Combine excel sheet #2
######################################################

combineXlsxSheet2 <- function(filepath, filepattern="^[A-Za-z].*\\.xlsx"){
  oldwd <- getwd()
  setwd(filepath)

  xlsxFiles <- dir(path=filepath,
                   pattern=filepattern) # Get all xlsx files in the directory

  # print(xlsxFiles)
  # }

  outputdf <- data.frame(name=character(), # Initialize empty dataframe
                         compound=character(),
                         score=numeric(),
                         normScore=numeric(),
                         pvalue=numeric(),
                         adjP=numeric())

  for(file_i in 1:length(xlsxFiles)){
    name <- strsplit(xlsxFiles[file_i], split=".", fixed=TRUE)[[1]][1]
    sheet <- read.xlsx(xlsxFiles[file_i],sheetIndex=2, colIndex = 1:5)
    page <- cbind(name, sheet)
    outputdf <- rbind(outputdf, page)
    jgc() # garbage collection
  }
  return(outputdf)
  setwd(oldwd)
}

combinedSheet2 <- combineXlsxSheet2('/Users/kwj/Projects/GelbRotation/DrugScreen/cgea/xlsxFiles')
combinedSheet2 <- data.table(combinedSheet2)
setkey(combinedSheet2, name, feature)

sigFeatures <- combinedSheet2[adjPvalue<0.10,.(name, feature)]
sort(table(as.factor(as.character(sigFeatures$feature))),decreasing = TRUE)

kable(data.frame(table(sigFeatures$name)), col.names=c('Network', 'Number of Drug targets'), format='latex')

par(mfrow=c(2,4))
namelevels <- levels(combinedSheet2$name)
for(i in 1:8){
  scores <- combinedSheet2[name==namelevels[i]]$score
  hist(scores, main=paste(eval(namelevels[i])), cex=5)
}

namelevels <- levels(combinedSheet2$name)
for(i in 1:8){
  scores <- combinedSheet2[name==namelevels[i]]$adjPvalue
  hist(scores, main=paste(eval(namelevels[i])), cex=5,xlab="Adjusted P",xlim=c(0,0.5))
}
