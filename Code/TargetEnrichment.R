### Kipp Johnson
### kipp.johnson@icahn.mssm.edu

### Required packages
rm(list=ls())
library(data.table)
library(knitr) # kable() function

### Read in the required datasets
dname <-  read.csv("~/Projects/GelbRotation/DrugScreen/drug_name.csv")
targets <- read.csv("~/Projects/GelbRotation/DrugScreen/drug_target_uniprot_links.csv")
  #screen0 <- read.csv("~/Projects/GelbRotation/DrugScreen/FDA_Celine_070715.csv")
  screen0 <- read.csv("~/Projects/GelbRotation/DrugScreen/FDA_Celine_top5.csv")
  #screen0 <- read.csv("~/Projects/GelbRotation/DrugScreen/FDA_Celine_gt5.csv")
  #screen0 <- read.csv("~/Projects/GelbRotation/DrugScreen/FDA_Celine_gt4.csv")

### Merge datasets, remove unnecessary columns
dname$IUPAC.name <- NULL
dname$SMILES <- NULL
screen <- merge(dname, screen0, by.x='Number', by.y='Drug', all.x=TRUE)
screen$Sylvain <- NULL

# Create a vector 1 if drug was positive in screen, 0 if drug was negative in screen
screen$positive <- ifelse(is.na(screen$Adults), 0, 1)

### Take first word from drug name
# For screen:
medpattern <- '^.\\w*' # Matches the first word in the string
med_results <- regexpr(pattern=medpattern, screen$Name, perl=TRUE) # compute the reg. expressions
mednames <- rep(NA, length(med_results)) # create an empty vector of appropriate length to store matches
screen$Name1w <- regmatches(screen$Name, med_results) # store the matched names in vector

### Manual RegEXP fixes
screen$Name1w[45] <- "Timolol"
screen$Name1w[70] <- "Triprolidine"
screen$Name1w[85] <- "Isoproterenol"
screen$Name1w[92] <- "Tubocuraraine"
screen$Name1w[93] <- "Butaclamol"
screen$Name1w[94] <- "Apomorphine"
screen$Name1w[99] <- "Raclopride"
screen$Name1w[101] <- "Sulpiride"
screen$Name1w[108] <- "Epinephrine"
screen$Name1w[110] <- "Norepinephrine"
screen$Name1w[119] <- "Vitamin A"
screen$Name1w[174] <- "Deprenyl"
screen$Name1w[270] <- "Hydroxycamptothecin"
screen$Name1w[275] <- "Methylsalicylate" # Doesn't match anyway
screen$Name1w[355] <- "Deoxythymidine"
screen$Name1w[362] <- "Aminoglutethimide"
screen$Name1w[363] <- "Aminosalicylic"
screen$Name1w[364] <- "Aminosalicylic"
screen$Name1w[415] <- "Dideoxycytidine"
screen$Name1w[441] <- "Fluorouracil"
screen$Name1w[463] <- "Hydroxyprogesterone"
screen$Name1w[470] <- "Ketoprofen"
screen$Name1w[638] <- "Levothyroxine"

### For target:
med_results <- regexpr(pattern=medpattern, targets$Name, perl=TRUE) # compute the reg. expressions
mednames <- rep(NA, length(med_results)) # create an empty vector of appropriate length to store matches
targets$Name1w <- regmatches(targets$Name, med_results) # store the matched names in vector

# Manual fix
targets[targets$Name=='Vitamin A',]$Name1w <- 'Vitamin A'

### Get a list of all matching drug names in our dataset
all_matches <- vector()
for(i in seq_along(screen$Name1w)){
  best_match <- agrep(screen$Name[i], targets$Name, value=TRUE, max.distance=0.10)[1]

  if(is.na(best_match)){
    best_match <- agrep(screen$Name1w[i], targets$Name1w, value=TRUE, max.distance=0.10)[1]
  }

  #cat(c(as.character(screen$Name)[i], "\t\t", screen$Name1w[i], "\t",best_match, "\n"))
  all_matches[i] <- best_match
}

### How many drugs (in our list of 640) did not match?
length(which(is.na(all_matches)))
screen$matchname <- all_matches

screen <- screen[which(!is.na(screen$matchname)),]

### Merge targets and
targ <- merge(screen, targets, by.x="matchname", by.y="Name")
targ <- merge(screen, targets, by.x="matchname", by.y="Name1w")

targ <- as.data.table(targ)
setkey(targ, matchname, UniProt.Name)
targ <- unique(targ)

target_matches <- subset(targ, positive==1)
target_matches$UPName <- as.character(target_matches$UniProt.Name)

### Sub-table of non-matched drugs in our dataset
target_nonmatches <- subset(targ, positive==0)
target_nonmatches$UPName <- as.character(target_nonmatches$UniProt.Name)

### Construct the following 2x2 table for enrichment
#############################################################
#                                                           #
#                  Drug screen +  Drug screen -             #
#        Target +       a              b                    #
#        Target -       c              d                    #
#                                                           #
#                    OR = (a/c)/(b/d)                       #
#############################################################

### List of match frequencies
tab1 <- sort(table(target_matches$UPName), decreasing=TRUE)
tab1 <- as.data.frame(tab1)
names(tab1) <- c('A')

### Compute box B
tab2 <- sort(table(target_nonmatches$UPName), decreasing=TRUE)
tab2 <- as.data.frame(tab2)
names(tab2) <- c('B')

### Merge the two tables and re-ordeer columns
tab1 <- cbind(tab1, targets=row.names(tab1));
tab2 <- cbind(tab2, targets=row.names(tab2));
tab <- merge(tab1,tab2, all=TRUE);
#tab <- tab[c('targets','A','B','C','D')];

tab$A <- ifelse(is.na(tab$A), 0, tab$A)
tab$B <- ifelse(is.na(tab$B), 0, tab$B)

### nScreen Matches
nMatches <- length(unique(target_matches$matchname))
tab$C <- nMatches - tab$A

### Compute box D
nNonMatches <- length(unique(target_nonmatches$matchname))
tab$D <- nNonMatches - tab$B

### Verify all cells sum to same amount
tab$AC <- tab$A+tab$C;
tab$BD <- tab$B+tab$D;
tab$AB <- tab$A+tab$B;
tab$CD <- tab$C+tab$D;
tab$N <- tab$A+tab$B+tab$C+tab$D

#####
##### Calculate the ORs for every target effect
#####
for(i in 1:dim(tab)[1]){ # for every target in our table of targets:

  # Here, we are constructing the 2x2 matrix for target i. The A,B,C,D columns were calculated earlier
  targetmat <- matrix(data = c(tab$A[i], #
                             tab$B[i], #
                             tab$C[i], #
                             tab$D[i]) #
                    ,nrow=2, byrow=TRUE) # fill A then B then C then D (default = A..C..B..D)

  #print(i) ;print(enrichmat) # progress printed in console

  # perform a two-way fisher test on the 2x2 table
  fishertest <- fisher.test(targetmat, or=1, alternative = "two.sided", conf.int=TRUE, B=10000)

  # extract the values from our fisher test result
  tab$OR[i] <- fishertest$estimate
  tab$pval[i] <- fishertest$p.value
  tab$LCI[i] <- fishertest$conf.int[[1]] # lower confidence interval
  tab$UCI[i] <- fishertest$conf.int[[2]] # upper confidence interval
}

tab$FDR.pval <- p.adjust(tab$pval, method="BH")
tab <- as.data.table(tab)

tab[,Sig:=ifelse(pval<0.05,1,0)] # is the effect significant at p<0.05?
tab[,FDR.sig:=ifelse(FDR.pval<0.05,1,0)] # is the effect significant at p<0.05?

tab <- tab[order(-FDR.sig, pval, FDR.pval)]
kable(tab[Sig==1])


# ###
# ### Poisson regression
# library(glmnet)
#
# targ$Adults0 <- ifelse(is.na(targ$Adults),0,targ$Adults)
# y <- targ$Adults0
#
# uniquetargets <- unique(targ$UniProt.ID)
#
# x <- list()
# adultvec <- unique(targ$Adults0)
#
# for(i in seq_along(unique(targ$Adults0))){
#   x[[i]] <- as.character(targ[Adults0==adultvec[i]]$UniProt.ID)
# }
#
# for(i in seq_along(uniquetargets)){
#
#   targ[, eval(paste0('exp_',uniquetargets[i])):=ifelse(uniquetargets[i] %in% sapply(x[Adults0], "["),1,0), by=.(Adults0)]
# }
#
#
# # for(i in seq_along(SigDrugs)){
# #   admit.drug.ml[,eval(paste0("Exp_",SigDrugs[i])):=ifelse(SigDrugs[i] %in% med_name_lc, 1, 0), by=.(mrn)]
# # }
#
# setkey(targ,DrugBank.ID)
# targ <- unique(targ)
#
# targ <- data.frame(targ)
# exposures <- targ[,14:590]
#
# cv.fit <- cv.glmnet(x=as.matrix(exposures),y=targ$Adults0, family="poisson")
# coef(cv.fit, s = "lambda.min")
#
# which(coef(cv.fit, s = "lambda.min")[,1]>0)
#
# mod_coef <- names(which(coef(cv.fit, s = "lambda.min")[,1]>0))
# mod_coef <-
#   strsplit(mod_coef, "_")[,2]

kable(tab[Sig==1, .(targets,A,B,C,D,N,OR,pval,FDR.pval,Sig,FDR.sig)])

tabc <- copy(tab)
tabc$Sig <- ifelse(tabc$Sig==1, tabc$Sig<-"Y",tabc$Sig<-"N")
tabc$FDR.sig <- ifelse(tabc$FDR.sig==1, tabc$FDR.sig<-"Y",tabc$FDR.sig<-"N")

xtable(tabc[Sig=="Y", .(targets,A,N,OR,pval,FDR.pval,Sig,FDR.sig)],digits=c(0,0,0,0,0,2,2,1,1))



















