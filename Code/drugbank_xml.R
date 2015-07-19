library(XML)

rm(list=ls())

db <- xmlTreeParse("~/Projects/GelbRotation/DrugScreen/drugbank.xml")
dbtop = xmlRoot(db)
db2 <- xmlChildren(dbtop)

xmlName(db2)
#names(dbtop) <- paste("drug",1:length(names(dbtop)), sep="_")

dbdf <- xmlSApply(dbtop, function(x) xmlSApply(x, xmlValue))
dbdf <- xmlSApply(db2, function(x) xmlSApply(x, xmlValue))

