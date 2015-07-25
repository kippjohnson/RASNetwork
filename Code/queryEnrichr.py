################################################################################
## Kipp Johnson
## kipp.johnson@icahn.mssm.edu
## July 24th, 2015
################################################################################

### HTTP Post code taken from RNA-seq analysis script:
### http://amp.pharm.mssm.edu/Enrichr/#help
### See bottom of page

import cookielib, urllib2, urllib
import poster
import sys
import os
import time
import argparse

## Python script to automatically query EnrichR database and retrieve result tables
## usage: python queryEnrichr.py -g genelist -d database
##        genelist should be a list of Gene symbols separated by newlines (1 per line)
##        database should be one of the EnrichR databases
##        --> I've only typed in the database names for those which are printed with the
##            -h or --help option thus far
##
##        Resulting tables are printed to stdout
##
##

### Recommend that you use this via a bash loop like the folllowing:
# for f in */*.nodes;
# do
#   echo "doing $f..."
#   python ../../Code/queryEnrichr.py -g $f -d KEGG_2015 > KEGG_2015_table.txt;
#   python ../../Code/queryEnrichr.py -g $f -d GO_Biological_Process > GO_biological_process.txt;
#   python ../../Code/queryEnrichr.py -g $f -d MSigDB_Computational > MSigDB_computational.txt ;
#   python ../../Code/queryEnrichr.py -g $f -d Drug_Perturbations_From_GEO > Drug_Perturbations_From_GEO.txt;
#   python ../../Code/queryEnrichr.py -g $f -d CMAP_up > CMAP_up.txt;
#   python ../../Code/queryEnrichr.py -g $f -d CMAP_down > CMAP_down.txt;
#   echo "...done"
# done;

parser = argparse.ArgumentParser(description="""

     Query the EnrichR database
     Usage: python queryEnrichr.py -g genelist -d database

     Currently supported database arguments:

        KEGG_2015
        WikiPathways_2015
        Reactome_2015
        BioCarta_2015
        PPI_Hub_Proteins
        KEA
        Panther
        GO_Biological_Process
        Go_Cellular_Component
        Go_Molecular_Function
        Human_Phenotype_Ontology
        CMAP_up
        CMAP_down
        GeneSigDB
        OMIM_Disease
        OMIM_Expanded
        MSigDB_Computational
        MSigDB_Oncogenic_Signatures
        Drug_Perturbations_From_GEO
        Disease_Signatures_from_GEO_up
        Disease_Signatures_from_GEO_down

    Results printed to STDOUT
     """,formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--genelist', '-g', required=True, help='Location of input genelist')
parser.add_argument('--database', '-d', required=True, help='Name of EnrichR database to query')
args = parser.parse_args()

genes = open(args.genelist, "r") # list of genes, 1 per line
db = str(args.database)

possible_dbs = ["KEGG_2015",
                "WikiPathways_2015",
                "Reactome_2015",
                "BioCarta_2015",
                "PPI_Hub_Proteins",
                "KEA",
                "Panther",
                "GO_Biological_Process",
                "Go_Cellular_Component",
                "Go_Molecular_Function",
                "Human_Phenotype_Ontology",
                "CMAP_up",
                'CMAP_down',
                "GeneSigDB",
                "OMIM_Disease",
                "OMIM_Expanded",
                "MSigDB_Computational",
                "MSigDB_Oncogenic_Signatures",
                "Drug_Perturbations_From_GEO",
                "Disease_Signatures_from_GEO_up",
                "Disease_Signatures_from_GEO_down"]

if db not in possible_dbs:
    print "\nList of Queryable databases:"
    for db in possible_dbs:
        print "\t", db
    print "\n",
    raise ValueError('Input DB is not in the list of queryable databases:')

geneString = ''
for g in genes:
    geneString = geneString+str(g)

def Enrichr(geneslist,optionalinfo='',querydb=db):
    #post a gene list to enrichr server and get the link.
    opener = poster.streaminghttp.register_openers()
    opener.add_handler(urllib2.HTTPCookieProcessor(cookielib.CookieJar()))
    params = {'list':geneslist,'description':optionalinfo}
    datagen, headers = poster.encode.multipart_encode(params)
    url = "http://amp.pharm.mssm.edu/Enrichr/enrich"
    request = urllib2.Request(url, datagen,headers)
    urllib2.urlopen(request)
    time.sleep(1) # absolutely critical for some reason...it will fail without it!
    if db == 'KEGG_2015':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=KEGG_2015_table&backgroundType=KEGG_2015').read()
    elif db =='WikiPathways_2015':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=WikiPathways_2015_table&backgroundType=WikiPathways_2015').read()
    elif db =='Reactome_2015':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=Reactome_2015.txt&backgroundType=Reactome_2015').read()
    elif db =='BioCarta_2015':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=BioCarta_2015.txt&backgroundType=BioCarta_2015').read()
    elif db =='PPI_Hub_Proteins':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=PPI_Hub_Proteins.txt&backgroundType=PPI_Hub_Proteins').read()
    elif db =='KEA':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=KEA.txt&backgroundType=KEA').read()
    elif db =='Panther':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=Panther.txt&backgroundType=Panther').read()
    elif db =='GO_Biological_Process':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=GO_Biological_Process.txt&backgroundType=GO_Biological_Process').read()
    elif db =='Go_Cellular_Component':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=Go_Cellular_Component.txt&backgroundType=Go_Cellular_Component').read()
    elif db =='Go_Molecular_Function':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=Go_Molecular_Function.txt&backgroundType=Go_Molecular_Function').read()
    elif db =='Human_Phenotype_Ontology':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=Human_Phenotype_Ontology.txt&backgroundType=Human_Phenotype_Ontology').read()
    elif db =='CMAP_up':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=CMAP_up.txt&backgroundType=CMAP_up').read()
    elif db =='CMAP_down':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=CMAP_down.txt&backgroundType=CMAP_down').read()
    elif db =='OMIM_Disease':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=OMIM_Disease.txt&backgroundType=OMIM_Disease').read()
    elif db =='OMIM_Expanded':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=OMIM_Expanded.txt&backgroundType=OMIM_Expanded').read()
    elif db =='MSigDB_Computational':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=MSigDB_Computational.txt&backgroundType=MSigDB_Computational').read()
    elif db =='Drug_Perturbations_From_GEO':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=Drug_Perturbations_From_GEO.txt&backgroundType=Drug_Perturbations_From_GEO').read()
    elif db =='Disease_Signatures_from_GEO_up':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=Disease_Signatures_from_GEO_up.txt&backgroundType=Disease_Signatures_from_GEO_up').read()
    elif db =='Disease_Signatures_from_GEO_down':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=Disease_Signatures_from_GEO_down.txt&backgroundType=Disease_Signatures_from_GEO_down').read()
    elif db =='GeneSigDB':
        y=urllib2.urlopen('http://amp.pharm.mssm.edu/Enrichr/export?filename=GeneSigDB.txt&backgroundType=GeneSigDB').read()
    else:
        y='Something happened here'
        raise ValueError('db not in list, should not be here...')
    return y

returned_table = Enrichr(geneString)
print returned_table # could easily have it write the table instead here
