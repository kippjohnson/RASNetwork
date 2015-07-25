###  Kipp Johnson
###  kipp.johnson@icahn.mssm.edu
###
###  Usage: python -n network -g geneset -l nlayer
###  Run with -h/--help for additional information
###
###  Input Network should be of form:
###  source target
###  RAS MAPK
###
###  (i.e. RAS\sMAPK)
###  Do not include source/target row; code doesn't deal with headers currently
###
###  Input gene set of interest sould be a file with just 1 gene symbol per line
###
###  Could easily be modified to only expand upstream or downstream
###  --> just don't call find_upstream_nodes() or find_downstream_nodes()
###      in get_layer()

import argparse
import re

parser = argparse.ArgumentParser('Finds subnetwork within layer vicinity of geneset')
parser.add_argument('--network', '-n', required=True, help='File path to input network')
parser.add_argument('--geneset', '-g', required=True, help='File path to geneset to expand around')
parser.add_argument('--layer', '-l', required=True, help='Layers to expand subnetwork around geneset')
args = parser.parse_args()

##
inputnetworkfile = open(args.network)
inputgenes = open(args.geneset)
nLayers = int(args.layer)

if nLayers < 1:
    raise ValueError('Number of layers must be greater than or equal to 1')

## Read the genes into list called geneset
def parse_geneset(inputgeneset):
    geneset = []
    for gene in inputgeneset:
        geneset.append(gene.rstrip("\n "))
    return geneset

def parse_networkfile(inputnetwork):
    networkfile = []
    for line in inputnetwork:
        networkfile.append(line.rstrip("\n"))
    return networkfile

def find_upstream_nodes(gene, network):
    upstreamNode = []
    for line in network:
        line = re.split("\s", line)
        if line[1] == gene:
            upstreamNode.append(str(line[0])+"\t"+str(line[1]))
    return upstreamNode

def find_downstream_nodes(gene, network):
    downstreamNode = []
    for line in network:
        line = re.split("\s", line)
        if line[0] == gene:
            downstreamNode.append(str(line[0])+'\t'+str(line[1]))
    return downstreamNode

def get_layer(inGene, inNetwork):
    layer_interaction = []

    upstreamNodes=find_upstream_nodes(gene=inGene, network=inNetwork)
    for item in upstreamNodes:
        if len(item) > 0:
            layer_interaction.append(item)

    downstreamNodes=find_downstream_nodes(gene=inGene, network=inNetwork)
    for item in downstreamNodes:
        if len(item) > 0:
            layer_interaction.append(item)

    return layer_interaction

def runlayer(inputgenes=inputgenes, inputnetworkfile=inputnetworkfile):
    networkLayer = []
    for gene in geneset:
        genelayer = get_layer(inGene=gene, inNetwork=networkfile)
        if len(genelayer) > 0:
            networkLayer.append(genelayer)

    interactionLayer = []
    for layer in networkLayer:
        for interaction in layer:
            interactionLayer.append(interaction)
    return interactionLayer

def collateGenes(inputlist):
    collatedGenes = []
    for sourcetarget in inputlist:
        sourcetarget = re.split("\t", sourcetarget)
        collatedGenes.append(sourcetarget[0])
        collatedGenes.append(sourcetarget[1])
    return collatedGenes

if __name__== "__main__":
    # main()
    geneset = parse_geneset(inputgenes)
    networkfile = parse_networkfile(inputnetworkfile)

    for i in range(nLayers):
        outgenes = runlayer(inputgenes=geneset, inputnetworkfile=networkfile)
        geneset = collateGenes(outgenes)
        geneset = list(set(geneset))

    for item in outgenes:
        print item
