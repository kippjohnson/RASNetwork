#!/usr/local/bin/python

import sys
import re

in_network = open(sys.argv[1], "r")
out_nodes = open(sys.argv[1].rstrip(".txt")+'.nodes.txt', "w")
genes = open("/Users/kwj/Projects/GelbRotation/GeneLists/rasopathy_genes.txt")

i=0
gene_list = [None]*15

for entry in genes:
    gene_list[i] = str(entry).rstrip(" \n")
    i += 1

netgenes = []
for line in in_network:
    line = line.rstrip() #strip trailing whitespace

    linesplit = re.split("\s", line)
    # Print len(linesplit), "\t", linesplit[0:5]
    if len(linesplit)==3:
        netgenes.append(linesplit[0])
        netgenes.append(linesplit[2])

    # Columns: source, target
    if len(linesplit)<=2:
        netgenes.append(linesplit[0])
        netgenes.append(linesplit[1])

netgenes = list(set(netgenes))

out_nodes.write("Gene\tStatus\n")
for i in range(0, len(netgenes)):
    if netgenes[i] in gene_list:
        out_nodes.write(str(netgenes[i])+"\tRASopathy_gene\n")
    else:
        out_nodes.write(str(netgenes[i])+"\tnon_indicated\n")

in_network.close()
out_nodes.close()
