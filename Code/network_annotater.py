#!/usr/local/bin/python

import sys
import re

in_network = open(sys.argv[1], "r")
out_network = open(sys.argv[1].rstrip(".txt")+'.anno.txt', "w")
genes = open("/Users/kwj/Projects/GelbRotation/GeneLists/rasopathy_genes.txt")

i=0
gene_list = [None]*15

for entry in genes:
    gene_list[i] = str(entry).rstrip(" \n")
    i += 1

#print(gene_list)

for line in in_network:
    line = line.rstrip() #strip trailing whitespace

    if line.startswith("source"): # header lines
        out_network.write(str(line)+"\t"+"status"+"\n")

    else:
        linesplit = re.split("\s", line)
        # Columns: source, interaction type, target
        # print len(linesplit), "\t", linesplit[0:5]
        if len(linesplit)==3:
            if linesplit[0] in gene_list:
                out_network.write(line.rstrip("\n")+'\t'+"RASopathy_gene"+"\n")
            elif linesplit[1] in gene_list:
                out_network.write(line.rstrip("\n")+'\t'+"RASopathy_gene"+"\n")
            elif linesplit[2] in gene_list:
                out_network.write(line.rstrip("\n")+'\t'+"RASopathy_gene"+"\n")
            else:
                out_network.write(line.rstrip("\n")+"\t"+"non_indicated"+"\n")

        # Columns: source, target
        if len(linesplit)<=2:
            if linesplit[0] in gene_list:
                out_network.write(line.rstrip("\n")+'\t'+"RASopathy_gene"+"\n")
            elif linesplit[1] in gene_list:
                out_network.write(line.rstrip("\n")+'\t'+"RASopathy_gene"+"\n")
            else:
                out_network.write(line.rstrip("\n")+"\t"+"non_indicated"+"\n")

in_network.close()
out_network.close()
