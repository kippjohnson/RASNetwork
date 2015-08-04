import sys
import re

nodeStore = []

in_network = open(sys.argv[1], "r")

for line in in_network:
    line = re.split("\s", line)
    nodeStore.append(line[0])
    nodeStore.append(line[2])

nodeStore = set(nodeStore)

for node in nodeStore:
    print node
