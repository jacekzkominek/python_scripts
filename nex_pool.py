#!/usr/bin/env python

from Bio import Phylo
import sys

trees = []
for t in sys.argv[1:]:
	for tr in Phylo.parse(t, "nexus"):
		trees.append(tr)
	#trees.append(Phylo.read(t, "nexus"))

Phylo.write(trees, "combined.nex", "nexus")
