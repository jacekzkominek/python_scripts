#!/usr/bin/env python

from Bio import Phylo
import sys, os, math, argparse

parser = argparse.ArgumentParser(description="Convert Internode Certainty into BS values.")
parser.add_argument('tree', default="", help="Input tree file (newick format)")
args = parser.parse_args()

filename = args.tree
tree = Phylo.read(filename, "newick")

nodes = tree.get_nonterminals()
bs_total = 0
bs_count = 0
for node in nodes:
	if (node.confidence and node.confidence != None):
		bs_total += node.confidence
		bs_count += 1
		
print "\n","Average bootstrap per branch:",float(bs_total/bs_count),"\n"

