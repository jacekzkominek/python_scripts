#!/usr/bin/env python

from ete2 import Tree

t = Tree("/home/jacek/RAxML_bipartitions.Cterm_B100", format=0)
for node in t.traverse("postorder"):
  # Do some analysis on node
  node.support = float(node.support)/100


t.write(outfile="/home/jacek/RAxML_bipartitions.Cterm_B100_b", format=2)
