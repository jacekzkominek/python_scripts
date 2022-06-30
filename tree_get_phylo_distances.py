#!/usr/bin/env python

import sys, argparse, os, subprocess, shutil
from collections import defaultdict
import re 
from ete3 import *

parser = argparse.ArgumentParser(description='')
parser.add_argument("--sp_tree", default="", help="Species tree")
parser.add_argument("--gene_tree", default="", help="Gene tree (or dir with trees)")
parser.add_argument("--name1", default="", help="Taxon name (\"all\" for all taxa")
#~ parser.add_argument("name2", default="", help="Input FASTA file")
args = parser.parse_args()

sp_dist = defaultdict(lambda: defaultdict(float))

files = []
if os.path.exists(args.gene_tree):
	if os.path.isdir(args.gene_tree) == True:
		files = os.listdir(args.gene_tree)
	else:
		files = args.gene_tree

t1 = Tree(args.sp_tree, format=1)
targets  = []
for node1 in t1.traverse("postorder"):
	if len(node1) == 1:
		n1 = str(node1)[3:]
		targets.append(n1)
		if n1.find(args.name1) == -1:
			sp_dist[args.sp_tree][n1] = t1.get_distance(args.name1,n1)
		else:
			sp_dist[args.sp_tree][n1] = 0

print args.name1+"\t",
for target in sorted(sp_dist[args.sp_tree].keys()):
	print target+"\t",
print ""
print args.sp_tree+"\t",
for target in sorted(targets):
	print str(sp_dist[args.sp_tree][target])+"\t",

sys.stdout.flush()
for f in sorted(files):
	t = Tree(args.gene_tree+"/"+f, format=1)	
	if str(t).find(args.name1) != -1:
		for node1 in t.traverse("postorder"):
			if len(node1) == 1:
				n1 = str(node1)[3:]
				if n1.find(args.name1) == -1:
					sp_dist[f][n1] = t.get_distance(args.name1,n1)
				else:
					sp_dist[f][n1] = 0

		print ""
		print f+"\t",
		sys.stdout.flush()
		for target in sorted(targets):
			print str(sp_dist[f][target])+"\t",
		sys.stdout.flush()
