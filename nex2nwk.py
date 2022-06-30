#!/usr/bin/env python

import sys, os, argparse
from ete2 import Tree
from collections import defaultdict
from natsort import natsorted
import dendropy

parser = argparse.ArgumentParser(description="")
parser.add_argument("infile", help="Query FASTA")
args = parser.parse_args()

post_trees = dendropy.TreeList()
post_trees.read(file=open(args.infile,"r"), schema="nexus")
post_trees.write(path=args.infile+".nwk",schema="newick")


#~ tree = ""
#~ with open(args.infile) as tf:
	#~ tipmap = False
	#~ tip = defaultdict(str)
	#~ for l in tf:
		#~ l = l.strip()
		#~ #print l
		#~ if tipmap == True:
			#~ tip[l.split()[0]] = l.split()[1][:-1]
			#~ #print l
		#~ if l.find("translate") != -1:
			#~ tipmap = True
		#~ if tipmap  == True and l.find(";") != -1:
			#~ #print l
			#~ break
	#~ print tip
	
	#~ for l in tf:
		#~ l = l.strip()
		
		#~ if l.find("(") != -1:
			#~ l = l[l.find("("):]
			#~ for t in natsorted(tip.keys(), reverse=True):
				#~ #print t
				#~ l = l.replace(t+":",tip[t]+":")
			#~ tree = l
			#~ t = Tree(l, format=0)
			#~ t.write(outfile=args.infile+".nwk", format=5)
