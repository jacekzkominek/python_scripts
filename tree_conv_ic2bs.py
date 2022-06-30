#!/usr/bin/env python

import sys, os, math, argparse
from collections import defaultdict
import re 
from ete3 import *

parser = argparse.ArgumentParser(description="Convert IC scores into BS.")
parser.add_argument('tree', default="", help="Input tree file (newick format)")
parser.add_argument('--ica', default=False, action="store_true", help="Print out ICA instead of IC scores")
parser.add_argument('--ref', default="", help="Ref tree")
parser.add_argument('--root', default="", help="Ref root")
args = parser.parse_args()


with open(args.tree) as t:
	for l in t:
		t_out = ""
		chars2 = ",)["
		bl = False
		#STRIP BRANCH LENGTHS
		for c in l:
			if c == ":":
				bl = True
			elif bl == True and c in chars2:
				t_out += c
				bl = False
			elif bl == False:
				t_out += c
		
		#FIX IC SCORES AS BRANCH LABELS
		if args.ica == False:
			t_out2 = re.sub(r"\[([0-9\.\-]+)\,([0-9\.\-]+)\]", r"\1",t_out)
		else:
			t_out2 = re.sub(r"\[([0-9\.\-]+)\,([0-9\.\-]+)\]", r"\2",t_out)
		
		if args.ref != "" and args.root != "":
			t_ic = Tree(t_out2, format=0)
			t_ic.set_outgroup(args.root)
			t_ref = Tree(args.ref, format=1)
			t_ref.set_outgroup(args.root)
			node_ic = {}
			for node in t_ic.traverse("postorder"):
				node_ic["_".join(node.get_leaf_names())] = node.support
			for node in t_ref.traverse("postorder"):
				if "_".join(node.get_leaf_names()) in node_ic:
					node.support = node_ic["_".join(node.get_leaf_names())]
			print t_ref.write(format=0)
			print t_out2
		else:
			t_ic = Tree(t_out2, format=0)
			t_ic.set_outgroup(args.root)
			print t_ic.write(format=0, outfile="outtree1.nwk")
		
		
