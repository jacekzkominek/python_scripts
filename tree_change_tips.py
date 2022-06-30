#!/usr/bin/env python

import sys, argparse, os, subprocess, shutil
from collections import defaultdict
from ete3 import Tree

parser = argparse.ArgumentParser(description="")
parser.add_argument("tree", default="", help="Tree")
parser.add_argument("--postfix", default="", help="")
parser.add_argument("--prefix", default="", help="")
parser.add_argument("--root", default="", help="")
args = parser.parse_args()

t1 = Tree(args.tree, format=1)
for node in t1.traverse("postorder"):
	if node.is_leaf() == True:
		if args.root != "" and node.name == args.root:
			t1.set_outgroup(node)
		node.name = args.prefix+node.name+args.postfix
		
t1.write(format=0, outfile=args.tree+"_new.nwk")


