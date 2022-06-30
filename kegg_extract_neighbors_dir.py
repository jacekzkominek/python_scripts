#!/usr/bin/env python

import os, sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description='Compare two kegg annotations')
parser.add_argument("dir", help="KEGG annotation file.")
parser.add_argument("--kegg", help="KEGG marker.")
parser.add_argument("--n", default="10", help="Neighbors.")
args = parser.parse_args()

for d in sorted(os.listdir(args.dir)):
	if d.find(".ko") != -1:
		print d
		redir = open(args.dir+"/"+d.replace(".ko",".out"),"w")
		subprocess.call(["kegg_extract_neighbors.py",args.dir+"/"+d,"--kegg",args.kegg,"--n",args.n], stdout=redir)
		redir.close()
