#!/usr/bin/env python

import os, sys
import argparse
import subprocess
from collections import defaultdict

parser = argparse.ArgumentParser(description='Compare two kegg annotations')
parser.add_argument("dir", help="KEGG annotation file.")
parser.add_argument("--kegg", help="KEGG marker.")
args = parser.parse_args()

keggs = defaultdict(int)
for d in sorted(os.listdir(args.dir)):
	if d.find(".out") != -1:
		#print d
		with open(args.dir+"/"+d) as f:
			for l in f:
				if len(l.split()) == 2:
					keggs[l.split()[1]] += 1

for k in sorted(keggs, key=lambda x: keggs[x], reverse=True):
	print k,keggs[k]
