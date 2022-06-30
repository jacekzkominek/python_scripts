#!/usr/bin/env python

import os, sys
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='Compare two kegg annotations')
parser.add_argument("file1", help="Files with input data.")
parser.add_argument("file2", help="Files with input data.")
#parser.add_argument("--transpose", action="store_true", default=False, help="Transpose the KEGG matrix")
args = parser.parse_args()

kegg1 = set()#defaultdict(str)
kegg2 = set()#defaultdict(str)

with open(args.file1) as f1:
	for l in f1:
		if len(l.split()) == 2:
			#~ kegg1[l.split()[0]] = l.split()[1]
			kegg1.add(l.split()[0]+"_"+l.split()[1])

with open(args.file2) as f2:
	for l in f2:
		if len(l.split()) == 2:
			#~ kegg2[l.split()[0]] = l.split()[1]
			kegg2.add(l.split()[0]+"_"+l.split()[1])

print args.file1, len(kegg1), len(kegg2), len(kegg1.intersection(kegg2))


