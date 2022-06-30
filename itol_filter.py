#!/usr/bin/env python

import warnings, argparse, os
from collections import defaultdict

parser = argparse.ArgumentParser(description="Extract proteinortho clusters")
parser.add_argument("ref_file", help="Input file1 (ref)")
parser.add_argument("filter_file", help="Input file2 (filter for ref)")
args = parser.parse_args()

is_data = False
data = defaultdict(list)
with open(args.ref_file) as f1:
	for l in f1:
		l = l.strip()
		if l.find("#") == -1:
			if is_data == True and len(l.split()) > 1:
				data[l.split()[0]] = l.split()[1:]
			if is_data == False and l == "DATA":
				is_data = True

is_data = False			
with open(args.filter_file) as f2:
	for l in f2:
		l = l.strip()
		if is_data == True and len(l.split(",")) > 1:
			if l.split(",")[0] in data.keys():
				print l
		if is_data == False and l == "DATA":
			is_data = True

