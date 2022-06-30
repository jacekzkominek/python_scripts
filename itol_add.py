#!/usr/bin/env python

import warnings, argparse, os
from collections import defaultdict

parser = argparse.ArgumentParser(description="Extract proteinortho clusters")
parser.add_argument("original_file", help="Input file1 (ref)")
parser.add_argument("addfile", help="Adding data")
args = parser.parse_args()

data = defaultdict(list)
with open(args.addfile) as f1:
	for l in f1:
		l = l.strip()
		data[l.split()[0].lower()] = l.split()[1]

is_data = False			
printed = []
with open(args.original_file) as f2:
	for l in f2:
		l = l.strip()
		if is_data == True and len(l.split(",")) > 1:
			val = -1
			for k in data.keys():
				if l.split(",")[0].lower().find(k) != -1:
					val = data[k]
					printed.append(k)
					break
			print l+","+str(val)
		if is_data == False and l == "DATA":
			is_data = True

for k in data.keys():
	if k not in printed:
		print k,data[k]
