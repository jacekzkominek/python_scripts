#!/usr/bin/env python3

import os, sys, argparse
import pandas as pd
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description="")
parser.add_argument("f1", help="Input TSV file with MI coefs per line")
parser.add_argument("--min_mi", default="-1", help="Lower limit of MI to filter")
parser.add_argument("--sort", default=False, action="store_true", help="Sort by MI")
# ~ parser.add_argument("--network", default=False, action="store_true", help="Print out a .sif file format")
parser.add_argument("--bins", default=[], nargs=1, help="Simplify MIs into discrete number of bins")
args = parser.parse_args()

mis = []
with open(args.f1) as f:
	for l in f:
		k1 = l.split()[0]
		k2 = l.split()[1]
		mi = l.split()[2]
		jaccard_index = l.split()[3]
		if float(mi) > float(args.min_mi):
			mis.append([k1,k2,mi,jaccard_index])

if args.sort == True:
	mis = sorted(mis, key=lambda x:float(x[2]), reverse=True)

bins = []
if args.bins != []:
	bins = args.bins[0].split("_")

for mi in mis:
	if bins == []:
		print("\t".join(mi))
	elif args.bins != []:
		edge = ""
		for i,b in enumerate(bins):
			if float(b) >= float(mi[2]):
				edge = str(i)
				break
		if edge == "":
			edge = len(bins)
		print("\t".join([mi[0],mi[1],mi[2],mi[3],str(edge)]))


# ~ arr = numpy.array(lst)
# ~ pd.cut(df['ext price'], bins=4)
