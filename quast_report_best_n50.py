#!/usr/bin/env python

import os, sys
import argparse
import shutil
from collections import defaultdict

parser = argparse.ArgumentParser(description="")
parser.add_argument("report", default = "", help="Transposed_report.tsv")
parser.add_argument("assemblies", default = "", help="Directory with assemblies")
parser.add_argument("target", default = "", help="Target directory to output best assemblies")
parser.add_argument("--noexp", default = False, action = "store_true", help="Don't export")
parser.add_argument("--extend", default = False, action = "store_true", help="Print extensive QUAST report")
args = parser.parse_args()

stats_n50 = defaultdict(list)
stats_gc = defaultdict(float)
stats_size = defaultdict(int)
with open(args.report) as f:
	for l in f:
		if l.find("#") == -1:
			name = "_".join(l.split()[0].split("_")[:-1] )
			n50 = l.split()[-5]
			stats_n50[name].append(str(l.split()[0]+"="+n50))
			gc = float(l.split()[-6])
			stats_gc[str(l.split()[0])] = gc
			
			size = int(l.split()[-7])
			stats_size[str(l.split()[0])] = size
			
if os.path.exists(args.target) == False:
	os.makedirs(args.target)
	
for n in stats_n50.keys():
	best_n50 = 0
	best_size = 0
	best_a_n50 = ""
	best_a_size = ""
	for a in stats_n50[n]:
		n50 = int(a.split("=")[1])
		if n50 > best_n50:
			best_n50 = n50
			best_a_n50 = a.split("=")[0]
		size = stats_size[a.split("=")[0]]
		if size > best_size:
			best_size = size
			best_a_size = a.split("=")[0]
	if stats_size[best_a_n50] >= 0.9*best_size:
		print best_a+"\t"+str(best_n50),
		if args.extend == True:
			print str(stats_gc[best_a_n50])+"\t"+str(stats_size[best_a_n50])
		else:
			print
		if args.noexp == False:
			shutil.copy(os.path.join(args.assemblies,best_a_n50+".fas"),os.path.join(args.target,best_a_n50+".fas"))


