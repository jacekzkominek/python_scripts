#!/usr/bin/env python3

import os, sys, argparse
from collections import defaultdict
import numpy as np

parser = argparse.ArgumentParser(description='Obtain distribution of KEGG modules')
parser.add_argument("--file_mi", default="", help="Input file")
parser.add_argument("--min_mi", default="-1", help="")
parser.add_argument("--cross", default=False, action="store_true", help="")
parser.add_argument("--ji_high", default=False, action="store_true", help="")
parser.add_argument("--ji_low", default=False, action="store_true", help="")
parser.add_argument("--ji_filter", default=0.5, type=float, help="")
args = parser.parse_args()

pathways = defaultdict(list)
pathways_name = defaultdict(str)
kegg_to_pw = defaultdict(list)

min_mi = float
min_mi_dist = []
if args.min_mi.find("%") != -1:
	with open(args.file_mi) as f_mi:
		for l in f_mi:
			l = l.strip()
			mi = l.split()[2]
			min_mi_dist.append(float(mi))
	min_mi = np.quantile(np.array(min_mi_dist),float(args.min_mi.replace("%",""))/100)
	print("Filtering out all MI below "+str(min_mi)+" (bottom "+args.min_mi+")")
elif args.min_mi.find("%") == -1 and args.min_mi != "-1":
	min_mi = float(args.min_mi)
	print("Filtering out all MI below "+args.min_mi)

#READ IN ALL PATHWAYS
with open("/home/jkominek/scripts/python/data/keggs_2020.txt") as f:
	pw = ""
	for l in f:
		l = l.strip()
		if l[:2] == "ko":
			pw = l.split()[0][2:]
			pathways[pw] = []
			pathways_name[pw] = l.replace("ko"+pw+" ","").replace(" ","_").strip()
		if pw != "" and len(l) > 1 and l[0] == "K":
			k = l.split()[0]
			pathways[pw].append(k)
			kegg_to_pw[k].append(pw)

pathways_mi = defaultdict(set)
with open(args.file_mi) as f_mi:
	for l in f_mi:
		l = l.strip()
		k1 = l.split()[0]
		k2 = l.split()[1]
		mi = float(l.split()[2])
		ji = float(l.split()[3])
		if ji < args.ji_filter and args.ji_low == False:
			continue
		if ji >= args.ji_filter and args.ji_high == False:
			continue
		if mi > min_mi:
			pw1 = sorted(kegg_to_pw[k1])
			pw2 = sorted(kegg_to_pw[k2])
			if pw1 == [] or pw2 == []:
				continue
			for p1 in pw1:
				if p1 != "":
					if p1 in pw2:
						pathways_mi[p1].add(k1)
						pathways_mi[p1].add(k2)
					elif args.cross == True and p1 not in pw2:
						if k1+"x"+k2 in pathways_mi[p1] or k2+"x"+k1 in pathways_mi[p1]:
							continue
						else:
							pathways_mi[p1].add(k1+"x"+k2)
						# ~ pathways_mi[p1].add(k2+"x"+k1)
							

for pw in sorted(pathways_mi.keys()):
	if pw != "":
		print("\t".join([pw,pathways_name[pw]]+list(pathways_mi[pw])))



