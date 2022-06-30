#!/usr/bin/env python

import os, sys, argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="")
parser.add_argument("f1", help="Input TSV file")
parser.add_argument("--min_mi", default="0.1", help="Input TSV file")
parser.add_argument("--print_all", default=False, action="store_true", help="Input TSV file")
parser.add_argument("--bin_file", default="", help="Input TSV file")
parser.add_argument("--bin_neg", default=False, action="store_true", help="Only negative associations")
parser.add_argument("--bin_pos", default=False, action="store_true", help="Only positive associations")
args = parser.parse_args()

ks = []
mis = []

with open(args.f1) as f:
	for l in f:
		if l.count("K") > 1:
			ks = l.strip().split()
		else:
			k = l.split()[0]
			for i,v in enumerate(l.split()):
				if i > 0:
					v = float(v)
					if k != ks[i-1]:
						mis.append(k+"@"+ks[i-1]+"@"+str(v))

bins = defaultdict(list)
if args.bin_file != "":
	with open(args.bin_file) as f1:
		for l in f1.readlines():
			l = l.strip() 
			s = l.split()
			if len(s) >= 8:
				k1 = s[0]
				k2 = s[2]
				a = s[3]
				b = s[4]
				c = s[5]
				d = s[6]
				jac = s[7]
				bins[k1+"_"+k2] = [a,b,c,d,jac]
				bins[k2+"_"+k1] = [a,b,c,d,jac]
				# ~ print k1+"_"+k2,bins[k1+"_"+k2]

# ~ for k in bins:
	# ~ print k,bins[k]

recorded = set()
for pair in sorted(mis, reverse=True, key=lambda x: float(x.split("@")[2])):
	s = pair.split("@")
	k1 = s[0]
	k2 = s[1]
	v = s[2]
	# ~ print k1,v,k2
	# ~ sys.stdout.flush()
	if args.print_all == True or str(k1+"_"+k2) not in recorded:
		if float(v) >= float(args.min_mi):
			if args.bin_file == "":
				print k1,v,k2
				sys.stdout.flush()
				recorded.add(str(k1+"_"+k2))
				recorded.add(str(k2+"_"+k1))
			else:
				if str(k1.strip("\"")+"_"+k2.strip("\"")) in bins:
					# ~ #if float(bins[(k1+"_"+k2)][4]) > 0.2:
					print k1,v,k2," ".join(bins[(k1.strip("\"")+"_"+k2.strip("\""))])
					# ~ sys.stdout.flush()
					recorded.add(str(k1+"_"+k2))
					recorded.add(str(k2+"_"+k1))
		elif float(v) < float(args.min_mi):
			break

