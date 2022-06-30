#!/usr/bin/env python

import warnings
import sys, argparse, os, subprocess
from collections import defaultdict

parser = argparse.ArgumentParser(description="")
parser.add_argument("dir", help="Target dir")
parser.add_argument("--parse", default=False, action="store_true", help="")
args = parser.parse_args()


first = True
for f in sorted(os.listdir(args.dir)):
	if args.parse == False:
		if f.find(".fas") != -1:
			print f
			sys.stdout.flush()
			prog = "/home/jacek/software/gene_analysis/fLPS/bin/linux/fLPS"
			with open(args.dir+"/"+f.replace(".fas",".txt"),"w") as outf:
				subprocess.call([prog,args.dir+"/"+f], stdout=outf)
	else:
		sigs = defaultdict(int)
		if f.find(".txt") != -1:
			with open(args.dir+"/"+f,"r") as f1:
				for l in f1:
					seq = l.split()[0]
					typ = l.split()[1]
					num = l.split()[2]
					start = l.split()[3]
					end = l.split()[4]
					count = l.split()[5]
					prob = float(l.split()[6])
					sig = l.split()[7]
					if prob <= 0.001 and len(sig) == 3 and sig != "{X}":
						sigs[sig] += 1
			if first == True:
				for s in sorted(sigs.keys()):
					print s,
				print
				first = False
			print f,
			for s in sorted(sigs.keys()):
				print sigs[s],
			print
			sys.stdout.flush()
