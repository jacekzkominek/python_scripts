#!/usr/bin/env python

import os, sys, argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Coverage")
parser.add_argument("f1", help="Input COV file")
parser.add_argument("f2", help="Input COV file")
parser.add_argument("start", help="")
parser.add_argument("end", help="")
args = parser.parse_args()

cov_table = defaultdict(list)
for p in range(int(args.start),int(args.end)+1):
	cov_table[p] = [0,0,0,0,0,0,0,0]

with open(args.f1,"r") as f1:
	for l in f1:
		l = l.strip()
		if len(l.split()) == 5 and l.split()[0] != "POS":
			pos = int(l.split()[0])
			cov = l.split()[1]
			span = l.split()[2]
			cov_span = l.split()[3]
			prop = l.split()[4]
			
			#~ print l 
			if pos < int(args.end) and pos >= int(args.start) :
				cov_table[pos][0] = int(cov)
				cov_table[pos][1] = int(span)
				cov_table[pos][2] = int(cov_span)
				cov_table[pos][3] = float(prop)

with open(args.f2,"r") as f2:
	for l in f2:
		if len(l.split()) == 5 and l.split()[0] != "POS":
			l = l.strip()
			pos = int(l.split()[0])
			cov = l.split()[1]
			span = l.split()[2]
			cov_span = l.split()[3]
			prop = l.split()[4]
			
			#~ print l
			if pos < int(args.end) and pos >= int(args.start) :
				cov_table[pos][4] = -1*int(cov)
				cov_table[pos][5] = -1*int(span)
				cov_table[pos][6] = -1*int(cov_span)
				cov_table[pos][7] = -1*float(prop)

print "POS","READ","SPAN","SUM","PROP","REV_READ","REV_SPAN","REV_SUM","REV_PROP"
for p in sorted(cov_table, key=lambda x: int(x)):
	print p," ".join(str(x) for x in cov_table[p])

