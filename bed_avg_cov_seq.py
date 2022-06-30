#!/usr/bin/env python

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Calculate average coverage across features")
parser.add_argument("f1", help="Input BED file")
args = parser.parse_args()

cov = defaultdict(list)

#READ COVERAGE DATA
with open(args.f1) as f:
	for l in f:
		cov[l.split()[0]].append([int(l.split()[1]),int(l.split()[2])])

#CALCULATE MEAN COVERAGE
cov_max = ["",0]
for g in cov:
	avg = 0
	total = 0
	for h in cov[g]:
		avg += h[0]*h[1]
		total += h[1]
	#cov_mean.append([g,float(avg)/float(total)])
	if float(avg)/float(total) > cov_max[1]:
		cov_max = [g,float(avg)/float(total)]
	#print g, str(float(avg)/float(total))
print cov_max[0],cov_max[1]
