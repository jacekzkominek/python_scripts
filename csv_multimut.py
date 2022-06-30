#!/usr/bin/env python

import os
import argparse
import csv
from collections import defaultdict
import scipy.stats
import statsmodels.sandbox.stats.multicomp as sssm

parser = argparse.ArgumentParser(description='')
parser.add_argument('csv_file', help="CSV file with input data.")
parser.add_argument('--only_multif', action='store_true', default=False, help="Only print genes hit in multiple different reactors")
args = parser.parse_args()

#LOAD GENES AND MUTATIONS
genes = defaultdict(list)
with open(args.csv_file, "rb") as f:
	for l in f:
		l1 = l.split(",")
		if l1[1].find("chr") != -1:
			genes[l1[5]].append(l1[2]+"_"+l1[0])

for k in genes.keys():
	if args.only_multif == False:
		print k, len(genes[k]), ",".join(genes[k])
	else:
		fs = set()
		for p in genes[k]:
			fs.add(p.split("_")[1])
			if len(fs) > 1:
				print k, len(genes[k]), ",".join(genes[k]),
				break

##DSK2
#print scipy.stats.hypergeom.pmf(6,9080922,1122,817)
#print scipy.stats.binom.pmf(6,817,1122/9080922)
##FMP27
#print scipy.stats.hypergeom.pmf(5,9080922,7887,817)
#print scipy.stats.binom.pmf(5,817,7887/9080922)
##ASG1
#print scipy.stats.hypergeom.pmf(5,9080922,2895,817)
##LAA1
#print scipy.stats.hypergeom.pmf(4,9080922,6045,817)
##HEM3
#print scipy.stats.hypergeom.pmf(4,9080922,984,817)
##ACE2
#print scipy.stats.hypergeom.pmf(3,9080922,2313,817)
##POL3
#print scipy.stats.hypergeom.pmf(3,9080922,3294,817)
##RAP1
#print scipy.stats.hypergeom.pmf(3,9080922,2484,817)
##PUF4
#print scipy.stats.hypergeom.pmf(2,9080922,2667,817)
##YJL070C
#print scipy.stats.hypergeom.pmf(2,9080922,2667,817)
##CSE1
#print scipy.stats.hypergeom.pmf(2,9080922,2883,817)
##RGC1
#print scipy.stats.hypergeom.pmf(2,9080922,3252,817)
##EPS1
#print scipy.stats.hypergeom.pmf(2,9080922,2106,817)
##GFA1
#print scipy.stats.hypergeom.pmf(2,9080922,2154,817)
##UTH1
#print scipy.stats.hypergeom.pmf(2,9080922,1098,817)
##GRR1
#print scipy.stats.hypergeom.pmf(2,9080922,3456,817)
##YER156C
#print scipy.stats.hypergeom.pmf(2,9080922,1017,817)
##JID1
#print scipy.stats.hypergeom.pmf(2,9080922,906,817)
##YJL182C
#print scipy.stats.hypergeom.pmf(2,9080922,318,817)
##RIM15
#print scipy.stats.hypergeom.pmf(2,9080922,5313,817)
##IRA2
#print scipy.stats.hypergeom.pmf(2,9080922,9240,817)
##ATG11
#print scipy.stats.hypergeom.pmf(2,9080922,3537,817)
##TRS130
#print scipy.stats.hypergeom.pmf(2,9080922,3309,817)
##HEM12
#print scipy.stats.hypergeom.pmf(2,9080922,1098,817)
##HRK1
#print scipy.stats.hypergeom.pmf(2,9080922,2280,817)
##DPB11
#print scipy.stats.hypergeom.pmf(2,9080922,2295,817)
##UBR1
#print scipy.stats.hypergeom.pmf(2,9080922,5853,817)
#$print "aaaaaaaa%.10f%%" % (scipy.stats.binom.pmf(2,817,5853/9080922))
##VPS74
#print scipy.stats.hypergeom.pmf(2,9080922,1038,817)
