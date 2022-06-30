#!/usr/bin/env python

import os
import argparse
import csv
from collections import defaultdict
import scipy.stats
import statsmodels.sandbox.stats.multicomp as sssm

parser = argparse.ArgumentParser(description='Check for mutations that drop by 10% from their max frequency')
parser.add_argument('csv_file', help="CSV file with input data.")
args = parser.parse_args()

#LOAD GENES AND MUTATIONS
with open(args.csv_file, "rb") as f:
	for l in f:
		if l.strip().split(",")[0].find("Fermentor") == -1:
			freqs = [float(x) for x in l.strip().split(",")[6:]]
			#print freqs
			if freqs.count(0) < len(freqs)-1:
				maxfreq = 0
				#FIND MAXFREQ
				for f in freqs:
					if f > maxfreq:
						maxfreq = f
				#CHECK FOR FREQ DROP
				max_f = False
				for f in freqs:
					if f == maxfreq:
						max_f = True
					if max_f == True and f < (maxfreq-0.1):
						print l.strip()
						break
