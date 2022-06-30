#!/usr/bin/env python

import warnings
import sys, argparse, os

parser = argparse.ArgumentParser(description="Inspect results of InterProScan.")
parser.add_argument('file', help="Input file")
parser.add_argument('--bitscore', action='store_true', default=False, help="Extract raw bitscores instead of E-values")
args = parser.parse_args()

f = open(args.file)
mcl_tabs = []
g1=""
g2=""
score=""
score_col = 2
if args.bitscore == True:
	score_col = 1

for l in f:
	if l.find("Query=") != -1:
		g1 = l.split("=")[1].strip()
	if len(l.split()) == 3:
		g2 = l.split()[0]
		score = l.split()[score_col]
		print g1+"\t"+g2+"\t"+score
