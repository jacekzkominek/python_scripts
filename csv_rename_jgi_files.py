#!/usr/bin/env python

import os
import argparse
import csv

parser = argparse.ArgumentParser(description='Check for mutations that drop by 10% from their max frequency')
parser.add_argument('csv_file', help="CSV file with input data.")
args = parser.parse_args()

names = {}
with open(args.csv_file, "rb") as f:
	for l in f:
		if len(l.split(",")) > 1 and l.split(",")[0].split() > 1 and len(l.split(",")[0].split()) >= 2:
			id1 = l.split(",")[1].strip("\"")
			full = l.split(",")[0].strip("\"")
			#names[id1.strip()] = full.split()[0]+"_"+full.split()[1]
			names[id1.strip()[0:5].strip("123")] = full.split()[0]+"_"+full.split()[1]

for f in os.listdir(os.getcwd()):
	if len(f.split("_")) >= 2:
		#short = f.split("_")[0]
		short = f.split("_")[0].strip("123")
		#print short
		if names.get(short) != None:
			#print f,names[short]+"_"+"_".join(f.split("_")[1:])
			os.rename(f,names[short]+"_"+"_".join(f.split("_")[1:]))