#!/usr/bin/env python

import sys, argparse, os
import csv
from collections import *

parser = argparse.ArgumentParser(description="Compare CNVs.")
parser.add_argument('file1', help="CSV file with CNV data")
parser.add_argument('file2', help="Input CSV")
args = parser.parse_args()

yeast_orf = defaultdict(str)

with open("/home/jacek/stuff/scripts/python/data/yeast_orfs.txt", 'rb') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:	
		yeast_orf[row[1]] = row[0]

homo_loss = []
hemi_loss = []
single_gain = []
multi_gain = []
with open(args.file1, 'rb') as f:
    reader = csv.reader(f)
    for row in reader:
		tp = row[0]
		if row[6] != "Call PValue" and float(row[6]) < 0.000000001:
			for gene in row[8].split(","):
				if row[2] == "Homozygous Copy Loss":
					homo_loss.append(gene+"+"+tp)
				elif row[2] == "CN Loss":
					hemi_loss.append(gene+"+"+tp)
				elif row[2] == "CN Gain":
					single_gain.append(gene+"+"+tp)
				elif row[2] == "High Copy Gain":
					multi_gain.append(gene+"+"+tp)


search_genes = []
with open(args.file2, 'rb') as f:
	for l in f:
		search_genes.append(l.strip())

for gene in search_genes:
	for g in homo_loss:
		if gene == g.split("+")[0]:
			print "Homozygous Copy Loss\t"+gene+"\t"+yeast_orf[gene]+"\t"+g.split("+")[1].split("_")[0]+"\t"+g.split("+")[1].split("_")[1]
	for g in hemi_loss:
		if gene == g.split("+")[0]:
			print "Hemizygous Copy Loss\t"+gene+"\t"+yeast_orf[gene]+"\t"+g.split("+")[1].split("_")[0]+"\t"+g.split("+")[1].split("_")[1]
	for g in single_gain:
		if gene == g.split("+")[0]:
			print "Single Copy Gain\t"+gene+"\t"+yeast_orf[gene]+"\t"+g.split("+")[1].split("_")[0]+"\t"+g.split("+")[1].split("_")[1]
	for g in multi_gain:
		if gene == g.split("+")[0]:
			print "Multi Copy Gain\t"+gene+"\t"+yeast_orf[gene]+"\t"+g.split("+")[1].split("_")[0]+"\t"+g.split("+")[1].split("_")[1]
