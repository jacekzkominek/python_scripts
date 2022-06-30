#!/usr/bin/env python

import os, sys
import argparse
from natsort import natsorted, realsorted
from collections import defaultdict

parser = argparse.ArgumentParser(description='Compare two kegg annotations')
parser.add_argument("file1", help="KEGG annotation file.")
#parser.add_argument("file2", help="FASTA sequence file.")
parser.add_argument("--kegg", help="KEGG marker.")
parser.add_argument("--n", default="10", help="Neighbors.")
#parser.add_argument("--transpose", action="store_true", default=False, help="Transpose the KEGG matrix")
args = parser.parse_args()

seq_dict = defaultdict(list)
gene_list = []
with open(args.file1) as f1:
	mark = ""
	for l in f1:
		gene = l.split()[0]
		kegg = ""
		if len(l.split()) == 2:
			kegg = l.split()[1]
		contig = ""
		locator = ""
		if l.find("genemark") == -1 and l.find("maker") == -1:
			contig = l.split("gene=")[1].split("masked-")[1].split("-gene-")[0]
			locator = l.split("gene=")[1].split("masked-")[1].split("-gene-")[1].split("_CDS")[0]
		elif l.find("genemark") != -1:
			contig = l.split("gene=")[1].split("genemark-")[1].split("-gene-")[0]
			locator = l.split("gene=")[1].split("genemark-")[1].split("-gene-")[1].split("_CDS")[0]
		elif l.find("maker") != -1:
			contig = l.split("gene=")[1].split("maker-")[1].split("-gene-")[0]
			locator = l.split("gene=")[1].split("maker-")[1].split("-gene-")[1].split("_CDS")[0]
		locator1 = locator.split(".")[0]
		locator2 = locator.split(".")[1]
		gene_list.append([gene,contig.split("-")[0],locator1,locator2,kegg])
		if len(l.split()) == 2 and l.find(args.kegg) != -1:
			mark = gene

prev_list = []
loc = ""
prev_loc = ""
next_list = []
markpass = False
for sort_gene in natsorted(gene_list, key=lambda x: (x[1], x[2], x[3])):
	if markpass == True and len(next_list) < int(args.n)+1:
		next_list.append([sort_gene[0],sort_gene[4]])
	loc = sort_gene[1]+"_"+sort_gene[2]
	if markpass == False and prev_loc == loc:
		prev_list.append([sort_gene[0],sort_gene[4]])
	prev_loc = loc
	#print sort_gene[0],sort_gene[1],sort_gene[2],sort_gene[3],
	if sort_gene[0] == mark:
		markpass = True

if markpass == True:
	print "#PREV"
	prev2 = []
	for i,p in enumerate(reversed(prev_list)):
		if i < int(args.n)+1:
			prev2.append(p)
	for p in reversed(prev2):	
		print p[0],p[1]
	print "#NEXT"
	for n in next_list:
		print n[0],n[1]
