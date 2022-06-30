#!/usr/bin/env python

import os
import fnmatch
from Bio import SeqIO
import argparse
import csv
from collections import defaultdict
import urllib  
  
parser = argparse.ArgumentParser(description='')
parser.add_argument('csv_file', default="", nargs='?', help="CSV file with input data.")
#parser.add_argument('--codon_pos', default="10000", type=int, nargs='?', help="Only include SNPs at positions before or at this codon.")
#parser.add_argument('--log_file', default="", help="Name of the log file to ouput.")
parser.add_argument('--skip_ok', action='store_true', default=False, help="")
args = parser.parse_args()

sgd_file = open("/home/jacek/scripts/python/sgd_gff.tsv","r")
sgd_genes = defaultdict(list)
sgd_ids = defaultdict(list)
for l in sgd_file:
	l1 = l.split("\t")
	if l1[2] == "gene":
			if sgd_genes[urllib.unquote(l1[10][5:])] == []:
				sgd_genes[urllib.unquote(l1[10][5:])]=[int(l1[3]),int(l1[4]),l1[6]]
			if sgd_genes[urllib.unquote(l1[9][5:])] == []:
				sgd_genes[urllib.unquote(l1[9][5:])]=[int(l1[3]),int(l1[4]),l1[6]]

sgd_file.seek(0)				
for l in sgd_file:
	l1 = l.split("\t")
	if l1[2] == "gene":
		if l1[11].find("Alias=") != -1:
			for alias in l1[11][6:].split(","):	
				if sgd_genes[urllib.unquote(alias)] == []:
					sgd_genes[alias]=[int(l1[3]),int(l1[4]),l1[6]]
sgd_file.close()

f = open(args.csv_file, "rb")
gene_pos = []
for l in f:
	l1 = l.split(",")
	if l1[6] != "" and (l1[5] == "\"Synonynous\"" or l1[5] == "\"Non-synonymous\""):
		gene_pos.append([str.replace(l1[6][1:-1],'_',','),int(l1[1])])
f.close()

#print gene_pos

for gene in gene_pos:
	if sgd_genes[gene[0]] == []:
		print gene[0]+"\t==NOT FOUND=="
	elif gene[1] < sgd_genes[gene[0]][0] or gene[1] > sgd_genes[gene[0]][1]:
		print gene[0]+"\t"+str(sgd_genes[gene[0]])+"\t"+str(gene[1])
	elif args.skip_ok == False:
		print gene[0]+"\t==OK=="
