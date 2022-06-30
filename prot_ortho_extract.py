#!/usr/bin/env python

import warnings, argparse, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close
from collections import defaultdict

parser = argparse.ArgumentParser(description="Extract proteinortho clusters")
parser.add_argument("f1", help="Input file")
parser.add_argument("--ac", default=0.95, help="Input file")
parser.add_argument("--species", default=0.90, help="Input file")
parser.add_argument("--genes", default=2, help="Input file")
args = parser.parse_args()

orthoclust = []
max_species = 0
species_list = []

with open(args.f1) as f:
	for l in f:
		if l.find("#") == -1:
			max_species = max(int(l.split("\t")[0]),max_species)
			
with open(args.f1) as f:
	for l in f:
		if l.find("#") != -1:
			species_list = l.strip().split("\t")[3:]
		else:
			species = int(l.split("\t")[0])
			genes = int(l.split("\t")[1])
			ac = float(l.split("\t")[2])
			
			if species >= max_species*float(args.species) and genes <= max_species*float(args.genes) and ac >= float(args.ac):
				orthoclust.append(l.strip().split("\t")[3:])

print len(orthoclust),"clusters identified"
#~ for c in orthoclust:
	#~ print c
#~ print species_list

seq_clust = defaultdict(list)
for j,s in enumerate(species_list):
	f_s = open(s)
	print s
	for seq in SeqIO.parse(f_s, "fasta"):
		for i,clust in enumerate(orthoclust):
			clust_seq = clust[j]
			if seq.name.strip() == clust_seq.strip():
				seq.id = seq.id+"_"+s
				seq.name = ""
				seq.description = ""
				seq_clust[i].append(seq)
				out = True
				break

out_dir = os.getcwd()+"/clusters_"+str(len(orthoclust))+"_ac"+str(args.ac)+"_s"+str(args.species)+"_g"+str(args.genes)+"/"
os.makedirs(out_dir)
for key in seq_clust.keys():
	SeqIO.write(seq_clust[key], out_dir+"/cluster"+str(key+1)+".fas", "fasta")
	
