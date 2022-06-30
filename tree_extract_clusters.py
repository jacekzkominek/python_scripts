#!/usr/bin/env python

import os
import argparse
import itertools
from Bio import SeqIO
from collections import defaultdict

parser = argparse.ArgumentParser(description="Run clustering pipeline for sequence data.")
parser.add_argument("f", help="FASTA file with source sequences")
parser.add_argument("c", help="File with cluster data")
parser.add_argument("--species", default="", help="File with species names to use for subfamily annotation")
args = parser.parse_args()

#EXTRACT THE CLUSTERS
clusters = defaultdict(list)
with open(args.c, "r") as c_file:
	flag = False
	for l in c_file:
		if flag == False and l.find("#Fully resolved clusters") != -1:
			flag = True
			continue
		if flag == True and len(l.strip().split("\t")) > 0:
			clusters[len(clusters.keys())] = l.strip().split("\t")


#EXTRACT THE SEQUENCES AND EXPORT CLUSTER FASTA FILES
if args.f != "" and os.path.exists(args.f):
	out_seqs = defaultdict(list)
	source_seqs = SeqIO.parse(args.f, "fasta")
	no_clus = []
	for seq in source_seqs:
		found = False
		for c in clusters.keys():
			if seq.id in clusters[c]:
				out_seqs[c].append(seq)
				found = True
				break
		if found == False:
			no_clus.append(seq)

	#SAVE CLUSTERS (AND NO CLUSTER)
	for c in out_seqs.keys():
		SeqIO.write(out_seqs[c], args.f+"_cluster"+str(c+1)+".fas", "fasta")
	SeqIO.write(no_clus, args.f+"_clusterX.fas", "fasta")


#EXTRACT SPECIES NAMES AND EXPORT SUBFAMILY ANNOTATIONS
if args.species != "":
	species = []
	with open(args.species, "r") as species_file:
		for l in species_file:
			species.append(l.strip())

	#SPECIES-BASED
	species_list = []
	for i,s in enumerate(species):
		s_clus = [s]
		for cluster_id in sorted(clusters.keys()):
			c = 0
			for seq in clusters[cluster_id]:
				if seq.find(s) != -1:
					c += 1
			s_clus.append(str(c))
		species_list.append(s_clus)

	for i in species_list:
		print "\t".join(i)

	#CLUSTER-BASED
	#species_dict = defaultdict(list)
	#for cluster_id in sorted(clusters.keys()):
		#print clusters[cluster_id]
		#for s in species:
			#c = "0"
			#for seq in clusters[cluster_id]:
				#if seq.find(s) != -1:
					#c = "1"
					#break
			#species_dict[cluster_id].append(c)

		#print species_dict[cluster_id]