#!/usr/bin/env python

import os, sys, argparse
from collections import defaultdict
import numpy

parser = argparse.ArgumentParser(description='')
parser.add_argument('infile', help="Infile with input data.")
args = parser.parse_args()

contig_genes = defaultdict(list)

with open(args.infile, "rb") as f:
	for l in f:
		if l[0] != "#" and len(l.split()) >= 5:
			l = l.strip()
			feat = l.split()[2]
			if feat == "similarity":
				contig = l.split()[0]
				start = l.split()[3]
				end = l.split()[4]
				contig_genes[contig].append(start+"_"+end)

global_dist = []
for c in contig_genes:
	if len(contig_genes[c]) >= 10:
		contig_genes_sorted = sorted(contig_genes[c], key= lambda x: int(x.split("_")[0]))
		dist = []
		prev_start = 0
		prev_end = 0
		#print contig_genes_sorted
		for i,g in enumerate(contig_genes_sorted):
			start = int(g.split("_")[0])
			end = int(g.split("_")[1])
			if i >= 1:
				dist.append(start-prev_end)
			prev_start = start
			prev_end = end
			##i += 1
			if i > len(contig_genes_sorted):
				break
		if len(dist) > 1:
			global_dist = global_dist + dist
			# ~ print c,dist
			# ~ for d in dist:
				# ~ print d
			print c,numpy.mean(dist),numpy.std(dist), numpy.median(dist), numpy.subtract(*numpy.percentile(dist, [75, 25]))

# ~ print global_dist
print numpy.mean(global_dist), numpy.std(dist), numpy.median(global_dist), numpy.subtract(*numpy.percentile(dist, [75, 25]))

# seqname source feature start end score strand frame attributes
#
#gi|1002316101|dbj|BCJV01000001.1|	exonerate:protein2genome:local	gene	301060	302619	2762	+	.	gene_id 1 ; sequence augustus_masked-gi|1002316101|dbj|BCJV01000001.1|-processed-gene-3.9-mRNA-1 ; gene_orientation .
#gi|1002316101|dbj|BCJV01000001.1|	exonerate:protein2genome:local	cds	301060	302619	.	+	.	
#gi|1002316101|dbj|BCJV01000001.1|	exonerate:protein2genome:local	exon	301060	302619	.	+	.	insertions 0 ; deletions 0
#gi|1002316101|dbj|BCJV01000001.1|	exonerate:protein2genome:local	similarity	301060	302619	2762	+	.	alignment_id 1 ; Query augustus_masked-gi|1002316101|dbj|BCJV01000001.1|-processed-gene-3.9-mRNA-1 ; Align 301060 1 1560
# --- END OF GFF DUMP ---
#
# --- START OF GFF DUMP ---
#
#
