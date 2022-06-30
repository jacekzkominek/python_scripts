#!/usr/bin/env python

from Bio.Blast.Applications import *
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SearchIO
import sys, argparse, os, subprocess, shutil
from collections import defaultdict
import copy, string

parser = argparse.ArgumentParser(description='Run BLAST over the local database.')
parser.add_argument('--blast_output', default=".", help="Format 6 BLAST results")
parser.add_argument('--genome_file', help="Genome FASTA file")
parser.add_argument('--flank', default=1000)
args = parser.parse_args()

hits = []
with open(args.blast_output,"r") as bf:
	for l in bf:
		sp = l.strip().split()
		query = sp[0]
		target = sp[1]
		target_start = sp[8]
		target_end = sp[9]
		hits.append([query,target,target_start,target_end])

prev_hit = ""
inc = 2
for h in hits:
	newseqs = []
	query = h[0]
	target = h[1]
	start = int(h[2])
	end = int(h[3])
	for inseq in SeqIO.parse(args.genome_file, "fasta"):
		if inseq.description.find(target) != -1:
			newseq = copy.deepcopy(inseq)
			if start < end:
				newseq.seq = inseq.seq[max(0,start-int(args.flank)):end+int(args.flank)]
			elif start > end:
				newseq.seq = inseq.seq[max(0,end-int(args.flank)):start+int(args.flank)].reverse_complement()
			print(len(newseq.seq))
			if query != prev_hit:
				newseq.id = query
				inc = 2
			else:
				newseq.id = query+"_"+str(inc)
				inc += 1
			newseq.title = ""
			newseq.description = ""
			newseqs.append(newseq)
			prev_hit = query
			break
	outfile = open(query+".fas","a+")
	SeqIO.write([newseq],outfile,"fasta")
	outfile.close()
			

