#!/usr/bin/env python

import os, sys, argparse, string
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='')
parser.add_argument("gtf", help="GTF file")
parser.add_argument("fasta", help="FASTA file")
parser.add_argument("--igs", default="300", help="")
parser.add_argument("--out", default="upstream.fas", help="")
args = parser.parse_args()

contig_genes = defaultdict(list)
with open(args.gtf, "rb") as f:
	start2 = -1
	stop2 = -1
	for l in f:
		l = l.strip()
		contig = l.split()[0]
		cds = l.split()[2]
		start = l.split()[3]
		end = l.split()[4]
		strand = l.split()[6]
		gene = l.split()[-1].replace("\"","").replace(";","")
		if cds == "start_codon":
			if strand == "+":
				start2 = start
			if strand == "-":
				start2 = end
		if cds == "stop_codon":
			contig_genes[contig].append([gene,start2,strand])
			start2 = -1

			
	
seqs = list(SeqIO.parse(args.fasta, "fasta"))
igs_seqs = []
for seq in seqs:
	title = seq.description
	for contig in contig_genes:
		if contig in title:
			for gene in contig_genes[contig]:
				name = gene[0]
				start = gene[1]
				strand = gene[2]
				s = ""
				if strand == "+":
					s = seq.seq[int(start)-int(args.igs)-1:int(start)-1]
				elif strand == "-":
					s = seq.seq[int(start):int(start)+int(args.igs)]
					trans = string.maketrans("ATGC","TACG")
					s = str(s).translate(trans)[::-1]
				if len(s) == int(args.igs):
					igs_seqs.append([name,s])

with open(args.out,"w") as out_f:		
	for s in igs_seqs:
		out_f.writelines(">"+s[0]+"\n"+s[1]+"\n")
		
