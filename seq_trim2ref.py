#!/usr/bin/env python

import warnings
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
import sys, argparse, os

parser = argparse.ArgumentParser(description="Trim an alignment to a specific reference sequence.")
parser.add_argument('f1', help="Input FASTA file")
parser.add_argument('--ref_id', nargs=1, type=int, default=1, type=int, help="Ordinal number of the reference sequence (first seq by default)")
args = parser.parse_args()

f = args.f1

seqs = list(SeqIO.parse(f, "fasta"))
refseq = seqs[args.ref_id-1]
refseq_name = seqs[args.ref_id-1].id
gap_pos = set()
for x in range (0, len(refseq)):
	if refseq[x].find("-") != -1:
		gap_pos.add(x)

seq_pos = set()
for x in range (0, len(refseq)):
	seq_pos.add(x)

trim_pos = sorted(seq_pos.difference(gap_pos))	
degapped = []
for seq in seqs:
	seq1 = ""
	for x in trim_pos:
		seq1 += seq[x]
	degapped.append(seq1)

new_f = open(f.replace(".fas", "_reftrim.fas"),'w')

for i, seq in enumerate(seqs):
	new_f.write(">"+seq.id+"\n")
	new_f.write(degapped[i]+"\n")
new_f.close()

print "\nTrimmed the alignment in",f,"to",refseq_name,"\n"

