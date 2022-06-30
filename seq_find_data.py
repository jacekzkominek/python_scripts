#!/usr/bin/env python

import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys, argparse, re
from collections import defaultdict

parser = argparse.ArgumentParser(description="Find a motif in a multi-sequence file and print sequences that don't contain it.")
parser.add_argument('file', help="Input file")
parser.add_argument('motif', help="Sequence motif to find")
parser.add_argument('--extract', action='store_true', help="Extract sequences with and without the searched motif")
parser.add_argument('--extract_terminus', action='store_true', help="Extract sequences with the searched motif based on the N and C terminus")
parser.add_argument('--full_info', action='store_true', help="Print positional information of the found motifs.")
parser.add_argument('--seq_sum', action='store_true', help="Print amount of motifs in a sequence.")
parser.add_argument('--degap', action='store_true', help="Degap sequences before processing.")
args = parser.parse_args()

f1 = open(args.file)
mot = args.motif
seqs1 = list(SeqIO.parse(f1, "fasta"))
f1.close()
seqs_out1 = []
seqs_out2 = []
motifs = defaultdict(list)
seqs = []

for seq in seqs1:
	title = seq.description
	seqs.append(title)
	hit = False
	s = seq.seq
	
	if args.degap == True:
		s = str(seq.seq).replace("-","")
	
	for m in re.finditer(str(mot), str(s)):
		if args.full_info == True:
			print mot+'\tfound at position\t'+str(m.start()+1)+"\t"+title
			motifs[title].append(str(m.start()))
		if args.extract == True and hit == False:
			seqs_out1.append(seq)  
		if args.extract_terminus == True and hit == False:
			if m.start() <= len(s)/2:
				seqs_out1.append(seq)
			elif m.start() > len(s)/2:
				seqs_out2.append(seq)
		hit = True
	if hit == False:
		print mot+"\tnot found in \t"+title
		if args.extract == True:
			seqs_out2.append(seq)

if args.extract == True:
	SeqIO.write(seqs_out1, args.file.replace(".fas","_motif_yes.fas"), "fasta")
	SeqIO.write(seqs_out2, args.file.replace(".fas","_motif_no.fas"), "fasta")

if args.extract_terminus == True:
	SeqIO.write(seqs_out1, args.file.replace(".fas","_Nterm.fas"), "fasta")
	SeqIO.write(seqs_out2, args.file.replace(".fas","_Cterm.fas"), "fasta")

if args.seq_sum == True:
	print "\nTotal number of motifs:"
	for seq in sorted(seqs):
		print seq+"\t"+str(len(motifs[seq]))
