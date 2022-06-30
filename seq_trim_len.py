#!/usr/bin/env python

from Bio import SeqIO
import os, sys, argparse

parser = argparse.ArgumentParser(description="Trim sequences shorter than specific length.")
parser.add_argument('--f1', default="", help="Input FASTA file")
parser.add_argument('--dir', default="", help="Input dir with FASTA files")
parser.add_argument('--length', type=int, default=200, help="Min length")
args = parser.parse_args()

dirname = ""
filenames = []
if args.dir == "":
	filenames.append(args.f1)
	dirname = os.getcwd()
else:
	for f in sorted(os.listdir(args.dir)):
		if f.find(".fas") != -1:
			filenames.append(f)
	dirname = args.dir
cut = int(args.length)

for filename in filenames:
	seqs = []
	for seq_record in SeqIO.parse(dirname+"/"+filename, "fasta"):
		if len(seq_record) >= cut:
			seqs.append(seq_record)
	print(filename)
	SeqIO.write(seqs, dirname+"/"+filename+"_cut"+str(cut)+".fas", "fasta")
