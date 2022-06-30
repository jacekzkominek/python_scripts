#!/usr/bin/env python

from Bio import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys
import argparse
import string 

parser = argparse.ArgumentParser(description="Filter stop codons (\"*\") from protein sequences.")
parser.add_argument('file1', help="Input file")
parser.add_argument('--conv', action='store_true', default=False, help="Convert ambiguity characters to Ns")
args = parser.parse_args()

f = args.file1
ambig = ["Y","R","W","S","K","M","D","V","H","B","X"]
out_seqs = []
for seq in SeqIO.parse(f, "fasta"):
	for a in ambig:
		if str(seq.seq).find(a) != -1:
			print f+"\t"+a+"\t"+seq.id
	if args.conv == True:
		seq2_seq = Seq(str(seq.seq).translate(string.maketrans("".join(ambig),"NNNNNNNNNNN")),IUPAC.unambiguous_dna)
		seq2 = SeqRecord(seq2_seq, id=seq.id, description=seq.description, name=seq.name)
		#print seq2
		out_seqs.append(seq2)

if args.conv == True:
	SeqIO.write(out_seqs,f.replace(".fas","_N.fas"),"fasta")
	print f,"converted ambiguous bases to N"
