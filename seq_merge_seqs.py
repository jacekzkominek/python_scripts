#!/usr/bin/env python

import os, sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Merge sequences from multiple FASTA files into one.')
parser.add_argument('--dir', default=".", help="Dir with FASTA files")
#parser.add_argument('--files', nargs='*', help="Subject FASTA files")
parser.add_argument('--MB', default="275", help="Max filesize")

args = parser.parse_args()

out_seqs = []
filesize = 0
i = 1
j = 1
for f in sorted(os.listdir(args.dir)):
	if f.find(".fas") != -1:
		if filesize+os.path.getsize(f) >= int(args.MB)*1000000:
			SeqIO.write(out_seqs,os.getcwd()+"/all_merged"+str(i)+".fas", "fasta")
			out_seqs = []
			i += 1
			filesize = 0
		seqs = SeqIO.parse(f,"fasta")
		filesize += os.path.getsize(f)
		print str(j),f,filesize/1000000
		j += 1
		sys.stdout.flush()
		for s in seqs:
			tmp = s.description
			s.description = ""
			s.name = ""
			s.id = ""
			s.id = f+"_"+tmp
			out_seqs.append(s)
		
SeqIO.write(out_seqs,os.getcwd()+"/all_merged"+str(i)+".fas", "fasta")


		
