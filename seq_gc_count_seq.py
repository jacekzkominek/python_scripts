#!/usr/bin/env python

import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, argparse, os
from collections import defaultdict

parser = argparse.ArgumentParser(description="")
parser.add_argument("--infile", help="Input file")
parser.add_argument("--dir", default="", help="")
parser.add_argument("--total", default=False, action="store_true", help="")
args = parser.parse_args()


if args.dir != "":
	for f in sorted(os.listdir(args.dir)):
		print f
		seqs = list(SeqIO.parse(args.dir+"/"+f,"fasta"))
		total_glob = 0
		gc_count_glob = 0
		for seq in seqs:
			gc_count = 0
			total = 0
			for pos in seq.seq:
				total += 1
				total_glob += 1
				if pos.upper() == "G" or pos.upper() == "C":
					gc_count += 1
					gc_count_glob += 1
			if args.total == False:
				print seq.description,gc_count,total, "{:.2f}".format(float(float(gc_count)/float(total))*100)
		if total_glob > 0:
			print "Total",str(gc_count_glob),str(total_glob),"{:.2f}".format(float(float(gc_count_glob)/float(total_glob))*100)+"\n"
			sys.stdout.flush()

	
else:
	seqs = list(SeqIO.parse(args.infile,"fasta"))
	total_glob = 0
	gc_count_glob = 0
	for seq in seqs:
		gc_count = 0
		total = 0
		for pos in seq.seq:
			total += 1
			total_glob += 1
			if pos.upper() == "G" or pos.upper() == "C":
				gc_count += 1
				gc_count_glob += 1
		print seq.description,gc_count,total, "{:.2f}".format(float(float(gc_count)/float(total))*100)
	print "Total",gc_count_glob,total_glob
