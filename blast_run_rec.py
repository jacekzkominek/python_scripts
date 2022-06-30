#!/usr/bin/env python

from Bio.Blast.Applications import *
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import sys, argparse, os

parser = argparse.ArgumentParser(description='Run reciprocal BLAST into S.cerevisiae database.')
parser.add_argument('file', help="Input FASTA file")
#parser.add_argument('species', help="Short species identifier")
args = parser.parse_args()

cmdline = NcbiblastpCommandline(query=args.file, db="fungi_asco/Scer_aa", evalue=1, outfmt=5, out=args.file+"blast.xml")
stdout, stderr = cmdline()

spc = args.file.split("_")[0][0].upper()+args.file.split("_")[1][0:3]

seqs = list(SeqIO.parse(args.file, "fasta"))
out_seqs = []

for blast_record in NCBIXML.parse(open(args.file+"blast.xml")):
	if blast_record.alignments:
		hit_count = 1
		for alignment in blast_record.alignments:
			for seq in seqs:
				if seq.description == blast_record.query:
					seq.id = spc+alignment.hit_def.split(" ")[1][0].upper()+alignment.hit_def.split(" ")[1][1:].lower()+"|"+seq.id
					seq.description = seq.id
					out_seqs.append(seq)
					break
			
			break

SeqIO.write(out_seqs, args.file+"_out.fas", "fasta")

