#!/usr/bin/env python

import warnings
import sys, argparse, os, subprocess

parser = argparse.ArgumentParser(description="")
parser.add_argument("query", help="Query FASTA")
parser.add_argument("dir", help="Target dir")
parser.add_argument("--tblastn", default=False, action="store_true", help="Use TBLASTN instead of BLASTP")
parser.add_argument("--blastn", default=False, action="store_true", help="Use BLASTN instead of BLASTP")
parser.add_argument("--pi", default=0, help="Percent identity threshold")
parser.add_argument("--evalue", default="1E-5", help="Evalue")
args = parser.parse_args()

for f in sorted(os.listdir(args.dir)):
	if f.find(".fas") != -1 and f.find(".tsv") == -1 and f.find(".phr") == -1 and f.find(".pin") == -1 and f.find(".psq") == -1:
		print f
		sys.stdout.flush()
		prog = "/home/jkominek/software/alignment/ncbi-blast-2.11.0+/bin/blastp"
		if args.tblastn == True:
			prog = "/home/jkominek/software/alignment/ncbi-blast-2.11.0+/bin/tblastn"
		if args.blastn == True:
			prog = "/home/jkominek/software/alignment/ncbi-blast-2.11.0+/bin/blastn"
		subprocess.call([prog,"-query",args.query,"-evalue",args.evalue,"-outfmt","6","-subject",args.dir+"/"+f])
		#subprocess.call([prog,"-query",args.query,"-perc_identity",str(args.pi),"-evalue",args.evalue,"-outfmt","6","-subject",args.dir+"/"+f])

