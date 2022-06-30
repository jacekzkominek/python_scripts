#!/usr/bin/env python

import warnings
import sys, argparse, os, subprocess

parser = argparse.ArgumentParser(description="")
parser.add_argument("query", help="Query FASTA")
parser.add_argument("dir", help="Target dir")
parser.add_argument("--cpu", default="1", help="Target dir")
args = parser.parse_args()

prog = "/home/jacek/software/gene_analysis/proteinortho_v5.16/proteinortho5.pl"
blast = "/home/jacek/software/alignment/ncbi-blast-2.6.0+/bin/"
for f in sorted(os.listdir(args.dir)):
	if f.find(".fas") != -1 and f.find(".tsv") == -1 and f.find(".fas.pin") == -1 and f.find(".fas.psq") == -1 and f.find(".fas.phr") == -1:
		print f
		sys.stdout.flush()
		subprocess.call([prog,args.query,args.dir+"/"+f,"-cpus="+args.cpu,"-project="+args.query+"_"+f,"-keep","-clean","-singles","-e=1e-5","-sim=0.50","-cov=25","-p=blastp+","-blastpath="+blast])
