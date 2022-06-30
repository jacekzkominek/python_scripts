#!/usr/bin/env python

import warnings, argparse, os, sys, subprocess, shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

parser = argparse.ArgumentParser(description="")
parser.add_argument('dir_fasta', help="Dir with input FASTA files")
parser.add_argument('dir_out', help="")
args = parser.parse_args()

startdir = os.getcwd()
if os.path.exists("/home/jacek/software/genome_analysis/CHROnicle/"+args.dir_out) == False:
	os.makedirs("/home/jacek/software/genome_analysis/CHROnicle/"+args.dir_out)
	os.makedirs("/home/jacek/software/genome_analysis/CHROnicle/"+args.dir_out+"/00rawGenom/")
for f in sorted(os.listdir(args.dir_fasta)):
	if os.path.exists("/home/jacek/software/genome_analysis/CHROnicle/"+args.dir_out+"/"+f.replace(".fasta","")) == False:
		os.makedirs("/home/jacek/software/genome_analysis/CHROnicle/"+args.dir_out+"/00rawGenom/"+f.replace(".fasta",""))
		shutil.copyfile(startdir+"/"+args.dir_fasta+"/"+f,"/home/jacek/software/genome_analysis/CHROnicle/"+args.dir_out+"/00rawGenom/"+f.replace(".fasta","")+"/"+f)
	os.chdir("/home/jacek/software/genome_analysis/CHROnicle/Programs/0Convert2InputF/")
	subprocess.call(["./ConvertFasta.py",args.dir_out,f.replace(".fasta",""),f[0]+"".join(f.split("_")[1:3]),str(1),str(2),str(3),str(4),str(5),str(0)])
	
	
	


