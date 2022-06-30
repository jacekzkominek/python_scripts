#!/usr/bin/env python

import warnings
import sys, argparse, os, subprocess, shutil

parser = argparse.ArgumentParser(description="")
parser.add_argument("query", help="Query FASTA")
parser.add_argument("dir", help="Target dir")
parser.add_argument("out_dir", help="Target dir")
parser.add_argument("--cpu", default="1", help="Target dir")
args = parser.parse_args()

prog = "/home/jacek/software/gene_analysis/OrthoFinder/orthofinder/orthofinder.py"
if os.path.exists(args.out_dir) == False:
	os.makedirs(args.out_dir)
for f in sorted(os.listdir(args.dir)):
	if f.find(".fas") != -1:
		if os.path.exists(args.out_dir+"/"+f) == False:
			os.makedirs(args.out_dir+"/"+f)	
			shutil.copyfile(args.query,args.out_dir+"/"+f+"/"+args.query.split("/")[-1])
			shutil.copyfile(args.dir+"/"+f,args.out_dir+"/"+f+"/"+f)
		print f
		sys.stdout.flush()
		subprocess.call([prog,"-f",args.out_dir+"/"+f+"/","-t",args.cpu,"-a",args.cpu,"-S","diamond","-og","-I","2"])
		
		
