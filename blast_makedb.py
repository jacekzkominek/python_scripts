#!/usr/bin/env python

import os
import sys
import argparse
import fnmatch
import shutil

parser = argparse.ArgumentParser(description='Check CSV file for synonymous SNP data and analyze them.')
parser.add_argument('fasta_file', default="", nargs='?', help="Fasta file with sequence data. (all for all files in dir)")
parser.add_argument('--prot', default=False, action='store_true', help="Input is a protein sequence.")
parser.add_argument('--skip_ext', default=False, action='store_true', help="Do not add the sequence type extension (aa/nn)")
parser.add_argument('--keep_name', default=False, action='store_true', help="DB name equal to filename")
args = parser.parse_args()
path = os.getcwd()

if args.fasta_file != "":
	dbtype = "nucl"
	ext = "_nn"
	if args.prot == True:
		dbtype = "prot"	
		ext = "_aa"
	if args.skip_ext == True:
		ext = ""
	if args.fasta_file == "all":
		for fname in os.listdir("."):
			if os.path.isfile(fname):
				if fnmatch.fnmatch (fname, "*.fas"):
					title = fname[:-4]
					if args.keep_name == False:
						title = fname.split("_")[0][0]+fname.split("_")[1][:3]
					os.system("makeblastdb -in "+fname+" -title "+title+ext+" -out "+title+ext+" -dbtype "+dbtype)
					shutil.copyfile(fname,title+ext+".fas")
	elif os.path.isfile(args.fasta_file) == True:
		title = args.fasta_file[:-4]
		if args.keep_name == False:
			title = args.fasta_file.split("_")[0][0]+args.fasta_file.split("_")[1][:3]
		os.system("makeblastdb -in "+args.fasta_file+" -title "+title+ext+" -out "+title+ext+" -dbtype "+dbtype)
		shutil.copyfile(args.fasta_file,title+ext+".fas")
		

		
			
