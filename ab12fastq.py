#!/usr/bin/env python

import os
import fnmatch
from Bio import SeqIO
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='Convert .ab1 files to .fastq files.')
parser.add_argument('input_pattern', default="", nargs='?', help="Input file or pattern.")
args = parser.parse_args()

if args.input_pattern != "":
		pattern = args.input_pattern
		print "Finding all files with",pattern,"in their filename"
		for fname in os.listdir("."):
			if os.path.isfile(fname): 
				if fnmatch.fnmatch (fname, "*"+pattern+"*"):
					print fname
					SeqIO.convert(fname,"abi",fname[:-4]+".fq", "fastq")
		
			
			
