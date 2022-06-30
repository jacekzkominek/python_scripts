#!/usr/bin/env python

import os, sys, argparse
import pysam
from collections import defaultdict

parser = argparse.ArgumentParser(description="Coverage")
parser.add_argument("f1", help="Input BAM file")
parser.add_argument("f2", help="Reads")
parser.add_argument("f3", help="Outfile")
parser.add_argument("--start", default=330000, type=int, help="Start")
parser.add_argument("--end", default=350000, type=int, help="End")
args = parser.parse_args()

#CVER
#~ seq = "NODE_5_length_936291_cov_13.3144_ID_5280" 
#~ regions = [["entH_entA",336789,336815,0,0],["entA_entE",337573,337693,0,0],["entE_entC",339306,339502,0,0],["entC_entF",340689,340877,0,0],["entF_entD",344800,344983,0,0],["entD_entB",345702,345938,0,0]]
#~ genes = [["entH",336382,336789,0],["entA",336815,337573,0],["entE",337693,339306,0],["entC",339502,340689,0],["entF",340877,344800,0],["entD",344983,345702,0],["entB",345938,346780,0]]


read_list = []
with open(args.f2) as f:
	for l in f:
		read_list.append(l.strip())

bamf = pysam.AlignmentFile(args.f1, "rb")
pairedreads = pysam.AlignmentFile(args.f3, "wb", template=bamf)

seq = "NODE_5_length_936291_cov_13.3144_ID_5280"
reads = bamf.fetch(seq, args.start, args.end)
passed_reads = defaultdict(int)
insert_limit = 10000
read_list2 = defaultdict()
for read in reads:
	read_data = str(read).split()
	name = read_data[0]
	flags = int(read_data[1])
	pos1 = int(read_data[3])
	pos2 = int(read_data[7])
	if pos2-pos1 > insert_limit or read.reference_length > 1000:
		continue
	if name in read_list:
		if (flags & 64) != 0 and name not in passed_reads:
			passed_reads[name] = pos1
			read_list2[name] = read
			
		if  (flags & 128) != 0 and name in passed_reads:
			pairedreads.write(read_list2[name])
			pairedreads.write(read)
			
		

