#!/usr/bin/env python

import os, sys, argparse
import pysam 
from collections import defaultdict

parser = argparse.ArgumentParser(description="Coverage")
parser.add_argument("f1", help="Input BAM file")
parser.add_argument("--start", default=368000, help="Start")
parser.add_argument("--end", default=381500, help="End")
parser.add_argument("--print_cov", default=False, action="store_true", help="End")
parser.add_argument("--print_reads", default=False, action="store_true", help="End")
args = parser.parse_args()

bamf = pysam.AlignmentFile(args.f1, "rb")
reads = bamf.fetch()
all_reads = []
all_reads_names = []
softclip_reads = []
softclip_reads_names = []
for read in reads:
	read_data = str(read).split()
	name = read_data[0]
	flags = int(read_data[1])
	cigar = read_data[5]
	seq1 = read_data[9]
	if name not in all_reads_names:
		all_reads.append(read)	
		all_reads_names.append(name)	
		#print seq1[:5]
		if (seq1[:5] == "AAAAA" or seq1[:5] == "TTTTT"):
			if cigar.find("S") != -1 and unicode(cigar.split("S")[0]).isnumeric() == True:
				softclip = int(cigar.split("S")[0])
				if softclip >= 5:
					if args.print_reads == True:
						print read
						sys.stdout.flush()
					softclip_reads.append(read)
		elif (seq1[-5:] == "AAAAA" or seq1[-5:] == "TTTTT"):
			if cigar.find("S") != -1:
				softclip = 0
				if unicode(cigar[cigar.find("S")-3:cigar.find("S")]).isnumeric() == True:
					softclip = int(cigar[cigar.find("S")-3:cigar.find("S")])
				elif unicode(cigar[cigar.find("S")-2:cigar.find("S")]).isnumeric() == True:
					softclip = int(cigar[cigar.find("S")-2:cigar.find("S")])
				elif unicode(cigar[cigar.find("S")-1:cigar.find("S")]).isnumeric() == True:
					softclip = int(cigar[cigar.find("S")-1:cigar.find("S")])
				if softclip >= 5:
					softclip_reads.append(read)
					if args.print_reads == True:
						print read
						sys.stdout.flush()

print "TOTAL READS",len(all_reads)
print "SOFTCLIP READS",len(softclip_reads)
