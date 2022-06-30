#!/usr/bin/env python

import os, sys
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description='')
parser.add_argument("dir", help="Directory with BUSCO subdirs.")
parser.add_argument("--print_only", default=False, action="store_true", help=".")
parser.add_argument("--out", default="out", help="Output directory.")
parser.add_argument("--min", default="0", help=".")
args = parser.parse_args()

busco = defaultdict(list)mi
total = 0
global_set = set()
for l in sorted(os.listdir(args.dir)):
	if os.path.isdir(args.dir+"/"+l) == True and os.path.isdir(args.dir+"/"+l+"/single_copy_busco_sequences"):
		print l,
		ct = 0
		seq_set = set()
		for f in os.listdir(args.dir+"/"+l+"/single_copy_busco_sequences"):
			if f.find("faa") != -1:
				seq_set.add(f.replace(".faa",""))
				seqs = list(SeqIO.parse(args.dir+"/"+l+"/single_copy_busco_sequences/"+f, "fasta"))
				if len(seqs) == 1:
					for s in seqs:
						#busco[f.split(".")[0]].append([f.split(".")[0]+l,str(s.seq)])
						busco[f.split(".")[0]].append([l,str(s.seq)])
						ct += 1
		diff = 0
		if len(global_set) == 0 and total == 0:
			global_set = global_set.union(seq_set)
		else:
			diff = len(global_set)-len(global_set.intersection(seq_set))
			global_set = global_set.intersection(seq_set)
		print ct,len(global_set),diff
		sys.stdout.flush()
		total += 1
print "Total",total
sys.stdout.flush()
print "Shared by",total,"taxa",len(global_set)
sys.stdout.flush()

if args.print_only == False:
	count_min = total
	if args.min != "0":
		if args.min.find("%") != -1:
			count_min = int(float(float(args.min.replace("%",""))/float(100))*float(total))
		elif args.min.find("-") != -1:
			count_min = int(total)-int(args.min.replace("-",""))
		else:
			count_min = int(args.min)
	print "Extracting shared by",int(count_min),"taxa",
	if os.path.exists(args.out) == False:
		os.makedirs(args.out)
	min_ct = 0
	for b in busco.keys():
		if len(busco[b]) >= count_min:
			min_ct += 1
		#if len(busco[b]) == total:
			with open(args.out+"/"+str(len(busco[b]))+"_"+b+".fas","w") as busco_out:
				for f in busco[b]:
					busco_out.write(">"+f[0]+"\n"+f[1]+"\n")
	print min_ct
