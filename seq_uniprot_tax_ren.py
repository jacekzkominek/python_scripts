#!/usr/bin/env python

from Bio import SeqIO
import sys, argparse
from collections import defaultdict 

parser = argparse.ArgumentParser(description="Crossref Uniprot sequences and output only those present in all species")
parser.add_argument('f1', help="Input FASTA file")
parser.add_argument('name', default="")
parser.add_argument('--tax', action="store_true", default=False, help="")
args = parser.parse_args()

species = set()
multi = list()
s = "Tax="
e = "RepID="
out = []
names = []
for seq in SeqIO.parse(args.f1, "fasta"):
	d = seq.description
	if d.find(s) != -1:
		sp = d[d.find(s)+len(s):d.find(e)]
		if len(sp.split()) >= 2 and sp.find(" sp.") == -1:
			id_num = d.split()[0].split("_")[1]
			seq.title = ""
			seq.id = ""
			seq.description = ""
			seq.id = sp.strip().replace(" ","_").replace("=","_").replace("-","_")+"_"+args.name+"_"+id_num
			if sp.find("(strain") != -1:
				seq.id = sp.strip()[:sp.find("(strain")-1].replace(" ","_").replace("=","_").replace("-","_")+"_"+args.name+"_"+id_num
				names.append(sp.strip()[:sp.find("(strain")])
			else:
				names.append(sp.strip())
			out.append(seq)
				
SeqIO.write(out, args.f1.replace(".fas","_ren.fas"), "fasta")

if args.tax == True:
	for n in sorted(names):
		print n
