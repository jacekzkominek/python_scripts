#!/usr/bin/env python

import warnings, argparse, os, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

parser = argparse.ArgumentParser(description="")
parser.add_argument('dir_fasta', help="Dir with FASTA files")
parser.add_argument('dir_gtf', help="Dir with GTF files")
parser.add_argument('dir_out', help="Output dir")
args = parser.parse_args()

if os.path.exists(args.dir_out) == False:
	os.makedirs(args.dir_out)
	
for f,g in zip(sorted(os.listdir(args.dir_fasta)),sorted(os.listdir(args.dir_gtf))):
	gtf = defaultdict(list)
	with open(args.dir_gtf+"/"+g) as gtf_file:
		colName = ""
		colChr = ""
		colStrand = ""
		colStart = ""
		colEnd = ""
		for l in gtf_file:
			if l.split()[2] == "start_codon":
				colName = l.split()[9].replace("\"","").replace(";","")
				colChr = l.split()[0]
				colStrand = l.split()[6]
				colStart = l.split()[3]
			if l.split()[2] == "stop_codon":
				colEnd = l.split()[3]
				gtf[colName.replace("-","_")] = [colChr,colStrand,colStart,colEnd]
				colName = ""
				colChr = ""
				colStrand = ""
				colStart = ""
				colEnd = ""

	seqs = list(SeqIO.parse(args.dir_fasta+"/"+f, "fasta"))
	out_seqs = []
	i = 1
	out_seqs_map = []
	for s in seqs:
		name = s.description
		name = name.split("gene=")[1].split()[0].replace("-","_")
		if name in gtf.keys():
			#print name
			colChr = gtf[name][0]
			colStrand = gtf[name][1]
			colStart = gtf[name][2]
			colEnd = gtf[name][3]
			seq2 = s
			seq2.id=""
			seq2.name=""
			seq2.description=""
			out_seqs_map.append("gene"+str(i)+"\t"+name)
			seq2.id="gene"+str(i)+" "+colChr+" "+colStart+" "+colEnd+" "+colStrand
			i+=1
			out_seqs.append(seq2)
	print f
	sys.stdout.flush()
	with open(args.dir_out+"/"+f.replace("max.pep","fasta")+"_map","w") as out_map:
		for s in out_seqs_map:
			out_map.writelines(s+"\n")
	SeqIO.write(out_seqs, args.dir_out+"/"+f.replace("max.pep","fasta"), "fasta")
	


