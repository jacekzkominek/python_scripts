#!/usr/bin/env python

import warnings
import sys, argparse, os, subprocess
from collections import defaultdict 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="")
parser.add_argument("query", help="Query FASTA")
parser.add_argument("dir", help="Target dir")
parser.add_argument("--fasta_dir", default="", help="Dir with FASTA files")
args = parser.parse_args()

queries = []
with open(args.query) as q:
	for l in q:
		if l.find(">") != -1:
			queries.append(l.split()[0][1:])

print "Input\t"+"\t".join(sorted(queries))+"\tTotal"			

q_dist = defaultdict(lambda : defaultdict(list))
out_seq_names = defaultdict(lambda : defaultdict(list))
for f in sorted(os.listdir(args.dir)):
	if f.find("proteinortho-graph") != -1:
		for q in queries:
			q_dist[f][q] = []
		with open(args.dir+"/"+f) as f1:
			t_best = defaultdict(float)
			t_best_pass = defaultdict(bool)
			for l in f1:
				if l.find("#") == -1:
					q = l.split()[0]
					t = l.split()[1]
					ab = float(l.split()[3])
					ba = float(l.split()[5])
					q_dist[f][q].append(t+"___"+str(ab+ba))
					t_best_pass[t] = False
					if ab+ba > t_best[t]:
						t_best[t] = ab+ba
			print f,
			total = 0
			for q in sorted(q_dist[f]):
				if len(q_dist[f][q]) == 0:
					print "\t0",
				else:
					ct = 0
					for t in q_dist[f][q]:
						t1 = t.split("___")[0]
						t1_abba = float(t.split("___")[1])		
						#print q,t,t_best[t1],t1_abba
						if t1_abba == t_best[t1] and t_best_pass[t1] == False:
							ct += 1
							total += 1
							t_best_pass[t1] = True
							#out_seq_names[f.split(".fas_")[1].split(".proteinortho-graph")[0]][q].append(t.split("___")[0])
							out_seq_names[f.split(".fas_")[1].split(".proteinortho-graph")[0]][q].append(t)
					print "\t"+str(ct),
				
			print "\t"+str(total)

if args.fasta_dir != "":
	out_seqs = defaultdict(list)
	for f in sorted(os.listdir(args.fasta_dir)):
		if os.path.isdir(args.fasta_dir+"/"+f) == False:
			for seq in SeqIO.parse(args.fasta_dir+"/"+f, "fasta"):
				for g in out_seq_names[f]:
					for n in out_seq_names[f][g]:
						name = n.split("___")[0]
						bit = n.split("___")[1]
						if seq.id == name:
					#if seq.id in out_seq_names[f][g]:
							s2 = seq
							s2.description = ""
							s2.name = ""
							s2.id = f+"_"+seq.id+"_"+bit
							out_seqs[g].append(s2)
	
	for g in out_seqs:
		SeqIO.write(out_seqs[g], g+"_out.fas", "fasta")