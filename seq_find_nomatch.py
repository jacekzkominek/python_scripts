#!/usr/bin/env python

import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from os import sys, close
from collections import defaultdict

csv = open(sys.argv[1])
seqs = list(SeqIO.parse(sys.argv[2], "fasta"))
nomatch_list = []
nomatch_seqs = []

for a in csv:
	if a.split(",")[4] == "":
		nomatch_list.append(a.split(",")[2][1:-1])

print len(nomatch_list)	

a = 1
nf = 0
for l in nomatch_list:
	print a,l,
	found = False
	for seq in seqs:
		found = False
		if seq.id.find(l) != -1:
			print "\t"+seq.id
			nomatch_seqs.append(seq)
			found = True
			break
	if found == False:
		print "\t"+"---===SEQ NOT FOUND===---"
		nf+=1	
	a+=1
		
#print nf

SeqIO.write(nomatch_seqs, sys.argv[1]+"_nom.fas", "fasta")

csv.close()


