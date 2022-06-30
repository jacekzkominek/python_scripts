#!/usr/bin/env python

import warnings, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from os import sys, close, getcwd
import os
import numpy
from collections import defaultdict

parser = argparse.ArgumentParser(description="Calculate identity and similarity statistics for an alignment.")
parser.add_argument('f1', help="Input FASTA file")
parser.add_argument('--blosum', default="62", help="BLOSUM matrix to use")
parser.add_argument('--id', action="store_true", default=False, help="Only print identity information")
#parser.add_argument('--ref_id', default=0, help="Sequence to use as reference (first by default, 0-based)")
#parser.add_argument('--print_pos_id', type=int, default=0, help="Print positions with identity level.")
args = parser.parse_args()

path = os.path.dirname(os.path.realpath(__file__))
blosum_f = open(path+"/blosum/blosum"+args.blosum+".blo")
blosum = defaultdict(int)
blosum_ref = ""
for line in blosum_f:
	if line.find('#') != -1:
		continue
	else:
		row = line.split()
		if len(row) == 24:
			blosum_ref = row
		elif len(row) == 25:
			for i, score in enumerate(row):
				if i > 0:
					blosum[blosum_ref[i-1]+row[0]] = int(score)
					#print blosum_ref[i-1],row[0],score
blosum_f.close()

f = open(args.f1)
seqs = list(SeqIO.parse(f, "fasta"))
seq_count = float(len(seqs))
seq_len = len(seqs[0])

id1_scores = []
id2_scores = []
id3_scores = []
sim1_scores = []
sim2_scores = []
sim3_scores = []
for i, seq in enumerate(seqs):
	seq1 = seqs[i].seq
	for j in range(i+1,int(seq_count)):
		seq2 = seqs[j].seq
		id_count = 0
		sim_count = 0
		total_len = 0
		aln_len = 0
		for x in range (0, seq_len):
			seq1_res = seq1[x]
			seq2_res = seq2[x]
			if seq1_res == "-" and seq2_res == "-":
				continue
			elif seq2_res != "-" and seq2_res != "-":
				aln_len += 1
				total_len += 1
				if seq1_res == seq2_res:
					id_count += 1
				if blosum[seq1_res+seq2_res] > 0:
					sim_count += 1
			elif (seq1_res != "-" and seq2_res == "-") or (seq1_res == "-" and seq2_res != "-"):
				total_len += 1
				
		id1_scores.append(100*(float(id_count)/float(total_len)))
		id2_scores.append(100*(float(id_count)/float(aln_len)))
		id3_scores.append(100*(float(id_count)/float((total_len+aln_len)/2)))
		sim1_scores.append(100*(float(sim_count)/float(total_len)))
		sim2_scores.append(100*(float(sim_count)/float(aln_len)))
		sim3_scores.append(100*(float(sim_count)/float((total_len+aln_len)/2)))
		#print i,j

print "Pairwise sequence identity for "+args.f1+": "+"%.2f%% +/- %.2f%% (over total length)" % (numpy.mean(id1_scores),numpy.std(id1_scores))
print "Pairwise sequence identity for "+args.f1+": "+"%.2f%% +/- %.2f%% (over aligned length)" % (numpy.mean(id2_scores),numpy.std(id2_scores))
print "Pairwise sequence identity for "+args.f1+": "+"%.2f%% +/- %.2f%% (over mean length)" % (numpy.mean(id3_scores),numpy.std(id3_scores))
if args.id == False:
	print "Pairwise sequence similarity for "+args.f1+": "+"%.2f%% +/- %.2f%% (over total length)" % (numpy.mean(sim1_scores),numpy.std(sim1_scores))
	print "Pairwise sequence similarity for "+args.f1+": "+"%.2f%% +/- %.2f%% (over aligned length)" % (numpy.mean(sim2_scores),numpy.std(sim2_scores))
	print "Pairwise sequence similarity for "+args.f1+": "+"%.2f%% +/- %.2f%% (over mean length)" % (numpy.mean(sim3_scores),numpy.std(sim3_scores))


	
			
#pos_id = []
#pos_sim = []

#f.seek(0)
#for x in range (0, seq_len):
	#first_res = ""
	#current_res = ""
	#non_sim_count = 0
	#for i, seq in enumerate(list(SeqIO.parse(f, "fasta"))):
		#if i == 0:
			#first_res = seq[x]
		#else:
			#current_res = seq[x]
			#if first_res == current_res and first_res != "-":
					#id_score
			
			#col = blosum[0].index(first_res)+1
			#sim_score = 0
			
			##RELAXED GAP PROCESSING
			#if current_res.find ('-') != -1:
				#sim_score = -1000
			#else:
				#for y in range (1,len(blosum)):
					#if blosum[y][0].find(current_res) != -1:
						#sim_score = int(blosum[y][col])
						#break
			#if sim_score <= 0:
				#non_sim_count += 1
				#if non_sim_count > threshold:
					#f.seek(0)				
					#break
			#if i == seq_count-1:
				#sim_pos_list.append(x+1)
				#f.seek(0)
#f.close()

#print
#print "Positions with similarity level within allowed threshold (%d/%d, %.2f%%):" % (int(threshold), int(seq_count), 100*(float(seq_count)-float(threshold))/float(seq_count))
#print str(sim_pos_list).strip('[]')
#print

