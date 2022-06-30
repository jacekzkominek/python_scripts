#!/usr/bin/env python

import os, sys, argparse
from Bio import SeqIO
from collections import defaultdict
from natsort import natsorted

parser = argparse.ArgumentParser(description="Change blast to tab format.")
parser.add_argument('blast_file', help="Input file")
parser.add_argument('seq_file', help="Input seq file")
parser.add_argument('--eval', default="1E-50", help="Threshold evalue")
parser.add_argument('--bit', default=False, action="store_true", help="Threshold evalue")
parser.add_argument('--bitval', default=100, help="Threshold evalue")
args = parser.parse_args()

seq_ids = []
with open(args.seq_file) as f:
	for s in SeqIO.parse(f,"fasta"):
		seq_ids.append(s.id)

#gene_list = ["Anid_SidA","Anid_SidC","Anid_SidD","Anid_SidF","Anid_SidG","Anid_SidL","Ecol_entA","Ecol_entB","Ecol_entC","Ecol_entD","Ecol_entE","Ecol_entF","Ecol_fur","Ecol_fepE","Ecol_fepA","Ecol_fepC","Ecol_fepD","Ecol_fepG","Ecol_fepB","Ecol_entH","Ecol_entS","Ecol_fes","Scer_Arn1","Scer_Arn2","Scer_Arn3_Sit1","Scer_Arn4_Enb1","Scer_Fit1","Scer_Fit2","Scer_Fit3","Scer_Fet4","Scer_Smf1","Scer_Smf3","Scer_Ccc1","Scer_Fre6","Scer_Fet5","Scer_Fth1","Scer_Ftr1","Scer_Fet3","Scer_Fre1","Scer_Fre2","Scer_Fre3","Scer_Fre4","Scer_Fre5","Scer_Fre7","Scer_Fre8","Scer_Hmx1","Scer_Aft1","Scer_Aft2","Scer_Cth2","Umay_Urb1","Pchry_SreP","Anid_SreA","Ncra_Sre","Calb_Sfu1","Spom_Fep1"]	
gene_list = ["ScerGAL4","ScerGAL80","ScerRGT2","ScerMIG1","ScerRTG1"]


blast = []
with open(args.blast_file) as f:
	source = ""
	genes = dict((s,[]) for s in seq_ids)
	target_brh = defaultdict(list)
	#print "\t",
	#~ for g in natsorted(genes.keys()):
	for g in gene_list:
		print "\t"+g,
	print "\tTotal"
	for l in f:
		if len(l.split()) == 1:
			if source != "":
				total = 0
				print source,
				#~ for g in natsorted(genes.keys()):
				for g in gene_list:
					print " ",
					if len(genes[g]) == 0:
						print "\t0",
					else:
						count = 0
						for t in genes[g]:
							if target_brh[t[0]][0] == g:
								count += 1
								
						print "\t"+str(count),
						total += count
				print "\t"+str(total)
			source = l.strip()
			genes = dict((s,[]) for s in seq_ids)
			target_brh = defaultdict(list)
		elif len(l.split()) > 1:
			query = l.split()[0]
			target = l.split()[1]
			ev = float(l.split()[10])
			bit = float(l.split()[11])
			if args.bit == False:
				if ev < float(args.eval):
					genes[query].append([target, ev])
					if target_brh[target] == []:
						target_brh[target] = [query, ev]
					elif ev < float(target_brh[target][1]):
						target_brh[target] = [query, ev]
			else:
				if bit > float(args.bitval):
					genes[query].append([target, bit])
					if target_brh[target] == []:
						target_brh[target] = [query, bit]
					elif bit > float(target_brh[target][1]):
						target_brh[target] = [query, bit]
	
	print source,
	total = 0

	
	for g in gene_list:
	#~ for g in natsorted(genes.keys()):
		print " ",
		if len(genes[g]) == 0:
			print "\t0",
		else:
			count = 0
			for t in genes[g]:
				if target_brh[t[0]][0] == g:
					count += 1
			print "\t"+str(count),
			total += count
	print "\t"+str(total)
		
			
	#~ print source," ".join([str(genes[k]) for k in sorted(genes.keys())])

	#~ for l in f:
		#~ if len(l.split()) == 1:
			#~ if source != "":
				#~ print source," ".join([str(genes[k]) for k in sorted(genes.keys())])
			#~ source = l.strip()
			#~ genes = dict((s,0) for s in seq_ids)
		#~ elif len(l.split()) > 1:
			#~ for g in genes.keys():
				#~ if l.split()[0] == g and float(l.split()[10]) <= float(args.eval):
					#~ genes[g] += 1
			
	#~ print source," ".join([str(genes[k]) for k in sorted(genes.keys())])
