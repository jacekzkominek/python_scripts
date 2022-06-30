#!/usr/bin/env python

import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, argparse, os
from collections import defaultdict

import numpy as np
import random
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
	
parser = argparse.ArgumentParser(description="")
parser.add_argument('dir', help="Input dir")
parser.add_argument('--plot', action='store_true', default=False, help="")
parser.add_argument('--by_pw', action='store_true', default=False, help="")
parser.add_argument('--min_ct', default=1, help="Minimal amount of sequences in a file")
parser.add_argument('--keggfile', default="", help="")
args = parser.parse_args()

pp = None
if args.plot == True:
	pp = PdfPages("AA_plot.pdf")

if args.keggfile == "":
	aa_dict = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"X":0,"Y":0}
	print "species"+"\t"+"\t".join(sorted(aa_dict.keys()))

	for f in sorted(os.listdir(args.dir)):
		aa_dict = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"X":0,"Y":0}
		if f.find(".fas") != -1 and f.find(".tsv") == -1 and f[-3:] == "fas":
			seqs = list(SeqIO.parse(args.dir+"/"+f,"fasta"))
			if len(seqs) >= int(args.min_ct):
				print f,
				total = 0
				prot_aa_dist = defaultdict(list)
				for s in seqs:
					prot_aa_dict = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"X":0,"Y":0}
					prot_aa_dict_rel = defaultdict(float)
					prot_total = 0
					for r in s.seq:
						if r in aa_dict.keys():
							aa_dict[r] += 1
							total += 1
							prot_aa_dict[r] += 1
							prot_total += 1
					for aa in prot_aa_dict:
						prot_aa_dict_rel[aa] = float(prot_aa_dict[aa])/float(prot_total)
						prot_aa_dist[aa].append(prot_aa_dict_rel[aa])
							
				if args.plot == True:
					for aa in sorted(prot_aa_dist):
						bins = np.linspace(0,0.2,21)
						plt.xlim([0,0.2])
						ax = plt.gca()
						ax.xaxis.set_ticks_position("bottom")
						ax.xaxis.set_tick_params(width=2)
						plt.xticks(np.linspace(0,0.2,11))
						plt.title(f+"_"+aa)
						#print prot_aa_dist[aa]
						plt.hist(prot_aa_dist[aa], bins=bins, alpha=0.5, color="red")
						plt.legend(loc="upper left")
						plt.savefig(pp,format="pdf")
						plt.close()
						
				for a in sorted(aa_dict.keys()):
					print "\t"+str(float(aa_dict[a])/float(total)),
				print
				sys.stdout.flush()

	if args.plot == True:
		pp.close()


if args.keggfile != "":
	ref_pw = defaultdict(list)
	ref_pw_names = defaultdict(str)
	with open("/home/jacek/scripts/python/data/ko00001.keg") as kegg_ref:
		c = ""
		for l in kegg_ref:
			if l.split()[0] == "C":
				c = l.split()[1]
				ref_pw_names[c] = "_".join(l.split()[2:])
			if l.split()[0] == "D":
				ref_pw[c].append(l.split()[1])
	print "Reference KEGG loaded"
	sys.stdout.flush()
	
	kegg_pw = defaultdict(list)
	with open(args.keggfile) as keggf:
		for l in keggf:
			if len(l.split()) >= 2:
				l = l.strip()
				seq = l.split()[0]
				kegg = l.split()[1]
				for p in ref_pw:
					if kegg in ref_pw[p]:
						kegg_pw[p].append(seq)
						break
	print "User KEGG loaded"
	sys.stdout.flush()
	
	aa_dict = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"X":0,"Y":0}
	print "pathway"+"\t"+"\t".join(sorted(aa_dict.keys()))
	sys.stdout.flush()
	
	if args.by_pw == True:
		#BY PATHWAY
		for f in sorted(os.listdir(args.dir)):
			seqs = list(SeqIO.parse(args.dir+"/"+f,"fasta"))
			for pw in sorted(kegg_pw):
				total = 0
				aa_dict = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"X":0,"Y":0}
				for s in seqs:
					if s.id in kegg_pw[pw]:
						for r in s.seq:
							aa_dict[r] += 1
							total += 1
				
				print pw+"_"+ref_pw_names[pw]+"\t",
				for a in sorted(aa_dict.keys()):
					print "\t"+str(float(aa_dict[a])/float(total)),
				print
				sys.stdout.flush()
	else:
		#BY SEQ
		pw_aa = defaultdict(lambda: defaultdict)
		pw_aa_total = defaultdict(int)
		for pw in sorted(kegg_pw):
			pw_aa[pw] = {"A":0,"C":0,"D":0,"E":0,"F":0,"G":0,"H":0,"I":0,"K":0,"L":0,"M":0,"N":0,"P":0,"Q":0,"R":0,"S":0,"T":0,"V":0,"W":0,"X":0,"Y":0}
			pw_aa_total[pw] = 0
			
		for f in sorted(os.listdir(args.dir)):
			#print f
			seqs = list(SeqIO.parse(args.dir+"/"+f,"fasta"))
			#sn = 1
			for s in seqs:
			#	print str(sn)+"/"+str(len(seqs))+"\r",
			#	sn += 1
				for pw in sorted(kegg_pw):
					if s.id in kegg_pw[pw]:
						for r in s.seq:
							pw_aa[pw][r] += 1
							pw_aa_total[pw] += 1
						break

			for pw in sorted(kegg_pw):
				print pw+"_"+ref_pw_names[pw]+"\t",
				for a in sorted(pw_aa[pw]):
					print "\t"+str(float(pw_aa[pw][a])/float(pw_aa_total[pw])),
				print
				sys.stdout.flush()
