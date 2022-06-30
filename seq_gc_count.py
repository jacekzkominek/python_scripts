#!/usr/bin/env python

import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, argparse, os
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import random
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from PyPDF2 import PdfFileMerger, PdfFileReader

def chunks(seq, win, step):
	seqlen = len(seq)
	for i in range(0,seqlen,step):
		j = seqlen if i+win>seqlen else i+win
		yield seq[i:j]
		if j==seqlen: break

def gc_count(f):
	nuc_glob = 0
	gc_glob = 0
	if f.find(".fas") != -1 and f.find(".tsv") == -1 and f[-3:] == "fas":
		seqs = list(SeqIO.parse(args.dir+"/"+f,"fasta"))
		seqs1 = sorted(seqs, key=lambda x: len(x), reverse=True)
		
		total_size = 0
		total_size_all = 0
		for s in seqs1:
			total_size_all += len(s.seq)
			if len(s.seq) >= int(args.min_contig):
				total_size += len(s.seq)
		if total_size >= int(args.total):
			n50 = 0
			l50 = 0
			small_size = 0
			for s in seqs1:
				# ~ if len(s.seq) >= int(args.min_contig):
				small_size += len(s.seq)
				n50 = len(s.seq)
				l50 += 1
				if float(small_size) >= (float(float(total_size_all)/float(2))):
					break
			if n50 > int(args.n50):
				outfile = open(os.getcwd()+"/"+args.outdir+"/"+f+"_gc.txt","w")
				for win in sorted(args.window, key=lambda x: int(x)):
					
					slide = args.window_slide
					if args.window_slide == "":
						slide = win
										
					print f,"calculating GC content in",win+"bp","windows sliding",slide+"bp in contigs of minimum",args.min_contig,"bp for assemblies of at least",args.n50,"N50 and",args.total,"bp in total size"
					sys.stdout.flush()
					seqs2 = sorted(seqs, key=lambda x: len(x), reverse=True)
					
					for s in seqs2:					
						if len(s.seq) >= int(args.min_contig):
							w_gc = []
							contig_windows = chunks(s.seq,int(win),int(slide))				
							for w in contig_windows:
								gc = 0
								nuc = 0
								if args.window_full_only == False or len(w) == int(win):
									for nn in w:
										if nn.upper() != "N":
											nuc += 1
											if nn.upper() == "G" or nn.upper() == "C":
												gc += 1
									if nuc > 0:
										gc_pc = float(float(gc)/float(nuc))
										w_gc.append(gc_pc)
										#gc_stats[f].append(gc_pc)
										nuc_glob += nuc
										gc_glob += gc
							gc_stats[f+win].append(w_gc)
							
					if nuc_glob > 0 and gc_glob > 0:
						if args.window_print == False:
							all_gc = []
							for w_gc in gc_stats[f+win]:
								for w in w_gc:
									all_gc.append(w)
							outfile.write(f+"\t"+win+"\t"+slide+"\t"+args.min_contig+"\t"+str(n50)+"\t"+"\t".join([str(x) for x in all_gc])+"\n")
						else:
							outfile.write(f+"\t"+win+"\t"+slide+"\t"+args.min_contig+"\t"+str(n50)+"\t")
							for w_gc in gc_stats[f+win]:
								outfile.write("\t"+"_".join([str(x) for x in w_gc]))
							outfile.write("\n")
				outfile.close()

parser = argparse.ArgumentParser(description="")
parser.add_argument("dir", help="Input dir")
parser.add_argument("outdir", default="", help="")
parser.add_argument("--window", nargs="+", default=["1000"], help="")
parser.add_argument("--window_full_only", action="store_true", default=False, help="")
parser.add_argument("--window_print", action="store_true", default=False, help="")
parser.add_argument("--window_slide", default="", help="")
parser.add_argument("--min_contig", default="1000", help="")
parser.add_argument("--n50", default="20000", help="")
parser.add_argument("--total", default="1", help="")
parser.add_argument("--cpu", default="18", help="")
args = parser.parse_args()

gc_stats = defaultdict(list)
if os.path.exists(os.getcwd()+"/"+args.outdir) == False:
	os.makedirs(os.getcwd()+"/"+args.outdir)
	
files = sorted(os.listdir(args.dir))

p = Pool(processes=int(args.cpu))
data = []
data = p.map(gc_count, files)
p.close()
