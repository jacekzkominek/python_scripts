#!/usr/bin/env python

from Bio.Blast.Applications import *
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import SearchIO
import sys, argparse, os, subprocess, shutil
from collections import defaultdict
import copy, string

parser = argparse.ArgumentParser(description='Run BLAST over the local database.')
parser.add_argument('file', help="Input FASTA file")
parser.add_argument('dir', default=".", help="Directory with fasta files.")
parser.add_argument('--outdir', default="", help="Output subdirectory.")
parser.add_argument('--blasthits', default=1)
parser.add_argument('--igs', default=500)
parser.add_argument('--fileout', default=False, action="store_true")
parser.add_argument('--evalue', default=1e-5)
parser.add_argument('--numal', default=1)
parser.add_argument('--short', default=False, action="store_true")
parser.add_argument('--blastp', default=False, action="store_true")
parser.add_argument('--output_all', default=False, action="store_true")

args = parser.parse_args()
if args.outdir != "" and not os.path.exists(args.outdir):
	os.makedirs(args.outdir)
if args.outdir == "":
	args.outdir = os.getcwd()

fnames = []
once = defaultdict(bool)
for dirpath,dirnames,filenames in os.walk(args.dir):
	for f in sorted(filenames):
		fnames.append(f)
		if args.blastp == False:
			cmdline = NcbitblastnCommandline(query=args.file, subject=os.path.join(dirpath,f), evalue=float(args.evalue), outfmt=5, num_alignments=int(args.numal))#, out=args.outdir+"/"+f+"_out")
		else:
			cmdline = NcbiblastpCommandline(query=args.file, subject=os.path.join(dirpath,f), evalue=float(args.evalue), outfmt=5, num_alignments=int(args.numal))#, out=args.outdir+"/"+f+"_out")
		
		#cmdline = NcbiblastpCommandline(query=args.file, subject=os.path.join(dirpath,f), evalue=float(args.evalue), outfmt=5, num_alignments=int(args.numal), out=args.outdir+"/"+f+"_out")
		outfile = subprocess.PIPE
		#outfile = None
		#if args.fileout == True:
		#	outfile = open(args.outdir+"/"+f+"_out","w")
			#print cmdline
			#print str(cmdline).split()
		#child = subprocess.Popen(str(cmdline),stdout=None,stderr=None, stdin=None, shell=True)
		child = subprocess.Popen(str(cmdline),stdout=outfile,stderr=subprocess.PIPE,universal_newlines=True,shell=(sys.platform!="win32"))
		
		#output = None
		output = child.stdout
		#stdout, stderr = cmdline()
		#if args.fileout == True:
		#	outfile.close()
		#	output = open(args.outdir+"/"+f+"_out")	
		
		#print output
		print f
		out_seqs=[]
		out_seqs_hits=set()
		blast_hits=0
		# ~ sp_set = set()
		q_set = set()
		seq_len = len(list(SeqIO.parse(args.file, "fasta")))
		
		for qresult in SearchIO.parse(output, 'blast-xml'):
			for hit in qresult.hits:
				h = hit.description
				q = qresult.id
				start = 0
				end = 0
				ev = 0
				# ~ print hit
				for hsp in hit.hsps:
					if hsp.hit_strand == 1:
						start = hsp.hit_start
						end = hsp.hit_end
					elif hsp.hit_strand == -1:
						start = hsp.hit_end
						end = hsp.hit_start
					ev = hsp.evalue
					break
				print q,h,start,end,ev
				sys.stdout.flush()
				if len(out_seqs) >= 1 and q == out_seqs[-1].split("___")[0] and float(ev) < float(out_seqs[-1].split("___")[5]):
					out_seqs = out_seqs[:-1]
				sp = f.split("_")[0]+"_"+f.split("_")[1]
				if f.split("_")[0].find("yHMPu") != -1 or f.split("_")[0].find("yHAB") != -1:
					sp = f.split("_")[0]+"_"+f.split("_")[1]+"_"+f.split("_")[2]
				out_seqs.append(q+"___"+sp+"___"+h+"___"+str(start)+"___"+str(end)+"___"+str(ev))
				out_seqs_hits.add(h)
				blast_hits+=1
				q_set.add(q)
		
		
		for q in q_set:
			once[q] = True
					
		# ~ print out_seqs
		# ~ print blast_hits
		sys.stdout.flush()
		
		
		#EXTRACT SEQUENCES
		for inseq in SeqIO.parse(os.path.join(dirpath,f), "fasta"):
			total = 0
			for s in out_seqs:
				q = s.split("___")[0]
				sp = s.split("___")[1]
				n = s.split("___")[2]
				start = int(s.split("___")[3])
				end = int(s.split("___")[4])
				
				if inseq.description == n:
				# ~ if inseq.description.find(n) != -1:
					desc = inseq.description
					newseq = copy.deepcopy(inseq)
					atg = 0
					stop = 0
					if args.blastp == True:
						atg = 1
						stop = len(str(newseq.seq))
					else:
						if start < end:
							#newseq.seq = inseq.seq[start:end]
							inc = -5
							while (atg == 0):
								cod = inseq.seq[max(0,start-(inc*3)):start-(inc*3)+3].lower()
								if cod == "atg":
									atg = start-(inc*3)
								if start-(inc*3) < 1:
									atg = 1
								inc += 1
							
							inc = -5
							while (stop == 0):
								cod = inseq.seq[end+(inc*3):end+(inc*3)+3].lower()
								if cod == "taa" or cod == "tga" or cod == "tag":
									stop = end+(inc*3)
								if end+(inc*3)+3 > len(inseq.seq):
									stop = len(inseq.seq)
								inc += 1
							
							newseq.seq = inseq.seq[atg:stop+3]
							
						elif start > end:
							inc = -5
							while (atg == 0):
								cod = inseq.seq[start+(inc*3):start+(inc*3)+3].reverse_complement().lower()
								if cod == "atg":
									atg = start+(inc*3)
								if start+(inc*3) > len(inseq.seq):
									atg = len(inseq.seq)
								inc += 1
		
							inc = -5
							while (stop == 0):
								cod = inseq.seq[end-(inc*3):end-(inc*3)+3].reverse_complement().lower()
								if cod == "taa" or cod == "tga" or cod == "tag":
									stop = end-(inc*3)
								if end-(inc*3) < 1:
									stop = 1
								inc += 1
							
							newseq.seq = inseq.seq[stop:atg+3].reverse_complement()
					sp = sp.replace(".fas","")	
					if args.short == True:
						newseq.id = q+"_"+sp
					else:
						newseq.id = q+"_"+sp+"_"+n+"_"+str(atg)+"_"+str(stop)	
					newseq.title = ""
					newseq.description = ""
					
					newseq_igs = copy.deepcopy(inseq)
					if start < end:
						newseq_igs.seq = inseq.seq[atg-int(args.igs):atg]
					elif start > end:
						newseq_igs.seq = inseq.seq[atg:atg+int(args.igs)].reverse_complement()
					
					if args.short == True:
						newseq_igs.id = q+"_"+sp+"_IGS"+str(args.igs)
					else:
						newseq_igs.id = q+"_"+sp+"_IGS"+str(args.igs)+"_"+n+"_"+str(atg)+"_"+str(stop)
					newseq_igs.title = ""
					newseq_igs.description = ""
					
					outfile = None
					#if once[q] == True:
					#	outfile = open(args.outdir+"/"+q+".fas","w")
					#	outfile.close()
					#	outfile = open(args.outdir+"/"+q+".fas","w")
					#	once[q] = False
					#else:	
					#	outfile = open(args.outdir+"/"+q+".fas","a+")
					if args.output_all == False and (len(newseq) == 0 or len(newseq_igs) == 0):
						continue
						
					outfile = open(args.outdir+"/"+q+".fas","a+")
					SeqIO.write([newseq],outfile,"fasta")
					outfile.close()
					
					if len(newseq_igs.seq) > 0:
						outfile = open(args.outdir+"/"+q+"_igs"+str(args.igs)+".fas","a+")
						SeqIO.write([newseq_igs],outfile,"fasta")
						outfile.close()
						outfile2 = open(args.outdir+"/all_igs"+str(args.igs)+".fas","a+")
						SeqIO.write([newseq_igs],outfile2,"fasta")
						outfile2.close()
					total+=1
				#if total == blast_hits:
					#break
		
		sys.stdout.flush()

