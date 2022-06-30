#!/usr/bin/env python

import os
import sys
import argparse
import urllib
import subprocess
import random
import string
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from subprocess import Popen, PIPE, STDOUT

parser = argparse.ArgumentParser(description='Check CSV file for synonymous SNP data and analyze them.')
parser.add_argument('csv_file', default="", nargs='?', help="CSV file with input data.")
parser.add_argument('--codon_pos', default="100000", type=int, nargs='?', help="Only include SNPs at positions before or at this codon.")
parser.add_argument('--nonsyn', default=False, action='store_true', help="Analyse nonsynonymous SNPs as well.")
parser.add_argument('--fold_ss', nargs='?', default = "", help="Fold mRNA secondary structure and calculate dG (mfold,vienna,unafold,rnastruct).")
args = parser.parse_args()
path = os.getcwd()


#LOAD CAI VALUES FOR CODONS
codon_file = open("/home/jacek/scripts/python/codons_w.csv","r")
codon = defaultdict(list)
for l in codon_file:
	l1 = l.split(",")
	codon[l1[0]]=[l1[1],l1[2].rstrip()]
codon_file.close()

#LOAD GENE START/END DATA FROM SGD, FIRST USING PRIMARY GENE NAMES, THEN USING ALIASES
sgd_file = open("/home/jacek/scripts/python/sgd_gff.tsv","r")
sgd_genes = defaultdict(list)

alias = False
for l in sgd_file:
	l1 = l.split("\t")

	if l1 == "###" and alias == False:
		alias = True
	elif l1 == "###" and alias == True:
		break

	if len(l1) >= 8 and l1[2] == "gene":

		gene_names = []
		gene_start = int(l1[3])
		gene_end = int(l1[4])
		gene_strand = l1[6]

		if alias == False:
			gene_names.append(urllib.unquote(l1[10][5:]))
			gene_names.append(urllib.unquote(l1[9][5:]))
		elif alias == True:
			if l1[11].find("Alias=") != -1:
				for alias in l1[11][6:].split(","):
					gene_names.append(alias)

		exon_struct = []
		while True:
			l2 = sgd_file.next()
			if l2.split("\t")[2] == "CDS":
				exon_start = int(l2.split("\t")[3])
				exon_end = int(l2.split("\t")[4])
				exon_struct.append([exon_start,exon_end])
			else:
				break

		for name in gene_names:
			if sgd_genes[name] == []:
				sgd_genes[name]=[gene_strand,exon_struct]

sgd_file.close()

f = open(args.csv_file, "rb")
snps = []
sgd_seqs = None
seqs = []

if args.fold_ss == "":
	print "#\tChr\tPos\tStrand\tGene\tCodon\tRef-AA\tRef-cod\tRef-CAI\tMut-AA\tMut-cod\tMut-CAI"
else:
	print "#\tChr\tPos\tStrand\tGene\tCodon\tRef-AA\tRef-cod\tRef-CAI\tMut-AA\tMut-cod\tMut-CAI\tRef-dG\tMut-dG\tdG(R-M)"
	sgd_seqs = open("/home/jacek/scripts/python/sgd_coding.fas","r")
	tmp_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
	d = os.path.dirname(path+"/tmp-"+args.fold_ss+"-"+tmp_string+"/")
	if not os.path.exists(d):
		os.makedirs(d)
	os.chdir(d)

i = 1
for l in f:
	if l != "":
		l1 = l.split(",")

		if l1[5] == "\"Synonynous\"" or (args.nonsyn == True and l1[5] == "\"Non-synonymous\""):
			ref = l1[7][1:-1].upper().replace("T","U")
			mut = l1[8][1:-1].upper().replace("T","U")
			gene_name = l1[6][1:-1]
			snp_codon = 0

			if sgd_genes[gene_name] != []:
				gene_start = 0
				snp_pos = int(l1[1])

				#PARSE EXON DATA TO IDENTIFY THE CODON OF THE SNP
				if sgd_genes[gene_name][0] == "+":
					length = 0
					for exon in sgd_genes[gene_name][1]:
						if snp_pos < exon[1]:
							snp_codon = length+1+(snp_pos-exon[0])/3
							break
						else:
							length = length+1+(exon[1]-exon[0])/3
				elif sgd_genes[gene_name][0] == "-":
					length = 0
					for exon in reversed(sgd_genes[gene_name][1]):
						if snp_pos > exon[0]:
							snp_codon = length+1+(exon[1]-snp_pos)/3
							break
						else:
							length = length+1+(exon[1]-exon[0])/3
			else:
				snp_codon = 0

			if snp_codon <= args.codon_pos:
				#PRINT BASIC OUTPUT FIRST, IN CASE WE NEED TO WAIT FOR MRNA FOLDING
				gene_chr = l1[0][1:-1]
				snp_pos = l1[1]
				print str(i)+"\t"+gene_chr+"\t"+snp_pos+"\t"+sgd_genes[gene_name][0]+"\t"+gene_name+"\t"+str(snp_codon)+"\t"+codon[ref][0]+"\t"+ref+"\t"+codon[ref][1]+"\t"+codon[mut][0]+"\t"+mut+"\t"+codon[mut][1],
				sys.stdout.flush()
				if args.fold_ss == "":
					print
				else:
				#PERFORM MRNA FOLDING PREDICTIONS, CALCULATE ENERGY DIFFERENCES BETWEEN REFERENCE AND MUTANT, AND TRY PACKAGE-SPECIFIC COMPARISON MEASURES
					print "\t",
					sys.stdout.flush()
					sgd_seqs.seek(0)
					seqs = SeqIO.parse(sgd_seqs,"fasta")

					for seq in seqs:
						if seq.description.split()[1] == gene_name:
							mrna_seq = str(seq.seq).replace("T","U")
							dG = []
							if mrna_seq[(snp_codon-1)*3:((snp_codon-1)*3)+3] == ref:

								mrna_ref_seq_rec = SeqRecord(Seq(mrna_seq,IUPAC.unambiguous_rna),id=gene_name,description="")
								SeqIO.write(mrna_ref_seq_rec, gene_name+".fas", "fasta")

								mrna_mut_seq = mrna_seq[:snp_codon*3]+mut+mrna_seq[(snp_codon+1)*3:]
								mrna_mut_seq_rec = SeqRecord(Seq(mrna_mut_seq,IUPAC.unambiguous_rna),id=gene_name+"_mut",description="")
								SeqIO.write(mrna_mut_seq_rec, gene_name+"_mut.fas", "fasta")

								dG = 0
								dG_mut = 0
								#USE MFOLD V3.6
								if args.fold_ss == "mfold":
									subprocess.call(["mfold","SEQ="+gene_name+".fas","RUN_TYPE=txt","MAX=1"],stdout=open(gene_name+".out","w+"))
									out = open(gene_name+".out","r")
									for l in out:
										if l.find("Minimum folding energy is") != -1:
											dG = float(l.split(" ")[4])
											break

									subprocess.call(["mfold","SEQ="+gene_name+"_mut.fas","RUN_TYPE=txt","MAX=1"],stdout=open(gene_name+"_mut.out","w+"))
									out = open(gene_name+"_mut.out","r")
									for l in out:
										if l.find("Minimum folding energy is") != -1:
											dG_mut = float(l.split(" ")[4])
											break

									subprocess.call(["ct_compare",gene_name+".fas.ct",gene_name+"_mut.fas.ct"],stdout=open(gene_name+"_cmp.out","w+"))
									out = open(gene_name+"_cmp.out","r")
									out_hx_bp = ""
									for l in out:
										if len(l.split()) == 14:
											hx_procent = float(l.split()[3])/float(l.split()[6])
											bp_procent = float(l.split()[8])/float(l.split()[11])
											out_hx_bp = "hx: %.2f%% bp: %.2f%%" % (hx_procent*100, bp_procent*100)
											
											break


									print str(dG)+"\t"+str(dG_mut)+"\t"+str(dG-dG_mut)+" ("+out_hx_bp+")"
								
								#USE VIENNA RNA PACKAGE
								elif args.fold_ss == "vienna":
									subprocess.call(["RNAfold"],stdin=open(gene_name+".fas","r"), stdout=open(gene_name+".out","w+"))
									out = open(gene_name+".out","r")

									for l in out:
										if l.find(")\n") != -1:
											dG = float(l[l.rfind("(")+1:-2])
											rna_ss = l.split(" ")[0]
											break

									subprocess.call(["RNAfold"],stdin=open(gene_name+"_mut.fas","r"), stdout=open(gene_name+"_mut.out","w+"))
									out = open(gene_name+"_mut.out","r")

									for l in out:
										if l.find(")\n") != -1:
											dG_mut = float(l[l.rfind("(")+1:-2])
											rna_mut_ss = l.split(" ")[0]
											break

									p = Popen(["RNAdistance","-DfhwcFHWC"],stdin=PIPE, stdout=PIPE, stderr=PIPE)
									out = p.communicate(input=rna_ss+"\n"+rna_mut_ss)[0]

									print str(dG)+"\t"+str(dG_mut)+"\t"+str(dG-dG_mut),"("+out.strip()+")"
								
								#USE UNAFOLD V3.8
								elif args.fold_ss == "unafold":
									subprocess.call(["hybrid-ss-min",gene_name+".fas"],stdout=open(gene_name+".out","w+"))
									out = open(gene_name+".fas.ct","r")
									for l in out:
										if l.find("=") != -1:
											dG = float(l.split("\t")[1].split("=")[1].strip())
											break

									subprocess.call(["hybrid-ss-min",gene_name+"_mut.fas"],stdout=open(gene_name+"_mut.out","w+"))
									out = open(gene_name+"_mut.fas.ct","r")
									for l in out:
										if l.find("=") != -1:
											dG_mut = float(l.split("\t")[1].split("=")[1].strip())
											break

									print str(dG)+"\t"+str(dG_mut)+"\t"+str(dG-dG_mut)

								#USE RNASTRUCTURE
								elif args.fold_ss == "rnastruct":
									subprocess.call(["/home/jacek/software/seq/RNAstructure-src/exe/Fold",gene_name+".fas",gene_name+".ct","-m","1"], stdout=open(gene_name+".out","w+"))
									out = open(gene_name+".ct","r")
									for l in out:
										if l.find("ENERGY") != -1:
											dG = float(l.split()[3])
											break
											
									subprocess.call(["/home/jacek/software/seq/RNAstructure-src/exe/Fold",gene_name+"_mut.fas",gene_name+"_mut.ct","-m","1"], stdout=open(gene_name+"_mut.out","w+"))
									out = open(gene_name+"_mut.ct","r")
									for l in out:
										if l.find("ENERGY") != -1:
											dG_mut = float(l.split()[3])
											break
									
									subprocess.call(["/home/jacek/software/seq/RNAstructure-src/exe/scorer",gene_name+".ct",gene_name+"_mut.ct",gene_name+"_cmp.out"],stdout=open(gene_name+"_cmp_stdout.out","w+"))
									out = open(gene_name+"_cmp.out","r")
									out_sens_ppv = ""
									for l in out:
										if l.find("Sensitivity") != -1:
											out_sens_ppv = "Sens:"+l.split()[5].rstrip()
											next
										if l.find("PPV") != -1:
											out_sens_ppv += " PPV:"+l.split()[5].rstrip()
											
									print str(dG)+"\t"+str(dG_mut)+"\t"+str(dG-dG_mut),"("+out_sens_ppv+")"
									
								break
							else:
								print gene_name,"SNP does not match reference sequence!"
								
				i+=1
f.close()

