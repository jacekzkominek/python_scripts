#!/usr/bin/env python2

from Bio import SeqIO
import os, sys, argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Trim sequences shorter than specific length.")
parser.add_argument('--f1', help="Genomes biosamples to filename map")
parser.add_argument('--genome_dir', help="Folder with genome files")
parser.add_argument('--cont_dir', help="Folder with contamination results")
parser.add_argument('--out_dir', default="outdir", help="Folder with contamination results")
args = parser.parse_args()

if os.path.exists(args.out_dir) == False:
	os.makedirs(args.out_dir)
biosamples = {}
with open(args.f1,"r") as inf1:
	for l in inf1:
		l = l.strip()
		biosamples[l.split("\t")[0]] = l.split("\t")[-1]

sample_list = set()
trim_list = defaultdict(list)
exclude_list = defaultdict(list)
dup_list = defaultdict(list)
for f in sorted(os.listdir(args.cont_dir)):
	with open(args.cont_dir+"/"+f,"r") as inf2:
		sample = ""
			
		exclude_flag = False
		exclude = []
		trim_flag = False
		trim = []
		dup_flag = False
		dup = []
		for l in inf2:
			l = l.strip()
			if l.find("PRJNA") != -1:
				sample = l.split()[2]
			if l.find("Exclude:") != -1:
				exclude_flag = True
			if exclude_flag == True and l.find("Sequence name") == -1 and len(l.split()) >= 1:
				exclude.append(l.split()[0])
			if l.find("Trim:") != -1:
				trim_flag = True
				exclude_flag = False
			if trim_flag == True and l.find("Sequence name") == -1 and len(l.split()) >= 3:
				if l.split()[2].find(",") == -1:
					trim.append(l.split()[0]+"_"+l.split()[2])
				else:
					ts = l.split()[2].split(",")
					for t in reversed(ts):
						trim.append(l.split()[0]+"_"+t)
			if l.find("Duplicated:") != -1:
				dup_flag = True
				trim_flag = False
			if dup_flag == True and l.find("Sequence names") == -1 and len(l.split()) >= 3:
				for d in l.split()[1:-2]:
					dup.append(d.replace("lcl|",""))
		
		sample_list.add(sample)
		trim_list[sample] = trim
		exclude_list[sample] = exclude
		dup_list[sample] = dup
		#print exclude
		#print 
		print trim
		#exit()
		#print biosamples[sample],len(exclude),len(trim),len(dup)
#print sample_list
#print sorted(sample_list)
for s in sorted(sample_list):
	seqs = []
	print biosamples[s]
	for seq_record in SeqIO.parse(args.genome_dir+"/"+biosamples[s], "fasta"):
		if seq_record.id not in exclude_list[s] and seq_record.id not in dup_list[s] and len(seq_record.seq) >= 200:
			intrim = False
			for t in trim_list[s]:
				if t.find(seq_record.id) != -1:
					intrim = True
			if intrim == False:
				seqs.append(seq_record)
			else:
				for t in trim_list[s]:
					if t.find(seq_record.id) != -1:
						seq_id = t.split("_")[0]
						# ~ print t
						trim_start = int(t.split("_")[-1].split("..")[0])-1
						trim_end = int(t.split("_")[-1].split("..")[1])
						seq_tmp = seq_record.seq
						print(len(seq_tmp))
						seq_tmp = seq_tmp[:trim_start]+seq_tmp[trim_end:]
						print(len(seq_tmp))
						seq_record.seq = seq_tmp
				seqs.append(seq_record)
	SeqIO.write(seqs, args.out_dir+"/"+biosamples[s], "fasta")
	#exit()
