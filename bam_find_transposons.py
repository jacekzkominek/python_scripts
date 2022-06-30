#!/usr/bin/env python

import pysam
import argparse
import sys
import os

parser = argparse.ArgumentParser(description='Find read pairs mapped to a transposon and to another region..')
#parser.add_argument("bam_file", help="Name-sorted BAM file to open")
args = parser.parse_args()

transposons = []
#t=[]
with open("/home/jacek/stuff/scripts/python/data/saccharomyces_cerevisiae_R64-1-1_20110208_ann.gff","r") as f:
	for l in f:
		#if l.find("Retrotransposon") != -1:
		if l.split()[2] == "LTR_retrotransposon":
			#transposons.append([l.split()[0],int(l.split()[3]),int(l.split()[4]),l.split()[8][l.split()[8].find("=")+1:l.split()[8].find(";")]])
			print l.strip()#t.append(l.strip())
			transposons.append([
				l.split()[8][l.split()[8].find("=")+1:l.split()[8].find(";")],
				l.split()[0],
				int(l.split()[3]),
				int(l.split()[4]),
				l.split()[6]])

#for tr in t:
#	print t
#for t in transposons:
#	print t[0],t[1],t[2],t[3],t[4]



##print transposons
#logfile = open("translog.txt","w")
#for d in sorted(os.listdir("/media/jacek/Data/reactors_sequencing/raw_bam_files_analysis_3.1/")):
	#if d.find("recal.bam") != -1 and d[-3:] == "bam":

		##infile = pysam.AlignmentFile(args.bam_file, "rb")
		#infile = pysam.AlignmentFile(d, "rb")
		#for trans in transposons:
			#trans_reads = infile.fetch(region=trans[1]+":"+str(trans[2])+":"+str(trans[3]))
			#for read in trans_reads:
				#r1_start = read.reference_start
				##r1_end = read.reference_end-1
				#try:
					#r1_mate = infile.mate(read)
				#except ValueError:
					#continue
				#else:
					#r2_start = r1_mate.reference_start
					##r2_end = mate(read).reference_end-1
					##if r1_start-r2_start >= 1000:
					#print >>logfile,d,trans[0],str(trans[1]),str(trans[2]),trans[3],read.query_name,r1_start,r2_start,r1_start-r2_start
					#sys.stdout.flush()

		#infile.close()
#logfile.close()
