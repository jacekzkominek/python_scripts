#!/usr/bin/env python

import os, sys, argparse
import pysam 
from collections import defaultdict

parser = argparse.ArgumentParser(description="Coverage")
parser.add_argument("f1", help="Input BAM file")
parser.add_argument("--start", default=332500, type=int, help="Start")
parser.add_argument("--end", default=348000, type=int, help="End")
parser.add_argument("--print_cov", default=False, action="store_true", help="End")
parser.add_argument("--print_regions", default=False, action="store_true", help="End")
parser.add_argument("--print_genes", default=False, action="store_true", help="End")
parser.add_argument("--print_reads", default=False, action="store_true", help="End")
parser.add_argument("--forward", default=False, action="store_true", help="End")
parser.add_argument("--f2r1", default=False, action="store_true", help="End")
args = parser.parse_args()

bamf = pysam.AlignmentFile(args.f1, "rb")

#CVER
#seq = "NODE_5_length_936291_cov_13.3144_ID_5280" 
#regions = []
#genes = [["entB",345938,346780,0]]

#CVER
#~ seq = "NODE_5_length_936291_cov_13.3144_ID_5280" 
#~ regions = [["entH_entA",336789,336815,0,0],["entA_entE",337573,337693,0,0],["entE_entC",339306,339502,0,0],["entC_entF",340689,340877,0,0],["entF_entD",344800,344983,0,0],["entD_entB",345702,345938,0,0]]
#~ genes = [["entH",336382,336789,0],["entA",336815,337573,0],["entE",337693,339306,0],["entC",339502,340689,0],["entF",340877,344800,0],["entD",344983,345702,0],["entB",345938,346780,0]]

#SBOM
#~ seq = "NODE_25_length_94690_cov_24.2565_ID_3467" 
#~ regions = [["entB_TM",33138,34080,0,0],["TM_Fre",35543,35607,0,0],["Fre_entA",37559,37823,0,0],["entA_entE",38575,38812,0,0],["entE_entC",40401,40877,0,0],["entC_entF",42022,42979,0,0],["entF_entD",46899,47350,0,0]]
#~ regions = [["entB_Fre",33138,35607,0,0],["Fre_entA",37559,37823,0,0],["entA_entE",38575,38812,0,0],["entE_entC",40401,40877,0,0],["entC_entF",42022,42979,0,0],["entF_entD",46899,47350,0,0]]
#~ genes = [["entB",32302,33138,0],["TM",34080,35543,0],["Fre",35607,37559,0],["entA",37823,38575,0],["entE",38812,40401,0],["entC",40877,42022,0],["entF",42979,46899,0],["entD",47350,48084,0]]

#CAPI
#~ seq = "LBNK01000008.1" 
#~ regions = [["entD_entF",512488,512836,0,0],["entF_entC",516879,519352,0,0],["entC_entE",520467,521051,0,0]]
#~ regions = [["entD_entF",512488,512836,0,0],["entF_SNZ",516879,517227,0,0],["SNZ_SNO",518081,518177,0,0],["SNO_entC",518806,519352,0,0],["entC_entE",520467,521051,0,0]]
#~ genes = [["entD",511757,512488,0],["entF",512836,516879,0],["SNZ",517227,518081,0],["SNO",518177,518806,0],["entC",519352,520467,0],["entE",521051,522658,0]]

#~ seq = "LBNK01000003.1" 
#~ regions = [["entA_Fre",943019,943463,0,0],["Fre_TM",945457,945509,0,0],["TM_entB",947014,947614,0,0]]
#~ genes = [["entA",942270,943019,0],["Fre",943463,945457,0],["TM",945509,947014,0],["entB",947614,948474,0]]

#COLI
seq = "NC_000913.3" 
#regions = [["entD_fepA",610079,610257,0,0],["fepA_fes",612494,612737,0,0],["fes_entF",613936,614157,0,0],["entF_fepE",618035,618254,0,0],["fepE_fepC",619384,619387,0,0],["fepC_fepG",620199,620199,0,0],["fepG_fepD",621188,621188,0,0],["fepD_entS",622189,622300,0,0],["entS_fepB",623547,623557,0,0],["fepB_entC",624510,624885,0,0],["entC_entE",626057,626070,0,0],["entE_entB",627677,627694,0,0],["entB_entA",628548,628551,0,0],["entA_entH",629294,629300,0,0]]
regions = [["entD_fepA",610079,610257,0,0],["fepA_fes",612494,612737,0,0],["fes_ybdZ",613936,613942,0,0],["ybdZ_entF",614157,614157,0,0],["entF_fepE",618035,618254,0,0],["fepE_fepC",619384,619387,0,0],["fepC_fepG",620199,620199,0,0],["fepG_fepD",621188,621188,0,0],["fepD_entS",622189,622300,0,0],["entS_fepB",623547,623557,0,0],["fepB_entC",624510,624885,0,0],["entC_entE",626057,626070,0,0],["entE_entB",627677,627694,0,0],["entB_entA",628548,628551,0,0],["entA_entH",629294,629300,0,0]]
genes = [["entD",609462,610079,0],["fepA",610257,612494,0],["fes",612737,613936,0],["ybdZ",613942,614157,0],["entF",614157,618035,0],["fepE",618254,619384,0],["fepC",619387,620199,0],["fepG",620199,621188,0],["fepD",621188,622189,0],["entS",622300,623547,0],["fepB",623557,624510,0],["entC",624885,626057,0],["entE",626070,627677,0],["entB",627694,628548,0],["entA",628551,629294,0],["entH",629300,629710,0],["fur",710203,710646,0]]

fr_flag = 64
rw_flag = 128
if args.f2r1 == True:
	fr_flag = 128
	rw_flag = 64



reads = bamf.fetch(seq, args.start, args.end)
read_cov = defaultdict(int)
span_cov = defaultdict(int)
read_names = defaultdict(list)
read_names2 = defaultdict(list)
passed_reads = defaultdict(int)
insert_limit = 1000
i = 1
for read in reads:
	read_data = str(read).split()
	name = read_data[0]
	flags = int(read_data[1])
	pos1 = int(read_data[3])
	pos2 = int(read_data[7])
	len1 = int(read_data[8])
	if pos2-pos1 > insert_limit:
		continue	
	#f1r2
	if args.forward == False:
		#if (flags & 64) != 0 and name not in passed_reads:
		if (flags & fr_flag) != 0 and name not in passed_reads:
			for p in range(pos1,pos1+len1):
				read_cov[p] += 1
			for p in range(pos1+len1,pos2):
				span_cov[p] += 1	
			for r in regions:
				if pos1 < r[1] and pos2 > r[2]:
					r[3] += 1
					read_names[r[0]].append(name)
					#read_names[r[0]].append(read)
				if (pos1 < r[1] and pos2 > r[1]) or (pos1 < r[2] and pos2 > r[2]):
					r[4] += 1
			for g in genes:	
				if pos1+len1 > g[1] and pos2 < g[2]:
					g[3] += 1
					#read_names2[g[0]].append(name)
		#elif (flags & 128) != 0 and name in passed_reads:
		elif (flags & rw_flag) != 0 and name in passed_reads:
			for p in range(pos1,pos1+len1):
				read_cov[p] += 1
			#if pos1 < r[2] and pos2 > r[1]:
			#		r[4] += 1
			for r in regions:
				if name not in read_names[r[0]]:
					if passed_reads[name] < r[1] and pos1+len1 > r[2]:
						r[3] += 1
						read_names[r[0]].append(name)
						#read_names[r[0]].append(read)
			#for g in genes:
			#	if name not in read_names2[g[0]]:
			#		if passed_reads[name] < g[1] and pos1+len1 > g[2]:
			#			r[3] += 1
			#			read_names[r[0]].append(name)
			#			#read_names[r[0]].append(read)
	elif args.forward == True:
		#if (flags & 128) != 0 and name not in passed_reads:
		if (flags & rw_flag) != 0 and name not in passed_reads:
			for p in range(pos1,pos1+len1):
				read_cov[p] += 1
			for p in range(pos1+len1,pos2):
				span_cov[p] += 1	
			for r in regions:
				if pos1 < r[1] and pos2 > r[2]:
					r[3] += 1
					read_names[r[0]].append(name)
					#read_names[r[0]].append(read)
				if (pos1 < r[1] and pos2 > r[1]) or (pos1 < r[2] and pos2 > r[2]):
					r[4] += 1
			for g in genes:	
				if pos1+len1 > g[1] and pos2 < g[2]:
					g[3] += 1
					#read_names2[g[0]].append(name)
		#elif (flags & 64) != 0 and name in passed_reads:
		elif (flags & fr_flag) != 0 and name in passed_reads:
			for p in range(pos1,pos1+len1):
				read_cov[p] += 1
			#if pos1 < r[2] and pos2 > r[1]:
			#		r[4] += 1
			for r in regions:
				if name not in read_names[r[0]]:
					if passed_reads[name] < r[1] and pos1+len1 > r[2]:
						r[3] += 1
						read_names[r[0]].append(name)
						#read_names[r[0]].append(read)
			#for g in genes:
			#	if name not in read_names2[g[0]]:
			#		if passed_reads[name] < g[1] and pos1+len1 > g[2]:
			#			r[3] += 1
			#			read_names[r[0]].append(name)
			#			#read_names[r[0]].append(read)
	passed_reads[name] = pos1

if args.print_regions == True:
	print "REGIONS COV" 
	for r in regions:
		print r[0],r[3],r[4]
	print 

if args.print_genes == True:
	print "GENES COV" 
	for g in genes:
		print g[0],g[1],g[2],g[2]-g[1],g[3]
	print 

if args.print_reads == True:
	print "REGIONS READS"
	for r in regions:
		print "\n"+r[0]
		for read in read_names[r[0]]:
			print read
	print 

if args.print_cov == True:
	print "POS READ SPAN SUM PROP"
	for pos in range(min(min(read_cov.keys(),span_cov.keys())),max(max(read_cov.keys(),span_cov.keys()))):
		prop = ""
		if read_cov[pos] > 0 or span_cov[pos] > 0:
			prop = "{0:.2f}".format(float(read_cov[pos])/(float(read_cov[pos])+float(span_cov[pos])))
			print pos,read_cov[pos],span_cov[pos],read_cov[pos]+span_cov[pos],prop
		
