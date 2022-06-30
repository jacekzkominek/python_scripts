#!/usr/bin/env python

import os, sys, argparse
import pysam 
from collections import defaultdict

parser = argparse.ArgumentParser(description="Coverage")
parser.add_argument("f1", help="Input BAM file")
parser.add_argument("gff", help="Input GFF file")
parser.add_argument("--start", default=332500, type=int, help="Start")
parser.add_argument("--end", default=348000, type=int, help="End")
parser.add_argument("--print_cov", default=False, action="store_true", help="End")
parser.add_argument("--print_regions", default=False, action="store_true", help="End")
parser.add_argument("--print_genes", default=False, action="store_true", help="End")
parser.add_argument("--print_reads", default=False, action="store_true", help="End")
parser.add_argument("--forward", default=False, action="store_true", help="End")
args = parser.parse_args()

bamf = pysam.AlignmentFile(args.f1, "rb")

#~ #CVER
#~ seq = "NODE_5_length_936291_cov_13.3144_ID_5280" 
#~ regions = [["entH_entA",336789,336815,0,0],["entA_entE",337573,337693,0,0],["entE_entC",339306,339502,0,0],["entC_entF",340689,340877,0,0],["entF_entD",344800,344983,0,0],["entD_entB",345702,345938,0,0]]
#~ genes = [["entH",336382,336789,0],["entA",336815,337573,0],["entE",337693,339306,0],["entC",339502,340689,0],["entF",340877,344800,0],["entD",344983,345702,0],["entB",345938,346780,0]]

print "Query\tGene\tContig\tStart\tEnd\tStrand\tReadpairs_Total\tReadpairs_forward\tReadpairs_reverse"

with open(args.gff,"r") as f:
	for l in f:
		#~ gi|1002316097|dbj|BCJV01000005.1|	maker	CDS	49107	49982	.	+	0	gene_id "snap_masked-gi|1002316097|dbj|BCJV01000005.1|-processed-gene-0.11"; transcript_id "snap_masked-gi|1002316097|dbj|BCJV01000005.1|-processed-gene-0.11-mRNA-1";
		seq = l.split()[0]
		args.start = int(l.split()[3])
		args.end = int(l.split()[4])
		strand = l.split()[6]
		gene = l.strip().split()[9]
		query = l.strip().split()[12]
		
		genes = [[gene,args.start,args.end,0]]
		regions = [[gene,args.start,args.end,0,0]]
		
		print query,gene.replace("\"","").replace(";",""),seq,args.start,args.end,strand,
		sys.stdout.flush()
		
		k = 0
		plus_k = 0
		minus_k = 0
		#~ passed_reads = defaultdict(int)
		#~ passed_reads = []
		passed_reads = set()
		#~ k_all = 0
		#~ reads = bamf.fetch(seq, args.start, args.end,until_eof=True)
		#~ for read in reads:
			#~ flags = int(str(read).split()[1])
			#~ if (flags & 8) == 0:
				#~ k_all+=1
		#~ print k_all,
		#~ sys.stdout.flush()
		
		reads = bamf.fetch(seq, args.start, args.end)
		#~ print
		for read in reads:
			read_data = str(read).split()
			flags = int(read_data[1])
			name = read_data[0]
			#~ pos1 = int(read_data[3])
			if name not in passed_reads:
				if (flags & 8) == 0:
					k+=1
					if (flags & 64) != 0:
						plus_k+=1
					elif (flags & 128) != 0:
						minus_k+=1
						
			#~ passed_reads[name] = flags
			passed_reads.add(name)
			#~ print str(k)+"\r",
		#~ print k,
		print k,plus_k,minus_k
		sys.stdout.flush()
		
		#~ reads = bamf.fetch(seq, args.start, args.end)
		#~ passed_reads = defaultdict(int)
		#~ for read in reads:
			#~ read_data = str(read).split()
			#~ flags = int(read_data[1])
			#~ name = read_data[0]
			#~ pos1 = int(read_data[3])
			#~ pos2 = int(read_data[7])
			#~ if (flags & 8) == 0 and (flags & 64) != 0 and name not in passed_reads.keys():
				#~ plus_k+=1
			#~ passed_reads[name] = pos1
		#~ print plus_k,
		#~ sys.stdout.flush()
				
		#~ reads = bamf.fetch(seq, args.start, args.end)
		#~ passed_reads = defaultdict(int)
		#~ for read in reads:
			#~ read_data = str(read).split()
			#~ flags = int(read_data[1])
			#~ name = read_data[0]
			#~ pos1 = int(read_data[3])
			#~ pos2 = int(read_data[7])
			#~ if (flags & 8) == 0 and (flags & 128) != 0 and name not in passed_reads.keys():
				#~ minus_k+=1
			#~ passed_reads[name] = pos1
			
		#~ print minus_k
		#~ sys.stdout.flush()
		continue
		
		
		reads = bamf.fetch(seq, args.start, args.end)
		read_cov = defaultdict(int)
		span_cov = defaultdict(int)
		read_names = defaultdict(list)
		read_names2 = defaultdict(list)
		passed_reads = defaultdict(int)
		insert_limit = 20000
		i = 1
		for read in reads:
			#~ k+=1
			sys.stdout.write('%s\r' % k)
			read_data = str(read).split()
			name = read_data[0]
			flags = int(read_data[1])
			pos1 = int(read_data[3])
			pos2 = int(read_data[7])
			len1 = int(read_data[8])
			if pos2-pos1 > insert_limit:
				continue	
			
			if args.forward == False:
				if (flags & 64) != 0 and name not in passed_reads.keys():
					if args.print_cov == True:
						for p in range(pos1,pos1+len1):
							read_cov[p] += 1
						for p in range(pos1+len1,pos2):
							span_cov[p] += 1	
					if args.print_regions == True:
						for r in regions:
							if pos1 < r[1] and pos2 > r[2]:
								r[3] += 1
								read_names[r[0]].append(name)
								#read_names[r[0]].append(read)
							if (pos1 < r[1] and pos2 > r[1]) or (pos1 < r[2] and pos2 > r[2]):
								r[4] += 1
					if args.print_genes == True:
						for g in genes:	
							if pos1+len1 > g[1] and pos2 < g[2]:
								g[3] += 1
								#read_names2[g[0]].append(name)
				elif (flags & 128) != 0 and name in passed_reads.keys():
					if args.print_cov == True:
						for p in range(pos1,pos1+len1):
							read_cov[p] += 1
					#if pos1 < r[2] and pos2 > r[1]:
					#		r[4] += 1
					if args.print_regions == True:
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
				if (flags & 128) != 0 and name not in passed_reads.keys():
					if args.print_cov == True:
						for p in range(pos1,pos1+len1):
							read_cov[p] += 1
						for p in range(pos1+len1,pos2):
							span_cov[p] += 1	
					if args.print_regions == True:
						for r in regions:
							if pos1 < r[1] and pos2 > r[2]:
								r[3] += 1
								read_names[r[0]].append(name)
								#read_names[r[0]].append(read)
							if (pos1 < r[1] and pos2 > r[1]) or (pos1 < r[2] and pos2 > r[2]):
								r[4] += 1
					if args.print_genes == True:
						for g in genes:	
							if pos1+len1 > g[1] and pos2 < g[2]:
								g[3] += 1
								#read_names2[g[0]].append(name)
				elif (flags & 64) != 0 and name in passed_reads.keys():
					if args.print_cov == True:
						for p in range(pos1,pos1+len1):
							read_cov[p] += 1
					#if pos1 < r[2] and pos2 > r[1]:
					#		r[4] += 1
					if args.print_regions == True:
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
			print "\nREGIONS COV" 
			for r in regions:
				print r[0],r[3],r[4]

		if args.print_genes == True:
			#~ print "\nGENES COV" 
			for g in genes:
				print seq,g[0],g[1],g[2],g[2]-g[1],g[3]

		if args.print_reads == True:
			print "\nREGIONS READS"
			for r in regions:
				print "\n"+r[0]
				for read in read_names[r[0]]:
					print read

		if args.print_cov == True:
			print "POS READ SPAN SUM PROP"
			for pos in range(min(min(read_cov.keys(),span_cov.keys())),max(max(read_cov.keys(),span_cov.keys()))):
				prop = ""
				if read_cov[pos] > 0 or span_cov[pos] > 0:
					prop = "{0:.2f}".format(float(read_cov[pos])/(float(read_cov[pos])+float(span_cov[pos])))
					print pos,read_cov[pos],span_cov[pos],read_cov[pos]+span_cov[pos],prop
				
