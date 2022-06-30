#!/usr/bin/env python

import sys, argparse, os, itertools
from collections import defaultdict


def readfq(fp):#,barcode_length): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        #name, seqs, last = last[1:].partition(" ")[0], [], None
        name, seqs, last = last[1:], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
            #seqs.append(l[:barcode_length])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def readfq_full(f1):
	reads = []
	name = ""
	seq = ""
	qual = ""
	for i,l in enumerate(f1):
		l = l.strip()
		#if i % 4 == 2:# and l == "+":
		#	continue
		if i % 4 == 0:# and l == "@":
			name = l
		elif i % 4 == 1:
			seq = l
		elif i % 4 == 3:
			qual = l
			reads.append([name,seq,qual])
	return reads
					
    
parser = argparse.ArgumentParser(description="Demultiplex sequencing run.")
parser.add_argument('barcodes', help="Barcodes")
parser.add_argument('--r1', help="Read1")
parser.add_argument('--r2', help="Read2")
parser.add_argument('--out', default="out", help="Output dir")
parser.add_argument('--force_pair_bc', default=False, action="store_true", help="Force same barcodes on mates")
parser.add_argument('--no_out', default=False, action="store_true", help="Do not output reads")
parser.add_argument('--flush', default="all", help="Reads to flush at a time")
parser.add_argument('--pos', default=False, action="store_true", help="Check for positional effects")
args = parser.parse_args()

barcodes = {}
barcodes_list = []
bar_max = 0
bar_min = 10
with open(args.barcodes) as f:
	for l in f:
		l = l.strip()
		barcodes[l.split()[1]] = l.split()[0]
		barcodes_list.append(l.split()[1])
		bar_max=max(bar_max,len(l.split()[1]))
		bar_min=min(bar_min,len(l.split()[1]))


barcode_r1r2_r1 = defaultdict(list)
barcode_r1r2_r1_ct = defaultdict(int)
barcode_r1r2_r2 = defaultdict(list)
barcode_r1r2_r2_ct = defaultdict(int)
barcode_r1n_r1 = defaultdict(list)
barcode_r1n_r1_ct = defaultdict(int)
barcode_r1n_r2 = defaultdict(list)
barcode_r1n_r2_ct = defaultdict(int)
barcode_r2n_r1 = defaultdict(list)
barcode_r2n_r1_ct = defaultdict(int)
barcode_r2n_r2 = defaultdict(list)
barcode_r2n_r2_ct = defaultdict(int)
barcode_mix_r1 = []
barcode_mix_r1_ct = 0
barcode_mix_r2 = []
barcode_mix_r2_ct = 0
barcode_none_r1 = []
barcode_none_r1_ct = 0
barcode_none_r2 = []
barcode_none_r2_ct = 0

f1 = open(args.r1)
f2 = open(args.r2)

#PREP OUTPUT FILES
if args.no_out == False:
	if os.path.exists(os.getcwd()+"/"+args.out) == False:
		os.makedirs(os.getcwd()+"/"+args.out)
	open(os.getcwd()+"/"+args.out+"/none_R1.fq", 'w').close()
	open(os.getcwd()+"/"+args.out+"/none_R2.fq", 'w').close()
	open(os.getcwd()+"/"+args.out+"/mix_R1.fq", 'w').close()
	open(os.getcwd()+"/"+args.out+"/mix_R2.fq", 'w').close()
	for b in barcodes.keys():
		open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R1.fq", 'w').close()
		open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R2.fq", 'w').close()
		if args.force_pair_bc == False:
			open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R1nR1.fq", 'w').close()
			open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R1nR2.fq", 'w').close()
			open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R2nR1.fq", 'w').close()
			open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R2nR2.fq", 'w').close()
i = 0
if args.force_pair_bc == True:
	for read1,read2 in itertools.izip(readfq(f1),readfq(f2)):
	#~ for read1,read2 in itertools.izip(readfq_full(f1),readfq_full(f2)):
		seq1 = read1[1]
		seq2 = read2[1]
	
		r1_hit_bar = ""
		r2_hit_bar = "_"
		
		i += 1
		if args.flush != "all" and i % int(args.flush) == 0:
			for b in barcodes.keys():
				if args.no_out == False:
					with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R1.fq","a") as f1:
						f1.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r1r2_r1[b]]))
					with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R2.fq","a") as f2:
						f2.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r1r2_r2[b]]))	

					with open(os.getcwd()+"/"+args.out+"/none_R1.fq","a") as n1:
						n1.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_none_r1]))
					with open(os.getcwd()+"/"+args.out+"/none_R2.fq","a") as n2:
						n2.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_none_r2]))
													
					with open(os.getcwd()+"/"+args.out+"/mix_R1.fq","a") as n1:
						n1.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_mix_r1]))
					with open(os.getcwd()+"/"+args.out+"/mix_R2.fq","a") as n2:
						n2.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_mix_r2]))
				
				barcode_r1r2_r1_ct[b] += len(barcode_r1r2_r1[b])
				barcode_r1r2_r2_ct[b] += len(barcode_r1r2_r2[b])
				barcode_r1r2_r1[b] = []
				barcode_r1r2_r2[b] = []			
			
				barcode_none_r1_ct += len(barcode_none_r1)
				barcode_none_r2_ct += len(barcode_none_r2)
				barcode_none_r1 = []
				barcode_none_r2 = []

				barcode_mix_r1_ct += len(barcode_mix_r1)
				barcode_mix_r2_ct += len(barcode_mix_r2)
				barcode_mix_r1 = []
				barcode_mix_r2 = []


		if seq1[:bar_min] not in {seq2[:bar_min]:1}:
			barcode_none_r1.append(read1)
			barcode_none_r2.append(read2)
			continue
		else:
			for bartest in range(bar_min,bar_max+1):
				if seq1[:bartest] in barcodes:
					r1_hit_bar = seq1[:bartest]
					barcode_r1r2_r1[r1_hit_bar].append(read1)
					barcode_r1r2_r2[r1_hit_bar].append(read2)
					break
		
else:
	for read1,read2 in itertools.izip(readfq(f1),readfq(f2)):
	#~ for read1,read2 in itertools.izip(readfq_full(f1),readfq_full(f2)):
		name1 = read1[0]
		name2 = read2[0]
		
		tile1 = name1.split(":")[4]
		x1 = name1.split(":")[5]
		y1 = name1.split(":")[6]
		
		tile2 = name2.split(":")[4]
		x2 = name2.split(":")[5]
		y2 = name2.split(":")[6]
		
		seq1 = read1[1]
		seq2 = read2[1]

		r1_hit_bar = ""
		r2_hit_bar = ""
		
		i += 1
		#PRINT OUT READS
		if args.flush != "all" and i % int(args.flush) == 0:
			for b in barcodes.keys():
				if args.no_out == False:
					with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R1.fq","a") as f1:
						f1.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r1r2_r1[b]]))
					with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R2.fq","a") as f2:
						f2.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r1r2_r2[b]]))	
				
					with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R1nR1.fq","a") as f1:
						f1.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r1n_r1[b]]))
					with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R1nR2.fq","a") as f2:
						f2.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r1n_r2[b]]))
				
					with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R2nR1.fq","a") as f1:
						f1.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r2n_r1[b]]))
					with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R2nR2.fq","a") as f2:
						f2.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r2n_r2[b]]))

					with open(os.getcwd()+"/"+args.out+"/none_R1.fq","a") as n1:
						n1.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_none_r1]))
					with open(os.getcwd()+"/"+args.out+"/none_R2.fq","a") as n2:
						n2.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_none_r2]))
													
					with open(os.getcwd()+"/"+args.out+"/mix_R1.fq","a") as n1:
						n1.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_mix_r1]))
					with open(os.getcwd()+"/"+args.out+"/mix_R2.fq","a") as n2:
						n2.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_mix_r2]))
				
				barcode_r1r2_r1_ct[b] += len(barcode_r1r2_r1[b])
				barcode_r1r2_r2_ct[b] += len(barcode_r1r2_r2[b])
				barcode_r1r2_r1[b] = []
				barcode_r1r2_r2[b] = []			
				
				barcode_r1n_r1_ct[b] += len(barcode_r1n_r1[b])
				barcode_r1n_r2_ct[b] += len(barcode_r1n_r2[b])
				barcode_r1n_r1[b] = []
				barcode_r1n_r2[b] = []
				
				barcode_r2n_r1_ct[b] += len(barcode_r2n_r1[b])
				barcode_r2n_r2_ct[b] += len(barcode_r2n_r2[b])
				barcode_r2n_r1[b] = []
				barcode_r2n_r2[b] = []				
				
				barcode_none_r1_ct += len(barcode_none_r1)
				barcode_none_r2_ct += len(barcode_none_r2)
				barcode_none_r1 = []
				barcode_none_r2 = []

				barcode_mix_r1_ct += len(barcode_mix_r1)
				barcode_mix_r2_ct += len(barcode_mix_r2)
				barcode_mix_r1 = []
				barcode_mix_r2 = []
				
		#CHECK BARCODES
		for bartest in range(bar_min,bar_max+1):
			if seq1[:bartest] == seq2[:bartest] and seq1[:bartest] in barcodes:
				r1_hit_bar = seq1[:bartest]
				r2_hit_bar = seq1[:bartest]
				break
			elif seq1[:bartest] in barcodes:
				r1_hit_bar = seq1[:bartest]
			elif seq2[:bartest] in barcodes:
				r2_hit_bar = seq2[:bartest]
		
		if r1_hit_bar == r2_hit_bar and r1_hit_bar != "":	
			barcode_r1r2_r1[r1_hit_bar].append(read1)
			barcode_r1r2_r2[r1_hit_bar].append(read2)
		if r1_hit_bar != "" and r2_hit_bar == "":
			barcode_r1n_r1[r1_hit_bar].append(read1)
			barcode_r1n_r2[r1_hit_bar].append(read2)
		elif r1_hit_bar == "" and r2_hit_bar != "":
			barcode_r2n_r1[r2_hit_bar].append(read1)
			barcode_r2n_r2[r2_hit_bar].append(read2)
		elif r1_hit_bar != "" and r2_hit_bar != "" and r1_hit_bar != r2_hit_bar:
			barcode_mix_r1.append(read1)
			barcode_mix_r2.append(read2)
		elif r1_hit_bar == "" and r2_hit_bar == "":
			barcode_none_r1.append(read1)
			barcode_none_r2.append(read2)

		
f1.close()
f2.close()

total_r1r2 = 0
total_r1n = 0
total_r2n = 0
for b in sorted(barcodes):
	if args.flush != "all":
		if args.no_out == False:
			with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R1.fq","a") as f1:
				f1.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r1r2_r1[b]]))
			with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R2.fq","a") as f2:
				f2.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r1r2_r2[b]]))	
			with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R1nR1.fq","a") as f1:
				f1.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r1n_r1[b]]))
			with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R1nR2.fq","a") as f2:
				f2.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r1n_r2[b]]))	
			with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R2nR1.fq","a") as f1:
				f1.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r2n_r1[b]]))
			with open(os.getcwd()+"/"+args.out+"/"+barcodes[b]+"_R2nR2.fq","a") as f2:
				f2.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_r2n_r2[b]]))	
			
		barcode_r1r2_r1_ct[b] += len(barcode_r1r2_r1[b])
		barcode_r1r2_r2_ct[b] += len(barcode_r1r2_r2[b])
		if args.force_pair_bc == False:

			barcode_r1n_r2_ct[b] += len(barcode_r1n_r2[b])		
			barcode_r1n_r1_ct[b] += len(barcode_r1n_r1[b])
		
			barcode_r2n_r1_ct[b] += len(barcode_r2n_r1[b])
			barcode_r2n_r2_ct[b] += len(barcode_r2n_r2[b])
	
	total_r1r2 += barcode_r1r2_r1_ct[b]
	total_r1n += barcode_r1n_r1_ct[b]
	total_r2n += barcode_r2n_r1_ct[b]
	if args.force_pair_bc == False and args.flush != "all" and args.no_out == False:
		print barcodes[b],b,barcode_r1r2_r1_ct[b],barcode_r1n_r1_ct[b],barcode_r2n_r1_ct[b]
	else:
		print barcodes[b],b,barcode_r1r2_r1_ct[b],barcode_r1n_r1_ct[b],barcode_r2n_r1_ct[b]

if args.no_out == False:
	with open(os.getcwd()+"/"+args.out+"/none_R1.fq","a") as n1:
		n1.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_none_r1]))
	with open(os.getcwd()+"/"+args.out+"/none_R2.fq","a") as n2:
		n2.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_none_r2]))
		
	with open(os.getcwd()+"/"+args.out+"/mix_R1.fq","a") as n1:
		n1.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_mix_r1]))
	with open(os.getcwd()+"/"+args.out+"/mix_R2.fq","a") as n2:
		n2.write("\n".join(["\n".join(["@"+x[0],x[1],"+",x[2]]) for x in barcode_mix_r2]))
	barcode_none_r1_ct += len(barcode_none_r1)
	barcode_none_r2_ct += len(barcode_none_r2)
	barcode_mix_r1_ct += len(barcode_mix_r1)
	barcode_mix_r2_ct += len(barcode_mix_r2)

print "Total_R1R2",total_r1r2
print "Total_R1n",total_r1n
print "Total_R2n",total_r2n	
print "No_barcode",barcode_none_r1_ct
print "Mixed_barcodes",barcode_mix_r1_ct
print "Grand_total",str(total_r1r2+total_r1n+total_r2n+barcode_none_r1_ct+barcode_mix_r1_ct)
	

	
