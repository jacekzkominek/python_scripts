#!/usr/bin/env python

import sys, argparse, os, itertools
from collections import defaultdict
from natsort import natsorted

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
					
    
parser = argparse.ArgumentParser(description="Demultiplex sequencing run.")
parser.add_argument('r1', help="Read1")
parser.add_argument('r2', help="Read2")
#~ parser.add_argument('--flush', default="all", help="Reads to flush at a time")
#~ parser.add_argument('--pos', default=False, action="store_true", help="Check for positional effects")
args = parser.parse_args()

f1 = open(args.r1)
f2 = open(args.r2)


tile1 = defaultdict(int)
x1 = defaultdict(int)
y1 = defaultdict(int)
xy1 = defaultdict(int)

tile2 = defaultdict(int)
x2 = defaultdict(int)
y2 = defaultdict(int)
xy2 = defaultdict(int)

tile = ""
#for read1 in itertools.izip(readfq(f1)):#,readfq(f2)):
for read1 in readfq(f1):
	name1 = read1[0]
	#name2 = read2[0]
	if len(name1.split(":")) >= 7:# and len(name2.split(":")) >= 7:
		if tile != name1.split(":")[4] and len(xy1.keys()) > 0:
			if tile != "":
				for t,x,y in natsorted(xy1.keys()):
					print t,x,y
				sys.stdout.flush()
			xy1 = defaultdict(int)
			#xy2 = defaultdict(int)
			tile = name1.split(":")[4]
		xy1[(name1.split(":")[4],name1.split(":")[5],name1.split(":")[6].split()[0])] += 1
		#xy2[(name2.split(":")[4],name2.split(":")[5],name2.split(":")[6].split()[0])] += 1
			
		#~ tile1[name1.split(":")[4]] += 1
		#~ x1[name1.split(":")[5]] += 1
		#~ y1[name1.split(":")[6].split()[0]] += 1

		#~ tile2[name2.split(":")[4]] += 1
		#~ x2[name2.split(":")[5]] += 1
		#~ y2[name2.split(":")[6].split()[0]] += 1

for t,x,y in sorted(xy1.keys()):
	print t,x,y#,xy1[(x,y)]
	

	
