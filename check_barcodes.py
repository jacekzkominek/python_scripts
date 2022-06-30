#!/usr/bin/env python

import sys, argparse, os, itertools
from collections import defaultdict
import operator 

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

f1 = open(sys.argv[1],"r")
f2 = open(sys.argv[2],"r")

br1 = defaultdict(int)
br2 = defaultdict(int)
for read1,read2 in itertools.izip(readfq(f1),readfq(f2)):
	seq1 = read1[1]
	seq2 = read2[1]
	br1[seq1[:5]] += 1
	br2[seq2[:5]] += 1
	
print "BR1"
sorted_br1 = sorted(br1.items(), reverse=True, key=operator.itemgetter(1))
for b in sorted_br1:
	print b,br1[b]

sorted_br2 = sorted(br2.items(), reverse=True, key=operator.itemgetter(1))
print "BR2"
for b in sorted_br2:
	print b,br2[b]
