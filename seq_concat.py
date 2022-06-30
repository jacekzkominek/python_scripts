#!/usr/bin/env python

from Bio import SeqIO
from collections import defaultdict
import sys, argparse, os
from collections import defaultdict

parser = argparse.ArgumentParser(description='Run BLAST over the local database.')
parser.add_argument('dir', default = "")
parser.add_argument('--order', default=False, action="store_true")
parser.add_argument('--missing_taxa', default=False, action="store_true")
args = parser.parse_args()

def lcs(s1, s2):
    m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in xrange(1, 1 + len(s1)):
        for y in xrange(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]
    
d = args.dir

als = defaultdict(str)
names = defaultdict(list)
names2 = []
files = sorted(os.listdir(d), reverse=True)

for f in files:
	if f.find("mfal.fas") != -1:
		seqs = SeqIO.parse(open(f),"fasta")
		for s in seqs:
			if s.id not in names2:
				names2.append(s.id)

# ~ print names2
for f in files:
	if f.find("mfal.fas") != -1:
		print f
		seqs = SeqIO.parse(open(f),"fasta")
		len_base = 0
		seqs2 = []
		for s in seqs:
			len_base = len(s.seq)
			seqs2.append(s)
		ls = 0
		if args.missing_taxa == True and args.order == True:
			for n in names2:
				dum_seq = "X"*len_base
				for s in seqs2:
					if s.id == n:
						dum_seq = s.seq
						break
				als[n] = als[n]+dum_seq
				if ls == 0:
					ls = len(s.seq)
				if len(s.seq) != ls:
					print "ERROR IN LENGTH"
		else:
			if args.order == False:
				for s in seqs:
					#name = "_".join(s.id.split("_")[1:])
					name = s.id
					als[name] = als[name]+s.seq
					if ls == 0:
						ls = len(s.seq)
					if len(s.seq) != ls:
						print "ERROR IN LENGTH"
				# ~ han.close()

				lens = 0
				stop = False
				for a in als.keys():
					if lens == 0:
						lens = len(als[a])
						continue
					else:
						if len(als[a]) != lens:
							print "ERROR IN TOTAL LENGTH",f,a,lens,len(als[a])
							stop = True;
							break
				if stop == True:
					break
			else:
				for i,s in enumerate(seqs):
					names[i].append(s.id)
					als[i] = als[i]+s.seq
					if ls == 0:
						ls = len(s.seq)
					if len(s.seq) != ls:
						print "ERROR IN LENGTH"
				# ~ han.close()

				lens = 0
				stop = False
				for a in als.keys():
					if lens == 0:
						lens = len(als[a])
						continue
					else:
						if len(als[a]) != lens:
							print "ERROR IN TOTAL LENGTH",f,a,lens,len(als[a])
							stop = True;
							break
				if stop == True:
					break

with open("concat.fas","w") as out_f:		
	for a in sorted(als.keys()):
		if args.order == False:
			out_f.writelines("\n>"+str(a))
		elif args.order == True:
			out_f.writelines("\n>"+a)
			#out_f.writelines("\n>"+names[a][0])
		# ~ print a, als[a]
		# ~ break
		out_f.writelines("\n"+als[a])

