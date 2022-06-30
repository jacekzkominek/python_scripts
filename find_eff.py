#!/usr/bin/env python

import sys, os, re
from collections import defaultdict
from collections import Counter
import operator

filename = sys.argv[1]
f = open(filename)
muts = []

for l in f:
	ef = l.find("EFF=")
	if ef != -1:
		ef_string = l[ef+4:ef+l[ef:].find(")")+1]
		if ef_string.split("(")[0].find("CODING") != -1 and ef_string.find("DSK2") != -1:
			mut = ef_string.split("|")[3]
			first = mut[0]
			second = first
			if mut[-1].isalpha():
				second = mut[-1]
			muts.append((int(re.findall(r'\d+', mut)[0]),first,second))

freqs = defaultdict(int)
for a in sorted(muts, key=lambda pos: pos[0]):
	freqs[a] += 1
	
sorted_x = sorted(freqs.iteritems(), reverse=True,key=operator.itemgetter(1))
for a in sorted_x:
	print a

