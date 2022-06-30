#!/usr/bin/env python

import os, sys
from collections import defaultdict
assemble_stats = defaultdict(int)
list_total = []
for f in os.listdir(sys.argv[1]):
	
	with open(f,"r") as f1:
		for l in f1:
			if l.find("NODE") != -1:
				list_total.append(f+"_SPADES")
			elif l.find("scf") != -1:
				list_total.append(f+"_MASURCA")
			elif l.find("flattened_line") != -1:
				list_total.append(f+"_DISCOVAR")
			else:
				list_total.append(f+"_X")
			break

for l in list_total:
	print l
