#!/usr/bin/env python

import os, sys
#import argparse
from collections import defaultdict

file1 = sys.argv[1]

annots = defaultdict(list)
with open(file1) as f:
	for l in f:
		source = ""
		seq = ""
		if len(l.split(".fas_")[0]) == 2:
			source = l.split(".fas_")[0]+".fas"
			seq = l.split(".fas_")[1]
		if len(l.split(".fas_")[0]) > 2:
			source = l[:l.rfind(".fas_")]+".fas"
			seq = l[l.rfind(".fas_")+5:]
		
		annots[source].append(seq) 
		
for k in annots.keys():
	f = open(k.replace(".fas",".ko"),"w")
	f.writelines(annots[k])
	f.close()
		
