#!/usr/bin/env python

import os, sys
import argparse
import shutil 
from collections import defaultdict

parser = argparse.ArgumentParser(description='')
parser.add_argument("file1", help="QUAST Report")

args = parser.parse_args()

rep = args.file1
assems = defaultdict(int)
with open(rep) as f:
	for l in f:
		if l.find("#") == -1 and len(l.split()) >= 17:
			assems[l.split()[0]] = int(l.split()[17])

if os.path.exists(os.getcwd()+"/below100kb") == False:
	os.makedirs(os.getcwd()+"/below100kb")
for f in sorted(os.listdir(os.getcwd())):
	f1=f
	f=f.replace(".fas","") 
	if f in assems.keys() and assems[f] < 100000:
		print f,assems[f]
		shutil.move(os.getcwd()+"/"+f1,os.getcwd()+"/below100kb/"+f1)
		
