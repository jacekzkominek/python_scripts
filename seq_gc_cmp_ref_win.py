#!/usr/bin/env python

import warnings
import sys, argparse, os
from collections import defaultdict
from natsort import natsorted 

parser = argparse.ArgumentParser(description="")
parser.add_argument("f1", help="Input file")
parser.add_argument("f2", nargs="+", help="Input dir")
parser.add_argument("--strict", action="store_true", default=False, help="")
args = parser.parse_args()

ref = defaultdict(float)
with open(args.f1) as f:
	for l in f:
		if len(l.split()) >= 2:
			ref[l.split()[0].lower()] = float(l.split()[1])

print "Input",
if len(args.f2) == 1:
	with open(args.f2[0]) as f:
		wins = []
		for l in f:
			if len(l.split()) >= 2:
				if l.find("Input") != -1:
					for i,j in enumerate(l.split()):
						if (i+1)%2 == 0:
							wins.append(j)
							print args.f2[0]+"_"+j,
					print
				else:
					sp = l.split()[0].lower()
					for r in ref:
						if sp.find(r) != -1:
							print sp,
							for i,j in enumerate(l.split()):
								if (i+1)%2 == 0:
									#print ref[r]
									print abs(float(ref[r])-float(float(j)*100)),
							print
							break
else:
	data_dict = defaultdict(lambda: defaultdict(float))
	glob_keys = set()
	for file2 in args.f2:
		with open(file2) as f:
			wins = []
			for l in f:
				if len(l.split()) >= 2:
					if l.find("Input") == -1:
						sp = l.split()[0].lower()
						gc = l.split()[1]
						for r in ref:
							if sp.find(r) != -1:
								for i,j in enumerate(l.split()):
									if (i+1)%2 == 0:
										gc = abs(float(ref[r])-float(float(gc)*100))
								break
						data_dict[sp][file2] = gc
						glob_keys.add(file2)
	for g in natsorted(glob_keys):
		print g,
	print
	for d in natsorted(data_dict.keys()):
		
		hits = 0
		for g in natsorted(glob_keys):
			if g in natsorted(data_dict[d].keys()):
				hits += 1
		if args.strict == False or hits == len(glob_keys):
			print d,
			for g in natsorted(glob_keys):
				if g in natsorted(data_dict[d].keys()):
					print "{:.5f}".format(float(data_dict[d][g])),
					#print data_dict[d][g],
				else:
					print "",
			print
		
	#for file2 in f2:	
		#with open(file2) as f:
			#wins = []
			#for l in f:
				#if len(l.split()) >= 2:
					#if l.find("Input") != -1:
						#for i,j in enumerate(l.split()):
							#if (i+1)%2 == 0:
								#wins.append(j)
								#print args.f2+"_"+j,
						#print
					#else:
						#sp = l.split()[0].lower()
						#for r in ref:
							#if sp.find(r) != -1:
								#print sp,
								#for i,j in enumerate(l.split()):
									#if (i+1)%2 == 0:
										##print ref[r]
										#print abs(float(ref[r])-float(j)),
								#print
								#break
				