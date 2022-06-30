#!/usr/bin/env python

import os
import argparse

parser = argparse.ArgumentParser(description="Read PREQC reports and extract Branch Var/Dup/Err data.")
parser.add_argument('file1', help="Input *.preqc files")
#parser.add_argument('--only_sig_pv', action='store_true', default=False, help="")
args = parser.parse_args()

all_dp = {}
if args.file1 == "all":
	for f in sorted(os.listdir(os.getcwd())):
		#print f
		if f.find(".preqc") != -1:
			br_clas = []
			with open(f) as fo:
				br_clas_found = False
				dp = {}
				for l in fo:
					if br_clas_found == False:
						if l.find("BranchClassification") != -1:
							br_clas_found = True
					elif br_clas_found == True:
						if l.find(":") != -1:
							dp[l.strip().split(":")[0].replace("\"","")] = l.strip().split(":")[1].replace(",","").strip()
						if l.find("}") != -1:
							br_clas.append(dp)
							dp = {}
						if l.find("]") != -1:
							break
			all_dp[f] = br_clas

#~ ks = ["variant_rate","repeat_rate"]
ks = ["num_variant_branches","num_repeat_branches","num_error_branches"]
for k in ks:
	dps = set()
	for f in sorted(all_dp.keys()):
		for p in all_dp[f]:
			dps.add(p["k"])
	
	print "\n"+k+"\t"+" ".join(sorted(dps))
	for f in sorted(all_dp.keys()):
		print f,
		data_matrix = []
		for dp in sorted(dps):
			data = 0
			for p in all_dp[f]:
				if dp == p["k"]:
					data = 0
					if float(p["num_kmers"]) != 0:
						data = float(p[k])/float(p["num_kmers"])
					break
			print data,
			data_matrix.append(data)
		print
		
