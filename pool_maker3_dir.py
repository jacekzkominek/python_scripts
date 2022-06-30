#!/usr/bin/env python

import os
import fileinput
import argparse
import shutil
from natsort import natsorted

parser = argparse.ArgumentParser(description="")
parser.add_argument("dir", default = "", help="Directory with maker subdirs")
#parser.add_argument("sub", default = "", help="Target directory with submit files")
#parser.add_argument("--nosub", default = False, action = "store_true", help="Do not submit")
args = parser.parse_args()

for f in sorted(os.listdir(args.dir)):
	if f.find("_out") != -1:
		print f,
		prot = []
		nuc = []
		gff = []
		for f2 in natsorted(os.listdir(args.dir+"/"+f)):
			if f2.find("proteins") != -1 and f2.find("proteins_maker") == -1:
				prot.append(f2)
			if f2.find("transcripts.fas") != -1:
				nuc.append(f2)
			if f2.find("transcripts.gff") != -1:
				gff.append(f2)
		if len(prot) >= 2 and len(nuc) >= 2 and len(gff) >= 2:
			print
			shutil.copyfile(args.dir+"/"+f+"/"+prot[-2],args.dir+"/"+prot[-2])
			shutil.copyfile(args.dir+"/"+f+"/"+nuc[-2],args.dir+"/"+nuc[-2])
			shutil.copyfile(args.dir+"/"+f+"/"+gff[-2],args.dir+"/"+gff[-2])
		else:
			print " - not found"