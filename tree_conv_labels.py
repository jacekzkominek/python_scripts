#!/usr/bin/env python

import sys, os, math, argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="Convert tree labels")
parser.add_argument('tree', default="", help="")
parser.add_argument('--conv', default="", help="")
args = parser.parse_args()

conv_list = defaultdict(str)
with open(args.conv) as f:
	for l in f:
		conv_list[l.strip().split()[0]] = l.strip().split()[1]

with open(args.tree) as f:
	for l in f:
		for c in conv_list:
			l = l.replace(c,c+"_"+conv_list[c])
		print l 
