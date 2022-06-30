#!/usr/bin/env python

import os
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('folder', nargs='?', help="Folder to search.")
parser.add_argument('--ext', default="fas", nargs='?', help="Folder to search.")
args = parser.parse_args()

files_all = []

for f in os.listdir(args.folder):
	if f.endswith(args.ext):
		files_all.append(f)

for f in os.listdir(args.folder):
	if f.endswith(".tsv") and f[:-4] in files_all:
		files_all.remove(f[:-4])

for f in sorted(files_all):
	print f
