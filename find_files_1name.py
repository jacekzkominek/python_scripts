#!/usr/bin/env python

import os

files_1name = set()

for f in os.listdir("/media/jacek/Data/bacteria/all_pep/"):
	if f.endswith(".fas"):
		files_1name.add(f.split("_")[0]+f.split("_")[1])

print len(files_1name)
