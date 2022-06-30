#!/usr/bin/env python

import sys
import fileinput
import string

table = string.maketrans(":;,[]'() ","_________")

for line in fileinput.input(sys.argv[1], inplace=True): 
	print line.translate(table),
