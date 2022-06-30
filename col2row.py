#!/usr/bin/env python
import csv
from os import sys, close

arr = []

f = open(sys.argv[1])

for l in f:
	arr.append(l.split())
f.close()

arr2 = zip(*arr)

with open("1.out", "wb") as f:
	writer = csv.writer(f,delimiter='\t')
	writer.writerows(arr2)
