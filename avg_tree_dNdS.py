#!/usr/bin/env python

import sys, os, math

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

filename = sys.argv[1]
f = open(filename)
total=0
count=0

for line in f:
	line2 = line.replace(')',"")
	line3 = line2.replace('(',"")
	line4 = line3.replace(',',":")
	line5 = line3.replace(',',";")
	line6 = line5.split(':')
	for a in line6:
		if is_number(a):
			total += float(a)
			count += 1
print			
print "Total\t\tBr.count\tTotal/count"
print total,"\t",count,"\t\t",total/count
print
