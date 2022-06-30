#!/usr/bin/env python

import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from os import sys, close
from collections import defaultdict
import sys, argparse, os

#CONVERT LIST TO STRING
def list2str(l1):
	out_str = ""
	for a in l1:
		out_str = out_str+a[0]+"-"+str(a[1])+" "
	return out_str.rstrip()

def dic2str(dic1):
	out_str = ""
	for key, value in dic1.items():
		out_str = out_str+key+"-"+str(value)+"\t"
	return out_str.rstrip()

def dicval2str(dic1):
	out_str = ""
	for key, value in dic1.items():
		out_str = out_str+"\t"+str(value)
	return out_str.rstrip()

def nonzeroval(dic1,tol):
	nzval = 0
	for key, value in dic1.items():
		if value > 0+tol:
			nzval += 1
	return nzval

def checkblo(dic1, blodic,blo_tol):

	blosite = True
	for key1, val1 in dic1.items():
		if val1 > 0+blo_tol:
			for key2, val2 in dic1.items():
				if val2 > 0+blo_tol:
					if int(blodic[key1+key2]) <= 0:
						blosite = False
						break
			if blosite == False:
				break
	return blosite


#RETURN EVOLUTIONARY INDEX FOR 2 RESIDUES
def check_ei_score(s):
	EI = ["ST","VI","SA","NS","DE","IL","NT","YF","EQ","LM","TA","RK","KQ","NH","GA","QP","SG","QH","VL","RH","AP","KN","RQ","SP","AV","DN","TM","TP","KT","VM","EA","SC","RS","RT","IM","QL","LW","PH","TI","LF","SL","KI","HY","DA","DH","LH","KM","RP","EG","VF","EK","DG","IF","SI","GV","RG","EV","SY","RI","RM","RL","GC","PL","RC","NY","SW","SF","DV","CF","NI","CW","CY","RW","GW","DY"]
	ei_rank = 0
	for i, res_pair in enumerate(EI):
		if (s == res_pair or reverse(s) == res_pair):
			ei_rank = i+1
			break
	return ei_rank


#START MAIN
parser = argparse.ArgumentParser(description="Calculate amino acid frequencies in a multi-sequence alignment.")
parser.add_argument('file', help="Input FASTA file")
parser.add_argument('--blosum', nargs='?', default="62", help="Specific BLOSUM matrix to use in similarity check (BLO62 by default)")
parser.add_argument('--print_sites', action='store_true', default=False, help="Print full frequencies for individual sites.")
parser.add_argument('--print_sites_simple', action='store_true', default=False, help="Print only count and BLO for individual sites.")
parser.add_argument('--print_blosites', action='store_true', default=False, help="Print a list of BLOsites.")
parser.add_argument('--skip', nargs='?', type=int, default="0", help="Skip N initial sequences.")
parser.add_argument('--tol', nargs='?', type=int, default=0, help="Tolerance level (default = 0).")
parser.add_argument('--blo_tol', nargs='?', type=int, default=-1, help="BLOsite tolerance level (default = 0, or equal to --tol if not specified).")
#parser.add_argument('--skip_gaps', action='store_true', default=False, help="Skip all positions containing gaps.")

args = parser.parse_args()

#READ IN THE BLOSUM MATRIX
pathname = os.path.dirname(sys.argv[0])
blosum_f = open(pathname+"/blosum/blosum"+args.blosum+".blo")
blosum = defaultdict(int)
blosum_ref = ""
for line in blosum_f:
	if line.find('#') != -1:
		continue
	else:
		row = line.split()
		if len(row) == 24:
			blosum_ref = row
		elif len(row) == 25:
			for i, score in enumerate(row):
				if i > 0:
					blosum[blosum_ref[i-1]+row[0]] = score
blosum_f.close()
blo_tol = 0
if args.blo_tol == -1:
	blo_tol = args.tol
else:
	blo_tol = args.blo_tol

#OPEN FILE, READ IN SEQUENCE DATA
f = open(args.file)
o = open(args.file+"_log.txt",'w')
seqs = []

for i, seq in enumerate(list(SeqIO.parse(f, "fasta"))):
	if i >= args.skip:
		seqs.append(seq)

if len(seqs) == 0:
	print "No sequences loaded! Maybe the input file is not FASTA format? Exiting."
	sys.exit(-1)

seq_count = len(seqs)-args.skip
seq_len = len(str(seqs[0].seq).replace("-",""))

print "\nAlignment in file", "\""+args.file+"\"", "contains", seq_count, "sequences and", seq_len, "positions"
print>>o, "\nAlignment in file", "\""+args.file+"\"", "contains", seq_count, "sequences and", seq_len, "positions"
print "\nCalculating residue frequencies using tolerance level "+str(args.tol)+" and BLOsite tolerance level "+str(blo_tol)
print>>o, "\nCalculating residue frequencies using tolerance level "+str(args.tol)+" and BLOsite tolerance level "+str(blo_tol)

aa_count = defaultdict(int)
aa_blocount = defaultdict(list)
aa_blocount_total = 0
if args.print_sites == True:
	print "Pos\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\tX\tCount\tBLOsite"
	print>>o, "Pos\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\tX\tCount\tBLOsite"
elif args.print_sites_simple == True:
	print "Pos\tCount\tBLOsite"
	print>>o, "Pos\tCount\tBLOsite"

for x in range (0, seq_len):
	gap = False
	res_freq = dict({'A': 0, 'C': 0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0, 'X': 0})
	for seq in seqs:
		if seq[x] == '-':
			gap = True
			break
		else:
			res_freq[seq[x]] += 1
	if gap == True:
		continue
	else:
		aa_count[nonzeroval(res_freq,args.tol)] += 1
		blomark = ""
		if checkblo(res_freq, blosum,blo_tol) == True:
			aa_blocount[nonzeroval(res_freq,args.tol)].append(x+1)
			blomark = "*"
			aa_blocount_total += 1
		if args.print_sites == True:
			print str(x+1)+dicval2str(res_freq)+"\t"+str(nonzeroval(res_freq,args.tol))+"\t"+blomark
			print>>o, str(x+1)+dicval2str(res_freq)+"\t"+str(nonzeroval(res_freq,args.tol))+"\t"+blomark
		elif args.print_sites_simple == True:
			print str(x+1)+"\t"+str(nonzeroval(res_freq,args.tol))+"\t"+blomark
			print>>o, str(x+1)+"\t"+str(nonzeroval(res_freq,args.tol))+"\t"+blomark

print
print>>o
print "AAs\tSites\tBLOsites"
print>>o, "AAs\tSites\tBLOsites"
real_pos_count = 0
for x in range(1, 21):
	print str(x)+"\t"+str(aa_count[x])+"\t"+str(len(aa_blocount[x]))
	print>>o, str(x)+"\t"+str(aa_count[x])+"\t"+str(len(aa_blocount[x]))
	real_pos_count += aa_count[x]
print "Total\t"+str(real_pos_count)+"\t"+str(aa_blocount_total)
print>>o, "Total\t"+str(real_pos_count)+"\t"+str(aa_blocount_total)

print
print>>o
if args.print_blosites == True:
	print "AAs\tBLOsites:"
	print>>o, "AAs\tBLOsites:"
	for k,v in aa_blocount.items():
		if (v != []):
			print str(k)+"\t"+(', '.join(str(x) for x in v))
			print>>o, str(k)+"\t"+(', '.join(str(x) for x in v))

