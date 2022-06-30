#!/usr/bin/env python

import sys, argparse, os
import csv
from collections import *
import mygene
import scipy.stats
import statsmodels.sandbox.stats.multicomp as sssm
mg = mygene.MyGeneInfo()

parser = argparse.ArgumentParser(description="Extract subnetwork from SIF file and identify core nodes.")
parser.add_argument('file1', help="SIF file with base network")
parser.add_argument('file2', help="List of genes in the subnetwork")
parser.add_argument('--only_sig_pv', action='store_true', default=False, help="Print only genes with significant enrichment (p-value)")
parser.add_argument('--only_sig_pvC', action='store_true', default=False, help="Print only genes with significant enrichment (FDR corrected p-value)")
parser.add_argument('--only_sig_pv2', action='store_true', default=False, help="Print only genes with significant enrichment (p-value of N or more)")
parser.add_argument('--only_sig_pvC2', action='store_true', default=False, help="Print only genes with significant enrichment (FDR corrected p-value of N or more)")
args = parser.parse_args()

if args.only_sig_pvC == True:
	args.only_sig_pv = True
if args.only_sig_pvC2 == True:
	args.only_sig_pv2 = True

#LOAD A MIRRORED BASE NETWORK
base_net = defaultdict(set)
base_total = 0
with open(args.file1, 'rb') as f1:
    reader = csv.reader(f1, delimiter='\t')
    for row in reader:
		node1 = row[0]
		node2 = row[2]
		base_net[node1].add(node2)
		base_net[node2].add(node1)
		base_total += 1


#LOAD LIST OF GENES IN THE SUBNETWORK
sub_net_list = set()
print "\n=== Reading",args.file2,"==="
with open(args.file2, 'rb') as f2:
	for l in f2:
		sub_net_list.add(l.strip())

#EXTRACT THE SUBNETWORK FROM THE BASE
sub_net = defaultdict(set)
missing = []
for g in sub_net_list:
	if g not in base_net.keys():
		#print g,"is missing from the base network!"
		missing.append(g)
	else:
		sub_net[g] = base_net[g].intersection(sub_net_list)

#CALCULATE THE TOTAL NUMBER OF INTERACTIONS IN THE SUBNETWORK
checked = list()
for g in sub_net.keys():
	for g2 in sub_net[g]:
		if (g+g2 in checked) or (g2+g in checked):
			next
		else:
			checked.append(g+g2)

total = len(checked)
print "Genes in the base network:",len(base_net)
print "Genes in the subnetwork:",len(sub_net_list)
print "Genes in the subnetwork (valid):",len(sub_net),"("+str(len(missing)),"missing",",".join(missing)+")"
print "Valid interactions:",total


sub_net_len = defaultdict(list)

#CALCULATE ENRICHMENT P-VALUES FOR ALL GENES IN THE SUBNETWORK
stats = []
stats2 = []
for g in sorted(sub_net, key=lambda s: len(sub_net[s])):
	pv = scipy.stats.hypergeom.pmf(len(sub_net[g]),len(base_net),len(base_net[g]),len(sub_net))
	pv2 = 1-scipy.stats.hypergeom.cdf(len(sub_net[g])-1,len(base_net),len(base_net[g]),len(sub_net))
	stats.append([g,pv])
	stats2.append([g,pv2])

#SORT P-VALUES
stats_sorted = []
pvs = []
for g in sorted(stats, key=lambda s: s[1]):
	stats_sorted.append(g)
	pvs.append(g[1])

stats2_sorted = []
pvs2 = []
for g in sorted(stats2, key=lambda s: s[1]):
	stats2_sorted.append(g)
	pvs2.append(g[1])


#CALCULATE FDR CORRECTION FOR THE P-VALUES
pvsC = sssm.fdrcorrection0(pvs, is_sorted=True)[1]
for (g, pvC) in zip(stats_sorted, pvsC):
	g.append(pvC)

pvsC2 = sssm.fdrcorrection0(pvs2, is_sorted=True)[1]
for (g, pvC2) in zip(stats2_sorted, pvsC2):
	g.append(pvC2)

#PRINT SUMMARY STATISTICS FOR ALL GENES IN THE SUBNETWORK
print
print "N".ljust(5)+"Gene".ljust(12)+"Base".rjust(12)+"BaseTotal%".rjust(12)+"Sub".rjust(12)+"SubTotal%".rjust(12)+"SubBase%".rjust(12)+"p-val".rjust(12)+"p-valC".rjust(12)+"p-val+".rjust(12)+"p-valC+".rjust(12)
n = 1
for g in sorted(sub_net, key=lambda s: len(sub_net[s])):
	pv = 0
	pvC = 0
	for s in stats_sorted:
		if s[0] == g:
			pv = s[1]
			pvC = s[2]
			break
	pv2 = 0
	pvC2 = 0
	for s2 in stats2_sorted:
		if s2[0] == g:
			pv2 = s2[1]
			pvC2 = s2[2]
			break

	if (args.only_sig_pvC2 == False or pvC2 < 0.05):
	#if (args.only_sig_p == False and args.only_sig_pC == False and args.only_sig_p+ == False and args.only_sig_pC+ == False) /
	#or (args.only_sig_pC2 == True and pvC2 <= 0.05)
	#or (args.only_sig_pC == True and pvC <= 0.05)
	#or (args.only_sig_p == True and pv <= 0.05):
		pv = "%.4f" % (pv)
		pvC = "%.4f" % (pvC)
		pv2 = "%.4f" % (pv2)
		pvC2 = "%.4f" % (pvC2)

		conv_id = mg.query(g, scopes='ensemblgene', fields='symbol', species="559292", as_dataframe=False)
		print str(n).ljust(4),
		if conv_id['hits'] == []:
			print g.ljust(12),
			sub_net_len[len(sub_net[g])].append(g)
		else:
			print conv_id['hits'][0]['symbol'].ljust(12),
			sub_net_len[len(sub_net[g])].append(conv_id['hits'][0]['symbol'])

		#bt = "%.2f%%" % (float(len(base_net[g]))*100/float(base_total))
		bt = "%.2f%%" % (float(len(base_net[g]))*100/float(len(base_net)))
		sys.stdout.write(str(len(base_net[g])).rjust(12)+bt.rjust(12)+str(len(sub_net[g])).rjust(12))

		if total == 0:
			a = "%.2f%%".rjust(12) % (0/1)
			sys.stdout.write(a.rjust(12))
		else:
			#a = "%.2f%%".rjust(12) % (float(len(sub_net[g]))*100/total)
			a = "%.2f%%".rjust(12) % (float(len(sub_net[g]))*100/len(sub_net))
			sys.stdout.write(a.rjust(12))

		if len(base_net[g]) == 0:
			b = "%.2f%%" % (0/1)
			print b.rjust(12),
		else:
			b = "%.2f%%" % (float(len(sub_net[g]))*100/float(len(base_net[g])))
			print b.rjust(12),


		sys.stdout.write(pv.rjust(12)+pvC.rjust(12)+pv2.rjust(12)+pvC2.rjust(12)+"\n")
		sys.stdout.flush()
		n += 1

print "\nSignificant genes:".ljust(21)+"<0.05".rjust(8)+"<0.01".rjust(8)+"<0.001".rjust(8)
print "At p-value:".ljust(20)+str(len([x for x in stats if x[1] <= 0.05])).rjust(8)+str(len([x for x in stats if x[1] <= 0.01])).rjust(8)+str(len([x for x in stats if x[1] <= 0.001])).rjust(8)
print "At p-valueC:".ljust(20)+str(len([x for x in stats if x[2] <= 0.05])).rjust(8)+str(len([x for x in stats if x[2] <= 0.01])).rjust(8)+str(len([x for x in stats if x[2] <= 0.001])).rjust(8)
print "At p-value+:".ljust(20)+str(len([x for x in stats2 if x[1] <= 0.05])).rjust(8)+str(len([x for x in stats2 if x[1] <= 0.01])).rjust(8)+str(len([x for x in stats2 if x[1] <= 0.001])).rjust(8)
print "At p-valueC+:".ljust(20)+str(len([x for x in stats2 if x[2] <= 0.05])).rjust(8)+str(len([x for x in stats2 if x[2] <= 0.01])).rjust(8)+str(len([x for x in stats2 if x[2] <= 0.001])).rjust(8)

print "\nInt\t#\tGenes"
for l in sorted(sub_net_len):
	print str(l)+"\t"+str(len(sub_net_len[l]))+"\t",
	for g in sub_net_len[l]:
		print g,
		sys.stdout.flush()
	print

print