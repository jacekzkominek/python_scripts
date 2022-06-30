#!/usr/bin/env python

import os, sys, argparse, scipy.stats
from collections import defaultdict

parser = argparse.ArgumentParser(description='')
parser.add_argument('gff_file', help="Dir.")
args = parser.parse_args()

with open(args.gff_file, "r") as f:
	cov_list = []
	# ~ gene_list = defaultdict(float)
	# ~ genes = ["entA","entB","entC","entD","entE","entF","entH"]
	for l in f:
		l = l.strip()
		if len(l.split()) > 4:
			contig = l.split()[0]
			cds = l.split()[2]
			start = l.split()[3]
			end = l.split()[4]
			gene = l.split()[-1]
			cov = l.split()[-5].replace("\"","").replace(";","")
			if cds == "transcript":
				# ~ for g in genes:
					# ~ if g in l:
						# ~ gene_list[g] = float(cov)
						# ~ break
				cov_list.append(float(cov))
	#CVER
	gene_list = {"entA":285.977,"entB":918.491,"entC":1549.001,"entD":926.413,"entE":954.094,"entF":1082.735,"entH":1472.135}
	#CAPI
	# ~ gene_list = {"entA":211.640,"entB":1305.492,"entC":87.806,"entD":183.955,"entE":140.243,"entF":28.833}
	#SBOM
	# ~ gene_list = {"entA":463.384,"entB":1752.682,"entC":240.183,"entD":109.934,"entE":1211.553,"entF":746.039}
	for g in sorted(gene_list):
		print g,str(gene_list[g])+" ("+str(int(scipy.stats.percentileofscore(cov_list,gene_list[g])))+")"
