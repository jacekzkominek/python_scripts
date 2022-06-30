#!/usr/bin/env python3

import os, sys, argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='Obtain distribution of KEGG modules')
parser.add_argument("f", help="Files with input data.")
parser.add_argument("--counts", default=1, type=int, help="")
parser.add_argument("--counts_only", action="store_true", default=False, help="Print detailed genes and KEGGs")
parser.add_argument("--print_genes", action="store_true", default=False, help="Print detailed genes and KEGGs")
parser.add_argument("--score", default=200, type=int, help="")
parser.add_argument("--kingdom", action="store_true", default=False, help="Level1 Taxonomy (Kingdom)")
parser.add_argument("--phylum", action="store_true", default=False, help="Level2 Taxonomy (Phylum)")
parser.add_argument("--genus", action="store_true", default=False, help="Level3 Taxonomy (Genus)")
args = parser.parse_args()

data = []
with open(args.f) as f:
	for l in f:
		l = l.strip()
		f = l.split("\t")
		if len(f) >= 7:
			gene = f[0]
			kegg = f[1]
			class1 = f[2]
			class2 = f[3]
			genus = f[4]
			kegg_gid = f[5]
			score = f[6]
			if float(score) > float(args.score) and class1 != "" and class2 != "" and genus != "":
				entry = defaultdict(str)
				entry["gene"] = gene
				entry["kegg"] = kegg
				entry["class1"] = class1
				entry["class2"] = class2
				entry["genus"] = genus
				entry["kegg_gid"] = kegg_gid
				entry["score"] = float(score)
				data.append(entry)


if args.kingdom == True:
	print("\n")
	current_entry = defaultdict(str)
	current_keggs = set()
	for g in sorted(data, key=lambda x: (x["class1"])):
		if current_entry["class1"] != None and current_entry["class1"] != g["class1"] and len(current_keggs) >= int(args.counts):
			print(current_entry["class1"],len(current_keggs), sep="\t", end="")
			if args.counts_only == False:
				print("\t", end="")
				print(" ".join(sorted(current_keggs)), end="\n")
			else:
				print("", end="\n")
			current_entry = g
			current_keggs = set()
		else:
			current_keggs.add(g["kegg"])

if args.phylum == True:
	print("\n")
	current_entry = defaultdict(str)
	current_keggs = set()
	current_genes = []
	for g in sorted(data, key=lambda x: (x["class1"], x["class2"])):
		if current_entry["class2"] != None and current_entry["class2"] != g["class2"] and len(current_keggs) >= int(args.counts):
			if args.print_genes == True:
				print("")
			print(current_entry["class1"],current_entry["class2"],len(current_genes), len(current_keggs), sep="\t", end="")
			if args.counts_only == False:
				print("\t", end="")
				print(" ".join(sorted(current_keggs)), end="\n")
			else:
				print("", end="\n")
			if args.print_genes == True:
				for k in sorted(current_keggs, reverse=True):
					print(k)
					for cg in sorted(current_genes, key=lambda x: x["kegg"]):
						if k == cg["kegg"]:
							print(cg["gene"])
			current_entry = g
			current_keggs = set()
			current_genes = []
		else:
			if g["kegg"] != "":
				current_keggs.add(g["kegg"])
			current_genes.append(g)

if args.genus == True:
	print("\n")
	current_entry = defaultdict(str)
	current_keggs = set()
	current_genes = []
	for g in sorted(data, key=lambda x: (x["class1"], x["class2"], x["genus"], -x["score"])):
		if current_entry["genus"] != None and current_entry["genus"] != g["genus"] and len(current_keggs) >= int(args.counts):
			if args.print_genes == True:
				print("")
			print(current_entry["class1"],current_entry["class2"],current_entry["genus"],len(current_genes), len(current_keggs), sep="\t", end="")
			if args.counts_only == False:
				print("\t", end="")
				print(" ".join(sorted(current_keggs, reverse=True)), end="\n")
			else:
				print("", end="\n")
			if args.print_genes == True:
				for k in sorted(current_keggs, reverse=True):
					print(k)
					for cg in sorted(current_genes, key=lambda x: x["kegg"]):
						if k == cg["kegg"]:
							print(cg["gene"])
			current_entry = g
			current_keggs = set()
			current_genes = []
		else:
			if g["kegg"] != "":
				current_keggs.add(g["kegg"])
			current_genes.append(g)


