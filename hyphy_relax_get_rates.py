#!/usr/bin/env python

import json
import os, sys, argparse

parser = argparse.ArgumentParser(description='Read Hyphy RELAX input.')
parser.add_argument('infile', default="", help="Input file")
parser.add_argument('--best', default=False, action='store_true', help="Output only best model.")
#parser.add_argument('--tree_k', default=False, action='store_true', help="Output only best model.")
args = parser.parse_args()
path = os.getcwd()

with open (args.infile, "r") as fh:
	relax_json = json.load(fh)
	print "p-value",relax_json["relaxation-test"]["p"]
	print "K",relax_json["fits"]["Alternative"]["K"]
	aic = 999999999
	for model in relax_json["fits"]:
		aic = min(relax_json["fits"][model]["AIC-c"],aic)
	for model in relax_json["fits"]:
		if args.best == False or relax_json["fits"][model]["AIC-c"] == aic:
			print "\nModel:",model
			print "AIC-c",relax_json["fits"][model]["AIC-c"]
			print "LnL",relax_json["fits"][model]["log-likelihood"]
			tree = relax_json["tree"]
			for cat_type in relax_json["fits"][model]["rate-distributions"]:
				print cat_type
				for cat in relax_json["fits"][model]["rate-distributions"][cat_type]:
					print str(cat[0])+"\t"+str(cat[1])
			branch_k = {}
			if relax_json["fits"][model].get("branch-annotations"):
				for b in relax_json["fits"][model]["branch-annotations"]:
					branch_k[b] = relax_json["fits"][model]["branch-annotations"][b]
				print "\nTree with",relax_json["fits"][model]["annotation-tag"]
				for b in branch_k:
					#~ print b,branch_k[b]
					tree = tree.replace(b+",",b+":"+str(branch_k[b])+",").replace(b+")",b+":"+str(branch_k[b])+")")
				print tree+";"
			else:
				print "###No tree data###"
		
	#~ if args.tree_k == True:
		#~ tree = relax_json["tree"]
		#~ branch_k = {}
		#~ for b in relax_json["fits"]["General Descriptive"]["branch-annotations"]:
			#~ branch_k[b] = relax_json["fits"]["General Descriptive"]["branch-annotations"][b]
		#~ print "\nTree K"
		#~ for b in branch_k:
			#~ tree = tree.replace(b+",",b+":"+str(branch_k[b])+",").replace(b+")",b+":"+str(branch_k[b])+")")
		#~ print tree+";"
	
