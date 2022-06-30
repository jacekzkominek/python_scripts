#!/usr/bin/env python

import sys, argparse, os, subprocess, shutil
from collections import defaultdict
from ete3 import Tree

parser = argparse.ArgumentParser(description='')
parser.add_argument("--sp_tree", default="", help="Species tree")
args = parser.parse_args()

#keep = ["babjeviella_inositovora","brettanomyces_anomalus","brettanomyces_bruxellensis","candida_albicans","candida_dubliniensis","candida_maltosa","candida_orthopsilosis","candida_parapsilosis","candida_sojae","candida_tanzawaensis","candida_tenuis","candida_tropicalis","clavispora_lusitaniae","geotrichum_candidum","hyphopichia_burtonii","kazachstania_africana","kluyveromyces_aestuarii","kluyveromyces_dobzhanskii","kluyveromyces_marxianus","kluyveromyces_wickerhamii","lachancea_lanzarotensis","lipomyces_starkeyi","lodderomyces_elongisporus","metschnikowia_bicuspidata","metschnikowia_fructicola","meyerozyma_caribbica","meyerozyma_guilliermondii","nadsonia_fulvescens","naumovozyma_dairenensis","pachysolen_tannophilus","saccharomyces_arboricola","saccharomyces_cerevisiae","saccharomyces_eubayanus","saccharomyces_mikatae","saccharomyces_paradoxus","saccharomyces_uvarum","saprochaete_clavata","scheffersomyces_stipitis","spathaspora_arborariae","spathaspora_passalidarum","starmerella_bombicola","tetrapisispora_phaffii","tortispora_caseinolytica","torulaspora_delbrueckii","vanderwaltozyma_polyspora","wickerhamomyces_anomalus","wickerhamomyces_ciferrii","yHMPu5000034635_nadsonia_fulvescens","yHMPu5000034965_brettanomyces_naardenensis","yHMPu5000034976_dekkera_anomala","yHMPu5000034998_cephaloascus_albidus","yHMPu5000034999_cephaloascus_fragrans","yHMPu5000035043_babjeviella_inositovora","yHMPu5000035269_wickerhamomyces_ciferrii","yHMPu5000035273_wickerhamomyces_anomalus","yHMPu5000035316_candida_homilentoma","yHMPu5000035682_hyphopichia_burtonii","yHMPu5000035720_candida_chilensis","yHMPu5000037913_candida_cylindracea","zygosaccharomyces_bailii"]
keep2 = ["babjeviella_inositovora","brettanomyces_anomalus","brettanomyces_bruxellensis","candida_albicans","candida_dubliniensis","candida_maltosa","candida_orthopsilosis","candida_parapsilosis","candida_sojae","candida_tanzawaensis","candida_tenuis","candida_tropicalis","clavispora_lusitaniae","debaryomyces_hansenii","geotrichum_candidum","hyphopichia_burtonii","kazachstania_africana","kazachstania_naganishii","kluyveromyces_aestuarii","kluyveromyces_dobzhanskii","kluyveromyces_lactis","kluyveromyces_marxianus","kluyveromyces_wickerhamii","lachancea_kluyveri","lachancea_lanzarotensis","lachancea_thermotolerans","lipomyces_starkeyi","lodderomyces_elongisporus","metschnikowia_bicuspidata","metschnikowia_fructicola","meyerozyma_caribbica","meyerozyma_guilliermondii","nadsonia_fulvescens","naumovozyma_castellii","naumovozyma_dairenensis","pachysolen_tannophilus","saccharomyces_arboricola","saccharomyces_cerevisiae","saccharomyces_eubayanus","saccharomyces_kudriavzevii","saccharomyces_mikatae","saccharomyces_paradoxus","saccharomyces_uvarum","saprochaete_clavata","scheffersomyces_stipitis","spathaspora_arborariae","spathaspora_passalidarum","starmerella_bombicola","tetrapisispora_blattae","tetrapisispora_phaffii","tortispora_caseinolytica","torulaspora_delbrueckii","vanderwaltozyma_polyspora","wickerhamomyces_anomalus","yarrowia_lipolytica","yHMPu5000034635_nadsonia_fulvescens","yHMPu5000034965_brettanomyces_naardenensis","yHMPu5000034976_dekkera_anomala","yHMPu5000034998_cephaloascus_albidus","yHMPu5000034999_cephaloascus_fragrans","yHMPu5000035269_wickerhamomyces_ciferrii","yHMPu5000035273_wickerhamomyces_anomalus","yHMPu5000035316_candida_homilentoma","yHMPu5000035682_hyphopichia_burtonii","yHMPu5000035720_candida_chilensis","yHMPu5000037913_candida_cylindracea","zygosaccharomyces_bailii","zygosaccharomyces_rouxii"]
keep = []
for k in keep2:
	# ~ keep.append(k.lower())
	keep.append(k)

nad = ["yHMPu5000034635_nadsonia_fulvescens","nadsonia_fulvescens"]
wic = ["yHMPu5000035269_wickerhamomyces_ciferrii","wickerhamomyces_ciferrii","wickerhamomyces_anomalus","yHMPu5000035273_wickerhamomyces_anomalus"]
bret = ["yHMPu5000034965_brettanomyces_naardenensis","yHMPu5000034976_dekkera_anomala","brettanomyces_anomalus","brettanomyces_bruxellensis"]
t1 = Tree(args.sp_tree, format=1)
t1.prune(keep)
t1.write(format=0, outfile=args.sp_tree+"_ref.nwk")


t1 = Tree(args.sp_tree, format=1)
t1.prune([item for item in keep if item not in nad])
t1.write(format=0, outfile=args.sp_tree+"_nad.nwk")

t1 = Tree(args.sp_tree, format=1)
t1.prune([item for item in keep if item not in wic])
t1.write(format=0, outfile=args.sp_tree+"_wic.nwk")

t1 = Tree(args.sp_tree, format=1)
t1.prune([item for item in keep if item not in bret])
t1.write(format=0, outfile=args.sp_tree+"_bret.nwk")

t1 = Tree(args.sp_tree, format=1)
t1.prune([item for item in keep if item not in nad+wic+bret])
t1.write(format=0, outfile=args.sp_tree+"_nad_wic_bret.nwk")

exit()
t1 = Tree(args.sp_tree, format=1)
t1.prune([item for item in keep if item not in wic+bret])
t1.write(format=0, outfile=args.sp_tree+"_nad.nwk")

t1 = Tree(args.sp_tree, format=1)
t1.prune([item for item in keep if item not in nad+bret])
t1.write(format=0, outfile=args.sp_tree+"_wic.nwk")

t1 = Tree(args.sp_tree, format=1)
t1.prune([item for item in keep if item not in nad+wic])
t1.write(format=0, outfile=args.sp_tree+"_bret.nwk")

t1 = Tree(args.sp_tree, format=1)
t1.prune([item for item in keep if item not in nad+wic+bret])
t1.write(format=0, outfile=args.sp_tree+"_nad_wic_bret.nwk")
