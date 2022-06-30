#!/usr/bin/env python

from Bio.Blast.Applications import *
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import sys, argparse, os

parser = argparse.ArgumentParser(description='Run BLAST over the local database.')# add_help=False)
parser.add_argument('file', help="Input FASTA file")
parser.add_argument('database', help="BLAST database(s) to search")
parser.add_argument('--dbtype', nargs='?', default="aa", help="Target database type [aa (default)|nn|gen]", dest="dbtype")
parser.add_argument('--subdir', nargs='?', default=".", help="Target subdirectory.")
parser.add_argument('--txt', action='store_true', help="Generate TXT output instead of XML (does not extract seqs)", dest="txt")
parser.add_argument('--best_only', '--best', action='store_true', help="Extract the best hits only", dest="best")
#parser.add_argument('--no_margin', action='store_true', help="Extract just the BLAST hit", dest="no_margin")
parser.add_argument('--margin', type=int, default=0, help="Positions around the hit to extract together with it", dest="margin")
#parser.add_argument('--no_extract', action='store_true', help="Do not extract any hits, only run the BLAST", dest="no_ext")
#parser.add_argument('--outfile', nargs='?', default="BLAST_output", help="Output file", dest="outfile")
parser.add_argument('--split_res', action='store_true', help="Split results from multiple queries into individual files", dest="split_res")
parser.add_argument('--add_E', action='store_true', help="Add E-values to the extracted hits", dest="add_E")
parser.add_argument('--rename', action='store_true', help="Change full hit names to short 4-letter names (e.g.Scer)", dest="rename")
parser.add_argument('--max_hits', nargs='?', type=int, default=500, help="Maximum number of hits")
parser.add_argument('--tab', action='store_true', help="Output as flat tab file")
args = parser.parse_args()

if (os.path.exists(args.file) == False):
	print "\n",args.file,"does not exist!\n"
	exit()

fmt = 5
if args.txt == True:
	fmt = 0
if args.tab == True:
	fmt = 6

best = ""
if args.best == True:
	best = "best_"

dblist=""
subdir=args.subdir
max_hits=10000
if args.max_hits != "":
	max_hits=args.max_hits

#SENSU STRICTO DATASET = 8 SPECIES
if args.database == "sac_sensu_stricto":
	dblist = "sarb sbay scer seub skud smik spar suva"
	subdir="sac_sensu_stricto"

#POST-WGD DATASET = 8 SPECIES
elif args.database == "sac_post_WGD":
	dblist = "cgla kafr knag ncas ndai tbla tpha vpol"
	subdir="sac_post_WGD"

#PRE-WGD DATASET = 8 SPECIES
elif args.database == "sac_pre_WGD":
	dblist = "ecym egos klac kthe lklu lwal tdel zrou"
	subdir="sac_pre_WGD"

#NO SPECIAL DATABASE SET - USE WHAT IS GIVEN IN THE PARAMETER
else:
	dblist=args.database
#else:
#	subdir = args.database.split("/")[0]
#	dblist = args.database.split("/")[1]

prog=""
if args.file.find("_nn") >= 0:
	if args.dbtype == "aa":
		prog = "blastx"
	else:
		prog = "blastn"
elif args.file.find("_aa") >= 0:
	if args.dbtype == "aa":
		prog = "blastp"
	else:
		prog = "tblastn"

outdir = "./blast_files"
if os.path.exists(outdir) == False:
	os.mkdir(outdir)

dbs = dblist.split();
for db in dbs:
	db1 = db
	if args.dbtype == "aa":
		db1 += "_aa"
	elif args.dbtype == "nn":
		db1 += "_nn"

	blastdb=os.environ['BLASTDB']+subdir+"/"+db1
	outfile = outdir+"/"+db1+"_"+best+"blast_hits"

	if fmt == 5:
		outfile += ".xml"
	elif fmt == 0 or fmt == 6:
		outfile += ".txt"
		print "Text output format, won't extract sequences!"

	cmdline = ""
	print ("BLASTing "+db1),
	if prog == "blastn":
		cmdline = NcbiblastnCommandline(query=args.file, db=blastdb, evalue=1, outfmt=fmt, out=outfile, max_target_seqs=max_hits)
	if prog == "blastx":
		cmdline = NcbiblastxCommandline(query=args.file, db=blastdb, evalue=1, outfmt=fmt, out=outfile, max_target_seqs=max_hits)
	if prog == "tblastn":
		cmdline = NcbitblastnCommandline(query=args.file, db=blastdb, evalue=1, outfmt=fmt, out=outfile, max_target_seqs=max_hits)
	if prog == "blastp":
		cmdline = NcbiblastpCommandline(query=args.file, db=blastdb, evalue=1, outfmt=fmt, out=outfile, max_target_seqs=max_hits)

	#print
	#print cmdline

	stdout, stderr = cmdline()
	print "- done"
	if fmt == 5:
		db_source = blastdb+".fas"
		db_dict = SeqIO.index(db_source, "fasta")
		blast_hits=""
		if (args.split_res == False):
			blast_hits = open(db1+"_"+best+"hits.fas", "w")
		for blast_record in NCBIXML.parse(open(outfile)):
			if (args.split_res == True):
				blast_hits = open(blast_record.query+"_"+db1+"_hits.fas", "w")
			if blast_record.alignments:
				hit_count = 1
				for alignment in blast_record.alignments:
					for hsp in alignment.hsps:
						db_seq=db_dict[alignment.hit_def.split(" ")[0]].seq

						hit_s_pos=min(hsp.sbjct_start,hsp.sbjct_end)
						hit_e_pos=max(hsp.sbjct_start,hsp.sbjct_end)
						hit_len = len(db_seq)
						if (args.margin > 0):
							hit_s_pos -= args.margin
							hit_e_pos += args.margin
							if (hit_s_pos <= 1):
								hit_s_pos = 0
							if (hit_e_pos >= hit_len):
								hit_e_pos = hit_len
							db_seq = db_seq[hit_s_pos:hit_e_pos]
							if (hit_e_pos >= hit_len):
								hit_e_pos = str(hit_len)+"E"
						else:
							db_seq = db_seq[hit_s_pos:hit_e_pos]

						Eval = ""
						if (args.add_E == True):
							Eval = "/"+str(hsp.expect)

						hit_title = ""
						if args.rename == False:
							hit_title = blast_record.query+"_hit_"+str(hit_count)+"|length:"+str(hsp.align_length)+"|score:"+str(hsp.bits)+Eval+"|"+alignment.hit_def+"["+str(hit_s_pos+1)+"-"+str(hit_e_pos)+"]("+str(hsp.frame[1])+")"
						else:
							hit_title = db+blast_record.query[4:]
						#print hit_title+"\n"+db_seq+"\n"
						blast_hits.write('> %s\n' % (hit_title))
						blast_hits.write('%s\n' % (db_seq))
						if (args.best == True):
							break
						hit_count += 1
					if (args.best == True):
						break
			if (args.split_res == True):
				blast_hits.close()
		if (args.split_res == False):
			blast_hits.close()




