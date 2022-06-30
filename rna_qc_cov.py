#!/usr/bin/env python

import argparse
from collections import defaultdict
import subprocess
import os

parser = argparse.ArgumentParser(description="Calculate average coverage across features")
parser.add_argument("f1", help="Input BED file")
parser.add_argument("--norna", default=False, action="store_true", help="Don't run SortMeRNA")
parser.add_argument("--nomap", default=False, action="store_true", help="Don't do mapping")
parser.add_argument("--nocov", default=False, action="store_true", help="Don't do coverage mapping")
parser.add_argument("--local", default=False, action="store_true", help="Local alignment with Bowtie2")
parser.add_argument("--bwa", default=False, action="store_true", help="Alignment with BWA")
parser.add_argument("--strain", default=False, action="store_true", help="Strain level")
parser.add_argument("--db", default="RDP", choices=["barcodes","RDP","genbank","unite","unite_big","clete","ena500","ena100"], help="RNA database to use")
parser.add_argument("--cpu", default="6", help="CPUs to use")
args = parser.parse_args()

#RUN SORTMERNA AND BLAST THE RNADB
cpu = args.cpu
read_file = args.f1
sortmerna_path = "/home/jacek/software/assembly_reads/sortmerna-2.0-linux-64/"
sortmerna_path = "/home/jacek/software/assembly_reads/sortmerna-2.1-linux-64/"

if args.norna == False:
	ret = subprocess.call([sortmerna_path+"sortmerna","--ref",sortmerna_path+"rRNA_databases/"+args.db+".fasta,"+sortmerna_path+"rRNA_databases/"+args.db+"_db","--reads",read_file,"--fastx","--aligned",read_file+"_rRNA","-m","8192","-a",cpu], stdout=None, stderr=subprocess.STDOUT)

if args.nomap == False:
	FNULL = open(os.devnull, 'w')
	loc="--end-to-end"
	if args.local == True:
		loc="--local "

	#print read_file
	if args.bwa == False:	
		if read_file.find(".fq") != -1 or read_file.find(".fastq") != -1:
			ret = subprocess.call(["/home/jacek/software/assembly_reads/bowtie2-2.2.6/bowtie2",loc,"--very-sensitive","-N","1","-x",sortmerna_path+"rRNA_databases/"+args.db+"_db","-U",read_file+"_rRNA."+read_file.split(".")[-1],"-S",read_file+"_rRNA.sam","-p",cpu], stdout=None, stderr=subprocess.STDOUT)
		elif read_file.find(".fa") != -1 or read_file.find(".fasta") != -1:
			ret = subprocess.call(["/home/jacek/software/assembly_reads/bowtie2-2.2.6/bowtie2",loc,"--very-sensitive","-N","1","-x",sortmerna_path+"rRNA_databases/RDP_db","-U",read_file+"_rRNA."+read_file.split(".")[-1],"-f","-S",read_file+"_rRNA.sam","-p",cpu], stdout=None, stderr=subprocess.STDOUT)
	else:
		out_f = open(read_file+"_rRNA.sai",'w')
		ret = subprocess.call(["/home/jacek/software/assembly_reads/bwa/bwa","aln",sortmerna_path+"rRNA_databases/RDP_db.fasta",read_file+"_rRNA.fasta"], stdout=out_f)
		out_f.close()
		out_f = open(read_file+"_rRNA.sam",'w')
		ret = subprocess.call(["/home/jacek/software/assembly_reads/bwa/bwa","samse",sortmerna_path+"rRNA_databases/RDP_db.fasta",read_file+"_rRNA.sai",read_file], stdout=out_f)
		out_f.close()

	ret = subprocess.call(["/home/jacek/software/assembly_reads/samtools-1.2/samtools","view","-Shbo",read_file+"_rRNA_temp1.bam",read_file+"_rRNA.sam"], stdout=None, stderr=subprocess.STDOUT)
	#ret = subprocess.call(["/home/jacek/software/samtools-1.2/samtools","view","-Shbq","1","-o",read_file+"_rRNA_temp1.bam",read_file+"_rRNA.sam"], stdout=None, stderr=subprocess.STDOUT)
	ret = subprocess.call(["/home/jacek/software/assembly_reads/samtools-1.2/samtools","sort","-@",cpu,read_file+"_rRNA_temp1.bam",read_file+"_rRNA_temp2"], stdout=None, stderr=subprocess.STDOUT)
	ret = subprocess.call(["/home/jacek/software/assembly_reads/samtools-1.2/samtools","rmdup",read_file+"_rRNA_temp2.bam",read_file+"_rRNA.bam"], stdout=FNULL, stderr=FNULL)#subprocess.STDOUT)
	ret = subprocess.call(["/home/jacek/software/assembly_reads/samtools-1.2/samtools","index",read_file+"_rRNA.bam"], stdout=None, stderr=subprocess.STDOUT)
	out_f = open(read_file+"_rRNA_cov.txt","w")
	ret = subprocess.call(["bedtools","genomecov","-ibam",read_file+"_rRNA.bam","-g",sortmerna_path+"rRNA_databases/RDP_db.fasta"], stdout=out_f, stderr=subprocess.STDOUT)
	out_f.close()


if args.nocov == False:
	cov = defaultdict(list)
	#READ COVERAGE DATA
	with open(read_file+"_rRNA_cov.txt","r") as f:
		for l in f:
			cov[l.split()[0]].append([int(l.split()[1]),int(l.split()[2])])

	if args.db.find("unite") == -1:
		#CALCULATE MEAN COVERAGE AND CUMULATIVE COVERAGE PER SPECIES
		cov_max = ["",0]
		final_cov = defaultdict(float)
		final_cov_species = defaultdict(list)
		for g in cov:
			if g != "genome":
				avg = 0
				total = 0
				for h in cov[g]:
					avg += h[0]*h[1]
					total += h[1]
				if float(avg)/float(total) > cov_max[1]:
					cov_max = [g,float(avg)/float(total)]
				final_cov[g]=float(avg)/float(total)
				if args.db == "RDP" or args.db == "genbank":
					final_cov_species["_".join([ g.split("_")[1],g.split("_")[2].rstrip(";")])].append(float(avg)/float(total))
				elif args.db.find("ena") != -1:
					final_cov_species["_".join(g.split("|")[2].split("_")[:2])].append(float(avg)/float(total))
				elif args.db == "clete":
					final_cov_species[g].append(float(avg)/float(total))
				elif args.db == "barcodes":
					final_cov_species[g.split("|")[1]].append(float(avg)/float(total))

		final_cov_species2 = defaultdict(float)
		for s in final_cov_species:
			#print s, final_cov_species[s]
			final_cov_species2[s] = sum(final_cov_species[s])#/len(final_cov_species[s])

		if args.strain == True:
			for k in sorted(final_cov, key=final_cov.get, reverse=True):
				print k,final_cov[k]
		else:
			for k in sorted(final_cov_species2, key=final_cov_species2.get, reverse=True):
				print k,final_cov_species2[k]	
	else:
		if args.db == "unite":
			with open(sortmerna_path+"rRNA_databases/unite_tax.txt") as tax_f:
				tax_db = {}
				for l in tax_f:
					unite_id = l.strip().split("\t")[0]
					unite_species = l.strip().split("\t")[1].split(";")[-1][3:]
					tax_db[unite_id] = unite_species
				
				cov_max = ["",0]
				final_cov = defaultdict(float)
				final_cov_species = defaultdict(list)
				for g in cov:
					if g != "genome":
						avg = 0
						total = 0
						for h in cov[g]:
							avg += h[0]*h[1]
							total += h[1]
						if float(avg)/float(total) > cov_max[1]:
							cov_max = [g,float(avg)/float(total)]
						final_cov[g]=float(avg)/float(total)
						final_cov_species[tax_db[g]].append(float(avg)/float(total))
				
				final_cov_species2 = defaultdict(float)
				for s in final_cov_species:
					final_cov_species2[s] = sum(final_cov_species[s])#/len(final_cov_species[s])
						
				for k in sorted(final_cov_species2, key=final_cov_species2.get, reverse=True):
					print k,final_cov_species2[k]	
		elif args.db == "unite_big":
			cov_max = ["",0]
			final_cov = defaultdict(float)
			final_cov_species = defaultdict(list)
			for g in cov:
				if g != "genome":
					avg = 0
					total = 0
					for h in cov[g]:
						avg += h[0]*h[1]
						total += h[1]
					if float(avg)/float(total) > cov_max[1]:
						cov_max = [g,float(avg)/float(total)]
					final_cov[g]=float(avg)/float(total)
					sp = g.strip().split(";")[-1][3:].split("|")[0]
					final_cov_species[sp].append(float(avg)/float(total))
			
			final_cov_species2 = defaultdict(float)
			for s in final_cov_species:
				final_cov_species2[s] = sum(final_cov_species[s])
					
			for k in sorted(final_cov_species2, key=final_cov_species2.get, reverse=True):
				print k,final_cov_species2[k]	
				

					
