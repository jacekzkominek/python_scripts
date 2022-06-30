#!/usr/bin/env python

import os, argparse

parser = argparse.ArgumentParser(description="Update local software repositories.")
args = parser.parse_args()

print "\nCheck RAxML source..."
os.chdir("/home/jacek/software/phylo/raxml-git")
ret = os.system("git pull")

print "\nCheck ExaBayes source..."
os.chdir("/home/jacek/software/phylo/exa-bayes/")
ret = os.system("git pull")

print "\nCheck MrBayes source..."
os.chdir("/home/jacek/software/phylo/mrbayes-svn")
ret = os.system("svn update")

#print "\nCheck Bcftools source..."
#os.chdir("/home/jacek/software/align/bcftools")
#ret = os.system("git pull")

#print "\nCheck Htslib source..."
#os.chdir("/home/jacek/software/align/htslib")
#ret = os.system("git pull")

#print "\nCheck Samtool source..."
#os.chdir("/home/jacek/software/align/samtools")
#ret = os.system("git pull")

#print "\nCheck Bowtie2 source..."
#os.chdir("/home/jacek/software/align/bowtie2-git")
#ret = os.system("git pull")

#print "\nCheck Guake source..."
#os.chdir("/home/jacek/software/other/guake")
#ret = os.system("git pull")

#print "\nCheck Yakuake source..."
#os.chdir("/home/jacek/software/other/yakuake")
#ret = os.system("git pull")

print "\nCheck bioseq source..."
os.chdir("/home/jacek/software/custom/bioseq-code")
ret = os.system("git pull")
