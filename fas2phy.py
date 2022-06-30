#!/usr/bin/env python

from Bio import SeqIO
import sys

SeqIO.write(SeqIO.parse(sys.argv[1], 'fasta'), sys.argv[1].replace("fas","phy"), "phylip-relaxed")
