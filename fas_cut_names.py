#!/usr/bin/env python

from Bio import Alphabet
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os

seqs = []
for s in SeqIO.parse(sys.argv[1], "fasta"):
	from Bio.SeqRecord import SeqRecord
	#~ print s.id
	#~ print s.name
	#~ print s.description
	s2 = SeqRecord(s.seq, id=s.id[0:90], description="", name="")
	seqs.append(s2)

SeqIO.write(seqs, sys.argv[1]+"_cut.fas", "fasta")

