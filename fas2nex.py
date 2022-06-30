#!/usr/bin/env python

from Bio import Alphabet
from Bio import SeqIO
import sys
import os

count = SeqIO.convert(sys.argv[1], "fasta", os.path.splitext(sys.argv[1])[0]+".nex", "nexus", alphabet=Alphabet.generic_protein)

