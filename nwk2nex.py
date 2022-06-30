#!/usr/bin/env python

from Bio import Phylo
import sys

Phylo.convert(sys.argv[1], 'newick', sys.argv[1]+".nex", "nexus")
