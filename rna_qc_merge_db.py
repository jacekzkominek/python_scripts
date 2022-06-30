#!/usr/bin/env python

import argparse
from collections import defaultdict
import subprocess
import os

parser = argparse.ArgumentParser(description="Calculate average coverage across features")
#~ parser.add_argument("f1", help="Input BED file")
#~ parser.add_argument("--norna", default=False, action="store_true", help="Don't run SortMeRNA")
#~ parser.add_argument("--local", default=False, action="store_true", help="Local alignment with Bowtie2")
#~ parser.add_argument("--bwa", default=False, action="store_true", help="Alignment with BWA")
#~ parser.add_argument("--strain", default=False, action="store_true", help="Strain level")
#~ parser.add_argument("--db", default="rdp", choices=["rdp","gb","unite","unite_big"], help="RNA database to use")
#~ parser.add_argument("--cpu", default="12", help="CPUs to use")
args = parser.parse_args()

