#!/usr/bin/env python

"""
	Prepare the sample_table.csv file for edgeR.
"""

import sys
from collections import Counter

count_genotype = Counter()

if len(sys.argv) == 2:
	files = sys.argv[1]
	print '"","samplePath","genotype"'
	for line in open(files):
		sample_type, sample_path = line.rstrip().split('\t')
		count_genotype[sample_type] += 1
		N = count_genotype[sample_type]

		print '"{0}{1}","{2}","{0}"'.format(sample_type, N, sample_path)