#!/usr/bin/env python

"""
	Traverse all folders of a cancer type and output
	the paths to some text file.

	cancer_type / case_id / data_format / sample_type / file_id
"""

import os
import sys

if len(sys.argv) == 2:
	path = sys.argv[1]
	for dirpath, dirs, files in os.walk(path):
		if len(dirs) == 0: 
			if all(["star.SJ.out.tab" in files, "star.Aligned.sortedByCoord.out.bam" in files, "kallisto_abundance.tsv" in files, "quant.sf" in files]):
				if 'Normal' in dirpath: print 'normal\t{}/'.format(dirpath)
				elif 'Tumor' in dirpath: print 'tumor\t{}/'.format(dirpath)
				else: sys.stderr.write("Something is wrong!!! {} is neither normal or tumor...".format(dirpath))
			else: print 'ERROR\t{}/'.format(dirpath)