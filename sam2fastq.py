#!/usr/bin/env python

""" sam2fastq.py """

import gzip
import argparse

def sam2single(sam_file, output_file):
	
	seen = set()

	fout = gzip.open(output_file, 'wt') if output_file.endswith('gz') else open(output_file, 'w')
	for line in open(sam_file):
		if not line.startswith('@'):
			cols = line.rstrip().split('\t')
			if cols[0] not in seen: 
				seen.add(cols[0])
				fout.write('@{}\n{}\n+\n{}\n'.format(cols[0], cols[9], cols[10]))
	fout.close()

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-s', '--sam', required=True, help="Input sam file.")
	parser.add_argument('-o', '--output', required=True, help="Output fastq file.")
	args = parser.parse_args()

	sam2single(args.sam, args.output)