#!/usr/bin/env python

"""
Re-format the CCDS text file to bed format.
"""

import argparse

def reformatCCDS(ccds_file):
	
	for line in open(ccds_file):

		c, nc_accesion, gene, gene_id, ccds_id, ccds_status, cds_strand, cds_from, cds_to, cds_locations, match_type = line.rstrip().split('\t')
		if ccds_status == "Public":
			try: 
				loc = cds_locations[1:-1].split(', ')
				for i, l in enumerate(loc):
					print '{}\t{}\t{}\t{}:{}-{}\t1000\t{}'.format(c, l.split('-')[0], l.split('-')[1], gene, ccds_id, i+1, cds_strand)
			except IndexError: pass

if __name__ == '__main__': 

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-c', '--ccds', required=True, help="CCDS text file.")
	args = parser.parse_args()

	reformatCCDS(args.ccds)