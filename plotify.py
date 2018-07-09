#!/usr/bin/env python

"""
	Plot the PSI values for the paired normal and 
	tumor tissues. 
"""

import argparse

def get_paired_values(files_file, PSI_files, exitron_ID):
	
	samples = {}
	f2c = {}
	for line in open(files_file):
		genotype, path = line.rstrip().split('\t')
		case_id = path.split('/')[5]
		file_id = path.split('/')[8]

		f2c[file_id] = { 'case_id': case_id, 'genotype': genotype }
		samples[case_id] = { 'normal': None, 'tumor': None }

	for f in PSI_files:
		file_id = f.split('/')[-1].split('.')[1]		
		for line in open(f):
			if exitron_ID in line:
				try: 
					PSI = float(line.rstrip().split('\t')[-1])
				except ValueError: 
					PSI = float('nan')

				case_id, genotype = [ f2c[file_id][x] for x in ['case_id', 'genotype'] ]
				samples[case_id][genotype] = PSI

	print 'SAMPLE\tGENOTYPE\tPSI'
	for i, x in enumerate(samples):
		print '{}\t{}\t{}'.format(i, samples[x]['normal'], samples[x]['tumor'])


if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-f', '--files', required=True, help="Text file containing the sample type (tumor or normal) in the first column and the sample full directory path in the second column.")
	parser.add_argument('--psi', required=True, nargs='+', help="PSI values")
	parser.add_argument('--exitron-ID', required=True, help="Exitron ID")

	args = parser.parse_args()

	get_paired_values(args.files, args.psi, args.exitron_ID)