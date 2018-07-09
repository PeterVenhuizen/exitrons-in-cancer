#!/usr/bin/env python

import os
import argparse
import subprocess
from natsort import natsorted
from fileParser import yield_salmon

suppa_path = "/home/venhuip8/Apps/SUPPA/suppa.py"

def add_slash(path):
	return path if path.endswith('/') else path + '/'

def calculate_PSI(work_dir, expr_file, ioe_files):
	
	condition = expr_file.split('/')[-1].split('.')[0]

	for ioe in ioe_files:

		event_type = ioe.split('_')[1]
		event_dir = '{}{}/'.format(work_dir, event_type)
		if not os.path.exists(event_dir): os.makedirs(event_dir)

		output_file = "{}{}".format(event_dir, condition)

		subprocess.call("python3.6 {} psiPerEvent --ioe {} --expression-file {} -o {}".format(suppa_path, ioe, expr_file, output_file), shell=True)

def combine_TPM_files(work_dir, files, condition="normal"):
	
	#expr = [ { record['transcript_id']: record['tpm'] for record in yield_salmon(line.rstrip().split('\t')[1]+"quant.sf") for f in files if line.split('\t')[0] == condition } ]
	expr, samples = [], []
	for line in open(files):
		sample_type, sample_path = line.rstrip().split('\t')
		file_id = sample_path.split('/')[-1] if sample_path.split('/')[-1] != '' else sample_path.split('/')[-2]

		if sample_type == condition:
			expr.append( { record['transcript_id']: record['tpm'] for record in yield_salmon(sample_path+"quant.sf") } )
			samples.append(file_id)

	with open("{}{}.expr".format(work_dir, condition), 'w') as fout:
		fout.write( '\t'.join(samples)+'\n' )
		for t_id in natsorted(expr[0]):
			fout.write( '{}\t{}\n'.format(t_id, '\t'.join([ str(x[t_id]) for x in expr ])) )

def calculate_dPSI(work_dir, ioe_files):
	
	for ioe in ioe_files:

		event_type = ioe.split('_')[1]
		event_dir = '{}{}/'.format(work_dir, event_type)

		subprocess.call("python3.6 {0} diffSplice --method empirical --psi {2}normal.psi {2}tumor.psi --tpm {1}normal.expr {1}tumor.expr --input {3} --paired --gene-correction --output {2}normal_vs_tumor".format(suppa_path, work_dir, event_dir, ioe), shell=True)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-w', '--work-dir', required=True, help="Output working directory.")

	subparsers = parser.add_subparsers(dest="command", help="sub-command help.")

	parser_a = subparsers.add_parser("combine-salmon", help="Combine a bunch of Salmon TPM files into tab-separated format.")
	parser_a.add_argument('-f', '--files', required=True, help="Text file containing the sample type (tumor or normal) in the first column and the sample full directory path in the second column.")
	parser_a.add_argument('-c', '--condition', choices=["normal", "tumor"], default="normal", help="Select the condition to work with.")

	parser_b = subparsers.add_parser("calculate-PSI", help="Calculate the PSI for condition X and the provided ioe files.")
	parser_b.add_argument('--expr', required=True, help="Expression (salmon/kallisto TPM) file (from combine-salmon).")
	parser_b.add_argument('--ioe', required=True, nargs='+', help="SUPPA event files to use for quantification.")

	parser_c = subparsers.add_parser("calculate-dPSI", help="Calculate the dPSI for the provided ioe files.")
	parser_c.add_argument('--ioe', required=True, nargs='+', help="SUPPA event files to use for quantification.")

	args = parser.parse_args()

	work_dir = add_slash(args.work_dir)
	if not os.path.exists(work_dir): os.makedirs(work_dir)

	if args.command == "combine-salmon":
		combine_TPM_files(work_dir, args.files, args.condition)

	elif args.command == "calculate-PSI":
		calculate_PSI(work_dir, args.expr, args.ioe)

	elif args.command == "calculate-dPSI":
		calculate_dPSI(work_dir, args.ioe)