#!/usr/bin/env python

"""
	Download Ido's processed data (STAR, kallisto/SALMON). 

	storage structure:
		-> cancer_type 
			-> case_id
				-> RNA-Seq
					-> Primary_Tumor
					-> Solid_Tissue_Normal
				-> WXS
					-> Primary_Tumor
					-> Solid_Tissue_Normal
					-> Blood_Derived_Normal
			...
			
		example: /mnt/kalyna/EI_in_cancer/TCGA-BRCA/9f6be944-83de-42ab-8738-f0022f475e61/RNA-Seq/Primary_Tumor/01ea694e-989b-4a35-9397-5e508656d1d8/*
"""

import os
import sys
import argparse
import subprocess

# Hard-coded parameters
storage = "/mnt/kalyna/EI_in_cancer/"
url = "https://ngs.vbcf.ac.at/courseRNA201705/exitron/"
GDC_info_file = "/home/venhuip8/EI_in_cancer/GDC/70M_files.txt"

# Ask for IDs
parser = argparse.ArgumentParser(description="__doc__")
parser.add_argument('--IDs', required=True, help="List of file_IDs. The file_ID must be in the first column.")
args = parser.parse_args()

# Read in the GDC info (cancer_type, case_id, data_type (RNA-Seq/WXS), sample_type (control, tumor))
info = {}
for line in open(GDC_info_file):
	if not line.startswith('#'):
		cancer_type, case_id, data_type, sample_type_id, sample_type, file_id = line.rstrip().split('\t')[:6]
		sample_type = sample_type.replace(' ', '_')
		info[file_id] = { 
			'cancer_type': cancer_type, 
			'case_id': case_id, 
			'data_type': data_type,
			'sample_type': sample_type,
			'storage_folder': '{}{}/{}/{}/{}/{}/'.format(storage, cancer_type, case_id, data_type, sample_type, file_id)
		}

# Read in already downloaded IDs
try: 
	done = set([ line.rstrip().split('\t')[0] for line in open(storage+"downloaded_ids.txt") ])
except IOError: 
	done = set()

download_log = open(storage+"downloaded_ids.txt", 'a')

# Iterate over IDs and check if already downloaded
for line in open(args.IDs):

	file_ID = line.rstrip().split('\t')[0]

	# If not, download everything in the folder
	# wget -r --no-parent -P <storage> <url>
	if file_ID not in done:
		if not os.path.exists(info[file_ID]['storage_folder']): os.makedirs(info[file_ID]['storage_folder'])
		subprocess.call( "wget -r --no-parent -nH --cut-dirs 3 -P {} {}{}/".format(info[file_ID]['storage_folder'], url, file_ID), shell=True)

		# Add file_ID to download_log
		download_log.write(file_ID+'\n')

download_log.close()