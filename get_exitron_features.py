#!/usr/bin/env python

from __future__ import division
from natsort import natsorted
from collections import Counter
import os
import re
import argparse
import subprocess
from STAR_exitrons import Sequence, yield_fasta

bedtools = "/usr/local/bin/bedtools"

def add_slash(path):
	return path if path.endswith('/') else path + '/'

def pretty_exitrons(work_dir, exitron_map):

	NTRE = re.compile('N(\d+);T(\d+)')

	seen = set()
	with open(work_dir+"exitrons_pretty.txt", 'w') as fout:
		fout.write("Gene ID\tExitron ID\tEI length (nt)\tEI length is x3\tSource\tSource label\tGene Name\n")
		for line in open(exitron_map):
			if not line.startswith('#'):
				cols = line.rstrip().split('\t')
				exitron_id = '{}:{}-{}:{}'.format(*cols[:4])

				if exitron_id not in seen:
					gene_id, gene_name = cols[5:7]
					EIx3 = "yes" if not int(cols[-1]) else "no"
					EIlen = (int(cols[2])-int(cols[1]))+1

					source = []
					N, T = map(int, re.search(NTRE, cols[9]).groups())
					if N >= 3 and T >= 3: 
						source.append("Normal")
						source.append("Tumour")
					elif N >= 3: source.append("Normal")
					elif T >= 3: source.append("Tumour")
					else: source.append("Gencode v27")

					fout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(gene_id, exitron_id, EIlen, EIx3, ','.join(source), cols[9], gene_name))
					seen.add(exitron_id)

def get_exitron_sequences(work_dir, exitron_map, genome_fasta, ccds_gtf, reference_proteins):

	from fileParser import parse_GTF

	# Create the exitron bed file
	ei, ref_trs, seen = {}, set(), set()
	with open('{}/exitron.bed'.format(work_dir), 'w') as fout:
		for line in open(exitron_map):
			if not line.startswith('#'):
				c, s, e, strand, ref_t, gene_id, gene_name = line.split('\t')[:7]
				EIx3 = line.rstrip().split('\t')[-1] == "0"
				exitron_id = '{}:{}-{}:{}'.format(c, s, e, strand)

				if EIx3 and exitron_id not in seen:
					fout.write( "{0}\t{1}\t{2}\t{3}|{4}|{5}\t1000\t{6}\n".format(c, int(s)-1, e, exitron_id, gene_id, gene_name, strand) )
					ei[exitron_id] = ref_t
					ref_trs.add(ref_t)
					seen.add(exitron_id)

	# Get the reference CCDS
	gtf = parse_GTF(ccds_gtf, select_feature="CDS", get_introns=False)

	# Get the exitrons DNA sequence
	subprocess.call("{0} getfasta -s -name -fi {1} -bed {2}/exitron.bed -fo {2}/exitron_dna.fa".format(bedtools, genome_fasta, work_dir), shell=True)

	# Read in the reference transcript proteins
	ref_prot = {}
	for record in yield_fasta(reference_proteins):
		ref_t = record.id.split('|')[1]
		if ref_t in ref_trs:
			ref_prot[ref_t] = record.seq

	with open("{}/exitron_ref_t_prot.fa".format(work_dir), 'w') as fout:
		for record in yield_fasta("{}/exitron_dna.fa".format(work_dir)):

			# Get the reference transcript sequence and translation
			full_id = record.id.split('::')[0]
			exitron_id, gene_id, gene_symbol = full_id.split('|')
			ref_t = ei[exitron_id]
			fout.write( ">{}|{}\n{}\n".format(ref_t, exitron_id, ref_prot[ref_t]) )

			# Get the exitron AA position
			c, coord, strand = exitron_id.split(':')
			ei_start, ei_end = map(int, coord.split('-'))

			ccds = gtf[ref_t]["exons"]
			ei_len = len(record.seq)

			track = 0
			if strand == '+':
				for x, y in ccds:
					if x < ei_start < y and x < ei_end < y:
						leftover = (ei_start-x) % 3
						ei_from = track + (((ei_start-x)-leftover) / 3) + 1
						ei_to = ei_from + (ei_len / 3) -1
					else:
						leftover = (y-x) % 3
						track += (((y-x)-leftover) / 3) + 1

			elif strand == '-':
				for x, y in ccds:
					if x < ei_start < y and x < ei_end < y:
						leftover = (y-ei_end) % 3
						ei_from = track + (((y-ei_end)-leftover) / 3) + 1
						ei_to = ei_from + (ei_len / 3) - 1
					else:
						leftover = (y-x) % 3
						track += (((y-x)-leftover) / 3) + 1

			print "{}\t{}\t{}\t{}\t{:.0f}\t{:.0f}".format(gene_id, gene_symbol, ref_t, exitron_id, ei_from, ei_to)

def get_protein_domains(work_dir, pfam_database, protein_fasta, exitron_AA_pos, max_e_value):
	
	#subprocess.call("hmmsearch --tblout {0}/exitron_hmmsearch_Pfam.txt -E 1e-2 --cpu 4 {1} {2} > {0}/exitron_hmmsearch_Pfam_alignments.txt".format(work_dir, pfam_database, protein_fasta), shell=True)

	# Parse exitron AA positions
	EI_pos = {}
	for line in open(exitron_AA_pos):
		gene_id, gene_symbol, ref_t, exitron_id, ei_from, ei_to = line.rstrip().split('\t')
		EI_pos[ref_t+'|'+exitron_id] = { 'gene_id': gene_id, 'gene_symbol': gene_symbol, 'exitron_id': exitron_id, 'ei_from': int(ei_from), 'ei_to': int(ei_to) }

	domains, domainCount, ei2g = {}, Counter(), {}
	domainDesc = {}

	parse, lineCount = False, 0
	domainRE = re.compile("Query:\s+(.*?)\s+\[M=(\d+)\]")
	for line in open(work_dir+"exitron_hmmsearch_Pfam_alignments.txt"):

		# Get the domain name and size
		if line.startswith('Query:'):
			domain, domain_size = re.search(domainRE, line.rstrip()).groups()

		# Get the domain description
		if line.startswith('Description:'):
			desc = line.rstrip().split('Description: ')[1]
			domainDesc[domain] = desc

		# Get the matching exitron and match position
		if line.startswith('>>'):
			full_id = line.rstrip().split('>> ')[1]
			ref_t, exitron_id = full_id.split('|')
			parse = True

		if parse and line.rstrip() == "":
			parse, lineCount = False, 0
		elif parse and lineCount >= 3:
			domain_match = re.sub('\s+', ';', line.lstrip().rstrip()).split(';')

			e_value, hmm_from, hmm_to, hmm_code, ali_from, ali_to = domain_match[5:11]
			env_from, env_to, env_code = domain_match[12:15]

			if float(e_value) <= max_e_value:

				# See if the exitron overlaps with the detected domain
				ei_from = EI_pos[full_id]['ei_from']
				ei_to = EI_pos[full_id]['ei_to']

				ei_range = set(range(ei_from, int(ei_to)+1))
				ali_range = set(range(int(ali_from), int(ali_to)+1))

				# Domain/exitron overlap label and percentage
				domain_ei_overlap = ali_range.intersection(ei_range)
				if ei_range.issubset(ali_range):
					overlap_label = "exitron_inside_domain"
					domain_ei_perc = 100
				elif ali_range.issubset(ei_range):
					overlap_label = "domain_inside_exitron"
					domain_ei_perc = 100
				else: 
					overlap_label = "part_domain_exitron"
					domain_ei_perc = (len(domain_ei_overlap) / int(domain_size)) * 100
					#if domain_ei_perc > 100: 
					#	print domain_ei_perc, domain, domain_size, ref_t, EI_pos[ref_t]["exitron_id"], env_from, env_to, ei_from, ei_to

				if len(domain_ei_overlap) and any([ overlap_label == "exitron_inside_domain", domain_ei_perc >= 20]): 

					gene_id, gene_symbol, exitron_id = [ EI_pos[full_id][x] for x in ["gene_id", "gene_symbol", "exitron_id"] ]

					output_line = "{}\t{}\t{}\t{}\t{}\t{}-{}\t{}\t{}-{}".format(gene_symbol, ref_t, exitron_id, domain, e_value, env_from, env_to, overlap_label, ei_from, ei_to)
					try: domains[exitron_id].append(output_line)
					except KeyError: domains[exitron_id] = [ output_line ]
					domainCount[domain] += 1

					ei2g[exitron_id] = gene_id
			lineCount += 1

		elif parse:
			lineCount += 1

	# Output the Pfam domain matches
	with open(work_dir+"exitron_Pfam_domains.txt", 'w') as fout:
		fout.write("Gene ID\tGene Symbol\tTranscript ID\tExitron ID\tPfam domain\tDomain E-value\tDomain location protein\tExitron location domain\tExitron location protein\n")
		for exitron_id in natsorted(domains):
			for line in domains[exitron_id]:
				fout.write( "{}\t{}\n".format(ei2g[exitron_id], line) )

	# Output the Pfam domain counts
	with open(work_dir+"exitron_Pfam_domains_counts.txt", 'w') as fout:
		fout.write("Pfam domain\tPfam domain counts\tPfam domain description\n")
		for i, d in enumerate(domainCount.most_common()):
			fout.write( "{}\t{}\t{}\n".format(d[0], d[1], domainDesc[d[0]]) )

def parse_iupred(f, ei_from, ei_to, col=2):

	# Parse disordered regions (at least 20 consecutive AA
	# with values above 0.5)

	disordered, region, start_pos = [], [], None
	lineCount = 0
	for line in open(f):
		if not line.startswith('#'):
			lineCount += 1
			pos, AA = line.split('\t')[:2]
			score = float(line.rstrip().split('\t')[col])

			if score > 0.5:
				region.append(AA)
				if start_pos == None: start_pos = int(pos)

			elif score <= 0.5 and len(region) >= 20:
				disordered.append([start_pos, int(pos)-1])
				region, start_pos = [], None

			elif score <= 0.5:
				region, start_pos = [], None

	if len(region) >= 20:
		disordered.append([start_pos, int(pos)])

	# Calculate the disordered percentage
	sumDisAA = sum([ abs(y-x) for x, y in disordered ])
	percDis = (sumDisAA/lineCount) * 100

	return disordered, percDis, any([ len(set(range(ei_from, ei_to)).intersection(set(range(x, y)))) >= 20 for x, y in disordered ])

def get_disordered_regions(work_dir, protein_fasta, exitron_AA_pos):

	iupred_path = "python3.6 ~/Apps/IUPred2A/iupred2a.py"
	vsl2b_path = "java -jar ~/Apps/VSL2/VSL2.jar"

	# Parse exitron AA positions
	EI_pos = {}
	for line in open(exitron_AA_pos):
		gene_id, gene_symbol, ref_t, exitron_id, ei_from, ei_to = line.rstrip().split('\t')
		EI_pos[ref_t+'|'+exitron_id] = { 'gene_id': gene_id, 'gene_symbol': gene_symbol, 'exitron_id': exitron_id, 'ei_from': int(ei_from), 'ei_to': int(ei_to) }

	iupred, anchor, vsl2b = {}, {}, {}

	# Iterate over the exitron fasta sequence
	for record in yield_fasta(protein_fasta):

		ei_from, ei_to = [ EI_pos[record.id][x] for x in ['ei_from', 'ei_to'] ]

		# Write individual sequence to file, since
		# IUPred2A only takes one sequence at a time
		with open(work_dir+'tmp.fa', 'w') as fout:
			fout.write( '{}\n'.format(record.seq) )

		### IUPred2A & ANCHOR2 ###
		# Run IUPred2A in short mode
		subprocess.call("{0} -a {1}tmp.fa short > {1}tmp_iupred_output.txt".format(iupred_path, work_dir), shell=True)

		# IUPred2A
		disRegions, disPerc, eiOverlap = parse_iupred(work_dir+"tmp_iupred_output.txt", ei_from, ei_to)
		if eiOverlap:
			iupred[record.id] = '{}\t{:.2f}%'.format(','.join([ '{}-{}'.format(x, y) for x, y in disRegions]), disPerc)

		# ANCHOR2
		bindRegions, bindPerc, eiOverlap = parse_iupred(work_dir+"tmp_iupred_output.txt", ei_from, ei_to, 3)
		if eiOverlap:
			anchor[record.id] = '{}\t{:.2f}%'.format(','.join([ '{}-{}'.format(x, y) for x, y in bindRegions]), bindPerc)

		### VSL2B ###
		# Run VSL in VSL2B mode
		subprocess.call("{0} -s:{1}tmp.fa > {1}tmp_vsl2b_output.txt".format(vsl2b_path, work_dir), shell=True)

		# Parse disordered regions (at least 20 
		# consecutive AA with values above 0.5)
		disRegions, region, start_pos, parse = [], [], None, False
		for line in open(work_dir+"tmp_vsl2b_output.txt"):
			if line.startswith('---'):
				parse = True

			elif parse and line.startswith('==='):
				parse = False

				# End of file check
				if len(region) >= 20:
					disRegions.append([start_pos, int(pos)])

			elif parse:

				pos, AA, score, disorder = line.rstrip().split('\t')
				if float(score) > 0.5:
					region.append(AA)
					if start_pos == None: start_pos = int(pos)

				elif float(score) <= 0.5 and len(region) >= 20:
					disRegions.append([start_pos, int(pos)-1])
					region, start_pos = [], None

				elif float(score) <= 0.5:
					region, start_pos = [], None

		sumDisAA = sum([ abs(y-x) for x, y in disRegions ])
		percDis = (sumDisAA/len(record.seq)) * 100
		if any([ len(set(range(ei_from, ei_to)).intersection(set(range(x, y)))) >= 20 for x, y in disRegions ]):
			vsl2b[record.id] = '{}\t{:.2f}%'.format(','.join([ '{}-{}'.format(x, y) for x, y in disRegions ]), percDis)

	# Remove IUPred2A and VSL2B temporary files
	subprocess.call("rm {}tmp*".format(work_dir), shell=True)

	# Output all disordered regions
	with open(work_dir+"exitron_disordered_regions.txt", 'w') as fout:
		fout.write("Gene ID\tGene Symbol\tTranscript ID\tExitron ID\tExitron location protein\tIUPred2A disordered regions\tIUPred2A disordered percentage\tANCHOR2 binding regions\tANCHOR2 binding percentage\tVSL2B disordered regions\tVSL2B disordered percentage\n")
		for full_id in natsorted(EI_pos):
			if any([ full_id in iupred, full_id in anchor, full_id in vsl2b ]):
				ref_t, exitron_id = full_id.split('|')
				fout.write( "{}\t{}\t{}\t{}\t{}-{}\t{}\t{}\t{}\n".format(EI_pos[full_id]['gene_id'], EI_pos[full_id]['gene_symbol'], ref_t, exitron_id, EI_pos[full_id]['ei_from'], EI_pos[full_id]['ei_to'], \
					iupred[full_id] if full_id in iupred else 'none_found\tNA', \
					anchor[full_id] if full_id in anchor else 'none_found\tNA',
					vsl2b[full_id] if full_id in vsl2b else 'none_found\tNA') )

def get_summary(work_dir, exitron_map, cgc_file):
	
	"""
		* Number of exitrons detected (annotated (yes/no), origin (normal/tumor/both/neither))
		* Size distribution (all/normal/tumor/both)
		* GC distribution (all/normal/tumor/both)
		* Cancer Gene Consensus (GSC) genes detected
		* Protein domains found
		* Disordered regions
	"""

	NTRE = re.compile('N(\d+);T(\d+)')

	# Assumes that all previous functions have been run 
	# and all files are in the working directory
	if all([ os.path.isfile(work_dir+f) for f in ['exitron_dna.fa', 'exitron_Pfam_domains.txt', 'exitron_disordered_regions.txt'] ]):
		
		# Number of exitrons detected (annotated => yes/no, origin => normal/tumor/both/neither)
		seen, EI, g2ei, origin = set(), {}, {}, { 'yes': Counter(), 'no': Counter() }
		for line in open(exitron_map):
			if not line.startswith('#'):
				cols = line.rstrip().split('\t')
				exitron_id = '{}:{}-{}:{}'.format(*cols[:4])

				if exitron_id not in seen:

					# Give the exitron some origin label based on the
					# number of occurences. An exitron is seen as "present"
					# if it is found in at least 3 replicates of a genotype
					N, T = map(int, re.search(NTRE, cols[9]).groups())

					if N >= 3 and T >= 3: label = "Normal & Tumour"
					elif N >= 3: label = "Normal"
					elif T >= 3: label = "Tumour"
					else: label = "gencode"
					origin[cols[10]][label] += 1

					EI[exitron_id] = { 'gene_name': cols[6], 'is_annotated': cols[10], 'origin': label,'pfam': [], 'disorder': [] }
					g2ei[cols[6]] = exitron_id

					seen.add(exitron_id)

		print "### Number of detected exitrons ###"
		print "# The 'yes' and 'no' indicate whether an exitron was annotated\n# as an intron in the gencode.v27 annotation. An exitron is assigned\n# to the normal, tumor, or both categories if junctions where\n# detected in at least 3 replicates for the respective category.The\n# exitrons in the 'gencode' column are those previously annotated for\n# which not sufficient evidence was detected in this sample.\n"
		print "\tnormal\ttumor\tnormal&tumor\tgencode\ttotal"
		print "yes\t{}\t{}\t{}\t{}\t".format(*[ origin["yes"][x] for x in ["Normal", "Tumour", "Normal & Tumour", "gencode"] ])+str(sum(origin["yes"].values()))
		print "no\t{}\t{}\t{}\t{}\t".format(*[ origin["no"][x] for x in ["Normal", "Tumour", "Normal & Tumour", "gencode"] ])+str(sum(origin["no"].values()))

		# Exitron size and GC distribution
		with open(work_dir+"exitron_length_gc.txt", 'w') as fout:
			fout.write( 'IS_ANNOTATED\tORIGIN\tLENGTH\tGC\n' )
			for record in yield_fasta(work_dir+"exitron_dna.fa"):
				exitron_id = record.id.split('|')[0]
				fout.write( '{}\t{}\t{}\t{}\t{:.3f}\n'.format(exitron_id, EI[exitron_id]['is_annotated'], EI[exitron_id]['origin'], len(record.seq), Sequence(record.seq).get_base_content('gc')) )

		# Cancer Gene Census
		CGC = {}
		print '\n### Cancer Gene Census ###'
		print '# Exitrons found in the CGC genes.'
		print '# Exitron ID\tIs annotated?\tFound in\tGene Symbol\tGene Name'
		for line in open(cgc_file):
			gene_symbol, gene_name = line.split('\t')[:2]
			#tumour_type = line.split('\t')[9]

			if gene_symbol in g2ei:
				exitron_id = g2ei[gene_symbol]
				print '{}\t{}\t{}\t{}\t{}'.format(exitron_id, EI[exitron_id]['is_annotated'], EI[exitron_id]['origin'], gene_symbol, gene_name)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-w', '--work-dir', required=True, help="Output working directory.")

	subparsers = parser.add_subparsers(dest="command", help="sub-command help.")

	parser_a = subparsers.add_parser("get-sequences", help="Get the exitron DNA sequences.")
	parser_a.add_argument('-e', '--exitron-map', required=True, help="Exitron CCDS map.")
	parser_a.add_argument('-g', '--genome-fasta', required=True, help="Genome fasta.")
	parser_a.add_argument('-c', '--ccds', required=True, help="CCDS gtf file.")
	parser_a.add_argument('-r', '--reference-proteins', required=True, help="Reference transcript translations.")

	parser_b = subparsers.add_parser("get-domains", help="Get Pfam protein domains in the exitrons.")
	parser_b.add_argument('-p', '--pfam', required=True, help="Pfam hmm database.")
	parser_b.add_argument('-f', '--fasta', required=True, help="Exitron protein fasta file.")
	parser_b.add_argument('--exitron-pos', required=True, help="Exitron position on reference transcript.")
	parser_b.add_argument('-e', '--e-value', default="0.001", type=float, help="Max e-value cut-off.")

	parser_c = subparsers.add_parser("get-disordered-regions", help="Identify disordered regions with IUPred2A and VSL2B.")
	parser_c.add_argument('-f', '--fasta', required=True, help="Exitron protein fasta file.")
	parser_c.add_argument('--exitron-pos', required=True, help="Exitron position on reference transcript.")

	parser_d = subparsers.add_parser("get-summary", help="Summarize all the obtained exitron information.")
	parser_d.add_argument('-e', '--exitron-map', required=True, help="Exitron CCDS map.")
	parser_d.add_argument('--CGC', required=True, help="CGC file in TSV format.")

	parser_e = subparsers.add_parser("pretty-exitrons", help="Pretty exitrons.")
	parser_e.add_argument('-e', '--exitron-map', required=True, help="Exitron CCDS map.")

	args = parser.parse_args()

	work_dir = add_slash(args.work_dir)
	if not os.path.exists(work_dir): os.makedirs(work_dir)

	if args.command == "get-sequences":
		get_exitron_sequences(work_dir, args.exitron_map, args.genome_fasta, args.ccds, args.reference_proteins)
	elif args.command == "get-domains":
		get_protein_domains(work_dir, args.pfam, args.fasta, args.exitron_pos, args.e_value)
	elif args.command == "get-disordered-regions":
		get_disordered_regions(work_dir, args.fasta, args.exitron_pos)
	elif args.command == "get-summary":
		get_summary(work_dir, args.exitron_map, args.CGC)
	elif args.command == "pretty-exitrons":
		pretty_exitrons(work_dir, args.exitron_map)