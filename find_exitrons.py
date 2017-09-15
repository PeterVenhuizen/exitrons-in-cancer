#!/usr/bin/env python

"""
	For human adopted version of exitron discovery pipeline.
	Requires samtools (tested with v1.2) and bedtools (tested with v2.26.0)
"""

from __future__ import division
from natsort import natsorted
import argparse
import subprocess
import sys
import os

from fileParser import yield_junctions, yield_fasta, yield_bed
from sequence import Sequence

def add_slash(path):
	return path if path.endswith('/') else path +'/'

def get_shared_junctions(junc_files, min_support=3):
	return set.intersection(*[ set([ (j['chr'], j['junc_start'], j['junc_end'], j['strand']) for j in yield_junctions(f) if j['depth'] >= min_support ]) for f in junc_files ])

def prepare_junction_lib(work_dir, ccds_gtf, junc_files, sample_names, exitron_bed=None, min_support=3):
	
	# Generate working directory
	if not os.path.exists(work_dir): os.makedirs(work_dir)

	# Get the shared sample junctions
	junction_sets = []
	uniq_names = set(sample_names)
	for n in natsorted(uniq_names):
		indices = [ i for i, x in enumerate(sample_names) if x == n ]
		junction_sets.append( get_shared_junctions([ junc_files[i] for i in indices ], min_support) )
	all_junctions = set.union(*junction_sets)

	# Add existing exitrons when supplied
	try: 
		if os.path.isfile(exitron_bed):
			existing_EI = set([ ( r['chr'], r['start'], r['end']-1, r['strand'] ) for r in yield_bed(exitron_bed) ])
			all_junctions = all_junctions.union(existing_EI)
	except TypeError: pass

	# Output shared junctions
	with open('{}all_junctions.min_{}_depth.bed'.format(work_dir, min_support), 'w') as fout:
		for j in natsorted(all_junctions):
			c, s, e, strand = j
			fout.write( '{}\t{}\t{}\tJUNCTION\t1000\t{}\n'.format(c, s, e, strand) )

	# intersectBed -s -f 1 -wa -wb -a all_junctions.bed -b ccds_gtf
	# -f 1 reports only 100% matches
	subprocess.call("intersectBed -s -f 1 -wa -wb -a {0}all_junctions.min_{1}_depth.bed -b {2} | awk -F\"\\t\" '{{ OFS=\"\t\"; print $1, $2, $3, $10, $11, $NF }}' > {0}junctions_CCDS_map.tmp".format(work_dir, min_support, ccds_gtf), shell=True)

	# Parse intersection
	with open('{}junctions_CCDS_map.txt'.format(work_dir), 'w') as fout:
		fout.write('#CHR\tJUNCTION_START\tJUNCTION_END\tTRANSCRIPT_ID\tCCDS_START\tCCDS_END\n')
		for line in open('{}junctions_CCDS_map.tmp'.format(work_dir)):
			c, j_start, j_end, c_start, c_end, attributes = line.rstrip().split('\t')
			attr = {}
			for a in attributes.split(';'):
				if len(a):
					attr_name, attr_value = a.split(' "')
					attr[attr_name.strip()] = attr_value.replace('\"', '')

			fout.write( '{}\t{}\t{}\t{}\t{}\t{}\n'.format(c, j_start, j_end, attr["transcript_id"], c_start, c_end) )
	subprocess.call("rm -f {0}junctions_CCDS_map.tmp".format(work_dir), shell=True)

def filter_exitrons(work_dir, ccds_gtf, introns_bed, junction_map, genome_fasta):
	
	# Get all annotated introns
	ann_introns = set([ (i['chr'], i['start'], i['end']) for i in yield_bed(introns_bed) ])

	# Get the CCDS annotation
	ccds = {}
	for line in open(ccds_gtf):
		if not line.startswith('#'):
			c, source, feature, start, end, score, strand, frame, attributes = line.rstrip().split('\t')
			attr = {}
			for a in attributes.split(';'):
				if len(a):
					attr_name, attr_value = a.split(' "')
					attr[attr_name.strip()] = attr_value.replace('\"', '')

			try: ccds[attr["transcript_id"]]["exons"].append([ int(start), int(end) ])
			except KeyError: ccds[attr["transcript_id"]] = { "chr": c, "strand": strand, "exons": [[ int(start), int(end) ]] }
	for t_id in ccds: ccds[t_id]['exons'].sort(key=lambda x: x[0])

	with open('{}exitrons_CCDS_map.txt'.format(work_dir), 'w') as map_out:
		map_out.write( '#CHR\tEXITRON_START\tEXITRON_END\tTRANSCRIPT_ID\tCCDS_START\tCCDS_END\tIS_ANNOTATED\tMOD3\n' )
		for line in open(junction_map):
			if not line.startswith('#'):
				c, j_start, j_end, t_id, c_start, c_end = line.rstrip().split('\t')
				j_start, j_end = int(j_start), int(j_end)
				mod3 = abs(j_end-j_start) % 3

				if (c, j_start, j_end) in ann_introns:
					
					# Generate the spliced in bed file
					with open("{}tmp_spliced_in.bed".format(work_dir), 'w') as fout:
						for s, e in ccds[t_id]["exons"]:
							fout.write( '{}\t{}\t{}\n'.format(c, s-1, e) )

					# Generate exitrons spliced out bed file
					with open("{}tmp_spliced_out.bed".format(work_dir), 'w') as fout:
						for s, e in ccds[t_id]["exons"]:
							if s <= j_start <= e and s <= j_end <= e:
								fout.write( '{}\t{}\t{}\n'.format(c, s-1, j_start) )
								fout.write( '{}\t{}\t{}\n'.format(c, j_end, e) )
							else:
								fout.write( '{}\t{}\t{}\n'.format(c, s-1, e) )

					# Get sequences
					subprocess.call("bedtools getfasta -fi {0} -bed {1}tmp_spliced_in.bed -fo {1}tmp_spliced_in.fa".format(genome_fasta, work_dir), shell=True)
					subprocess.call("bedtools getfasta -fi {0} -bed {1}tmp_spliced_out.bed -fo {1}tmp_spliced_out.fa".format(genome_fasta, work_dir), shell=True)
					i = ''.join([ record.seq for record in yield_fasta("{}tmp_spliced_in.fa".format(work_dir)) ])
					o = ''.join([ record.seq for record in yield_fasta("{}tmp_spliced_out.fa".format(work_dir)) ])

					# Get reverse complement if necessary
					if ccds[t_id]['strand'] == '-':
						i = Sequence(i).get_reverse_complement()
						o = Sequence(o).get_reverse_complement()

					# Get the protein sequences
					inProtein, in_leftover = Sequence(i).translate()
					inProtein = inProtein[:inProtein.index('*')] if '*' in inProtein else inProtein[::]

					outProtein, out_leftover = Sequence(o).translate()
					outProtein = outProtein[:outProtein.index('*')] if '*' in outProtein else outProtein[::]

					# Check if the exitron containing protein (i.e. spliced in variant)
					# is at least of the same length, or longer, than the spliced out
					# variant and that the stop codon is at the same, or on a later 
					# position as with the spliced out variant. 
					keep = False

					if len(inProtein) >= len(outProtein):
						if len(inProtein) == len(outProtein) + abs(j_end-j_start)/3: keep = True

					if keep:
						map_out.write( '{}\tyes\t{}\n'.format(line.rstrip(), mod3) )

				else: 
					map_out.write( '{}\tno\t{}\n'.format(line.rstrip(), mod3) )

	subprocess.call("rm -f {}tmp_spliced*".format(work_dir), shell=True)

def prepare_bam_files(work_dir, bam_file, handle, genome_fasta):
	""" Extract the unique reads and unique exonic reads from the 
	supplied bam file, index and output to the working directory """

	# Create work_dir if not exists
	if not os.path.exists(work_dir): os.makedirs(work_dir)

	# Extract unique reads
	uniq_reads_bam = '{}accepted_hits_uniq_reads.{}.bam'.format(work_dir, handle)
	cmd = "samtools view %s | grep -w \"NH:i:1\" | perl -n -e '@line=split(/\\t/,$_); if ($line[5]=~/N/){ print \"$_\"; }' > %s" % (bam_file, uniq_reads_bam.replace('.bam', '.sam'))
	subprocess.call(cmd, shell=True)

	# Convert sam to bam
	subprocess.call("samtools view -bT {} {} > {}".format(genome_fasta, uniq_reads_bam.replace('.bam', '.sam'), uniq_reads_bam), shell=True)

	# Index bam
	subprocess.call("samtools index {}".format(uniq_reads_bam), shell=True)

	# Remove sam
	subprocess.call("rm -f %s" % (uniq_reads_bam.replace('.bam', '.sam')), shell=True)

	# Extract unique exonic reads
	uniq_exon_reads_bam = '{}accepted_hits_uniq_exonic_reads.{}.bam'.format(work_dir, handle)
	cmd = "samtools view %s | grep -w \"NH:i:1\" | perl -n -e '@line=split(/\\t/,$_); if ($line[5]!~/N/){ print \"$_\"; }' > %s" % (bam_file, uniq_exon_reads_bam.replace('.bam', '.sam'))
	subprocess.call(cmd, shell=True)

	# Convert sam to bam
	subprocess.call("samtools view -bT {} {} > {}".format(genome_fasta, uniq_exon_reads_bam.replace('.bam', '.sam'), uniq_exon_reads_bam), shell=True)

	# Index bam
	subprocess.call("samtools index {}".format(uniq_exon_reads_bam), shell=True)

	# Remove sam
	subprocess.call("rm -f %s" % (uniq_exon_reads_bam.replace('.bam', '.sam')), shell=True)

def calculate_PSI(work_dir, exitron_map, uniq_reads, uniq_exonic, handle):
	''' Calculate exitron PSI values, based on the coverage of the 
	unique exonic reads. 
						
						A		   B		 C
					   ----		  ----		----
	EEEEEEEEEEEEEEEEEEEEEXXXXXXXXXXXXXXXXXXXXXEEEEEEEEEEEEEEEEEEEEE
					   --	                  --
						 \                   / 
						  \        D        / 							
                           -----------------
                           	
    E = exon, X = exitron
    A = Reads aligning from -10 to +10 around exitron 5'SS
    B = Reads aligning from -10 to +10 around exitron middle point
    C = Reads aligning from -10 to +10 around exitron 3'SS
    D = Reads supporting exitron splicing, output from get_junction_uniq_read_support()
    
    PSI = ( ((A+B+C)/3) / ((A+B+C)/3) + D ) * 100 '''

    import re

    # Get junction support
    rc, info = {}, {}
    for line in open(exitron_map):
    	if not line.startswith('#'):

    		# CHR, EXITRON_START, EXITRON_END, TRANSCRIPT_ID, CCDS_START, CCDS_END, IS_ANNOTATED, MOD3
    		c, s, e, strand, t_id, ccds_start, ccds_end, is_annotated, mod3 = line.rstrip().split('\t')
    		info['{}:{}-{}'.format(c, s, e)] = { 'strand': strand, 't_id': t_id, 'is_annotated': is_annotated, 'mod3': mod3 }

    		subprocess.call("samtools view {0}{1} {2}:{3}-{4} > {0}tmp.sam".format(work_dir, uniq_reads, c, s, e), shell=True)
    		N = "{}N".format(int(e)-int(s)) # Get the required read/junction gap signature

    		uniq_count = 0
    		for aln in open("{}tmp.sam".format(work_dir)):
    			qname = aln.split('\t')[0]
    			pos = int(aln.split('\t')[3])
    			cigar = aln.split('\t')[5]
    			start = (pos + int(re.search('^[0-9]+)M', cigar).group(1))) - 1

    			# Check if the junction is at the correct position
    			# and if the junction size is correct
    			if N in cigar and start == int(s): uniq_count += 1

    		rc['{}:{}-{}'.format(c, s, e)+1] = { 'A': 0, 'B': 0, 'C': 0, 'D': uniq_count }

    # Get A, B, C read support
    with open("{}tmp_coverageBed_input.bed", 'w') as fout:
    	for line in open(exitron_map):
    		if not line.startswith('#'):

				c, s, e, strand, t_id, ccds_start, ccds_end, is_annotated, mod3 = line.rstrip().split('\t')
				s, e = int(s), int(e)
				middle_point = int(s + ((e-s)/2))
				locus = '{}:{}-{}'.format(c, s, e+1)

				fout.write( '{}\t{}\t{}\t{}_A\n'.format(c, s-10, s+10, locus) )
				fout.write( '{}\t{}\t{}\t{}_B\n'.format(c, middle_point-10, middle_point+10, locus) )
				fout.write( '{}\t{}\t{}\t{}_C\n'.format(c, e-10, e+10, locus) )

	subprocess.call("sort -k1,1 -k2,2n {0}tmp_coveragebed_input.bed > {0}tmp_coverageBed_input.sorted.bed".format(work_dir), shell=True)
	subprocess.call("coverageBed -sorted -counts -a {0}tmp_coverageBed_input.sorted.bed -b {0}{1} > {0}tmp_coverageBed_output.bed".format(work_dir, uniq_exonic), shell=True)
	for line in open("{}tmp_coverageBed_output.bed".format(work_dir)):
		c, start, end, locus, coverage = line.rstrip().split('\t')
		locus, letter = locus.split('_')
		rc[locus][letter] = int(coverage)

	# Calculate PSI
	with open("{}{}.PSI".format(work_dir, handle), 'w') as fout: 
		for x in natsorted(rc):
			try: PSI = ( ( ( rc[x]['A'] + rc[x]['B'] + rc[x]['C'] ) / 3 ) / ( ( ( rc[x]['A'] + rc[x]['B'] + rc[x]['C'] ) / 3 ) + rc[x]['D'] ) ) * 100
			except ZeroDivisionError: PSI = 'NA'
			fout.write( '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(x, info[x]['t_id'], info[x]['strand'], info[x]['is_annotated'], info[x]['mod3'], rc[x]['A'], rc[x]['B'], rc[x]['C'], rc[x]['D'], PSI) )

	# Clean up
	subprocess.call( "rm -f {}tmp*".format(work_dir), shell=True )

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-w', '--work-dir', required=True, help="Output working directory.")

	subparsers = parser.add_subparsers(dest='command', help="sub-command help")

	parser_a = subparsers.add_parser('prepare-junctions', help="Select all junctions with a minimum read support (default is 3) and map to the CCDS annotation.")
	parser_a.add_argument('-c', '--ccds', required=True, help="CCDS gtf file.")
	parser_a.add_argument('-j', '--junctions', required=True, nargs='+', help="Junction bed files.")
	parser_a.add_argument('-n', '--names', required=True, nargs='+', help="List of condition names, e.g. N N N D D D.")
	parser_a.add_argument('-e', '--exitrons', default=None, help="Existing exitron bed file.")
	parser_a.add_argument('-m', '--min-support', type=int, default=3, help="Minimum required junction coverage depth.")

	parser_b = subparsers.add_parser('filter-exitrons', help="Filter out intron retention events from the potential exitrons, based on the coding potential.")
	parser_b.add_argument('-c', '--ccds', required=True, help="CCDS bed file.")
	parser_b.add_argument('-i', '--introns', required=True, help="Transcriptome introns bed file.")
	parser_b.add_argument('-j', '--junction-map', required=True, help="Junction map file (from prepare-junctions).")
	parser_b.add_argument('-g', '--genome-fasta', required=True, help="Genome fasta.")

	parser_c = subparsers.add_parser('prepare-bam', help="Extract the unique (exonic) reads from the TopHat2 mapping bam file.")
	parser_c.add_argument('-b', '--bam', required=True, help="TopHat2 bam mapping file.")
	parser_c.add_argument('-f', '--file-handle', required=True, help="Unique file handle. The output files will be [work_dir]/accepted_hits_uniq_reads.[file-handle].bam and [work-dir]/accepted_hits_uniq_exonic_reads.[file-handle].bam.")
	parser_c.add_argument('-g', '--genome-fasta', required=True, help="Genome fasta.")

	parser_d = subparsers.add_parser('calculate-PSI', help="Calculate the exitron PSI.")
	parser_d.add_argument('--exitron-map', required=True, help="Exitron mapping file (from filter-exitrons).")
	parser_d.add_argument('--uniq-reads', required=True, help="Unique reads bam file (from prepare-bam).")
	parser_d.add_argument('--uniq-exonic', required=True, help="Unique exonic reads bam file (from prepare-bam).")
	parser_d.add_argument('-f', '--file-handle', required=True, help="Unique file handle. The output files will be [work_dir]/accepted_hits_uniq_reads.[file-handle].bam and [work-dir]/accepted_hits_uniq_exonic_reads.[file-handle].bam.")

	args = parser.parse_args()

	# Create work_dir if not exists
	if not os.path.exists(args.work_dir): os.makedirs(args.work_dir)

	if args.command == "prepare-junctions":
		if len(args.junctions) == len(args.names):
			prepare_junction_lib(add_slash(args.work_dir), args.ccds, args.junctions, args.names, args.exitrons, args.min_support)
		else: 
			sys.stderr.write("The number of junction files (-j, --junctions) and sample names (-n, --names) should be the same! Quitting...")
	elif args.command == "filter-exitrons":
		filter_exitrons(add_slash(args.work_dir), args.ccds, args.introns, args.junction_map, args.genome_fasta)
	elif args.command == "prepare-bam":
		prepare_bam_files(add_slash(args.work_dir), args.bam, args.file_handle, args.genome_fasta)
	elif args.command == "calculate-PSI":
		if all([ os.path.exists(args.exitron_map), os.path.exists(args.uniq_reads), os.path.exists(args.uniq_exonic) ]):
			calculate_PSI(add_slash(args.work_dir), args.exitron_map, args.uniq_reads, args.uniq_exonic, args.file_handle)
		else: 
			sys.stderr.write("One (or multiple) input file could not be accessed.\n")