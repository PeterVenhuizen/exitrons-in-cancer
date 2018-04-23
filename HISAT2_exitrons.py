#!/usr/bin/env python

"""
	HISAT2 aligner exitron identification and quantification pipeline.
	Requires samtools (tested with v1.2) and bedtools (tested with v.2.26.0)
"""

from __future__ import division
from collections import Counter
from natsort import natsorted
import argparse
import subprocess
import sys
import os

samtools_path = "/home/venhuip8/Apps/samtools-1.5/bin/samtools"
bedtools_path = "/usr/local/bin/bedtools"

def add_slash(path):
	return path if path.endswith('/') else path + '/'

def yield_fasta(f):
	''' Simple fasta parser that yield's the identifier and sequence of each record '''
	
	class SeqRecord:
		def __init__(self, seq_id, seq):
			self.id = seq_id
			self.seq = seq

	seq = ''
	for line in open(f):
		if line.startswith('>'):
			if len(seq):
				yield(SeqRecord(identifier, seq))
				seq = ''
			identifier = line[1:].rstrip()
		else: seq += line.rstrip()
	yield(SeqRecord(identifier, seq))

def yield_junctions(f):
	""" Yield each line as a dictionary. """

	for line in open(f):
		c, s, e, strand = line.rstrip().split('\t')
		yield({
			'chr': c,
			'start': int(s)+2,
			'end': int(e),
			'strand': strand
			})

def yield_bed(f):

	# https://genome.ucsc.edu/FAQ/FAQformat#format1
	bed_fields = [ 'chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts' ]

	for line in open(f):
		cols = line.rstrip().split('\t')
		d = { bed_fields[i]: cols[i] if i in xrange(len(cols)) else None for i in xrange(len(bed_fields)) }
		try: 
			d['start'] = int(d['start'])
			d['end'] = int(d['end'])
		except KeyError:
			print "Missing start or end information. Aborting..."
			sys.exit()
		yield(d)

class Sequence(object):
	def __init__(self, seq, seq_id="", start=0, end=0):
		self.sequence = seq.upper()
		self.id = seq_id
		self.start = int(start)
		self.end = int(end)

	def get_length(self):
		''' Return Sequence object length. '''
		
		if len(self.sequence) > 0:
			length = len(self.sequence)
		else:
			length = (self.end-self.start)+1
		
		return length

	def get_reverse_complement(self):
		''' Get the reverse complement of the sequence. '''
		
		com_dict = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A', 
					'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 
					'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 
					'N': 'N' }
		return ''.join([com_dict[self.sequence[-i]] for i in xrange(1, len(self.sequence)+1)])

	def translate(self, frame=0):
		''' Translate to protein '''
		
		seq = self.sequence.replace('U', 'T')
		
		aa_dict = { ("GCT", "GCC", "GCA", "GCG"): 'A', # Alanine (Ala/A)
			("TGT", "TGC"): 'C', # Cysteine (Cys/C)
			("GAT", "GAC"): 'D', # Aspartic acid (Asp/D)
			("GAA", "GAG"): 'E', # Glumatic acid (Glu/E)
			("TTT", "TTC"): 'F', # Phenylalanine (Phe/F)
			("GGT", "GGC", "GGA", "GGG"): 'G', # Glycine (Gly/G)
			("CAT", "CAC"): 'H', # Histidine (His/H)
			("ATT", "ATC", "ATA"): 'I', # Isoleucine (Iso/I)
			("AAA", "AAG"): 'K', # Lysine (Lys/K)
			("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"): 'L', # Leucine (Leu/L)
			"ATG": 'M', # Methionine (Met/M) START
			("AAT", "AAC"): 'N', # Asparagine (Asn/N)
			("CCT", "CCC", "CCA", "CCG"): 'P', # Proline (Pro/P)
			("CAA", "CAG"): 'Q', # Glutamine (Gln/Q)
			("CGT", "CGA", "CGG", "CGC", "AGA", "AGG"): 'R', # Arginine (Arg/R)
			("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"): 'S', # Serine (Ser/S)
			("ACT", "ACC", "ACA", "ACG"): 'T', # Threonine (Thr/T)
			("GTT", "GTC", "GTA", "GTG"): 'V', # Valine (Val/V)
			"TGG": 'W', # Tryptophan (Trp/W)
			("TAT", "TAC"): 'Y', # Tyrosine (Tyr/Y)
			("TAA", "TAG", "TGA"): '*' # Stop
		}
		
		protein = []
		for i in xrange(frame, self.get_length(), 3):
			codon = seq[i:i+3]
			if len(codon) == 3:
				try: protein.append(next(v for k, v in aa_dict.items() if codon in k))
				except StopIteration: protein.append('X')
				codon = ''
		
		return ''.join(protein), codon

def get_shared_junctions(junc_files, min_N=3):

	# Count the number of times a junctions has been found
	jN = Counter({}) 
	for f in junc_files:
		for j in yield_junctions(f):
			if j['strand'] in ['-', '+']:
				jN[(j['chr'], j['start'], j['end'], j['strand'])] += 1

	return set([ x for x in jN if jN[x] >= min_N ])

def prepare_junction_lib(work_dir, ccds_gtf, junc_files, sample_names, exitron_bed=None, min_N=3):

	# Get the shared sample junctions
	junction_sets = []
	uniq_names = set(sample_names)
	for n in natsorted(uniq_names):
		indices = [ i for i, x in enumerate(sample_names) if x == n ]
		junction_sets.append( get_shared_junctions([ junc_files[i] for i in indices ], min_N) )
	all_junctions = set.union(*junction_sets)

	# Add existing exitrons when supplied
	try: 
		if os.path.isfile(exitron_bed):
			existing_EI = set([ ( r['chr'], r['start'], r['end']-1, r['strand'] ) for r in yield_bed(exitron_bed) ])
			all_junctions = all_junctions.union(existing_EI)
	except TypeError: pass

	# Output shared junctions
	with open('{}all_junctions.N{}.bed'.format(work_dir, min_N), 'w') as fout:
		for j in natsorted(all_junctions):
			c, s, e, strand = j

			# Skip junctions with an undefined strand
			if strand in ['+', '-']: fout.write( '{}\t{}\t{}\tJUNCTION\t1000\t{}\n'.format(c, s, e, strand) )

	# Match the junctions to CCDS exons
	# bedtools intersect -s -f 1 -wa -wb -a all_junctions.bed -b ccds_gtf
	# -f 1 reports only 100% matches
	subprocess.call("{0} intersect -s -f 1 -wa -wb -a {1}all_junctions.N{2}.bed -b {3} | awk -F\"\\t\" '{{ OFS=\"\t\"; print $1, $2, $3, $6, $10, $11, $NF }}' > {1}junctions_CCDS_map.tmp".format(bedtools_path, work_dir, min_N, ccds_gtf), shell=True)

	# Parse intersection
	with open("{}junctions_CCDS_map.txt".format(work_dir), 'w') as fout:
		fout.write("#CHR\tJUNCTION_START\tJUNCTION_END\tSTRAND\tTRANSCRIPT_ID\tCCDS_START\tCCDS_END\n")
		for line in open("{}junctions_CCDS_map.tmp".format(work_dir)):
			j_chr, j_start, j_end, j_strand, ccds_start, ccds_end, attributes = line.rstrip().split('\t')
			attr = {}
			for a in attributes.split(';'):
				if len(a):
					attr_name, attr_value = list(filter(None, a.split(' ')))
					attr[attr_name.strip()] = attr_value.replace('\"', '')
			fout.write( "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(j_chr, j_start, j_end, j_strand, attr["transcript_id"], ccds_start, ccds_end) )
	subprocess.call("rm -f {}junctions_CCDS_map.tmp".format(work_dir), shell=True)

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
					attr_name, attr_value = list(filter(None, a.split(' ')))
					attr[attr_name.strip()] = attr_value.replace('\"', '')

			try: ccds[attr["transcript_id"]]["exons"].append([ int(start), int(end) ])
			except KeyError: ccds[attr["transcript_id"]] = { "chr": c, "strand": strand, "exons": [[ int(start), int(end) ]] }
	for t_id in ccds: ccds[t_id]['exons'].sort(key=lambda x: x[0])

	with open(junction_map.replace('junctions', 'exitrons'), 'w') as map_out:
		map_out.write( '#CHR\tEXITRON_START\tEXITRON_END\tSTRAND\tTRANSCRIPT_ID\tCCDS_START\tCCDS_END\tIS_ANNOTATED\tMOD3\n' )
		for line in open(junction_map):
			if not line.startswith('#'):
				c, j_start, j_end, strand, t_id, c_start, c_end = line.rstrip().split('\t')
				j_start, j_end = int(j_start), int(j_end)
				mod3 = (abs(j_end-j_start)+1) % 3

				if (c, j_start-1, j_end) in ann_introns:
					
					# Generate the spliced in bed file
					with open("{}tmp_spliced_in.bed".format(work_dir), 'w') as fout:
						for s, e in ccds[t_id]["exons"]:
							fout.write( '{}\t{}\t{}\n'.format(c, s-1, e) )

					# Generate exitrons spliced out bed file
					with open("{}tmp_spliced_out.bed".format(work_dir), 'w') as fout:
						for s, e in ccds[t_id]["exons"]:
							if s <= j_start <= e and s <= j_end <= e:
								fout.write( '{}\t{}\t{}\n'.format(c, s-1, j_start-1) )
								fout.write( '{}\t{}\t{}\n'.format(c, j_end, e) )
							else:
								fout.write( '{}\t{}\t{}\n'.format(c, s-1, e) )

					# Get sequences
					subprocess.call("{0} getfasta -fi {1} -bed {2}tmp_spliced_in.bed -fo {2}tmp_spliced_in.fa".format(bedtools_path, genome_fasta, work_dir), shell=True)
					subprocess.call("{0} getfasta -fi {1} -bed {2}tmp_spliced_out.bed -fo {2}tmp_spliced_out.fa".format(bedtools_path, genome_fasta, work_dir), shell=True)
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
						if len(inProtein) == len(outProtein) + (abs(j_end-j_start)+1)/3 + len(out_leftover): keep = True

					if keep:
						map_out.write( '{}\tyes\t{}\n'.format(line.rstrip(), mod3) )

				else: 
					map_out.write( '{}\tno\t{}\n'.format(line.rstrip(), mod3) )

	subprocess.call("rm -f {}tmp_spliced*".format(work_dir), shell=True)

def prepare_bam_files(work_dir, bam_file, handle, genome_fasta):
	""" Extract the unique reads and unique exonic reads from the 
	supplied bam file, index and output to the working directory. """

	# Extract unique reads
	uniq_reads_bam = "{}uniq_reads.{}.bam".format(work_dir, handle)
	cmd = "%s view %s | grep -w \"NH:i:1\" | perl -n -e '@line=split(/\\t/,$_); if ($line[5]=~/N/){ print \"$_\"; }' > %s" % (samtools_path, bam_file, uniq_reads_bam.replace('.bam', '.sam'))
	subprocess.call(cmd, shell=True)

	# Convert sam to bam
	#subprocess.call("{0} view -bT {1} {2} | {0} sort -o {3} -".format(samtools_path, genome_fasta, uniq_reads_bam.replace('.bam', '.sam'), uniq_reads_bam), shell=True)
	subprocess.call("{0} view -bT {1} {2} > {3}".format(samtools_path, genome_fasta, uniq_reads_bam.replace('.bam', '.sam'), uniq_reads_bam), shell=True)

	# Index bam
	subprocess.call("{} index {}".format(samtools_path, uniq_reads_bam), shell=True)

	# Remove sam
	subprocess.call("rm -f {}".format(uniq_reads_bam.replace('.bam', '.sam')), shell=True)

	# Extract unique exonic reads
	uniq_exon_reads_bam = '{}uniq_exonic_reads.{}.bam'.format(work_dir, handle)
	cmd = "%s view %s | grep -w \"NH:i:1\" | perl -n -e '@line=split(/\\t/,$_); if ($line[5]!~/N/){ print \"$_\"; }' > %s" % (samtools_path, bam_file, uniq_exon_reads_bam.replace('.bam', '.sam'))
	subprocess.call(cmd, shell=True)

	# Convert sam to bam
	#subprocess.call("{0} view -bT {1} {2} | {0} sort -o {3} -".format(samtools_path, genome_fasta, uniq_exon_reads_bam.replace('.bam', '.sam'), uniq_exon_reads_bam), shell=True)
	subprocess.call("{0} view -bT {1} {2} > {3}".format(samtools_path, genome_fasta, uniq_exon_reads_bam.replace(".bam", ".sam"), uniq_exon_reads_bam), shell=True)

	# Index bam
	subprocess.call("{} index {}".format(samtools_path, uniq_exon_reads_bam), shell=True)

	# Remove sam
	subprocess.call("rm -f {}".format(uniq_exon_reads_bam.replace('.bam', '.sam')), shell=True)

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
			info['{}:{}-{}'.format(c, s, int(e)+1)] = { 'strand': strand, 't_id': t_id, 'is_annotated': is_annotated, 'mod3': mod3 }
			#c, s, e, j_id, score, strand = line.rstrip().split('\t')
			#info['{}:{}-{}'.format(c, s, int(e)+1)] = { 'strand': strand }

			subprocess.call("{} view {} {}:{}-{} > {}tmp.sam".format(samtools_path, uniq_reads, c, s, e, work_dir), shell=True)
			N = "{}N".format((int(e)-int(s))+1) # Get the required read/junction gap signature

			uniq_count = 0
			for aln in open("{}tmp.sam".format(work_dir)):
				qname = aln.split('\t')[0]
				pos = int(aln.split('\t')[3])
				cigar = aln.split('\t')[5]
				try: start = (pos + int(re.search('^([0-9]+)M', cigar).group(1)))
				except AttributeError: pass

				# Check if the junction is at the correct position
				# and if the junction size is correct
				if N in cigar and start == int(s): uniq_count += 1

			rc['{}:{}-{}'.format(c, s, int(e)+1)] = { 'A': 0, 'B': 0, 'C': 0, 'D': uniq_count }

	# Get A, B, C read support
	all_the_single_ladies = set()
	with open("{}tmp_coverageBed_input.bed".format(work_dir), 'w') as fout:
		for line in open(exitron_map):
			if not line.startswith('#'):

				c, s, e, strand, t_id, ccds_start, ccds_end, is_annotated, mod3 = line.rstrip().split('\t')

				if (c, s, e, strand) not in all_the_single_ladies:

					all_the_single_ladies.add((c, s, e, strand))

					s, e = int(s), int(e)
					middle_point = int(s + ((e-s)/2))
					locus = '{}:{}-{}'.format(c, s, e+1)

					fout.write( '{}\t{}\t{}\t{}_A\n'.format(c, s-10, s+10, locus) )
					fout.write( '{}\t{}\t{}\t{}_B\n'.format(c, middle_point-10, middle_point+10, locus) )
					fout.write( '{}\t{}\t{}\t{}_C\n'.format(c, e-10, e+10, locus) )

	subprocess.call("sort -k1,1V -k2,2n {0}tmp_coverageBed_input.bed > {0}tmp_coverageBed_input.sorted.bed".format(work_dir), shell=True)
	#subprocess.call("sort -k1,1 -k2,2n {0}tmp_coverageBed_input.bed > {0}tmp_coverageBed_input.sorted.bed".format(work_dir), shell=True)
	subprocess.call("{0} coverage -sorted -counts -a {1}tmp_coverageBed_input.sorted.bed -b {2} > {1}tmp_coverageBed_output.bed".format(bedtools_path, work_dir, uniq_exonic), shell=True)
	for line in open("{}tmp_coverageBed_output.bed".format(work_dir)):
		c, start, end, locus, coverage = line.rstrip().split('\t')
		locus, letter = locus.split('_')
		rc[locus][letter] = int(coverage)

	# Calculate PSI
	with open("{}{}.PSI".format(work_dir, handle), 'w') as fout:
		fout.write( "EXITRON_ID\tTRANSCRIPT_ID\tSTRAND\tIS_ANNOTATED\tMOD3\tA\tB\tC\tD\tPSI\n" )
		for x in natsorted(rc):
			try: PSI = ( ( ( rc[x]['A'] + rc[x]['B'] + rc[x]['C'] ) / 3 ) / ( ( ( rc[x]['A'] + rc[x]['B'] + rc[x]['C'] ) / 3 ) + rc[x]['D'] ) ) * 100
			except ZeroDivisionError: PSI = 'NA'
			fout.write( '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(x, info[x]['t_id'], info[x]['strand'], info[x]['is_annotated'], info[x]['mod3'], rc[x]['A'], rc[x]['B'], rc[x]['C'], rc[x]['D'], PSI) )

	# Clean up
	subprocess.call( "rm -f {}tmp*".format(work_dir), shell=True )

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-w', '--work-dir', required=True, help="Output working directory.")

	subparsers = parser.add_subparsers(dest='command', help="sub-command help.")

	parser_a = subparsers.add_parser('prepare-junctions', help="Select all junction with minimal N (default = 3) occurences and map to the CCDS annotation.")
	parser_a.add_argument('-c', '--ccds', required=True, help="CCDS gtf file.")
	parser_a.add_argument('-j', '--junctions', required=True, nargs='+', help="Junction bed (SJ.out.tab) files.")
	parser_a.add_argument('-n', '--names', required=True, nargs='+', help="List of condition names, e.g. N N N D D D.")
	parser_a.add_argument('-e', '--exitrons', default=None, help="Supplemental exitron bed file.")
	parser_a.add_argument('--min-support', type=int, default=3, help="Minimum number of samples per condition in which a junction must occur.")

	parser_b = subparsers.add_parser('filter-exitrons', help="Filter out intron retention events from the potential exitrons, based on the coding potential.")
	parser_b.add_argument('-c', '--ccds', required=True, help="CCDS bed file.")
	parser_b.add_argument('-i', '--introns', required=True, help="Transcriptome introns bed file.")
	parser_b.add_argument('-j', '--junction-map', required=True, help="Junction map file (from prepare-junctions).")
	parser_b.add_argument('-g', '--genome-fasta', required=True, help="Genome fasta.")

	parser_c = subparsers.add_parser('prepare-bam', help="Extract the unique (exonic) reads from the alignment bam file.")
	parser_c.add_argument('-b', '--bam', required=True, help="Bam alignment file.")
	parser_c.add_argument('-f', '--file-handle', required=True, help="Unique file handle. The output files will be [work_dir]/uniq_reads.[file-handle].bam and [work_dir]/uniq_exonic_reads.[file-handle].bam.")
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
			sys.stderr.write("The number of junction files (-j, --junctions) must be equal to the number of sample names (-n, --names)!")

	elif args.command == "filter-exitrons":
		filter_exitrons(add_slash(args.work_dir), args.ccds, args.introns, args.junction_map, args.genome_fasta)

	elif args.command == "prepare-bam":
		prepare_bam_files(add_slash(args.work_dir), args.bam, args.file_handle, args.genome_fasta)

	elif args.command == "calculate-PSI":
		if all([ os.path.exists(args.exitron_map), os.path.exists(args.uniq_reads), os.path.exists(args.uniq_exonic) ]):
			calculate_PSI(add_slash(args.work_dir), args.exitron_map, args.uniq_reads, args.uniq_exonic, args.file_handle)
		else: 
			sys.stderr.write("One (or multiple) input file could not be accessed.\n")