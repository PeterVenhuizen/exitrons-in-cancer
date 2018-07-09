#!/usr/bin/env python

import argparse
import subprocess
from fileParser import parse_GTF
from natsort import natsorted

def exitrons_from_GTF(full_gtf, ccds_gtf):

	# Extract all introns from the CCDS transcripts
	# and output to a bed file
	gtf = parse_GTF(full_gtf, select_feature="exon", get_introns=True)

	with open('tmp.introns.bed', 'w') as fout:
		for t in natsorted(gtf):

			for i, intron in enumerate(gtf[t]['introns']):
				s, e = intron
				size = abs(s-e)+1
				fout.write( '{0}\t{1}\t{2}\t{3}|intron-{4} | x-x | {0}:{5}-{2} {6} LENGTH={7}\t1000\t{8}\n'.format(gtf[t]['chr'], s+1, e, t, i+1, s, "FORWARD" if gtf[t]['strand'] == '+' else "REVERSE", size, gtf[t]['strand']) )

	# intersect with the GTF
	subprocess.call("/usr/local/bin/bedtools intersect -s -f 1 -wa -wb -a tmp.introns.bed -b {} | awk -F\"\\t\" '{{ OFS=\"\t\"; print $1, $2, $3, $6, $10, $11, $NF }}' > junctions_CCDS_map.tmp".format(ccds_gtf), shell=True)

	# Parse intersection
	with open("junctions_CCDS_map.txt", 'w') as fout:
		fout.write("#CHR\tJUNCTION_START\tJUNCTION_END\tSTRAND\tTRANSCRIPT_ID\tCCDS_START\tCCDS_END\n")
		for line in open("junctions_CCDS_map.tmp"):
			j_chr, j_start, j_end, j_strand, ccds_start, ccds_end, attributes = line.rstrip().split('\t')
			attr = {}
			for a in attributes.split(';'):
				if len(a):
					attr_name, attr_value = list(filter(None, a.split(' ')))
					attr[attr_name.strip()] = attr_value.replace('\"', '')
			fout.write( "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(j_chr, j_start, j_end, j_strand, attr["transcript_id"], ccds_start, ccds_end) )
	subprocess.call("rm -f junctions_CCDS_map.tmp", shell=True)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g', '--gtf', required=True, help="GTF file.")
	parser.add_argument('--ccds', required=True, help="CCDS GTF file.")
	args = parser.parse_args()

	exitrons_from_GTF(args.gtf, args.ccds)