#!/usr/bin/env python

"""
Generate exons and introns bed files from a GTF file.
"""

from natsort import natsorted
from fileParser import parse_GTF
import argparse
import re
import os

def gtf2exon_intron(gtf_file, output_dir, select_feature="exon", t_id_attr="transcript_id", attr_sep=' "', base="", for_FA=False):

	# Add dot after base name
	if len(base) and not base.endswith('.'): base += '.'
	if not os.path.exists(output_dir): os.makedirs(output_dir)

	# Parse GTF
	gtf = parse_GTF(gtf_file, select_feature, t_id_attr, attr_sep, True)

	# Output introns
	fintron = open( '{}{}introns.bed'.format(output_dir, base), 'w' )

	for t in natsorted(gtf):

		for i, intron in enumerate(gtf[t]['introns']):
			s, e = intron
			size = abs(s-e)+1
			if for_FA: fintron.write( '{0}\t{1}\t{2}\t{3}|intron-{4} | x-x | {0}:{5}-{2} {6} LENGTH={7}\t1000\t{8}\n'.format(gtf[t]['chr'], s, e, t, i+1, s, "FORWARD" if gtf[t]['strand'] == '+' else "REVERSE", size, gtf[t]['strand']) )
			else: fintron.write( '{0}\t{1}\t{2}\t{3}|intron-{4} | x-x | {0}:{5}-{2} {6} LENGTH={7}\t1000\t{8}\n'.format(gtf[t]['chr'], s, e, t, i+1, s, "FORWARD" if gtf[t]['strand'] == '+' else "REVERSE", size, gtf[t]['strand']) )

	fintron.close()

	# Output exons
	fexon = open( '{}{}exons.bed'.format(output_dir, base), 'w' )
	for t in natsorted(gtf):

		strand = gtf[t]["strand"]
		c = gtf[t]["chr"]
		isFwd = strand == '+'
		exons = [ [s, e] for s, e in natsorted(gtf[t]['exons']) ] if isFwd else [ [s, e] for s, e in natsorted(gtf[t]["exons"], reverse=True) ]

		for i, exon in enumerate(exons):

			s, e = exon
			size = abs(s-e)+1
			if for_FA: fexon.write( '{0}\t{1}\t{2}\t{3}|exon-{4} | x-x | {0}:{5}-{2} {6} LENGTH={7}\t1000\t{8}\n'.format(c, s-1, e, t, i+1, s, "FORWARD" if isFwd else "REVERSE", size, strand) )
			else: fexon.write( '{0}\t{1}\t{2}\t{3}|exon-{4} | x-x | {0}:{5}-{2} {6} LENGTH={7}\t1000\t{8}\n'.format(c, s, e, t, i+1, s, "FORWARD" if isFwd else "REVERSE", size, strand) )

	fexon.close()

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g', '--gtf', default='-', help="Input GTF file")
	parser.add_argument('-o', '--output-dir', help="Output directory", required=True)
	parser.add_argument('--select-feature', default="exon", help="Feature to select in the GTF file")
	parser.add_argument('--transcript-attr', default='transcript_id', help="The transcript identifier attribute name")
	parser.add_argument('--attr-sep', default=' "', help="GFF attribute and attribute value separator")
	parser.add_argument('--base', help="Additional file description. Files will be BASE.exons.bed and BASE.introns.bed.", default="")
	parser.add_argument('--for-FASTA', dest='for_fasta', action='store_true')
	parser.set_defaults(for_fasta=False)
	args = parser.parse_args()

	gtf2exon_intron(args.gtf, args.output_dir+'/', args.select_feature, args.transcript_attr, args.attr_sep, args.base, args.for_fasta)
