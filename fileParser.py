#!/usr/bin/env python

"""
Collection of often used file parsers
"""

from natsort import natsorted
import sys
import gzip

def yield_salmon(f):
	""" Yield salmon output per line """

	for line in open(f):
		try: 
			t_id, length, eff_length, tpm, num_reads = line.rstrip().split('\t')
			yield( { 'transcript_id': t_id,  \
				'length': int(length), \
				'eff_length': float(eff_length), \
				'tpm': float(tpm), \
				'num_reads': float(num_reads) } )
		except ValueError: pass

def parse_salmon(f, group="gene"):
	""" Return the TPM values and group per gene or transcript """

	if group not in ["gene", "transcript"]: group = "gene"
	
	salmon = {}
	for record in yield_salmon(f):
		t_id = record["transcript_id"]
		if group == "gene":
			g_id = t_id.split('.')[0]
			try: salmon[g_id][t_id] = float(record["tpm"])
			except KeyError: salmon[g_id] = { t_id: float(record["tpm"]) }
		else: 
			salmon[t_id] = float(record["tpm"])

	return salmon

def yield_kallisto(f):
	""" Yield kallisto output per line """

	for line in open(f):
		try: 
			t_id, length, eff_length, est_counts, tpm = line.rstrip().split('\t')
			yield( { 'transcript_id': t_id,  \
				'length': int(length), \
				'eff_length': float(eff_length), \
				'tpm': float(tpm), \
				'est_counts': float(num_reads) } )
		except ValueError: pass

	pass

def parse_kallisto(f, group="gene"):
	""" Return the TPM values and group per gene or transcript """

	if group not in ["gene", "transcript"]: group = "gene"
	
	kallisto = {}
	for record in yield_kallisto(f):
		t_id = record["transcript_id"]
		if group == "gene":
			g_id = t_id.split('.')[0]
			try: kallisto[g_id][t_id] = float(record["tpm"])
			except KeyError: kallisto[g_id] = { t_id: float(record["tpm"]) }
		else: 
			kallisto[t_id] = float(record["tpm"])

	return kallisto

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

class parseFastQ(object):
	""" Returns read-by-read fastQ records as tuple: (sequence identifier, raw sequence, quality identifier, quality sequence). """

	def __init__(self, fastq_file):
		if fastq_file.endswith('gz'): self.fastq_file = gzip.open(fastq_file)
		else: self.fastq_file = open(fastq_file)
		self.lineNumber = 0
	
	def __iter__(self):
		return self
		
	def next(self):
		
		# A fastQ record should consist of 4 lines
		record = []
		for i in xrange(4):
			line = self.fastq_file.readline()
			self.lineNumber += 1
			
			if line: record.append(line.rstrip())
			else: record.append(None)
			
		# Inspect the record
		nValues = [ bool(x) for x in record ].count(True)
		
		# End of file
		if record.count(None) == 4:
			raise StopIteration
		
		# Check for premature EOF or empty line
		assert nValues == 4, \
			"ERROR: Premature EOF or empty line found near line number %d" % (self.lineNumber)
		
		# Check for "@" header
		assert record[0].startswith('@'), \
			"ERROR: The 1st line of the record does not start with '@'. Please check near %d" % (self.lineNumber-3)

		# Check for "+" header
		assert record[2].startswith('+'), \
			"Error: The 3rd line of the record does not start with '+'. Please check near %d" % (self.lineNumber-1)
			
		# Compare sequence and quality length
		#assert len(record[1]) == len(record[3]), \
		#	"Error: The length of the sequence data and the quality data are not equal. Please check near %d" % (self.lineNumber-2)
			
		return tuple(record)

def parse_GTF(f, select_feature="exon", t_id_attr="transcript_id", attr_sep=' "', get_introns=True):
	""" Parse a GTF file. Select transcripts based on the t_id_attr. """

	gtf = {}
	for line in open(f):
		if not line.startswith('#'):
			c, source, feature, start, end, score, strand, frame, attributes = line.rstrip().split('\t')
			attr = {}
			for a in attributes.split(';'):
				if len(a):
					attr_name, attr_value = list(filter(None, a.split(' ')))
					attr[attr_name.strip()] = attr_value.replace('\"', '')

			if feature == select_feature:
				try: gtf[attr[t_id_attr]]["exons"].append([ int(start), int(end) ])
				except KeyError: gtf[attr[t_id_attr]] = { "chr": c, "strand": strand, "exons": [[ int(start), int(end) ]], "introns": [] }

	# Get the transcript intron coordinates
	if get_introns:

		import re
		intron_RE = re.compile('I+')

		for t_id in natsorted(gtf):

			c = gtf[t_id]['chr']
			isFwd = gtf[t_id]['strand'] == '+'
			gtf[t_id]['exons'] = [ [s, e] for s, e in natsorted(gtf[t_id]['exons']) ] if isFwd else [ [s, e] for s, e in natsorted(gtf[t_id]['exons'], reverse=True) ]

			# Get transcript start and end
			start = min([ gtf[t_id]['exons'][i][0] for i in xrange(len(gtf[t_id]['exons'])) ])
			end = max([ gtf[t_id]['exons'][i][1] for i in xrange(len(gtf[t_id]['exons'])) ])

			IE_list = list('I' * abs(end-start))
			for s, e in gtf[t_id]['exons']: IE_list[ s-start:e-start ] = list('E' * (e-s))
			if isFwd:
				for i, m in enumerate(intron_RE.finditer(''.join(IE_list))):
					gtf[t_id]['introns'].append([ m.start()+start, m.start()+start+len(m.group())-1 ])
			else:
				for i, m in enumerate(reversed(list(intron_RE.finditer(''.join(IE_list))))):
					gtf[t_id]['introns'].append([ m.start()+start, m.start()+start+len(m.group())-1 ])

	return gtf

def yield_junctions(f, source="TopHat2"):
	""" Yield each line as a dictionary. """

	if source == "TopHat2":
		for line in open(f):
			try: 
				cols = line.rstrip().split('\t')
				yield({
					'chr': cols[0],
					'start': cols[1],
					'end': cols[2],
					'junction_id': cols[3],
					'depth': int(cols[4]),
					'strand': cols[5],
					'rgb': cols[8],
					'block_size': cols[9],
					'blocks': cols[10],
					'junc_start': int(cols[1]) + int(cols[10].split(',')[0]),
					'junc_end': int(cols[2]) - int(cols[10].split(',')[1]),
					'misc': cols[11]
					})
			except IndexError: pass

	elif source == "HISAT2":
		for line in open(f):
			c, s, e, strand = line.rstrip().split('\t')
			yield({
				'chr': c,
				'start': int(s),
				'end': int(e),
				'strand': strand
				})

	elif source == "STAR":
		for line in open(f):
			cols = line.rstrip().split('\t')
			yield({
				'chr': cols[0],
				'start': cols[1],
				'end': cols[2],
				'strand': { '0': 'undefined', '1': '+', '2': '-' }[cols[3]],
				'intron_motif': { '0': 'non-canonical', '1': 'GT-AG', '2': 'CT-AC', '3': 'GC-AG', '4': 'CT-GC', '5': 'AT-AC', '6': 'GT-AT' }[cols[4]],
				'is_annotated': cols[5] == '1',
				'uniq_reads': int(cols[6]),
				'multimap_reads': int(cols[7]),
				'max_overhang': int(cols[8])
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

def parse_vcf(f):

	vcf = {}
	for line in open(f):
		if not line.startswith('#'):
			CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = line.rstrip().split('\t')[:8]
			vcf['{}:{}:{}-{}'.format(CHROM, POS, REF, ALT)] = { 'ID': ID, 'REF': REF, 'ALT': ALT, 'QUAL': float(QUAL), 'FILTER': FILTER, 'INFO': INFO }
	return vcf