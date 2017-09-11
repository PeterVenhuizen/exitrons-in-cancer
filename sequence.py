from __future__ import division

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
	
	def get_base_content(self, chars):
		''' Return the percentage of bases in the sequence. 
		If the chars variable is set to 'gc', than the 
		GC-content of the sequence is returned. '''
		
		return sum([self.sequence.count(char) for char in chars.upper()])/self.get_length()
    
	def get_base_frequency(self, to_percent=False):
		''' Return the base composition counts for the sequence.
		A dictionary is created containing the counts for each
		unique character in the sequence. '''
		
		if to_percent: return { b: self.sequence.count(b)/self.get_length() for b in set(self.sequence) }
		else: return { b: self.sequence.count(b) for b in set(self.sequence) }
	
	def get_kmer_frequency(self, kmer_size=3):
		''' Return the kmer frequencies for the sequence. A counter
		is created containing the counts for each unique kmer. '''
		
		from collections import Counter
		
		d = Counter({})
		seq_len = len(self.sequence)-kmer_size+1
		for i in xrange(0, seq_len):
			kmer = self.sequence[i:i+kmer_size]
			try: d[kmer] += 1
			except KeyError: d[kmer] = 1
		return d

	def get_ORF(self):
		''' Look for ORFs in the sequence. '''
		pass
	
	def explode_IUPAC(self):
		''' Get all the unique sequences from an IUPAC sequence. '''
		
		IUPAC_translation = { 'A': 'A',
			'B': 'CGT',
			'C': 'C',
			'D': 'AGT',
			'G': 'G',
			'H': 'ACT',
			'I': 'ACGT',
			'K': 'GT',
			'M': 'AC',
			'N': 'ACGT',
			'R': 'AG',
			'S': 'CG',
			'T': 'T',
			'U': 'U',
			'V': 'ACG',
			'W': 'AT',
			'Y': 'CT'
		}
	
		seqArray = []
	
		for base in self.sequence:
			bases = IUPAC_translation[base]
			nBases = len(bases)
		
			if nBases > 1:
				if len(seqArray) == 0: seqArray = [ b for b in bases ]
				else: seqArray = [ s+b for s in seqArray for b in bases ]
			elif len(seqArray) > 0: seqArray = [ s+bases for s in seqArray ]
			else: seqArray = [ ''.join(seqArray)+bases ]
	
		return seqArray
	
	def iupac2regex(self):
		''' Convert a IUPAC sequence (string) to a python regular expression. '''

		code = { 'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 
			'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]',
			'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]' }
		return ''.join([ code[b.upper()] for b in list(self.sequence) ])

	def get_PCR_product(self, fwd, rev, regex=False):
		''' Get the PCR product(s) based on a forward and reverse primer. '''
		
		# Translate IUPAC to regex
		fwd_regex = self.iupac2regex(fwd)
		rev_regex = self.iupac2regex(rev)
		
		# Get the pcr products
		import re
		p = re.compile('(?=(%s.*%s))' % (fwd_regex, rev_regex))
		for i, m in enumerate(p.finditer(self.sequence)):
			print '>%d) start=%d; end=%d; size=%d\n%s' % (m.start()+1, m.start()+len(m.group(1)), len(m.group(1)), m.group(1))	