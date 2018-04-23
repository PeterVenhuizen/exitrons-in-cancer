#!/usr/bin/env python

from natsort import natsorted

tx2gene = {}
for line in open("gencode.v27.annotation.gtf"):
    if not line.startswith('#'):
        feature = line.split('\t')[2]
        attributes = line.rstrip().split('\t')[-1]
        attr = {}
        for a in attributes.split(';'):
            if len(a):
                attr_name, attr_value = list(filter(None, a.split(' ')))
                attr[attr_name.strip()] = attr_value.replace('\"', '') 

        if feature == "transcript":
            try: tx2gene[attr["gene_id"]].add(attr["transcript_id"])
            except KeyError: tx2gene[attr["gene_id"]] = set([attr["transcript_id"]])

with open('tx2gene.csv', 'w') as fout: 
    for gene_id in natsorted(tx2gene): 
        for tx in natsorted(tx2gene[gene_id]):
            fout.write( '{}\t{}\n'.format(tx, gene_id) )