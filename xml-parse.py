#!/usr/bin/env python3

from Bio.Blast import NCBIXML

fh = open('all-blast-hits.xml')
for i, rec in enumerate(NCBIXML.parse(fh)):
    print('>>>>>>{}\n{}\n'.format(i, rec.description))
    for hit in rec.alignments:
        print(hit)
        for hsp in hit.hsps:
            print('>>>HSP')
            print(hsp)
