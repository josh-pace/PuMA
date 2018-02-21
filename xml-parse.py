#!/usr/bin/env python3

from Bio.Blast import NCBIXML

fh = open('all-blast-hits.xml')
for i, rec in enumerate(NCBIXML.parse(fh)):
    print('>>>>>>{}\n'.format(i, rec))
    print('QUERY "{}"\n'.format(rec.query))
    for hit in rec.alignments:
        print('>>> HIT');
        print('title ' + hit.title)
        print('hit_id ' + hit.hit_id)
        print('hit_def ' + hit.hit_def)
        for hsp in hit.hsps:
            print('>>> HSP')
            print(hsp)
    break

print('Done')
