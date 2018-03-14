#!/usr/bin/env python3
"""docstring"""

import re
import os
import sys

args = sys.argv[1:]

if len(args) != 2:
    print('Usage: {} BLASTOUT GENOME'.format(os.path.basename(sys.argv[0])))
    sys.exit(1)

blast_out, genome = args

print('blast_out is "{}"'.format(blast_out))

if not os.path.isfile(blast_out):
    print('"{}" is not a file'.format(blast_out))
    sys.exit(1)

flds = 'qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split()

regex = re.compile(r'^ORF\d+\((\d+)-(\d+)\)$')

L1_stop = -1;
starts = []
for line in open(blast_out):
    row = dict(zip(flds, line.rstrip().split('\t')))
    print(row)
    id = row['qaccver']
    match = regex.search(id)
    if not match:
        print("Something weird!")
        continue

    start, stop = match.groups()

    #print("id {} start {} stop {}".format(id, start, stop))

    if row['saccver'] == 'L1':
        L1_stop = stop
    else:
        starts.append(int(start))

starts.sort() # inplace! mutate! danger! (BAD)
print("L1 stop {} next start {}".format(L1_stop, starts[0]))
