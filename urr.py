#!/usr/bin/env python3
"""docstring"""

import re
import os
import sys
from Bio import SeqIO

# -----------------------------------------------------------------------------------------
def get_args():
    args = sys.argv[1:]

    if len(args) != 2:
        print('Usage: {} BLASTOUT GENOME'.format(os.path.basename(sys.argv[0])))
        sys.exit(1)
    #print(args)
    blast_out, virus = args

    return
# -----------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------
def main():

    for seq_record in SeqIO.parse("{}".format(virus), 'fasta'):
        genome = seq_record.seq
    # Opening file and storing the genome

    #print('blast_out is "{}"'.format(blast_out))

    if not os.path.isfile(blast_out):
        print('"{}" is not a file'.format(blast_out))
        sys.exit(1)

    flds = 'qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split()

    regex = re.compile(r'^ORF\d+\((\d+)-(\d+)\)$')

    L1_stop = -1;
    starts = []
    for line in open(blast_out):
        row = dict(zip(flds, line.rstrip().split('\t')))
        #print(row)
        id = row['qaccver']
        match = regex.search(id)
        if not match:
            #print("Something weird!")
            continue

        start, stop = match.groups()

        #print("id {} start {} stop {}".format(id, start, stop))

        if row['saccver'] == 'L1':
            L1_stop = stop
        else:
            starts.append(int(start))

    starts.sort() # inplace! mutate! danger! (BAD)
    #print("L1 stop {} next start {}".format(L1_stop, starts[0]))

    URRstart = int(L1_stop) + 1
    URRstop = int(starts[0]) - 1

    genomelen = int(len(genome))

    if URRstop == 0:
        URRstop = genomelen

    if URRstop > URRstart:
        URRfound = str(genome[URRstart-1:URRstop]).lower()
        #Finding the URR if it stops at or before the end of the genome
        print(">{}-{}\n{}".format(URRstart,URRstop,URRfound))

    else:
        URRfound = str(genome[URRstart-1:] + genome[:URRstop]).lower()
        #Finding the URR if it goes past the end of the genome
        print(">{}-{}-1-{}\n{}".format(URRstart,genomelen, URRstop, URRfound))

# -----------------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
