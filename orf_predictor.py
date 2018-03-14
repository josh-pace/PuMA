#!/usr/bin/env python3
"""docstring"""

import argparse
import os
import re
import sys
from Bio import SeqIO

# --------------------------------------------------
def get_args():
    """get args"""
    parser = argparse.ArgumentParser(description='Argparse Python script',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', metavar='FILE', help='Input sequence file')
    parser.add_argument('-f', '--format', help='Input file format',
                        metavar='STR', type=str, default='fasta')
    parser.add_argument('-o', '--outfile', help='Output file',
                        metavar='STR', type=str, default='orfs.fa')
    parser.add_argument('-m', '--min', help='Minimum length',
                        metavar='INT', type=int, default=50)
    parser.add_argument('-t', '--trans', help='Translation table',
                        metavar='INT', type=int, default=1)
    return parser.parse_args()

# --------------------------------------------------
def find_orfs(seq, trans_table, min_prot_len):
    """Find all ORFs in the sequence"""
    orfs_re = re.compile(r'(M[^*]+\*)?')
    orfs = []

    for frame in range(3):
        # ensure we get a segment evenly divisible by 3
        frame_end = len(seq)
        while (frame_end - frame) % 3 != 0:
            frame_end -= 1

        trans = str(seq[frame:frame_end].translate(trans_table))

        for match in orfs_re.finditer(trans):
            aa_seq = match.group(0)
            if len(aa_seq) >= min_prot_len:
                orfs.append({'seq': aa_seq,
                             'start': frame + (match.start() * 3),
                             'stop': frame + (match.end() * 3) - 1})

    return orfs

# --------------------------------------------------
def main():
    """main"""
    args = get_args()
    infile = args.file
    file_format = args.format
    min_len = args.min
    trans_table = args.trans
    out_file = args.outfile

    if not os.path.isfile(infile):
        print('"{}" is not a file'.format(infile))
        sys.exit(1)

    out_dir = os.path.abspath(os.path.dirname(out_file))
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    fh = open(out_file, 'wt')

    for rec in SeqIO.parse("{}".format(infile), file_format):
        orfs = sorted(find_orfs(rec.seq, trans_table, min_len),
                      key=lambda orf: orf['start'])
        print('{} has {}'.format(rec.id, len(orfs)))
        for j, orf in enumerate(orfs):
            fh.write('>ORF{}({}-{})\n{}\n'.format(j + 1,
                                                  orf['start'] + 1,
                                                  orf['stop'] + 1,
                                                  orf['seq']))

    print('Done, see output file "{}"'.format(out_file))

# --------------------------------------------------
if __name__ == '__main__':
    main()
