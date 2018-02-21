#!/usr/bin/env python3
"""Run PUMA"""

import re
import sys
import os
import argparse
import subprocess
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline as blastp

# --------------------------------------------------
def get_args():
    """get args"""
    args = sys.argv
    bin_path = os.path.abspath(os.path.dirname(args[0]))
    parser = argparse.ArgumentParser(description='Find HPV proteins',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', metavar='FILE',
                        help='HPV genome file', required=True)

    parser.add_argument('-f', '--format', metavar='FORMAT', default='fasta',
                        help='FASTA/GenBank')

    parser.add_argument('-b', '--blastdb_dir', metavar='DIR', type=str,
                        default=os.path.join(bin_path, 'blast_database'),
                        help='BLAST db directory')

    parser.add_argument('-e', '--evalue', metavar='FLOAT', type=float,
                        default=0.001, help='BLAST evalue')

    parser.add_argument('-m', '--min_prot_len', metavar='NUM', type=int,
                        default=25,
                        help='Minimum protein length')

    parser.add_argument('-o', '--outdir', metavar='DIR', type=str,
                        default='puma-out',
                        help='Output directory')

    sites = 'L1 L2 E1 E2 E4 E5 E6 E7 E10 E2BS E1BS URR ALL'
    parser.add_argument('-s', '--sites', metavar='STR', type=str, default='ALL',
                        help='Comma-separated string of ' + sites)

    return  parser.parse_known_args()

# --------------------------------------------------
def get_orfs(seq, trans_table, min_prot_len):
    """Find all ORFs in the sequence"""
    orfs = []
    has_m = re.compile('M')

    # Frame shift
    for frame in range(3):
        # ensure we get a segment evenly divisible by 3
        frame_end = len(seq)
        while (frame_end - frame) % 3 != 0:
            frame_end -= 1

        trans = str(seq[frame:frame_end].translate(trans_table))
        trans_len = len(trans)
        aa_start = 0
        aa_end = 0

        while aa_start < trans_len:
            # returns the lowest index where "*" is found in
            # aa_start because "*" represents a stop codon
            aa_end = trans.find("*", aa_start)

            # .find returns -1 if no "*" is found in aa_start
            if aa_end == -1:
                aa_end = trans_len

            # Only finding proteins over a certain length
            if aa_end - aa_start >= min_prot_len:
                # Multiplying by 3 to get start position of the nucleotides
                start = frame + aa_start * 3
                aa_seq = trans[aa_start:aa_end]

                if has_m.search(aa_seq):
                    orfs.append([aa_seq, start])

            aa_start = aa_end + 1

    return orfs

# --------------------------------------------------
def main():
    """main"""
    args, _ = get_args()
    sites = re.split(r'\s*,\s*', args.sites.upper())
    input_file = args.input
    out_dir = os.path.abspath(args.outdir)
    blast_dir = os.path.abspath(args.blastdb_dir)
    input_format = args.format.lower()
    min_prot_len = args.min_prot_len
    evalue = args.evalue

    valid_sites = set('L1 L2 E1 E2 E4 E5 E6 E7 E10 E2BS E1BS URR ALL'.split())
    if not sites:
        print('--sites is required')
        sys.exit(1)

    if not input_file:
        print("--input is required")
        sys.exit(1)

    if not os.path.isfile(input_file):
        print('--input "{}" is not a file.'.format(input_file))
        sys.exit(1)

    if not os.path.isdir(blast_dir):
        print('--blastdb_dir "{}" is not a directory.'.format(blast_dir))
        sys.exit(1)

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    bad_sites = list(filter(lambda s: s not in valid_sites, sites))

    if bad_sites:
        print('Invalid site(s): {}'.format(', '.join(bad_sites)))
        sys.exit(1)

    valid_format = set(['fasta', 'genbank'])
    if not input_format in valid_format:
        msg = 'Invalid format ({}), please choose from {}'
        print(msg.format(input_format, ', '.join(valid_format)))
        sys.exit(1)

    recs = list(SeqIO.parse(input_file, input_format))
    num_seqs = len(recs)
    if num_seqs != 1:
        print('Expected 1 sequence in "{}," got "{}"'.format(input_file,
                                                             num_seqs))
        sys.exit(1)

    genome = recs[0]
    orfs = get_orfs(genome.seq, 1, min_prot_len)

    if not orfs:
        print('No ORFs, must stop.')
        sys.exit(1)

    orfs_fa = os.path.join(out_dir, 'orfs.fa')
    orfs_fh = open(orfs_fa, 'wt')
    for i, orf in enumerate(orfs):
        orfs_fh.write('\n'.join(['>' + str(i + 1), orf[0], '']))
    orfs_fh.close()

    print('Wrote {} ORFs to "{}"'.format(len(orfs), orfs_fa))

    blast_db = os.path.join(blast_dir, 'blast_database.txt')
    blast_out = os.path.join(out_dir, 'blast_results.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)

    print('BLASTing')
    cmd = blastp(query=orfs_fa,
                 db=blast_db,
                 evalue=evalue,
                 outfmt=6,
                 out=blast_out)

    stdout, stderr = cmd()
    if stdout:
        print("STDOUT = ", stdout)
    if stderr:
        print("STDERR = ", stderr)

    if not os.path.isfile(blast_out) or not os.path.getsize(blast_out):
        print('No BLAST output "{}" (must have failed)'.format(blast_out))
        sys.exit(1)

    print('See BLAST out "{}"'.format(blast_out))

    # So at this point I have a tab-delimited file of all the BLAST hits
    # and I'm trying to figure out the L1 business and URR

    print('Done.')

# --------------------------------------------------
if __name__ == '__main__':
    main()
