#!/usr/bin/env python3
"""Run PUMA"""

from puma_functions import *

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
def Trans_ORF(seq, trans_table, min_protein_length):
    ORFs = {}
    has_m = re.compile('M')
    for frame in range(3):  # Accounting for all 3 different frames
        trans = str(seq[frame:].translate(trans_table))
        trans_len = len(trans)
        aa_start = 0
        aa_end = 0
        while aa_start < trans_len:
            aa_end = trans.find("*",
                                aa_start)  # returns the lowest index where "*" is found in aa_start because "*" represents a stop codon
            if aa_end == -1:  # .find returns -1 if no "*" is found in aa_start
                aa_end = trans_len
            if aa_end - aa_start >= min_protein_length:  # Only finding proteins over a certain length
                start = frame + aa_start * 3  # Multiplying by 3 to get the start position of the nucleotides
                aa_seq = trans[aa_start:aa_end]
                if has_m.search(aa_seq):
                    ORFs[aa_seq] = start
            aa_start = aa_end + 1

    return ORFs
# --------------------------------------------------
def blast_proteins(genome,min_prot_len,evalue,blast_dir, out_dir):
    protein_start = {}
    protein_seq = {}
    found_proteins = {}

    orfs_dict = Trans_ORF(genome.seq, 1, min_prot_len)
    if not orfs_dict:
        print('No ORFs, must stop.')
        sys.exit(1)

    orfs_fa_dict = os.path.join(out_dir, 'orfs_dict.fa')
    orfs_fh_dict = open(orfs_fa_dict, 'wt')

    for orf in orfs_dict:
        orfs_fh_dict.write('\n'.join(['>' + str(orfs_dict[orf]), orf, '']))
    orfs_fh_dict.close()

    print('Wrote {} ORFs to "{}"'.format(len(orfs_dict), orfs_fa_dict))

    blast_db = os.path.join(blast_dir, 'blast_database.txt')
    blast_out = os.path.join(out_dir, 'blast_results.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)

    print('BLASTing')
    cmd = blastp(query=orfs_fa_dict,
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

    with open(blast_out) as tab_file:
        for line in csv.reader(tab_file, delimiter="\t"):
            protein_start[line[1]] = int(line[0])

    with open(orfs_fa_dict) as fasta_file:
        for line in fasta_file:
            for num in protein_start:
                try:
                    start = int(line[1:])
                    if start == protein_start[num]:
                        seq = next(fasta_file)
                        protein_seq[num] = seq[:-1]
                        # print("YES")
                except:
                    pass

    for seq in protein_seq:
        for start in protein_start:
            if seq == start:
                M = re.search('M',protein_seq[seq])
                real_start = protein_start[start] + M.start() + M.start() + M.start()
                end = protein_start[start] + ((len(protein_seq[seq])+1) * 3)
                found_proteins[seq] = [int(real_start) + 1, int(end)]




    return found_proteins

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
    # orfs = get_orfs(genome.seq, 1, min_prot_len)
    # orfs_dict = Trans_ORF(genome.seq, 1, min_prot_len)
    # #print(orfs)
    # if not orfs:
    #     print('No ORFs, must stop.')
    #     sys.exit(1)
    #
    # orfs_fa = os.path.join(out_dir, 'orfs.fa')
    # orfs_fh = open(orfs_fa, 'wt')
    # for i, orf in enumerate(orfs):
    #     orfs_fh.write('\n'.join(['>' + str(i + 1), orf[0], '']))
    # orfs_fh.close()
    #
    # orfs_fa_dict = os.path.join(out_dir, 'orfs_dict.fa')
    # orfs_fh_dict = open(orfs_fa_dict, 'wt')
    # for orf in orfs_dict:
    #     orfs_fh_dict.write('\n'.join(['>' + str(orfs_dict[orf]),orf, '']))
    # orfs_fh_dict.close()
    #
    # print('Wrote {} ORFs to "{}"'.format(len(orfs), orfs_fa))
    # print('Wrote {} ORFs to "{}"'.format(len(orfs_dict), orfs_fa_dict))
    #
    # blast_db = os.path.join(blast_dir, 'blast_database.txt')
    # blast_out = os.path.join(out_dir, 'blast_results.tab')
    #
    # if os.path.isfile(blast_out):
    #     os.remove(blast_out)
    #
    # print('BLASTing')
    # cmd = blastp(query=orfs_fa_dict,
    #              db=blast_db,
    #              evalue=evalue,
    #              outfmt=6,
    #              out=blast_out)
    #
    # stdout, stderr = cmd()
    # if stdout:
    #     print("STDOUT = ", stdout)
    # if stderr:
    #     print("STDERR = ", stderr)
    #
    # if not os.path.isfile(blast_out) or not os.path.getsize(blast_out):
    #     print('No BLAST output "{}" (must have failed)'.format(blast_out))
    #     sys.exit(1)
    #
    # print('See BLAST out "{}"'.format(blast_out))
    #
    # # So at this point I have a tab-delimited file of all the BLAST hits
    # # and I'm trying to figure out the L1 business and URR
    #
    # print('Done.')

    found = blast_proteins(genome,min_prot_len,evalue,blast_dir,out_dir)
    print("Testing blast_proteins")
    print("E6 seq:{}".format(str(genome.seq[found['E6'][0]-1:found['E6'][1]]).lower()))
    found_E6 = str(genome.seq[found['E6'][0]-1:found['E6'][1]]).lower()
    found_L1 = str(genome.seq[found['L1'][0] - 1:found['L1'][1]]).lower()
    real_E6 = """atggacctgcaaagtttttccagaggcaatcctttctcaggattggcctgtgtttggtgcagggagcctctcacagaagttgatgcttttaggtgcatgataaaagactttcatgttgtataccgagatggtgtgaaatttggtgcatgtaccacttgtcttgagaactgcttagataaagaaagaagactgtggaaaggtgtgccagtaacaggtgaggaagctcaattattgcatggcaaatcccttgataggctttgcataagatgctgctactgtgggggaaaactaaccaaaaacgagaagcagcggcatgtgctttataatgagcctttttgcaaaacgagatctaacataataagaggacgctgctacgactgctgcagacatggttcaaggtccaactacccatag"""
    real_L1 = '''atggcgttgtggcaacaaggccaaaagctgtatctccctccaacccctgtaagcaaggtgctatgcagtgaaacctatgtgcaaagaaaaagcatattctatcatgcagaaacggaacgcctgttaactgtaggacatccatactaccaagtcactgtgggggacaaaactgttcccaaagtgtctgctaatcaatttagagtttttaaaatacagctccccgatcccaatcagtttgcattgcctgataggactgtgcacaatccaagcaaggagcgcctggtttgggctgtaataggggttcaagtatctcgtggccaaccactaggaggcacagttactgggcaccccacttttaatgctctgcttgatgcagaaaatgttaatagaaaagttactgcacaaacaacagatgacaggaagcaaacaggattagatgctaagcaacaacagattctgttgctgggctgtacccctgcagaaggggaatactggaccacagcccgtccatgtgttactgatagactagaaaatggtgcgtgtcctcctttagaattaaagaacaaacacatagaagatggagacatgatggaaatagggtttggtgctgctgactttaaaacactaaatgccagtaaatcagatctacctcttgacattcaaaatgaaatatgcctgtatccagactacctcaaaatggctgaagatgctgctggaaacagtatgttcttctttgcaagaaaagaacaagtgtatgtaaggcatatatggactcgggggggctctgaaaaagaagcacccagtaaagacttctacctcaaaaatggtagaggtgaagaaactctaaaaatacctagtgtgcactttggcagtcccagtggatccttggtgtccactgataatcaaatatttaacaggccttattggctattcagggctcagggcatgaacaatgggattgcatggaataatttattatttttaactgtaggggataacacacggggaactaaccttagtattagtgtagctgcagatggaaacgcattgtcagagtatgatactggcaaatttaacctataccataggcatatggaagaatataagctagcatttatattggagctgtgctctgttgagattactgcacaaacactgtcacatctgcaaggactgatgccctctgtgctacaaaactgggaaatcggggtgcaacctcctgcttcttctattttagaagatacttataggtacatagagtctcctgcaactaaatgtgcaagtaatgttataccacccaaagaagacccttatgcagggcttaagttttggagcatagacttaaaagaaaagctgtctttggacttagaccaatttcccttgggaagaagattcttagctcagcaaggggcaggatgttcaactgtgagaaagagagctgttgcaaccagaaattccagtaagcctgcaaaaagaaaaaaaatcaaagcttaa'''
    print("Real:{}".format(real_E6))
    if found_E6 ==real_E6:
        print("YES")
    else:
        print("NOOOO")
    if found_L1 ==real_L1:
        print("YES")
    else:
        print("NOOOO")
    print("found_proteins:{}".format(found))

# --------------------------------------------------
if __name__ == '__main__':
    main()
