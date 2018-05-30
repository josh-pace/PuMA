#!/usr/bin/env python3
# File that runs PuMA
# Koenraad Van Doorslaer, Ken Youens-Clark, Josh Pace
# New blast function and dictionary setup


from puma_functions import *  # Importing all functions for PuMA


# -----------------------------------------------------------------------------------------
def get_args():
    args = sys.argv
    bin = os.path.dirname(args[0])
    parser = argparse.ArgumentParser(description='Displays identified protein '
                                                 'information within a given papillomavirus '
                                                 'genome.')

    parser.add_argument('-i', '--input', metavar='FILE',
                        help='Path to a genbank file formatted file that ' +
                             'contains a papillomavirus genome.',
                        required=True)

    parser.add_argument('-f', '--format', metavar='FORMAT', default='fasta',
                        help='FASTA/GenBank')

    parser.add_argument('-b', '--blastdb_dir', metavar='DIR', type=str,
                        default=os.path.join(bin, 'blast_database'),
                        help='BLAST db directory')

    parser.add_argument('-o', '--outdir', metavar='DIR', type=str,
                        default='puma-out',
                        help='Output directory (default: %(default)s)')
    parser.add_argument('-e', '--evalue', metavar='FLOAT', type=float,
                        default=0.001, help='BLAST evalue')

    parser.add_argument('-m', '--min_prot_len', metavar='NUM', type=int,
                        default=25,
                        help='Minimum protein length')

    parser.add_argument('-s', '--sites', metavar='STR', type=str,
                        default='ALL',
                        help='Comma-separated string of L1 L2 E1 E2 E4 E5 E5_delta '
                             'E5_zeta E5_epsilon E6 E7 '
                             'E9 E10 '
                             'E2BS E1BS URR ALL')

    return parser.parse_known_args()


# -----------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------
def main():
    """main"""
    #
    # Suppressing Biopython warning
    #
    warnings.simplefilter('ignore', BiopythonWarning)

    args, _ = get_args()
    sites = re.split(r'\s*,\s*', args.sites.upper())
    input_file = args.input
    out_dir = args.outdir
    blast_dir = args.blastdb_dir
    input_format = args.format.lower()
    min_prot_len = args.min_prot_len
    evalue = args.evalue

    valid_sites = set('L1 L2 E1 E2 E4 E5 E5_delta E5_zeta E5_epsilon E6 E7 E9 E10 E2BS '
                      'E1BS '
                      'URR '
                      'ALL'.split())
    if not sites:
        print('--sites is required')
        sys.exit(1)

    if not input_file:
        print("--i is required")
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

    startStop = []
    # Stores the start position of the URR and then all the possible stop positions

    URR = {}
    # Stores the imformation for the URR

    virus = {}
    # Main dictionary that will store information about proteins, URR, E2BS etc
    blasted = {}
    #Stores the blasted results to find URR start and stop

    for seq_record in SeqIO.parse("{}".format(input_file), input_format):
        Origseq = seq_record.seq
        ID = seq_record.name  # Name is accesion number
        name = seq_record.description.split(",")[0]  # description is name
        # Opening file and storing the genome

    virus['name'] = name
    virus['accession'] = ID
    virus['genome'] = Origseq
    # Adding name etc to dictionary
    #ORF = Trans_ORF(Origseq, 1, 25)  # Calling transORF function to translate genome

    print("\nThis is the gene information for {}:".format(virus['name']))



    blasted.update(blast_proteins(Origseq,min_prot_len,evalue,blast_dir,out_dir))

    virus.update(blasted)
    for protein in blasted:
        if protein == 'L1':  # Finding URR start position
            startStop.append((blasted[protein][1]))
            URRstart = blasted[protein][1]
        else:  # Getting start postions to find URR stop
            startStop.append(blasted[protein][0])

    startStop = sorted(startStop)
    #print(startStop)
    # Putting postions in increasing order to find URR stop position
    for numbers in startStop:  # Finding URR stop position
        if numbers == URRstart:
            if numbers == startStop[-1]:
                URRstop = startStop[0]
            else:
                position = startStop.index(numbers)
                URRstop = startStop[position + 1]

    #print('E1:{}'.format(virus['E1']))

    URRstart = int(URRstart) + 1
    URRstop = int(URRstop) - 1

    genomelen = int(len(Origseq))

    if URRstop == 0:
        URRstop = genomelen

    if URRstop > URRstart:
        URRfound = str(Origseq[URRstart-1:URRstop]).lower()
    # Finding the URR if it stops at or before the end of the genome

    else:
        URRfound = str(Origseq[URRstart-1:] + Origseq[:URRstop]).lower()
    # Finding the URR if it goes past the end of the genome

    if URRstop > URRstart:
        URR['URR'] = [int(URRstart), int(URRstop),URRfound]
        '''CASE WHEN URR DOES NOT GO PAST THE LENGTH OF THE GENOME. Makes the URR a 
            dictionary 
            with the key being 'URR' and the value being a list with the order of URR 
            start 
            position in genome,  the start of the genome, and the end of the URR in the 
            genome'''

    else:
        URR['URR'] = [int(URRstart), int(genomelen), 1, int(URRstop),URRfound]
        '''CASE WHEN URR GOES PAST THE LENGTH OF THE GENOME. Makes the URR a dictionary 
        with the key being 'URR' and the value being a list with the order of URR start 
        position in genome, the last position in the genome, the start of the genome, and 
        the end of the URR in the genome'''


    virus.update(URR)  # Adding URR to main dictionary

    # Calling E2BS function
    E2BS = find_E2BS(Origseq, URRfound, URRstart, ID, out_dir)
    virus.update(E2BS)

    for protein in virus:  # Getting E2 nucleotide sequence for the E4 function
        if protein == "E2":
            E2Start = virus[protein][0]
            E2Stop = virus[protein][1]
            E2seq = Origseq[E2Start - 1:E2Stop]

    E4 = find_E4(E2seq, Origseq)  # Calling E4 function

    virus.update(E4)  # Adding E4 to main dictionary

    E1BS = find_E1BS(Origseq, URRfound, URRstart, ID, out_dir)  # Calling E1BS function

    virus.update(E1BS)  # Adding E1BS to main dictionary

    if sites[0] == 'ALL':
        sites = {}
        sites.update(virus)
        del sites['name']
        del sites['genome']
        del sites['accession']


    #print(virus['URR'])


    #print(virus)

    # for protein in virus:
    #     print(protein)

    # for name in sites:
    #     if name == 'E2BS':
    #         print("\n{} E2 binding sites found:".format(len(virus['E2BS'])))
    #         for i in range(0, len(virus['E2BS'])):
    #             print('\n{} start and stop position:\n{},{}\n'.format(name,
    #                 virus[name][i], virus[name][i] + 11))
    #             print('{} sequnce:\n{}\n'.format(name, str(virus['genome'][virus['E2BS'][i]
    #                                                 - 1:virus['E2BS'][i] + 11]).lower()))
    #     elif name == 'E1BS':
    #         if type(virus[name][2]) == int:
    #             print('\n{} start and stop position:\n{},{},{},{}\n'.format(name, virus
    #             [name][0], virus[name][1], virus[name][2], virus[name][3]))
    #             print('{} seqeunce:\n{}\n'.format(name, virus[name][4]))
    #         else:
    #             print('\n{} start and stop position:\n{},{}\n'.format(name, virus[name][0]
    #                                                                   , virus[name][1]))
    #             print('{} sequnce:\n{}\n'.format(name, virus[name][2]))
    #     else:
    #         try:
    #             if type(virus[name][3]) == int:
    #                 print('\n{} start and stop position:\n{},{},{},{}\n'.format(name, virus
    #                       [name][0], virus[name][1], virus[name][2], virus[name][3]))
    #                 print('{} seqeunce:\n{}\n'.format(name, virus[name][4]))
    #                 if name != 'URR':
    #                     print('{} translated sequnce:\n{}\n'.format(name, virus[name][
    #                         5][:-1]))
    #             else:
    #                 print('\n{} start and stop position:\n{},{}\n'.format(name, virus[name][0]
    #                                                                   , virus[name][1]))
    #                 print('{} sequnce:\n{}\n'.format(name, virus[name][2]))
    #                 if name != 'URR':
    #                     print('{} translated seqeunce:\n{}\n'.format(name, virus[name][
    #                         3][:-1]))
    #         except IndexError:
    #             print('\n{} start and stop position:\n{},{}\n'.format(name, virus[name][0]
    #                                                                   , virus[name][1]))
    #             print('{} sequnce:\n{}\n'.format(name, virus[name][2]))
    #             if name != 'URR':
    #                 print('{} translated seqeunce:\n{}\n'.format(name, virus[name][
    #                     3][:-1]))


    #to_gff3(virus,genomelen,out_dir)
    #print(virus['E1BS'])
    #result_E1BS(virus)
    #print(len(virus['E4'][3]))
    #print('E1: {}'.format(virus['E1']))
    #E1_E4 = find_E1E4(virus['E1'],virus['E2'],virus['E4'],ID,Origseq,out_dir)
    #E8_E2 = find_E8E2(virus['E1'], virus['E2'], ID, Origseq, out_dir, blast_dir)

    print("testing splice function:{}".format(find_splice_acceptor( virus['E2'], ID,
                                                               blast_dir, out_dir)))

    # if E1_E4['E1^E4'][4] == 'E8':
    #     results = os.path.join('puma_results')
    #     if not os.path.isdir(results):
    #         os.makedirs(results)
    #     wrong = os.path.join(results, 'E1^E4_problem.txt')
    #     with open(wrong,'a') as file_out:
    #         # all = virus['name']
    #         # name = re.search('\(([^)]+)', all).group(1)
    #         virus['name'] = name
    #         file_out.write(virus['name'] + '\n')
    #         print('Name added for E1^E4')
    #
    #
    # if E8_E2['E8^E2'][4] == 'E8':
    #
    #     results = os.path.join('puma_results')
    #     if not os.path.isdir(results):
    #         os.makedirs(results)
    #     wrong = os.path.join(results, 'E8^E2_problem.txt')
    #     with open(wrong,'a') as file_out:
    #         # all = virus['name']
    #         # name = re.search('\(([^)]+)', all).group(1)
    #         virus['name'] = name
    #         file_out.write(virus['name'] + '\n')
    #         print('Name added for E8^E2')


    #print(virus['URR'])

    #virus.update(E1_E4)
    #virus.update(E8_E2)

    #print('E8^E2: {}'.format(virus['E8^E2']))
    #print('E8: {}'.format(virus['E8^E2'][6]))
    #print('E1^E4:{}'.format(virus['E1^E4']))
    #print('E1 for E1^E4: {}'.format(virus['E1^E4'][6]))
    #print('E1: {}'.format(virus['E1']))
    #to_results(virus)
    #check_E8_E1(virus,'E8')
    #result_E1BS(virus)
    #export_to_csv(virus)
    #print('Good')





    return


# -----------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()


