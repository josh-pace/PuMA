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
                        help='Comma-separated string of L1 L2 E1 E2 E4 E5 E6 E7 E10 E2BS E1BS URR ALL')

    return parser.parse_known_args()


# -----------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------
def main():
    """main"""
    #
    # Supressing Biopython warning
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

    valid_sites = set('L1 L2 E1 E2 E4 E5 E6 E7 E10 E2BS E1BS URR ALL'.split())
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

    print("\nThis is the protein information for {}:".format(virus['name']))

    #
    # Calling Blast function for each open reading frame found
    #

    # for key in ORF:
    #     endOfSeq = ORF[key] + ((len(key) + 1) * 3)
    #     blasted = Blast(key, ORF[key], endOfSeq, Origseq, blast_dir, out_dir)
    #     '''key is protein sequence, ORF[key] (value of ORF) is start position, endOfSeq
    #     is calculated end position'''
    #     if blasted != {}:
    #         #        print blasted
    #         virus.update(blasted)
    #         for keys in blasted:
    #             if keys == 'L1':  # Finding URR start position
    #                 startStop.append((blasted[keys][1]))
    #                 URRstart = blasted[keys][1]
    #             else:  # Getting start postions to find URR stop
    #                 startStop.append(blasted[keys][0])


    blasted.update(blast_proteins(Origseq,25,0.001,blast_dir,out_dir))
    #print(virus)
    virus.update(blasted)
    for keys in blasted:
        if keys == 'L1':  # Finding URR start position
            startStop.append((blasted[keys][1]))
            URRstart = blasted[keys][1]
        else:  # Getting start postions to find URR stop
            startStop.append(blasted[keys][0])

    startStop = sorted(startStop)

    # Putting postions in increasing order to find URR stop position

    for numbers in startStop:  # Finding URR stop position
        if numbers == URRstart:
            if numbers == startStop[-1]:
                URRstop = startStop[0]
            else:
                URRstop = startStop[numbers + 1]

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

    for key in virus:  # Getting E2 nucleotide sequence for the E4 function
        if key == "E2":
            E2Start = virus[key][0]
            E2Stop = virus[key][1]
            E2seq = Origseq[E2Start - 1:E2Stop]

    E4 = find_E4(E2seq, Origseq)  # Calling E4 function
    # print E4

    virus.update(E4)  # Adding E4 to main dictionary

    E1BS = find_E1BS(Origseq, URRfound, URRstart, ID, out_dir)  # Calling E1BS function
    try:
        if E1BS['E1BS'][2]:  # Finding if E1BS wraps around
            E1BSaround = 'Yes'
    except IndexError:
        E1BSaround = 'No'

    virus.update(E1BS)  # Adding E1BS to main dictionary
    # print E1BS
    '''At this point in the code, everything should be found and everything below this 
    comment is working with MySQL'''

    if sites[0] == 'ALL':
        sites = {}
        sites.update(virus)
        del sites['name']
        del sites['genome']
        del sites['accession']
        del sites['E1']
        del sites['E6']
        del sites['E2']
        del sites['L2']
        del sites['E7']
        del sites['URR']
        del sites['E4']
        del sites['L1']
        del sites['E1BS']




    for name in sites:
        if name == 'E2BS':
            print("\n{} E2 binding sites found:".format(len(virus['E2BS'])))
            for i in range(0, len(virus['E2BS'])):
                print('\n{} start and stop position:\n{},{}\n'.format(name,
                                                                      virus[name][i]
                                                                , virus[name][i] + 11))
                print('{} sequnce:\n{}\n'.format(name, str(virus['genome'][virus['E2BS'][i]
                                                    - 1:virus['E2BS'][i] + 11]).lower()))
        elif name == 'E1BS':
            if type(virus[name][2]) == int:
                print('\n{} start and stop position:\n{},{},{},{}\n'.format(name, virus
                [name][0], virus[name][1], virus[name][2], virus[name][3]))
                print('{} seqeunce:\n{}\n'.format(name, virus[name][4]))
            else:
                print('\n{} start and stop position:\n{},{}\n'.format(name, virus[name][0]
                                                                      , virus[name][1]))
                print('{} sequnce:\n{}\n'.format(name, virus[name][2]))
        else:

                if type(virus[name][3]) == int:
                    print('\n{} start and stop position:\n{},{},{},{}\n'.format(name, virus
                          [name][0], virus[name][1], virus[name][2], virus[name][3]))
                    print('{} seqeunce:\n{}\n'.format(name, virus[name][4]))
                    if name != 'URR':
                        print('{} translated sequnce:\n{}\n'.format(name, virus[name][
                            5][:-1]))
                else:
                    print('\n{} start and stop position:\n{},{}\n'.format(name, virus[name][0]
                                                                      , virus[name][1]))
                    print('{} sequnce:\n{}\n'.format(name, virus[name][2]))
                    if name != 'URR':
                        print('{} translated seqeunce:\n{}\n'.format(name, virus[name][
                            3][:-1]))



# -----------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()


