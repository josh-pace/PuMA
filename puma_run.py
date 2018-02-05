#!/usr/local/bin/python3
#File that runs PuMA
#Koenraad Van Doorslaer, Ken Youens-Clark, Josh Pace
#Version with input from Ken Youens-Clark


from puma_functions_no_absolute_path import * #Importing all functions for PuMA
#-----------------------------------------------------------------------------------------
def get_args():

    parser = argparse.ArgumentParser(description='Displays identified protein '
                                             'information within a given papillomavirus '
                                             'genome.')

    parser.add_argument('-inputfile',help='Path to a genbank file formatted file that '
                                    'contains a '
                               'papillomavirus '
                               'genome.',required=True)

    parser.add_argument('-o', '--outfile', metavar='FILE', type=str,
                        default='puma.out',
                        help='Output file (default: %(default)s)')

    parser.add_argument('-l1', action='store_true',
                        help='L1 protein, enter l1 or L1')

    parser.add_argument('-l2', action='store_true',
                        help='L2 protein, enter L2 or L2')

    parser.add_argument('-e1', action='store_true',
                        help='E1 protein, enter e1 or E1')

    parser.add_argument('-e2', action='store_true',
                        help='E1 protein, enter e2 or E2')

    parser.add_argument('-e4', action='store_true',
                        help='E4 protein, enter e4 or E4')

    parser.add_argument('-e5', action='store_true',
                        help='E5 protein, enter e5 or E5')

    parser.add_argument('-e6', action='store_true',
                        help='E6 protein, enter e6 or E6')

    parser.add_argument('-e7', action='store_true',
                        help='E7 protein, enter e7 or E7')

    parser.add_argument('-e10', action='store_true',
                        help='E10 protein, enter e10 or E10')

    parser.add_argument('-e2bs', action='store_true',
                        help='E2 binding sites, enter e2bs or E2BS')

    parser.add_argument('-e1bs', action='store_true',
                        help='E1 binding sites, enter e1bs or E1BS')

    parser.add_argument('-urr', action='store_true',
                        help='upstream regulatory region, enter urr or URR')

    parser.add_argument('-all', action='store_true',
                        help='displays all protein information, enter all or ALL')










    return  parser.parse_known_args()


#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
def main():
    """main"""

    args, unkown = get_args()
    options = []
    inputfile = args.inputfile
    outfile = args.outfile
    show_l1 = args.l1
    show_l2 = args.l2
    show_e1 = args.e1
    show_e2 = args.e2
    show_e4 = args.e4
    show_e5 = args.e5
    show_e6 = args.e6
    show_e7 = args.e7
    show_e10 = args.e10
    show_e2bs = args.e2bs
    show_e1bs = args.e1bs
    show_urr = args.urr
    show_all = args.all


    if not inputfile:
        print("--inputfile is required")
        sys.exit(1)

    if not os.path.isfile(inputfile):
        print('"{}" is not a file.'.format(inputfile))
        sys.exit(1)

    if show_l1:
        options.append('L1')
    if show_l2:
        options.append('L2')
    if show_e1:
        options.append('E1')
    if show_e2:
        options.append('E2')
    if show_e4:
        options.append('E4')
    if show_e5:
        options.append('E5')
    if show_e6:
        options.append('E6')
    if show_e7:
        options.append('E7')
    if show_e10:
        options.append('E10')
    if show_urr:
        options.append('URR')
    if show_e1bs:
        options.append('E1BS')
    if show_e2bs:
        options.append('E2BS')
    if show_all:
        options.append('ALL')




    startStop = []
    #Stores the start position of the URR and then all the possible stop positions

    URR = {}
    #Stores the imformation for the URR

    virus = {}
    # Main dictionary that will store information about proteins, URR, E2BS etc

    for seq_record in SeqIO.parse("{}" .format(inputfile), "genbank"):
        Origseq = seq_record.seq
        ID = seq_record.name #Name is accesion number
        name = seq_record.description.split(",")[0] #description is name
        # Opening file and storing the genome

    virus['name'] =name
    virus['accession'] = ID
    virus['genome'] = Origseq
    # Adding name etc to dictionary
    ORF = Trans_ORF(Origseq, 1, 25)#Calling transORF function to translate genome

    print("\nThis is the protein information for {}:\n".format(virus['name']))


    for key in ORF:#Calling Blast function for each open reading frame found
        endOfSeq = ORF[key] + ((len(key) + 1) * 3)
        blasted = Blast(key, ORF[key], endOfSeq,Origseq)
        '''key is protein sequence, ORF[key] (value of ORF) is start position, endOfSeq 
        is calculated end position'''
        if blasted != {}:
    #        print blasted
            virus.update(blasted)
            for keys in blasted:
                if keys == 'L1':#Finding URR start position
                    startStop.append((blasted[keys][1]))
                    URRstart = blasted[keys][1]
                else:#Getting start postions to find URR stop
                    startStop.append(blasted[keys][0])


    startStop = sorted(startStop)
    # Putting postions in increasing order to find URR stop position

    for numbers in startStop:#Finding URR stop position
        if numbers == URRstart:
            if numbers == startStop[-1]:
                URRstop = startStop[0]
            else:
                URRstop = startStop[numbers + 1]

    URRstart = int(URRstart) +1
    URRstop = int(URRstop) - 1

    genomelen = len(Origseq)

    if URRstop == 0:
        URRstop = genomelen


    if URRstop > URRstart:
        URRfound = Origseq[URRstart:URRstop - 1]
    #Finding the URR if it stops at or before the end of the genome

    else:
        URRfound = Origseq[URRstart:] + Origseq[:URRstop - 1]
    #Finding the URR if it goes past the end of the genome

    if URRstop > URRstart:
        URR['URR'] = [URRstart,URRstop]
        '''CASE WHEN URR DOES NOT GO PAST THE LENGTH OF THE GENOME. Makes the URR a 
            dictionary 
            with the key being 'URR' and the value being a list with the order of URR 
            start 
            position in genome,  the start of the genome, and the end of the URR in the 
            genome'''

    else:
        URR['URR'] = [URRstart,genomelen,1,URRstop]
        '''CASE WHEN URR GOES PAST THE LENGTH OF THE GENOME. Makes the URR a dictionary 
        with the key being 'URR' and the value being a list with the order of URR start 
        position in genome, the last position in the genome, the start of the genome, and 
        the end of the URR in the genome'''

    try:
        if URR['URR'][2]:
            URRaround = 'Yes'
    except IndexError:
        URRaround = 'No'
    # Finding if URR wraps around genome



    virus.update(URR)#Adding URR to main dictionary

    E2BS = find_E2BS(Origseq,URRfound, URRstart,ID)#Calling E2BS function
    #print E2BS
    virus.update(E2BS)

    for key in virus:#Getting E2 nucleotide sequence for the E4 function
        if key == "E2":
            E2Start = virus[key][0]
            E2Stop = virus[key][1]
            E2seq = Origseq[E2Start-1:E2Stop]

    E4 = find_E4(E2seq, Origseq)#Calling E4 function
    #print E4

    virus.update(E4)#Adding E4 to main dictionary

    E1BS = find_E1BS(Origseq,URRfound,URRstart,ID)#Calling E1BS function
    try:
        if E1BS['E1BS'][2]:#Finding if E1BS wraps around
            E1BSaround = 'Yes'
    except IndexError:
        E1BSaround = 'No'

    virus.update(E1BS)#Adding E1BS to main dictionary
    #print E1BS
    '''At this point in the code, everything should be found and everything below this 
    comment is working with MySQL'''

    if options:
        for name in options:
            if name == 'ALL':
                print("Working on this option")
                sys.exit(1)
            if name =='E2BS':
                print("{} E2 binding sites found:".format(len(virus['E2BS'])))
                for i in range(0,len(virus['E2BS']),1):
                    print('\n{} start and stop position:\n{},{}\n'.format(name,virus[name]
                                                                   [i],virus[name][i]+11))
                    print('{} sequnce:\n{}\n'.format(name, str(virus['genome'][virus[name]
                                                    [i] - 1:virus[name][i]+11]).lower()))
            try:
                virus[name][3]

                print('\n{} start and stop position:\n{},{},{},{}\n'.format(name,virus
                                 [name][0], virus[name][1],virus[name][2],virus[name][3]))
                print('{} sequnce:\n{}\n'.format(name, str(virus['genome'][virus[name]
                   [0]-1:] +virus['genome'][virus[name][2]-1:virus[name][ 3]-1]).lower()))
                print('{} translated sequnce:\n{}'.format(name,str(virus['genome'][
                    virus[name][0]-1:] + virus[name][virus[name][2]-1:virus[name][3]-1])))


            except IndexError:
                print('\n{} start and stop position:\n{},{}\n'.format(name, virus[name][0]
                                                                        , virus[name][1]))
                print('{} sequnce:\n{}\n'.format(name, str(virus['genome'][virus[name]
                                                        [0] - 1:virus[name][1]]).lower()))
                print('{} translated sequnce:\n{}'.format(name, Seq(str(virus['genome']
                                 [virus[name][0] - 1:virus[name][1]])).translate(1)[:-1]))




    #export_to_mysql(virus,E1BSaround,URRaround)


    #accessionList=get_accession()


    # for value in accessionList:
    #     export_to_csv(get_genome_info(value))
    #

    #print export_from_mysql('E2BS','BPV1REF')

    #result = get_genome_info('BPV1REF')

    #
    #export_to_csv(get_genome_info('AaPV1REF'))
  
    # for values in result['E2BS']:
    #     print values
#-----------------------------------------------------------------------------------------


if __name__ == '__main__':
    main()


