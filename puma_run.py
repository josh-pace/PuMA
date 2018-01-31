#!/usr/bin/python2.7
#File that runs PuMA
#Koenraad Van Doorslaer, Ken Youens-Clark, Josh Pace
#Version with input from Ken Youens-Clark


from puma_functions import * #Importing all functions for PuMA
#-----------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='This program displays identified protein '
                                             'information within a given papillomavirus '
                                             'genome.')

parser.add_argument('-inputfile',help='Path to a genbank file formatted file that '
                                    'contains a '
                               'papillomavirus '
                               'genome.',required=True)

args = parser.parse_args()


#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
def main():
    startStop = []
    #Stores the start position of the URR and then all the possible stop positions

    URR = {}
    #Stores the imformation for the URR

    virus = {}
    # Main dictionary that will store information about proteins, URR, E2BS etc

    for seq_record in SeqIO.parse("/Users/joshpace/Downloads/BPV1_REF.txt", "genbank"):
        Origseq = seq_record.seq
        ID = seq_record.name #Name is accesion number
        name = seq_record.description.split(",")[0] #description is name
        # Opening file and storing the genome

    virus['name'] =name
    virus['accession'] = ID
    virus['genome'] = Origseq
    # Adding name etc to dictionary
    ORF = Trans_ORF(Origseq, 1, 25)#Calling transORF function to translate genome



    for key in ORF:#Calling Blast function for each open reading frame found
        endOfSeq = ORF[key] + ((len(key) + 1) * 3)
        blasted = Blast(key, ORF[key], endOfSeq,Origseq)
        '''key is protein sequence, ORF[key] (value of ORF) is start position, endOfSeq 
        is calculated end position'''
        if blasted != {}:
            print blasted
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

    virus.update(E2BS)

    for key in virus:#Getting E2 nucleotide sequence for the E4 function
        if key == "E2":
            E2Start = virus[key][0]
            E2Stop = virus[key][1]
            E2seq = Origseq[E2Start-1:E2Stop]

    E4 = find_E4(E2seq, Origseq)#Calling E4 function


    virus.update(E4)#Adding E4 to main dictionary

    E1BS = find_E1BS(Origseq,URRfound,URRstart,ID)#Calling E1BS function
    try:
        if E1BS['E1BS'][2]:#Finding if E1BS wraps around
            E1BSaround = 'Yes'
    except IndexError:
        E1BSaround = 'No'

    virus.update(E1BS)#Adding E1BS to main dictionary

    '''At this point in the code, everything should be found and everything below this 
    comment is working with MySQL'''


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


