#All functions needed for PuMA
#Koenraad Van Doorslaer, Ken Youens-Clark, Josh Pace
#Version with input from Ken Youens-Clark

from Bio import SeqIO, GenBank, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
import os, glob, re, csv, time, operator, argparse, sys
import pymysql

#
# This function will translate the genome in all 3 frames and 
# split the ORF on the stop codon (*)
#
def Trans_ORF(seq, trans_table, min_protein_length):
    ORFs = {}
    for frame in range(3):  # Accounting for all 3 different frames
        trans = str(seq[frame:].translate(trans_table))
        trans_len = len(trans)
        aa_start = 0
        aa_end = 0
        while aa_start < trans_len:
            aa_end = trans.find("*", aa_start  )# returns the lowest index where "*" is found in aa_start because "*" represents a stop codon
            if aa_end == -1:# .find returns -1 if no "*" is found in aa_start
                aa_end = trans_len
            if aa_end- aa_start >= min_protein_length:  # Only finding proteins over a certain length
                start = frame + aa_start * 3  # Multiplying by 3 to get the start position of the nucleotides
                ORFs[trans[aa_start:aa_end]] = start
            aa_start = aa_end + 1

    return ORFs


#
# This function uses blast to find the different proteins in 
# the translated genome
#
def Blast(protein_sequence, start, end, genomic_sequence, blast_dir, out_dir):
    result = {}
    M = re.search('M', protein_sequence)
    if M:
        query = protein_sequence[M.start():]
        query = query + "*"
        with open(os.path.join(out_dir, "tempORF.txt"), "w") as tempfile:
            tempfile.write('>Blasting\n')
            tempfile.write(query)
            # print >> tempfile, '>Blasting'
            # print >> tempfile, query
        tempfile = "tempORF.txt"
        blast_db = os.path.join(blast_dir, "blast_database.txt")

        if not os.path.isfile(blast_db):
            print('BLAST db "{}" does not exist'.format(blast_db))
            return

        blast_out = os.path.join(out_dir, "Blasted.xml")
        if os.path.isfile(blast_out):
            os.remove(blast_out)

        blasting = blastp(query=tempfile, 
                          db=blast_db,
                          evalue=0.001,
                          outfmt=5,
                          out=blast_out)
        blasting()

        if not os.path.isfile(blast_out):
            print('No BLAST output "{}" (must have failed)'.format(blast_out))
            return

        if not os.path.getsize(blast_out):
            print('No hits!')
            return

        print('blast_out "{}"'.format(blast_out))
        blastedfile = open(blast_out, 'rt')
        blasted = str(blastedfile.read())
        print(blasted)
        DEF = re.search("<Hit_def>((.*))</Hit_def>", blasted)
        if DEF:
            if DEF.group(1) == 'L1':
                real_start = start + M.start() + M.start() + M.start()
                result['L1'] = []
                L1_pre = genomic_sequence[(start + 3 * M.start()):int(end)]
                splice = '(C|T)(C|T)(A|C|G|T)(C|T)AG(A)TG'
                spliced = re.search(splice, str(L1_pre))
                if spliced:
                    start_L1 = int(spliced.start()) + 6
                    if start_L1 % 3 == 0:
                        if start_L1 > 600:
                            L1_post = L1_pre
                            result['L1'] = [int(start_L1), int(end)]  #, str(L1_post),
                            # Seq(str(L1_post)).translate()]
                        else:
                            L1_post = L1_pre[start_L1:]
                            result['L1'] = [int(real_start) + 1 + int(start_L1),
                                int(end)] #, str(L1_post), Seq(str(L1_post)).translate()]
                    else:
                        L1_post = L1_pre
                        result['L1'] = [int(real_start ) +1, int(end)] #, L1_post,
                        # Seq(str(
                        # L1_post)).translate()]
                else:
                    L1_post = L1_pre
                    result['L1'] = [int(real_start) +1,int(end)] #, L1_post, Seq(str(
                    # L1_post)).translate()]
            else:
                real_start = start + M.start() + M.start() + M.start()
                name = re.search(r'random', DEF.group(1), re.I)
                if name:
                    pass
                else:
                    result[DEF.group(1)] = [int(real_start) + 1, int(end)] #,
                    # genomic_sequence[int(real_start):int(end)], query]


    return result


# --------------------------------------------------
# This function finds the E2BS in a genome using the URR
#
def find_E2BS(genome, URR, URRstart, ID, out_dir):
    genomeLength = len(genome) #Getting length of genome
    startListURR = []  # Storing the nucleotide start positions in URR of the E2BS
    startListGenome= []  # Storing the nucleotide start positions in genome of the E2BS
    E2BS = {} #Storing all E2BS

    # Writting URR to a file so FIMO can be used
    tmp = os.path.join(out_dir, "PuMA_URR_tempfile.fa")
    with open(tmp, "w") as tempfile:
        # print >> tempfile, '>URR for %s' %ID
        # print >> tempfile, URR
        tempfile.write('>URR for {}\n'.format(ID))
        tempfile.write(str(URR))
        
    # Executing FIMO, Using URR
    fimo_dir = os.path.join(out_dir, 'E2BS')
    fimo_cmd = 'fimo --oc {} --norc --verbosity 1 --thresh 1.0E-3 {} {}'
    cline = (fimo_cmd.format(fimo_dir, 'meme_3000_TOTAL.txt', tmp))

    os.system(str(cline))  # Executing FIMO

    fimo_out = os.path.join(fimo_dir, 'fimo.txt')

    if not os.path.isfile(fimo_out):
        print('Failed to create fimo out "{}"'.format(fimo_out))
        return

    # Getting nucleotide start positions from FIMO output file
    for column in csv.reader(open(fimo_out, "rU"), delimiter='\t'):  
        if column[3] == 'start':
            startListURR = []
        else:
            startListURR.append(column[3])

    startListURR = list(map(int, set(startListURR)))#Making the positions an integer
    # value to use as index and only using every unique start position
    startListGenome = list(map(int, startListGenome))#Making the positions an integer
    # value


    for i in range(0,len(startListURR),1): #Finding the positions of E2BS in genome
        genomestart = startListURR[i]
        genomestart = (genomestart + URRstart)
        if genomestart > genomeLength:
            genomestart = genomestart - genomeLength
            startListGenome.append(genomestart)
        else:
            startListGenome.append(genomestart)


    # for value in startListGenome:#Accounting for the case where the sequence is longer then the genome length
    #     if ((value + 12) > genomeLength):
    #         stop = (value + 12) - genomeLength
    #         index = startListGenome.index(value)
    #         startListGenome[index] = (str(value) + '..' + str(stop))

    E2BS['E2BS'] = startListGenome#Putting values into dictionary

    return E2BS


def find_E4(E2, genomic_sequence):  # Finds E4
    E4 = {}  # Storing E4 information
    trans = E2[1:len(E2)].translate()#Translates E2 nucleotide sequence
    E4protein = max(trans.split("*"), key=len)#Splits sequences on the stop codon, takes longest sequence
    E4_start = re.search(str(E4protein), str(trans)).start()#Finding the start position of E4
    E4_end = re.search(str(E4protein), str(trans)).end()#Finding the end position of E4
    E4_nt = str(E2[(E4_start * 3) + 1:((E4_end + 1) * 3) + 1])#Getting the E4 nucleotide sequence
    E4_nt_start = re.search(E4_nt, str(genomic_sequence)).start()#Finding nuceotide start position of E4
    E4_nt_end = E4_nt_start + len(E4_nt)#Finding nucleotide end position of E4
    E4['E4'] = [E4_nt_start + 1, E4_nt_end]# Storing all information into a dictionary
    return E4

def find_E1BS(genome,URR,URRstart,ID):
    genomeLength = len(genome)
    E1BS = {} #Storing E1BS
    startURR = 0

    with open("PuMA_URR_tempfile.fa", "w") as tempfile:#Writting URR to a file so FIMO can be used
        tempfile.write('>URR for {}\n'.format(ID))
        tempfile.write(str(URR))
        # print >> tempfile, '>URR for %s' %ID
        # print >> tempfile, URR
    cline = ("fimo --oc E1BS --norc --verbosity 1 --thresh 1.0E-4 --bgfile background_model_E1BS.txt meme_E1BS_1motif_18_21.txt PuMA_URR_tempfile.fa")#Executing FIMO, Using URR
    os.system(str(cline))  # Executing FIMO

    for column in csv.reader(open("E1BS/fimo.txt", "rU"), delimiter='\t'):  # Getting nucleotide start positions from FIMO output file
        if column[3] == 'start':
            startListURR = []
        else:
            startListURR.append(column[3])

    startURR = int(startListURR[0])#Finding the start position of URR
    genomestart = (startURR + URRstart)

    if genomestart > genomeLength:#Case where the start of E1BS is after the end of the genome
        genomestart = genomestart - genomeLength

    genomestop = genomestart +19

    if genomestop > genomeLength:#For printing all info
        genomestop = genomestop - genomeLength
        E1BS['E1BS'] = [genomestart, genomeLength, 1, genomestop]

    else:
        E1BS['E1BS'] = [genomestart, genomestop]



    # if genomestop > genomeLength:#For printing all info
    #     genomestop = genomestop - genomeLength
    #     #print genomestop
    #     E1BS[genomestart] = genome[genomestart-1:] + genome[:genomestop]
    # else:
    #     E1BS[genomestart] = genome[genomestart - 1:genomestop]


    return E1BS

def export_to_mysql(dict,E1BSlong,URRlong):#Exports start and stop positions for each protein etc to a MySQL database
    connection = pymysql.connect(host='localhost', user='root', password='HobbesSadie96', db='Practice')

    # with connection.cursor() as cursor:#Without genome
    #     sql = "INSERT INTO PuMA_Test (accession,name,URR,L1,L2,E1,E2,E4,E5,E6,E7,E1BS,E2BS) VALUE (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
    #     if URRlong == 'No' and E1BSlong == 'No':#Easy case (URR and E1BS don't wrap around)
    #         cursor.execute(sql, (str(dict['accession'][0]),str(dict['name'][0]),(str(dict['URR'][0]) + ".." + str(dict['URR'][1])),(str(dict['L1'][0]) + ".." + str(dict['L1'][1])),(str(dict['L2'][0]) + ".." + str(dict['L2'][1])),(str(dict['E1'][0]) + ".." + str(dict['E1'][1])),(str(dict['E2'][0]) + ".." + str(dict['E2'][1])),(str(dict['E4'][0]) + ".." + str(dict['E4'][1])),(str(dict['E5'][0]) + ".." + str(dict['E5'][1])),(str(dict['E6'][0]) + ".." + str(dict['E6'][1])),(str(dict['E7'][0]) + ".." + str(dict['E7'][1])),(str(dict['E1BS'][0]) + ".." + str(dict['E1BS'][1])),str(dict['E2BS'])))
    #
    #     elif URRlong == 'Yes' and E1BSlong == 'Yes':#Case where both URR and E1BS wrap around
    #         cursor.execute(sql, (str(dict['accession'][0]),str(dict['name'][0]),(str(dict['URR'][0]) + ".." + str(dict['URR'][1]) + "," + str(dict['URR'][2])+ ".." + str(dict['URR'][3])),(str(dict['L1'][0]) + ".." + str(dict['L1'][1])),(str(dict['L2'][0]) + ".." + str(dict['L2'][1])),(str(dict['E1'][0]) + ".." + str(dict['E1'][1])),(str(dict['E2'][0]) + ".." + str(dict['E2'][1])),(str(dict['E4'][0]) + ".." + str(dict['E4'][1])),(str(dict['E5'][0]) + ".." + str(dict['E5'][1])),(str(dict['E6'][0]) + ".." + str(dict['E6'][1])),(str(dict['E7'][0]) + ".." + str(dict['E7'][1])),(str(dict['E1BS'][0]) + ".." + str(dict['E1BS'][1])+","+str(dict['E1BS'][2])+".."+str(dict['E1BS'][3])),str(dict['E2BS'])))
    #         print 1
    #     elif URRlong =='No' and E1BSlong == 'Yes':#Case where only E1BS wrap around
    #         cursor.execute(sql, (str(dict['accession'][0]),str(dict['name'][0]),(str(dict['URR'][0]) + ".." + str(dict['URR'][1])),(str(dict['L1'][0]) + ".." + str(dict['L1'][1])),(str(dict['L2'][0]) + ".." + str(dict['L2'][1])),(str(dict['E1'][0]) + ".." + str(dict['E1'][1])),(str(dict['E2'][0]) + ".." + str(dict['E2'][1])),(str(dict['E4'][0]) + ".." + str(dict['E4'][1])),(str(dict['E5'][0]) + ".." + str(dict['E5'][1])),(str(dict['E6'][0]) + ".." + str(dict['E6'][1])),(str(dict['E7'][0]) + ".." + str(dict['E7'][1])),(str(dict['E1BS'][0]) + ".." + str(dict['E1BS'][1])+","+str(dict['E1BS'][2])+".."+str(dict['E1BS'][3])),str(dict['E2BS'])))
    #
    #
    #     elif URRlong == 'Yes' and E1BSlong == 'No':#Case where only URR wraps around
    #         cursor.execute(sql, (str(dict['accession'][0]),str(dict['name'][0]),(str(dict['URR'][0]) + ".." + str(dict['URR'][1]) + "," + str(dict['URR'][2])+ ".." + str(dict['URR'][3])),(str(dict['L1'][0]) + ".." + str(dict['L1'][1])),(str(dict['L2'][0]) + ".." + str(dict['L2'][1])),(str(dict['E1'][0]) + ".." + str(dict['E1'][1])),(str(dict['E2'][0]) + ".." + str(dict['E2'][1])),(str(dict['E4'][0]) + ".." + str(dict['E4'][1])),(str(dict['E5'][0]) + ".." + str(dict['E5'][1])),(str(dict['E6'][0]) + ".." + str(dict['E6'][1])),(str(dict['E7'][0]) + ".." + str(dict['E7'][1])),(str(dict['E1BS'][0]) + ".." + str(dict['E1BS'][1])),str(dict['E2BS'])))

    with connection.cursor() as cursor:#With genome
        sql = "INSERT INTO PuMA_Test (accession,name,genome,URR,L1,L2,E1,E2,E4,E5,E6,E7,E1BS,E2BS) VALUE (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
        if URRlong == 'No' and E1BSlong == 'No':#Easy case (URR and E1BS don't wrap around)
            cursor.execute(sql, (str(dict['accession'][0]),str(dict['name'][0]),str(dict['genome']),(str(dict['URR'][0]) + ".." + str(dict['URR'][1])),(str(dict['L1'][0]) + ".." + str(dict['L1'][1])),(str(dict['L2'][0]) + ".." + str(dict['L2'][1])),(str(dict['E1'][0]) + ".." + str(dict['E1'][1])),(str(dict['E2'][0]) + ".." + str(dict['E2'][1])),(str(dict['E4'][0]) + ".." + str(dict['E4'][1])),(str(dict['E5'][0]) + ".." + str(dict['E5'][1])),(str(dict['E6'][0]) + ".." + str(dict['E6'][1])),(str(dict['E7'][0]) + ".." + str(dict['E7'][1])),(str(dict['E1BS'][0]) + ".." + str(dict['E1BS'][1])),str(dict['E2BS'])))

        elif URRlong == 'Yes' and E1BSlong == 'Yes':#Case where both URR and E1BS wrap around
            cursor.execute(sql, (str(dict['accession'][0]),str(dict['name'][0]),str(dict['genome']),(str(dict['URR'][0]) + ".." + str(dict['URR'][1]) + "," + str(dict['URR'][2])+ ".." + str(dict['URR'][3])),(str(dict['L1'][0]) + ".." + str(dict['L1'][1])),(str(dict['L2'][0]) + ".." + str(dict['L2'][1])),(str(dict['E1'][0]) + ".." + str(dict['E1'][1])),(str(dict['E2'][0]) + ".." + str(dict['E2'][1])),(str(dict['E4'][0]) + ".." + str(dict['E4'][1])),(str(dict['E5'][0]) + ".." + str(dict['E5'][1])),(str(dict['E6'][0]) + ".." + str(dict['E6'][1])),(str(dict['E7'][0]) + ".." + str(dict['E7'][1])),(str(dict['E1BS'][0]) + ".." + str(dict['E1BS'][1])+","+str(dict['E1BS'][2])+".."+str(dict['E1BS'][3])),str(dict['E2BS'])))

        elif URRlong =='No' and E1BSlong == 'Yes':#Case where only E1BS wrap around
            cursor.execute(sql, (str(dict['accession'][0]),str(dict['name'][0]),str(dict['genome']),(str(dict['URR'][0]) + ".." + str(dict['URR'][1])),(str(dict['L1'][0]) + ".." + str(dict['L1'][1])),(str(dict['L2'][0]) + ".." + str(dict['L2'][1])),(str(dict['E1'][0]) + ".." + str(dict['E1'][1])),(str(dict['E2'][0]) + ".." + str(dict['E2'][1])),(str(dict['E4'][0]) + ".." + str(dict['E4'][1])),(str(dict['E5'][0]) + ".." + str(dict['E5'][1])),(str(dict['E6'][0]) + ".." + str(dict['E6'][1])),(str(dict['E7'][0]) + ".." + str(dict['E7'][1])),(str(dict['E1BS'][0]) + ".." + str(dict['E1BS'][1])+","+str(dict['E1BS'][2])+".."+str(dict['E1BS'][3])),str(dict['E2BS'])))


        elif URRlong == 'Yes' and E1BSlong == 'No':#Case where only URR wraps around
            cursor.execute(sql, (str(dict['accession'][0]),str(dict['name'][0]),str(dict['genome']),(str(dict['URR'][0]) + ".." + str(dict['URR'][1]) + "," + str(dict['URR'][2])+ ".." + str(dict['URR'][3])),(str(dict['L1'][0]) + ".." + str(dict['L1'][1])),(str(dict['L2'][0]) + ".." + str(dict['L2'][1])),(str(dict['E1'][0]) + ".." + str(dict['E1'][1])),(str(dict['E2'][0]) + ".." + str(dict['E2'][1])),(str(dict['E4'][0]) + ".." + str(dict['E4'][1])),(str(dict['E5'][0]) + ".." + str(dict['E5'][1])),(str(dict['E6'][0]) + ".." + str(dict['E6'][1])),(str(dict['E7'][0]) + ".." + str(dict['E7'][1])),(str(dict['E1BS'][0]) + ".." + str(dict['E1BS'][1])),str(dict['E2BS'])))

    connection.commit()
    connection.close()

    return

def get_accession():#Retrieves all accession numbers stored in MySQL database
    connection = pymysql.connect(host='localhost', user='root', password='HobbesSadie96', db='Practice')
    accession = []
    with connection.cursor() as cursor:
        sql = "SELECT accession FROM PuMA_TEST"
        cursor.execute(sql, accession)
        numbers = cursor.fetchall()
        for i in range(0,len(numbers),1):
            accession.append(numbers[i][0])




    connection.close()

    return accession





def export_from_mysql(choice,accession):#Exports start and stop positions for each protein etc from a MySQL database

    connection = pymysql.connect(host='localhost', user='root', password='HobbesSadie96', db='Practice')

    with connection.cursor() as cursor:
        if choice == 'name':
            name ={}
            sql = "SELECT name FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql, accession)
            Name = str(cursor.fetchall()).split("'")[1]
            name['name'] = Name
            return name

        elif choice == 'genome':
            genome = {}
            sql = "SELECT genome FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql, accession)
            Genome = str(cursor.fetchall()).split("'")[1]
            genome['genome'] = Genome
            return genome

        elif choice == 'URR':
            URRlist = []
            URR ={}
            sql = "SELECT URR FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql,accession)
            positions = str(str(str(str(cursor.fetchall()).split("'")[1]).split(",")).split("[")[1]).split("]")[0]
            one = int(str(positions.split("'")[1]).split("..")[0])#Finds position one
            two = int(str(positions.split("'")[1]).split("..")[1])#Finds position two
            URRlist.append(one)
            URRlist.append(two)
            try:#Checking to see if the sequence wraps around the genome, if it does, returns all four positions
                three = int(str(str(positions.split(",")).split("'")[3]).split("..")[0])
                four = int(str(str(positions.split(",")).split("'")[3]).split("..")[1])
                URRlist.append(three)
                URRlist.append(four)
                URR['URR'] = URRlist
                return URR
            except IndexError:#If sequence does not wrap around,returns the two positions found first
                URR['URR'] = URRlist
                return URR

        elif choice == 'L1':
            L1list =[]
            L1 ={}
            sql = "SELECT L1 FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql,accession)
            L1list = str(cursor.fetchall()).split("'")[1].split("..")
            L1list = map(int, L1list)
            L1['L1'] = L1list
            return L1

        elif choice == 'L2':
            L2list = []
            L2 ={}
            sql = "SELECT L2 FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql,accession)
            L2list = str(cursor.fetchall()).split("'")[1].split("..")
            L2list = map(int, L2list)
            L2['L2'] = L2list
            return L2

        elif choice == 'E1':
            E1list = []
            E1 ={}
            sql = "SELECT E1 FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql,accession)
            E1list = str(cursor.fetchall()).split("'")[1].split("..")
            E1list = map(int, E1list)
            E1['E1'] = E1list
            return E1

        elif choice == 'E2':
            E2list =[]
            E2 ={}
            sql = "SELECT E2 FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql,accession)
            E2list = str(cursor.fetchall()).split("'")[1].split("..")
            E2list = map(int, E2list)
            E2['E2'] = E2list
            return E2

        elif choice == 'E4':
            E4list = []
            E4 = {}
            sql = "SELECT E4 FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql,accession)
            E4list = str(cursor.fetchall()).split("'")[1].split("..")
            E4list = map(int, E4list)
            E4['E4']=E4list
            return E4

        elif choice == 'E5':
            E5list =[]
            E5 ={}
            sql = "SELECT E5 FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql,accession)
            E5list = str(cursor.fetchall()).split("'")[1].split("..")
            E5list = map(int, E5list)
            E5['E5']=E5list
            return E5

        elif choice == 'E6':
            E6list = []
            E6 = {}
            sql = "SELECT E6 FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql,accession)
            E6list = str(cursor.fetchall()).split("'")[1].split("..")
            E6list = map(int, E6list)
            E6['E6']=E6list
            return E6

        elif choice == 'E7':
            E7list = []
            E7 ={}
            sql = "SELECT E7 FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql,accession)
            E7list = str(cursor.fetchall()).split("'")[1].split("..")
            E7list = map(int, E7list)
            E7['E7'] = E7list
            return E7

        elif choice == 'E10':
            E10list = []
            E10 = {}
            sql = "SELECT E10 FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql,accession)
            positions = str(cursor.fetchall())
            if positions == '((None,),)':#Case where E10 does not exist in genome
                E10list = 'Null'
            else:#Case where E10 exists in genome
                one = int(str(positions.split("'")[1]).split("..")[0])
                two = int(str(positions.split("'")[1]).split("..")[1])
                E10list.append(one)
                E10list.append(two)
            E10['E10'] = E10list
            return E10

        elif choice == 'E1BS':
            E1BSlist =[]
            E1BS = {}
            sql = "SELECT E1BS FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql, accession)
            positions = str(cursor.fetchall())
            one = int(str(str(str(positions.split("'")[1]).split(",")).split("'")[1]).split("..")[0])#Finds position 1 and two
            two = int(str(str(str(positions.split("'")[1]).split(",")).split("'")[1]).split("..")[1])
            E1BSlist.append(one)
            E1BSlist.append(two)
            try:#Checking to see if the sequence wraps around the genome, if it does, returns all four values
                three = int(str(str(str(positions.split("'")[1]).split(",")).split("'")[3]).split("..")[0])
                four = int(str(str(str(positions.split("'")[1]).split(",")).split("'")[3]).split("..")[1])
                E1BSlist.append(three)
                E1BSlist.append(four)
                E1BS['E1BS']= E1BSlist
                return E1BS
            except IndexError:#If sequence does not wrap around,returns two values
                E1BS['E1BS'] = E1BSlist
                return E1BS

        elif choice == 'E2BS':#Figure out
            E2BSlist = []
            E2BS = {}
            sql = "SELECT E2BS FROM PuMA_TEST WHERE accession = %s"
            cursor.execute(sql,accession)
            sites = cursor.fetchall()
            for value in sites[0][0][1:-1].split(", "):
                E2BSlist.append(int(value))
            E2BS['E2BS'] = E2BSlist
            return E2BS



    connection.close()





def import_from_mysql(accession):#Function for testing retrevial methods from MySQL
    L1 = ''
    connection = pymysql.connect(host='localhost', user='root', password='HobbesSadie96', db='Practice')

    with connection.cursor() as cursor:
        sql = "SELECT E2BS FROM PuMA_TEST WHERE accession = %s"
        cursor.execute(sql, accession)
        L1 = cursor.fetchall()#.split("'")[1].split("..")
        #L1 = map(int,L1)
        #L1.append(str(cursor.fetchall()).split("'")[1].split("..")[1])



    connection.close()




    return L1

def get_genome_info(acccesion):#Gets all genome information from MySQL database given an accession number and puts all information into a dictionary
        sequences = {}
        options = ['name','genome','URR','L1','L2','E1','E2','E4','E5','E6','E7','E10','E1BS','E2BS']#NEED GENOME FOR AFTER TESTING
        sequences['accession'] = acccesion
        for values in options:
            sequences.update(export_from_mysql(values,acccesion))

        sequences['L1'] = [sequences['L1'][0],sequences['L1'][1],str(sequences['genome'][sequences['L1'][0]-1:sequences['L1'][1]]),Seq(str(sequences['genome'][sequences['L1'][0]-1:sequences['L1'][1]])).translate(1)]
        sequences['L2'] = [sequences['L2'][0],sequences['L2'][1],str(sequences['genome'][sequences['L2'][0]-1:sequences['L2'][1]]),Seq(str(sequences['genome'][sequences['L2'][0]-1:sequences['L2'][1]])).translate(1)]
        sequences['E1'] = [sequences['E1'][0],sequences['E1'][1],str(sequences['genome'][sequences['E1'][0]-1:sequences['E1'][1]]),Seq(str(sequences['genome'][sequences['E1'][0]-1:sequences['E1'][1]])).translate(1)]
        sequences['E2'] = [sequences['E2'][0],sequences['E2'][1],str(sequences['genome'][sequences['E2'][0]-1:sequences['E2'][1]]),Seq(str(sequences['genome'][sequences['E2'][0]-1:sequences['E2'][1]])).translate(1)]
        sequences['E4'] = [sequences['E4'][0],sequences['E4'][1],str(sequences['genome'][sequences['E4'][0]-1:sequences['E4'][1]]),Seq(str(sequences['genome'][sequences['E4'][0]-1:sequences['E4'][1]])).translate(1)]
        sequences['E5'] = [sequences['E5'][0],sequences['E5'][1],str(sequences['genome'][sequences['E5'][0]-1:sequences['E5'][1]]),Seq(str(sequences['genome'][sequences['E5'][0]-1:sequences['E5'][1]])).translate(1)]
        sequences['E6'] = [sequences['E6'][0],sequences['E6'][1],str(sequences['genome'][sequences['E6'][0]-1:sequences['E6'][1]]),Seq(str(sequences['genome'][sequences['E6'][0]-1:sequences['E6'][1]])).translate(1)]
        sequences['E7'] = [sequences['E7'][0],sequences['E7'][1],str(sequences['genome'][sequences['E7'][0]-1:sequences['E7'][1]]),Seq(str(sequences['genome'][sequences['E7'][0]-1:sequences['E7'][1]])).translate(1)]

        if sequences['E10'] == 'Null':#Seeing if E10 was found, if not it it deleted from sequence dictionary
            del sequences['E10']
        else:#If E10 exists then sequence and translated sequence are found
            sequences['E10'] = [sequences['E10'][0],sequences['E10'][1],str(sequences['genome'][sequences['E10'][0]-1:sequences['E10'][1]]),Seq(str(sequences['genome'][sequences['E10'][0]-1:sequences['E10'][1]])).translate(1)]
        try:#Case where URR wraps around,Need to look at the sequence slicing
            sequences['URR'] = [sequences['URR'][0], sequences['URR'][1],sequences['URR'][2],sequences['URR'][3], str(sequences['genome'][sequences['URR'][0] - 1:] + sequences['genome'][sequences['URR'][2] - 1:sequences['URR'][3]-1])]
        except IndexError:#Case where URR does not wrap around, Need to look at the sequence slicing
            sequences['URR'] = [sequences['URR'][0], sequences['URR'][1], str(sequences['genome'][sequences['URR'][0] - 1:sequences['URR'][1]-1])]
        try:#Case where E1BS wraps around,Need to look at the sequence slicing
            sequences['E1BS'] = [sequences['E1BS'][0], sequences['E1BS'][1],sequences['E1BS'][2],sequences['E1BS'][3], str(sequences['genome'][sequences['E1BS'][0]-2:] + sequences['genome'][sequences['E1BS'][2]-1:sequences['E1BS'][3]-1])]
        except IndexError:#Case where E1BS does not wrap around,Need to look at the sequence slicing
            sequences['E1BS'] = [sequences['E1BS'][0], sequences['E1BS'][1], str(sequences['genome'][sequences['E1BS'][0] - 1:sequences['E1BS'][1]-1])]

        return sequences

def export_to_csv(annotations):#Exports all information to a file
    out_file = csv.writer(open('PuMA_out_test.csv','a'))

    for value in annotations:
        if value == 'genome':
            out_file.writerows([[annotations['accession'],'CG',"",annotations[value]]])
        elif value == 'URR':
            try:
                out_file.writerows([[annotations['accession'], value, 'join(' + str(annotations[value][0]) + '..' + str(annotations[value][1]) + '+' + str(annotations[value][2]) + ".." + str(annotations[value][3]) +')', annotations[value][4]]])
            except IndexError:
                out_file.writerows([[annotations['accession'], value, str(annotations[value][0]) + '..' + str(annotations[value][1]), annotations[value][2], annotations[value][3]]])
        elif value == 'E1BS':
            try:
                out_file.writerows([[annotations['accession'], value, 'join(' + str(annotations[value][0]) + '..' + str(annotations[value][1]) + '+' + str(annotations[value][2]) + ".." + str(annotations[value][3]) + ')', annotations[value][4]]])
            except IndexError:
                out_file.writerows([[annotations['accession'], value, str(annotations[value][0]) + '..' + str(annotations[value][1]), annotations[value][2], annotations[value][3]]])
        elif value == 'name':
            pass

        elif value == 'accession':
            pass
        elif value == 'E2BS':
            for i in range(0,len(annotations[value]),1):
                out_file.writerows([[annotations['accession'], value, str(annotations[value][i]) + '..' + str(annotations[value][i] + 11), annotations['genome'][annotations[value][i]-1:annotations[value][i]+11]]])

        else:
            out_file.writerows([[annotations['accession'], value, str(annotations[value][0]) + '..' + str(annotations[value][1]), annotations[value][2], annotations[value][3]]])



    return
