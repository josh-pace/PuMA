# All functions needed for PuMA
# Koenraad Van Doorslaer, Ken Youens-Clark, Josh Pace
# New blast function and dictionary setup

from distutils.spawn import find_executable
from Bio import SeqIO, GenBank, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from Bio.Align.Applications import MuscleCommandline
import os, glob, re, csv, time, operator, argparse, sys
import warnings
from Bio import BiopythonWarning
import itertools
import shutil


# --------------------------------------------------
#
#This functions finds all the open reading frames and translates them
#
def trans_orf(seq, trans_table, min_protein_length):
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

# --------------------------------------------------
#
# This function uses blast to find the different proteins in
# the translated genome
#

def blast_proteins(genome,min_prot_len,evalue,blast_dir, out_dir):
    protein_start = {}
    protein_seq = {}
    found_proteins = {}

    orfs = trans_orf(genome, 1, min_prot_len)

    if not orfs:
        print('No ORFs, must stop.')
        sys.exit(1)

    orfs_fa = os.path.join(out_dir, 'orfs.fa')
    orfs_fh = open(orfs_fa, 'wt')

    for orf in orfs:
        orfs_fh.write('\n'.join(['>' + str(orfs[orf]), orf, '']))
    orfs_fh.close()


    blast_sub = os.path.join(blast_dir, 'conserved_test_E5.fa')
    blast_out = os.path.join(out_dir, 'blast_results.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)


    #print('BLASTing')
    cmd = blastp(query=orfs_fa,
                 subject=blast_sub,
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


    with open(blast_out) as tab_file:
        for line in csv.reader(tab_file, delimiter="\t"):
            protein_start[line[1]] = int(line[0])

    with open(orfs_fa) as fasta_file:
        for line in fasta_file:
            for num in protein_start:
                try:
                    start = int(line[1:])
                    if start == protein_start[num]:
                        seq = next(fasta_file)
                        protein_seq[num] = seq[:-1]
                except:
                    pass

    for seq in protein_seq:
        for start in protein_start:
            if seq == start:
                if seq == 'L1':
                    M = re.search('M', protein_seq[seq])
                    real_start = protein_start[start] + M.start() + M.start() + M.start()
                    end = protein_start[start] + ((len(protein_seq[seq])+1) * 3)
                    found_proteins['L1'] = []
                    L1_pre = genome[(protein_start[start] + 3 * M.start()):int(end)]
                    splice = '(C|T)(C|T)(A|C|G|T)(C|T)AG(A)TG'
                    spliced = re.search(splice, str(L1_pre))
                    if spliced:
                        start_L1 = int(spliced.start()) + 6
                        if start_L1 % 3 == 0:
                            if start_L1 > 600:
                                L1_post = L1_pre
                                found_proteins['L1'] = [int(start_L1),
                                    int(end),str(L1_post).lower(),Seq(str(
                                        L1_post)).translate()]
                            else:
                                L1_post = L1_pre[start_L1:]
                                found_proteins['L1'] = [int(real_start) + 1 +
                                       int(start_L1),int(end), str(L1_post).lower(),
                                    Seq(str(L1_post)).translate()]
                        else:
                            L1_post = L1_pre
                            found_proteins['L1'] = [int(real_start) + 1, int(end),
                                str(L1_post).lower(),
                            Seq(str(L1_post)).translate()]
                    else:
                        L1_post = L1_pre
                        found_proteins['L1'] = [int(real_start) + 1,
                            int(end), str(L1_post).lower(), Seq(str(L1_post)).translate()]
                else:
                    try:
                        M = re.search('M', protein_seq[seq])
                        real_start = protein_start[start] + M.start() + M.start() + M.start()
                        end = protein_start[start] + ((len(protein_seq[seq])+1) * 3)
                        sequence = str(genome[int(real_start):int(end)]).lower()
                        translated = Seq(sequence).translate()
                        found_proteins[seq] = [int(real_start) + 1, int(end), sequence,
                        translated]
                    except AttributeError:
                        pass




    return found_proteins

# --------------------------------------------------
#
#Finds the E1 binding site in the genome using the URR
#
# --------------------------------------------------

def find_E1BS(genome, URR, URRstart, ID, out_dir):
    genomeLength = len(genome)
    E1BS = {}  # Storing E1BS
    startURR = 0

    tmp = os.path.join(out_dir, "puma_urr.fa")
    with open(tmp, "w") as tempfile:  # Writting URR to a file so FIMO can be used
        tempfile.write('>URR for {}\n'.format(ID))
        tempfile.write(str(URR))
        # print >> tempfile, '>URR for %s' %ID
        # print >> tempfile, URR

    fimo_exe = find_executable('fimo')
    if not fimo_exe:
        print('Please install "fimo" into your $PATH')
        return

    fimo_dir = os.path.join(out_dir, 'E1BS')
    if not os.path.isdir(fimo_dir):
        os.makedirs(fimo_dir)



    fimo_cmd = '{} --oc {} --norc --verbosity 1 --thresh 1.0E-1 --bgfile {} {} {}'
    cline = (fimo_cmd.format(fimo_exe, fimo_dir, 'background_model_E1BS.txt',
                             'meme_E1BS_1motif_18_21.txt', tmp))

    os.system(str(cline))  # Executing FIMO

    fimo_out = os.path.join(fimo_dir, 'fimo.txt')

    if not os.path.isfile(fimo_out):
        print('Failed to create fimo out "{}"'.format(fimo_out))
        return

    for column in csv.reader(open(fimo_out, "rU"),
                             delimiter='\t'):  # Getting nucleotide start positions from FIMO output file
        if column[3] == 'start':
            startListURR = []
        else:
            startListURR.append(column[3])

    startURR = int(startListURR[0])  # Finding the start position of URR
    genomestart = (startURR + URRstart)



    if genomestart > genomeLength:  # Case where the start of E1BS is after the end of the genome
        genomestart = genomestart - genomeLength


    genomestop = genomestart + 19

    if genomestop > genomeLength:  # For printing all info
        genomestop = genomestop - genomeLength
        sequence = str(genome[int(genomestart)-2:] + genome[:genomestop]).lower()
        E1BS['E1BS'] = [int(genomestart), int(genomeLength), 1, int(genomestop),sequence]

    else:
        if genomestart == 1:
            sequence = str(genome[-1]).lower() + str(genome[int(genomestart - 1):int(
                genomestop)]).lower()
            E1BS['E1BS'] = [int(genomestart), int(genomestop), sequence]

        else:
            sequence = str(genome[int(genomestart-2):int(genomestop)]).lower()
            E1BS['E1BS'] = [int(genomestart), int(genomestop),sequence]

    return E1BS
# --------------------------------------------------

# --------------------------------------------------

#
# This function finds the E2BS in a genome using the URR
#
def find_E2BS(genome, URR, URRstart, ID, out_dir):
    genomeLength = len(genome)  # Getting length of genome
    startListURR = []  # Storing the nucleotide start positions in URR of the E2BS
    startListGenome = []  # Storing the nucleotide start positions in genome of the E2BS
    E2BS = {}  # Storing all E2BS

    # Writting URR to a file so FIMO can be used
    tmp = os.path.join(out_dir, "puma_urr.fa")
    with open(tmp, "w") as tempfile:
        tempfile.write('>URR for {}\n'.format(ID))
        tempfile.write(str(URR))

    # Executing FIMO, Using URR
    fimo_exe = find_executable('fimo')
    if not fimo_exe:
        print('Please install "fimo" into your $PATH')
        return

    fimo_dir = os.path.join(out_dir, 'E2BS')
    if not os.path.isdir(fimo_dir):
        os.makedirs(fimo_dir)

    fimo_cmd = '{} --oc {} --norc --verbosity 1 --thresh 1.0E-4 {} {}'
    cline = (fimo_cmd.format(fimo_exe, fimo_dir, 'meme_3000_TOTAL.txt', tmp))

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

    startListURR = list(map(int, set(startListURR)))  # Making the positions an integer
    # value to use as index and only using every unique start position
    startListGenome = list(map(int, startListGenome))  # Making the positions an integer
    # value

    for i in range(0, len(startListURR), 1):  # Finding the positions of E2BS in genome
        genomestart = startListURR[i]
        genomestart = (genomestart + URRstart)
        if genomestart > genomeLength:
            genomestart = genomestart - genomeLength
            startListGenome.append(genomestart-1)
        else:
            startListGenome.append(genomestart-1)

    E2BS['E2BS'] = startListGenome  # Putting values into dictionary

    return E2BS
# --------------------------------------------------
def find_E4(E2, genome):  # Finds E4
    E4 = {}  # Storing E4 information

    trans_E2 = Seq(E2[1:len(E2)]).translate()  # Translates E2 nucleotide sequence
    E4protein_long = str(max(trans_E2.split("*"),key=len))
    #print("E4 Before:{}".format(E4protein_long))
    # Splits sequences on the stop codon, takes longest sequence
    if 'M' in E4protein_long:
        M = re.search('M', E4protein_long)
        #print("M start:{}".format(M.start()))
        if M.start() > 41:
            E4protein = E4protein_long
        else:
            E4protein_tmp = E4protein_long.split('M',1)[1]
            E4protein = 'M' + str(E4protein_tmp)
    else:
        E4protein = E4protein_long
    #print("E4 After:{}".format(E4protein))
    E4_start = re.search(str(E4protein),str(trans_E2)).start()#Finding the start position of E4
    E4_end = re.search(str(E4protein), str(trans_E2)).end()#Finding the end position of E4
    E4_nt = str(E2[(E4_start * 3) + 1:((E4_end + 1) * 3) + 1])  # Getting the E4 nucleotide sequence
    E4_nt_start = re.search(E4_nt, str(genome)).start()  # Finding nuceotide start position of E4
    E4_nt_end = E4_nt_start + len(E4_nt)  # Finding nucleotide end position of E4
    sequence = str(genome[int(E4_nt_start):int(E4_nt_end)]).lower()
    translated = Seq(sequence).translate()
    E4['E4'] = [int(E4_nt_start) + 1, int(E4_nt_end) + 1]  # Storing start and stop
    # dictionary
    return E4
# --------------------------------------------------
#
#Finds E1^E4
#
def find_E1E4(E1_whole,E2_whole,ID,genome,start_E4_nt):
    E1_E4 = {}
    genome = str(genome).lower()
    startListE2 = []
    E2_seq = str(E2_whole[2])

    donor_options = ['aggta','aggtg', 'cggta','agagt']
    stopE1_options = []



    try:
        E1_seq = str(E1_whole[2])
        #print('E1_seq:{}'.format(E1_seq))
        for sites in donor_options:
            if sites in E1_seq:
                stopE1_options.append(re.search(sites, E1_seq).start())# Finding the

        stopE1_options = sorted(stopE1_options) #splice donor site list
        print('stopE1:{}'.format(stopE1_options))
        stopE1 = stopE1_options[0]
        stopE1 = stopE1 + 2 # Accounting for the ag that has to be apart of the sequence
        start_E1_nt = E1_whole[0]
        stop_E1_nt = (stopE1 - 2) + 1 + E1_whole[0]


        if start_E4_nt == False:
            E1_E4['E1^E4'] = False
            return E1_E4
        else:

            whole_E4 = find_E4(E2_seq, genome)
            stop_E4_nt = whole_E4['E4'][1] - 1

            E1_E4_seq = str(genome[start_E1_nt-1:stop_E1_nt]+ genome[
            start_E4_nt:stop_E4_nt])
            E1_E4_trans = Seq(E1_E4_seq).translate()[:-1]
            #start_seq = startListE2[0]# Finding the start position of E4
            # start_E4_nt = (re.search(start_seq[:-4], E2_seq).end() + E4_whole[0])
            # stop_E4_nt = E4_whole[1]



            E1_E4['E1^E4'] = [start_E1_nt,stop_E1_nt,start_E4_nt + 1,stop_E4_nt,E1_E4_seq,
                E1_E4_trans]
            return E1_E4
    except IndexError:
        E1_E4['E1^E4'] = [0, 0, 0, 0, 'Function Crashed','',genome[start_E1_nt:stop_E1_nt]]
        return E1_E4
# --------------------------------------------------
def find_E8E2(E1_whole, E2_whole, ID, genome,startE2_nt):
    E8_E2 = {}
    E1_seq = str(E1_whole[2])
    stopE8List = []
    genome = str(genome).lower()
    donor_options = ['aggtg','aggtt','agaga','aggga']

    if startE2_nt == False:

        E8_E2['E8^E2'] = False
        return E8_E2

    for stop in re.finditer('aggta', E1_seq):
        stopE8List.append(stop.start())
    print('Looking for aggta: {}'.format(stopE8List))
    actualSites = []
    outOfRange = 0
    for stop in stopE8List:
        if stop > 700 or stop < 20:
            outOfRange += 1

    print("StopList after first check:{}".format(stopE8List))

    print("Out of range:{}".format(outOfRange))

    if len(stopE8List) == outOfRange:#or len(stopE8List) == 0
        stopE8List.clear()
        for sites in donor_options:
            print("sites in donor_options:{}".format(sites))
            if sites in E1_seq:
                stopE8List.append(re.search(sites, E1_seq).start())
                print("Site found:{} and position:{}".format(sites,re.search(sites, E1_seq).start()))
        stopE8List = sorted(stopE8List)  # splice donor site list
    for stop in stopE8List:
        #print(stop)
        if stop > 325 and stop < 600:
            actualSites.append(stop)


    actualSites = list(sorted(set(actualSites)))

    print("Stop sites:{}".format(actualSites))

    if len(actualSites) == 0:
        E8_E2['E8^E2'] = False
        return E8_E2

        if startE2_nt != False:
            stopE2_nt = E2_whole[1]

    for position in actualSites:

        #try:
        stopE8 = position#actualSites[0]
        search_seq = E1_seq[:stopE8]
        i = len(search_seq)
        end_seq = i

        while(i>= 0 and i > end_seq - 70):
            #print("search_seq:{}".format(search_seq[i-3:i]))
            if search_seq[i-3:i] == 'atg':
                E8_seq = search_seq[i - 3:]
                startE8 = re.search(E8_seq, E1_seq).start()
                print("search_seq:{}".format(search_seq[i-3:i]))
                break
            else:
                print("in else")
                i = i - 3
                print("search_seq in else:{}".format(search_seq[-3:i]))
        continue

    print("startE8:{}".format(startE8 + E1_whole[0]))
    print("stopE8:{}".format((stopE8 + E1_whole[0])+1))

        #except

    startE8_nt = startE8 + E1_whole[0]
    stopE8_nt = (stopE8 + E1_whole[0]) + 1

    if startE2_nt != False:
        stopE2_nt = E2_whole[1]



    E8_E2_seq = Seq(genome[startE8_nt - 1:stopE8_nt] + genome[
    startE2_nt:stopE2_nt])
    E8_E2_trans = E8_E2_seq.translate()

    print("E8_seq:{}".format(E8_E2_seq))

    E8_E2['E8^E2'] = [startE8_nt, stopE8_nt, startE2_nt + 1, stopE2_nt, E8_E2_seq,
            E8_E2_trans]

    return E8_E2
# --------------------------------------------------
#
# --------------------------------------------------

#New version of find_splice_acceptor
#
def find_splice_acceptor( E2_whole, ID,genome, blast_dir, out_dir):

    E2_trans = str(E2_whole[3])
    E2_seq = E2_whole[2]

  
    blast_subject = os.path.join(blast_dir, 'accession_e2_trans.fa')
    #blast_db = os.path.join(blast_dir,'accession_e2_half.fa')
    blast_out = os.path.join(splice_acceptor_dir, 'blast_result.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)

    query_file = os.path.join(splice_acceptor_dir, 'query.fa')

    #print("ID:{} E2_Seq:{}".format(ID, E2_seq))
    with open(query_file, 'a') as query:
        query.write('>{}\n'.format(ID))
        query.write(E2_trans)

    cmd = blastp(query=query_file,
                 subject=blast_subject,
                 evalue=1e-10,
                 outfmt=6,
                 out=blast_out)

    stdout, stderr = cmd()
    if stdout:
        print("STDOUT = ", stdout)
    if stderr:
        print("STDERR = ", stderr)

    if not os.path.isfile(blast_out) or not os.path.getsize(blast_out):
        print('No BLAST output "{}" (must have failed)'.format(blast_out))
        startE2_nt = 15
        return startE2_nt

    blast_options = []
    with open(blast_out) as blast_file:
        blast_result = csv.reader(blast_file, delimiter='\t')
        for row in blast_result:
            blast_options.append(row[1])
    blast_options = blast_options[0:1]
    splice_sites = []
    aligned_starts = []
    for options in blast_options:
        query = options
        known_E2 = {}
        #print("Query:{}".format(query))
        csv_database = os.path.join(blast_dir, 'all_pave.csv')


        with open(csv_database, 'r') as csvfile:
            read = csv.DictReader(csvfile,('accession', 'gene', 'positions', 'seq'))
            for row in read:
                if row['accession'] == query and row['gene'] == 'E8^E2':
                    splice_acceptor_positions = row['positions']
                if row["accession"] == query and row["gene"] == 'E2':
                    known_E2[query] = str(row['seq']).lower()
                    known_E2_start = str(row["positions"])
                if row['accession'] == query and row['gene'] == 'CG':
                    known_CG = str(row['seq']).lower()

        try:
            splice_start_genome = splice_acceptor_positions.split('+')[1]

            splice_start_genome = int(str(splice_start_genome).split('..')[0])

            known_E2_start = int(known_E2_start.split('..')[0])

            splice_start_known = (splice_start_genome - known_E2_start)

            splice_sites.append(splice_start_known)

        except UnboundLocalError:
            print('Query does not have E1^E4 or E8^E2')
            startE2_nt = False
            return startE2_nt

        unaligned = os.path.join(splice_acceptor_dir, 'unaligned.fa')
        aligned = os.path.join(splice_acceptor_dir, 'aligned.fa')

        if os.path.isfile(unaligned):
            os.remove(unaligned)

        if os.path.isfile(aligned):
            os.remove(aligned)

        for key in known_E2:
            with open(unaligned, 'a') as sequence_file:
                sequence_file.write(">{}\n".format(ID))
                sequence_file.write("{}\n".format(E2_seq))
                sequence_file.write(">{}\n".format(key))
                sequence_file.write("{}\n".format(known_E2[key]))

        cline = MuscleCommandline(input=unaligned, out=aligned, verbose=False)

        stdout, stderr = cline()

        #How to output to stop printing



        # if stdout:
        #     print("STDOUT = ", stdout)
        # if stderr:
        #     print("STDERR = ", stderr)

        align_seq = []
        for aln in AlignIO.read(aligned, 'fasta'):
            align_seq.append(aln.seq)

        unknown_seq = str(align_seq[0]).lower()
        known_seq = str(align_seq[1]).lower()

        j = 0
        aligned_splice_start = 0

        for position in known_seq:
            aligned_splice_start = aligned_splice_start + 1
            if position.lower() in ['a','c','t','g']:
                j = j + 1
                if j == splice_start_known:
                    break

        aligned_starts.append(aligned_splice_start)

    aligned_splice_start = aligned_starts[-1]


    search_seq = unknown_seq[aligned_splice_start:aligned_splice_start + 50].replace('-','')

    startE2_nt = re.search(search_seq,str(genome).lower()).start()


    return startE2_nt


# --------------------------------------------------



# --------------------------------------------------
# def to_gff3(dict, genomelen, out_dir):
#     del dict['genome']
#     del dict['accession']
#     del dict['E1BS']
#     del dict['E2BS']
#     del dict['URR']
#     all = dict['name']
#     name = re.search('\(([^)]+)', all).group(1)
#     dict['name'] = name
#
#     gff3_out = os.path.join(out_dir, '{}.gff3'.format(dict['name']))
#
#     with open(gff3_out, 'a') as out_file:
#         out_file.write("##gff-version 3\n")
#         out_file.write("##sequence-region {} 1 {}\n".format(dict['name'],genomelen))
#
#     # seqname source feature start end score strand frame attribute
#     for protein in dict:
#         if protein == 'name':
#             pass
#         else:
#             with open(gff3_out,'a') as out_file:
#                 out_file.write("\t".join([
#                     dict['name'],
#                     'PuMA',
#                     'CDS',
#                     dict[protein][0],
#                     dict[protein][1],
#                     '.',
#                     '+',
#                     '.',
#                     ';'.join([
#                         'ID=' + protein,
#                         'Note=' + str(dict[protein][0]) + '-' + str(dict[protein][1])])
#                     ])
#                 )
#
#
#     return


# --------------------------------------------------

#
#Output to gff3 file


#NEED TO FIX FORMAT

def to_gff3(dict, genomelen, out_dir):
    del dict['genome']
    del dict['E1BS']
    del dict['E2BS']
    del dict['URR']

    name = dict['accession']
    dict['name'] = name
    gff3_out = os.path.join(out_dir, '{}.gff3'.format(dict['name']))

    with open(gff3_out, 'a') as out_file:
        out_file.write("##gff-version 3\n")
        out_file.write("##sequence-region {} 1 {}\n".format(dict['name'], genomelen))

    for protein in dict:
        if protein == 'name':
            pass
        elif "^" in protein:
            with open(gff3_out, 'a') as out_file:
                out_file.write(
            "{}\tPuMA\tCDS\t{}\t{}\t{}\t{}\t.\t+\t.\tID={};Note=[{}-{} + {}-{}]\n".format(
                 dict['name'], dict[protein][0], dict[protein][1],
                dict[protein][2],dict[protein][3],protein, dict[protein][ 0],
                dict[protein][1],dict[protein][2],dict[protein][3]))
        else:
            with open(gff3_out, 'a') as out_file:
                out_file.write(
            "{}\tPuMA\tCDS\t{}\t{}\t.\t+\t.\tID={};Note=[{}-{}]\n".format(
                 dict['name'], dict[protein][0], dict[protein][1],
                protein, dict[protein][ 0], dict[protein][1]))

    return
# --------------------------------------------------

# --------------------------------------------------

#
#Output of found sequences for verification against PaVE data
#

def to_results(dict):
    del dict['genome']
    del dict['accession']
    del dict['E1BS']
    del dict['E2BS']


    all = dict['name']

    short_name = re.search('\(([^)]+)', all).group(1)

    results_dir = os.path.join('puma_results')

    results = os.path.join(results_dir, 'puma_results_9_25_newBlast.fa')
    for protein in dict:
        if protein == 'name':
            pass

        elif protein == 'URR':
            try:
                if type(dict[protein][3]) == int:
                    with open(results,'a') as out_file:
                        out_file.write(">{}, {}\n".format(dict['name'],protein))
                        out_file.write("{}\n".format(dict[protein][4]))
                else:
                    with open(results,'a') as out_file:
                        out_file.write(">{}, {}\n".format(dict['name'],protein))
                        out_file.write("{}\n".format(dict[protein][2]))
            except IndexError:
                with open(results, 'a') as out_file:
                    out_file.write(">{}, {}\n".format(dict['name'], protein))
                    out_file.write("{}\n".format(dict[protein][2]))

        elif '^' in protein:
            with open(results, 'a') as out_file:
                out_file.write(">{}, {} gene\n".format(dict['name'], protein))
                out_file.write("{}\n".format(dict[protein][4]))



        else:
            with open(results,'a') as out_file:
                out_file.write(">{}, {} gene\n".format(dict['name'],protein))
                out_file.write("{}\n".format(dict[protein][2]))


    return

# --------------------------------------------------
#
#Output to csv file
#
def export_to_csv(annotations,out_dir):

    #csv_out = os.path.join(out_dir, 'puma_results.csv')

    results_dir = os.path.join('puma_results')
    csv_out = os.path.join(results_dir, 'puma_results.csv')


    with open(csv_out,'a') as out:
        out_file = csv.writer(out)

        for value in annotations:
            if value == 'genome':
                out_file.writerows([[annotations['accession'],
                                   'CG',"",str(annotations[value]).lower()]])
            elif value == 'URR':
                try:
                    out_file.writerows([[annotations['accession'],
                                       value, 'join(' + str(annotations[value][0])
                                       + '..' +
                                      str(annotations[value][1]) + '+' +
                                      str(annotations[value][2]) + ".." +
                                      str(annotations[value][3]) +')',
                                      annotations[value][4]]])
                except IndexError:
                    out_file.writerows([[annotations['accession'],
                                       value, str(annotations[value][0]) + '..' +
                                       str(annotations[value][1]),
                                        annotations[value][2]]])
            elif value == 'E1BS':
                try:
                    out_file.writerows([[annotations['accession'],
                                       value, 'join(' + str(annotations[value][0]) + '..'
                                       + str(annotations[value][1]) + '+' +
                                       str(annotations[value][2]) + ".." +
                                       str(annotations[value][3]) + ')',
                                       annotations[value][4]]])
                except IndexError:
                    out_file.writerows([[annotations['accession'],
                                        value, str(annotations[value][0]) + '..' +
                                        str(annotations[value][1]),
                                        annotations[value][2]]])
            elif value == 'E2BS':
                for i in range(0,len(annotations[value]),1):
                                  out_file.writerows([[annotations['accession'], value,
                                   str(annotations[value][i]) + '..' +
                                   str(annotations[value][i] +11),str(annotations['genome']
                                   [annotations[value][i]-1:annotations[value][
                                    i]+11]).lower()]])
            elif value == 'E1^E4':
                out_file.writerows([[annotations['accession'], value,
                                    'join(' + str(annotations[value][0]) + '..' +
                                     str(annotations[value][1]) + '+' +
                                    str(annotations[value][2]) + ".."+
                                    str(annotations[value][3]) + ')',
                                    annotations[value][4],annotations[value][5]]])
            elif value == 'E8^E2':
                out_file.writerows([[annotations['accession'], value,
                                      'join(' + str(annotations[value][0]) + '..' +
                                      str(annotations[value][1]) + '+' +
                                      str(annotations[value][2]) + ".."
                                      + str(annotations[value][3]) + ')'
                                      ,annotations[value][4],annotations[value][5]]])
            elif value == 'name':
                pass
            elif value == 'accession':
                pass
            else:
                out_file.writerows([[annotations['accession'],
                                    value, str(annotations[value][0]) + '..' +
                                    str(annotations[value][1]),
                                     annotations[value][2], annotations[value][3]]])

    return
# --------------------------------------------------
