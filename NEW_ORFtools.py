#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Authors: Koenraad Van Doorslaer, Josh Pace
University of Arizona
"""
from Bio import SeqIO, GenBank, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete2 import Tree
from collections import Counter
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from mechanize import Browser
import os, glob, re, csv, time, operator

def TransORF(seq, trans_table, min_protein_length):   #this function will translate the genome in all 3 frames and split the ORF on the stop codon *
    ORFs = {}
    seq_len = len(seq)
    for frame in range(3):#Accounting for all 3 different frames
        trans = str(seq[frame:].translate(trans_table))
        trans_len = len(trans)
        aa_start = 0
        aa_end = 0
        while aa_start < trans_len:
            aa_end = trans.find("*", aa_start)# returns the lowest index where "*" is found in aa_start because "*" represents a stop codon 
            if aa_end == -1:# .find returns -1 if no "*" is found in aa_start
               aa_end = trans_len 
            if aa_end-aa_start >= min_protein_length:#Only finding proteins over a certain length
                start = frame+aa_start*3# Multiplying by 3 to get the start position of the nucleotides
                ORFs[trans[aa_start:aa_end]] = start
            aa_start = aa_end+1 
    return ORFs
