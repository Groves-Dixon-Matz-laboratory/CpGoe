#!/usr/bin/env python
##getCpGoe.py
##written 6/26/14 by Groves Dixon
ProgramName = 'getCpGoe.py'
LastUpdated = '10/16/14'
By = 'Groves Dixon'
VersionNumber = '1.0'
print "\nRunning Program {}...".format(ProgramName)

VersionString = '{} version {} Last Updated {} by {}'.format(ProgramName, VersionNumber, LastUpdated, By)

Description = '''
Description:
This program reads though a set of nucleotide sequences and returns a table with CpGoe values
and other relevent data.
'''

AdditionalProgramInfo = '''
Additional Program Information:

'''

##Import Modules 

import time
import argparse
from sys import argv
from sys import exit
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
Start_time = time.time() ##keeps track of how long the script takes to run

##Set Up Argument Parsing
parser = argparse.ArgumentParser(description=Description, epilog=AdditionalProgramInfo) ##create argument parser that will automatically return help texts from global variables above
parser.add_argument('-i', required = False, dest = 'input', help = 'The the name of the file with the gene body sequences')
parser.add_argument('-sub', required = False, default = 'False', dest = 'subset', help = 'optional integer argument to take just the first part of genes. Designates how much to take.')
parser.add_argument('-output', '-o', required = True, dest = 'out', help = 'The desired name for the output file')
parser.add_argument('-cds', required = False, dest = 'cds', help = "Tell whether this is a CDS file or not, these are Misha's protein coding nucleotide sequences. Assign as 'y' for yes or leave blank")
args = parser.parse_args()

#Assign Arguments
inputFile = args.input
outfileName = args.out
subset = args.subset
CDS = args.cds


def read_genome(inputFile, outfileName, subset):
    with open(outfileName, 'w') as out:
        header = "EST\tC\tG\tT\tCpG\tGpC\tTpG\tlength"
        out.write(header)
        geneLengths = {}
        isogroupList = []
        dataDict = {}
        for seqRecord in SeqIO.parse(inputFile, "fasta"):
            if subset != 'False':
                subset = int(subset)
                seqRecord = seqRecord[0:subset]
            if CDS == 'y':
                geneName = seqRecord.description.split()[2]
                geneNumber = geneName.split('isogroup')[1]
                geneName = "isogroup=" + geneNumber
                print seqRecord.description
                print geneName
            else:
                geneName = seqRecord.id
            length = len(seqRecord.seq) 
            seq = seqRecord.seq
            seq = seq.upper() ##make sure all nucleotide letters are upper case
            Tcount = seq.count('T')
            Ccount = seq.count('C')
            Gcount = seq.count('G')
            Ncount = seq.count('N')
            CGcount = seq.count('CG')
            GCcount = seq.count('GC')
            TGcount = seq.count('TG')
            length = length - Ncount ##don't include Ns in the length measure for getting CpGoe because they had no opportunity to be either one
            if CDS == 'y':
                try:
                    if geneLengths[geneName] > length:
                        continue ##this way we skip repeats, selecting the longer of the two repeats for an isogroup
                except KeyError:
                    geneLengths[geneName] = length
            dataList = [geneName, Ccount, Gcount, Tcount, CGcount, GCcount, TGcount, length]
            stringList = []
            for i in dataList:
                stringList.append(str(i))
            dataString = '\t'.join(stringList)
            dataDict[geneName] = dataString #this will replace the shorter one for a duplicate if it was recorded, so we always keep the longer gene
        for gene in dataDict.keys():
            dataString = dataDict[gene]
            out.write('\n' + dataString)
            


read_genome(inputFile, outfileName, subset)
Time = time.time() - Start_time
print('\nTime took to run: {}'.format(Time))


