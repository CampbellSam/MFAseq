#!/usr/bin/env python

""" Simple Python script to generate wiggle files for the MFAseq data.
    Takes two bam files as input and generates a wiggle file output.
"""

from __future__ import division

__author__ = "Nick Dickens, Samantha Campbell"
__copyright__ = "Copyright 2016, Nicholas J. Dickens"
__credits__ = ["Nick Dickens", "Samantha Campbell"]
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Samantha Campbell"
__email__ = "s.campbell.1@research.gla.ac.uk"


import pysam
import sys
import os.path
import argparse


def get_args():
    parser = argparse.ArgumentParser(description='Generates a wiggle file from two bam files with the ration of ')
    parser.add_argument('--file1', required=True, help='the first DNAseq file (the ratio is first/second)')
    parser.add_argument('--file2', required=True, help='the second DNAseq file (the ratio is first/second)')
    parser.add_argument('--window', required=False, help='the window size in bases (tile) for binning the depth data, defaults to 2500')
    parser.add_argument('--wig', required=False, help='output wiggle file name, defaults to STDOUT')
    return parser.parse_args()

def wigOut (thisLine, fileHandle):
    fileHandle.write(str(thisLine) + "\n")
    return

def gbrowseHint ():
    thisText = """
Here are edits that you need to make to the Gbrowse config for these tracks,
change everything below the key = line in the existing config.
-------- 8< ----- 8< -------

autoscale = chromosome
description =

glyph           = wiggle_xyplot
graph_type      = line
height          = 100
color           = black
bgcolor         = mediumblue
fgcolour        = mediumblue
linewidth       = 2
max_score       = 2
min_score       = 0

-------- 8< ----- 8< -------
"""
    return thisText

def sortAndIndex (samfile):
    sys.stderr.write("WARNING: file %s was not indexed, sorting and indexing now..." % samfile)
    pysam.sort(samfile, samfile + '.srt')
    pysam.index(samfile + '.srt.bam')
    sys.stderr.write('done\n')
    return samfile + '.srt.bam'

# watch for low values/trap division by zero errors
# will set the ratio to 1000 if G2 is zero
def checkCount(readCount, conversionFactor):
    if readCount < 1:
        readCount = 0.001
    return readCount * conversionFactor



def main():
    args = get_args()

    fileE = args.file1
    fileG = args.file2

    windowSize = 2500

    if args.window:
        windowSize = int(args.window)


    fileHandle = sys.stdout

    if args.wig:
        fileHandle = open(args.wig,"w")

    # check the file is indexed
    if not os.path.exists(fileE + ".bai"):
        fileE = sortAndIndex(fileE)

    if not os.path.exists(fileG + ".bai"):
        fileG = sortAndIndex(samfileG)

    # connect to the bam files
    try:
        samfileE = pysam.Samfile(fileE, "rb" )
        samfileG = pysam.Samfile(fileG, "rb" )
    except:
        e = sys.exc_info()[0]
        sys.stderr.write('ERROR: %s\n' % e)
        sys.exit(1)

    # check that the references sequences in both files are the same list
    # and each chromosome is the same length
    if samfileE.nreferences != samfileG.nreferences:
        sys.exit('ERROR: There is a different number of sequences in the two file headers!\n')

    for chromE in samfileE.references:
        index=samfileE.references.index(chromE)
        chromG=samfileG.references[index]
        if chromE != chromG:
            sys.exit('ERROR: There is a mismatch in the names of the reference sequences (or their order)!\n')


    # downsize samfileG values to match total for samfileE
    try:
        gConversionFactor=samfileE.mapped/samfileG.mapped
    except ZeroDivisionError as detail:
        sys.stderr.write ('ERROR: there is a problem calculating gConversion factor:%s\n' % detail)
        sys.exit(1)


    # iterate through each chromosome
    for chromosome in samfileE.references:
        #get the length of the chromosome
        #print "fixedStep  chrom=%s  start=1  step=%d  span=%d" % (chromosome, windowSize, windowSize)
        headerLine = "fixedStep  chrom=%s  start=1  step=%d  span=%d" % (chromosome, windowSize, windowSize-1)
        wigOut(headerLine, fileHandle)
        chromosomeLength=samfileE.lengths[samfileE.references.index(chromosome)]
        #read along the length of each chromosome
        for windowStart in range(0, chromosomeLength-(windowSize), windowSize) :
            #count reads in E
            readsE = samfileE.count(chromosome, windowStart, windowStart+windowSize)
            readsE= checkCount(readsE , 1)
            #count reads in G
            readsG = samfileG.count(chromosome, windowStart, windowStart+windowSize)
            readsG = checkCount(readsG,gConversionFactor)
            #generate ratio using corrected counts
            try:
                ratio = readsE/readsG
            except ZeroDivisionError:
                sys.stderr.write('WARNING: %s:%d-%d had a strange ratio!' % (chromosome, windowStart, windowStart+windowSize))
                ratio = readsE/0.001
                sys.exit(1)
            #print "%4.3f" % ratio
            thisLine = "%4.3f" % ratio
            wigOut(thisLine, fileHandle)
            #print "%s %d %f" % (chromosome, windowStart, ratio)

    samfileE.close()
    samfileG.close()

    if args.wig:
        fileHandle.close()
        sys.stderr.write(gbrowseHint())

    '''
    The edits to make to the Gbrowse formatting for these tracks, change everything below the key = line
    -------- 8< ----- 8< -------
    autoscale = chromosome
    description =

    glyph           = wiggle_xyplot
    graph_type      = line
    height          = 100
    color           = black
    bgcolor         = mediumblue
    fgcolour        = mediumblue
    linewidth       = 2
    max_score       = 2
    min_score       = 0
    -------- 8< ----- 8< -------
    '''

if __name__ == '__main__':
    main()