#!/usr/bin/env bash

usage="$0 <earlyS.bam> <lateS.bam> > <outfile.wig>"
eS=$1
g2=$2

if [ -z $eS ]; then
echo "You must specify an early-S bam file! $usage"
exit
fi

if [ -z $g2 ]; then
echo "You must specify an early-S bam file! $usage"
exit
fi



./mfaseq.py --file1 $eS --file2 $g2