#!/usr/bin/python3

'''
    Author: Jose G. Macia-Vicente
    Date created: 2023-05-03
    Date last modified: 2023-05-03
    Function: pair-end assemble sanger sequences. Requires system-wide access 
      to PEAR (https://cme.h-its.org/exelixis/web/software/pear/).
'''

from Bio import SeqIO
import sys
import argparse
import os


# input arguments parameters
parser = argparse.ArgumentParser(
    description='Pair-end assemble Sanger sequences provided as ab1 files.')

parser.add_argument(
    "fwd",
    help = 'input fwd sequence')

parser.add_argument(
    "rev",
    help = 'input rev sequence')

parser.add_argument(
    "-o", help = 'prefix for output files',
    type = str,
    default = "output")

parser.add_argument(
    "-c", help = 'clean-up environment?',
    action=argparse.BooleanOptionalAction)

# parse arguments
args = parser.parse_args()
fwd = args.fwd
rev = args.rev
fq_fwd = fwd.replace(".ab1", ".fastq")
fq_rev = rev.replace(".ab1", ".fastq")

# convert ab1 files to fastq
SeqIO.convert(fwd, "abi-trim", fq_fwd, "fastq")
SeqIO.convert(rev, "abi-trim", fq_rev, "fastq")

# concatenate fwd and reverse reads
os.system("pear -f "+fq_fwd+" -r "+fq_rev+" -o "+args.o)

# save final output to fasta
os.system("sed -n '1~4s/^@/>/p;2~4p' "+args.o+".assembled.fastq > " \
  +args.o+".fasta")

# change fasta header
os.system("sed -i 's/>.*/>"+args.o+"/' " +args.o+".fasta")

# clean-up environment
if args.c:
    os.system("rm "+args.o+"*.fastq "+fq_fwd+" "+fq_rev)

# end
