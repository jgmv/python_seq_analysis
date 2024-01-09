#!/usr/bin/python3

'''
    Author: Jose G. Macia-Vicente
    Date created: 2024-01-09
    Date last modified: 2024-01-09
    Function: Extract and quality-trim sequences from ABI Sanger sequencing
	files. Includes the option for tabular format to include in custom database.
'''

from Bio import SeqIO
from glob import glob
import argparse
import numpy as np
import os
import sys

# input arguments parameters
parser = argparse.ArgumentParser(
    description='Pair-end assemble Sanger sequences provided as ab1 files.')

parser.add_argument(
	"-i",
	type = str,
	help = 'input idividual file name',
	default = None)

parser.add_argument(
	"-f",
	type = str,
	help = 'input folder name',
	default = None)

parser.add_argument(
	"-o",
	type = str,
	help = 'name of output file',
	default = "output.fasta")

parser.add_argument(
	"-l",
	type = int,
	help = 'minimum sequence lenght',
	default = 100)
	
parser.add_argument(
	"-u",
	help = 'upper-case sequence',
	action=argparse.BooleanOptionalAction)

parser.add_argument(
	"-db",
	help = 'format for database input',
	action=argparse.BooleanOptionalAction)

# parse arguments
args = parser.parse_args()

if args.i == None and args.f == None:
	sys.exit("Please provide an input file or folder.")

if args.i != None and args.f == None:
	path = [args.i]
	print("Reading file '"+args.i+"'")
else:
	path = os.listdir(args.f)
	path = sorted(glob(os.path.join(args.f, "*.ab1")), key=os.path.getctime)
	print("Reading all files in folder '"+args.f+"'")

seqOut = open(args.o, "w")

# extract sequences and write to file
count = 0
if args.db:
	seqOut.write("accession\tseq_length\tseq_qmean\tseq_q20perc\tsequence\n")

for filename in path:
	seq_record = SeqIO.read(filename, "abi-trim")
	seq_len = len(seq_record.seq)
	if args.u:
		seq = seq_record.seq.upper()
	else:
		seq = seq_record.seq.lower()
	if seq_len < args.l:
		print(seq_record.id+" too short. Excluded.\n")
		pass
	if args.db:
		seq_qual = np.array(seq_record.letter_annotations["phred_quality"])
		seq_q20 = round(100 * (seq_qual > 20).sum() / len(seq_record.seq), 1)
		seq_qual = round(seq_qual.mean(), 1)
		seqOut.write(seq_record.id+"\t"+str(seq_len)+"\t"+str(seq_qual)+"\t"+str(seq_q20)+"\t"+str(seq)+"\n")
		count +=1
	else:
		seqOut.write(">"+seq_record.id+" "+str(seq_len)+"\n"+str(seq)+"\n")
		count +=1

print(str(count)+" sequences extracted")
seqOut.close()

# end