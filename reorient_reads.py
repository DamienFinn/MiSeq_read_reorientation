# coding: utf-8
#reorient_reads.py
#Written by: Damien Finn, 20.11.20

#The purpose of this code is to sort Illumina MiSeq forward and reverse reads that have their orientations mixed up.
#For example, if a Forward reads file has some sequences starting with the forward Univ 16S primer, and others that 
#start with the reverse Univ 16S primer.

#Note that this script is dependent on having the following Python libraries installed: 'regex', 'gzip'

#To use this script, use the following commands:
#python reorient_reads.py -f <path to forward reads fastq.gz input file> -r <path to reverse reads fastq.gz file> 
#	-fp <Your forward primer> -rp <Your reverse primer> -o <path to output files>

import os.path
import sys
import regex
import gzip

forwardreads_file = sys.argv[sys.argv.index('-f')+1]
reversereads_file = sys.argv[sys.argv.index('-r')+1]
forwardprimer = sys.argv[sys.argv.index('-fp')+1]
reverseprimer = sys.argv[sys.argv.index('-rp')+1]

mycwd = os.getcwd()

#Define a function for chopping up the fastq inputs into lists of 4 elements (SampleID, sequence, '+', quality score)
def list_slice(S, n):
	for i in range(0, len(S), n):
		yield S[i: i + n]

#Create a list of degenerate primers

iupac = {'[':'[', ']':']','A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': '[AG]', 'Y': '[CT]',
            'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
            'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}

forward_primer = ''.join([iupac[symbol] for symbol in forwardprimer])
reverse_primer = ''.join([iupac[symbol] for symbol in reverseprimer])


#------------Importing reads from forward and reverse files------------#

myforwardlist = []
myreverselist = []

with gzip.open(forwardreads_file, 'rt') as unzipforward:
	for x in unzipforward:
		myforwardlist.append(x)
unzipforward.close()

choppedforward = list(list_slice(myforwardlist, 4))

print("Forward reads imported")

with gzip.open(reversereads_file, 'rt') as unzipreverse:
	for x in unzipreverse:
		myreverselist.append(x)
unzipreverse.close()

choppedreverse = list(list_slice(myreverselist, 4))

print("Reverse reads imported")
print("Sum of read pairs imported:" + '\t' + str(len(choppedreverse)))

myforwardlist = []
myreverselist = []

#Note that we have three levels of list: 1) all 'choppedreverse'; 2) the four elements for each sequence; 3) each element 
#Each element needs to be a list because there are commas in the quality scores, and python will break them up if it's not a list


#------------Sorting reads into the correct orientations---------------#

outputforward = []
outputreverse = []

count = 0

for x,y in zip(choppedforward,choppedreverse):
	w = x[1]
	query = w[0:60]
	querymatch = regex.search('('+forward_primer+'){e<=3}', query)
	if querymatch:
		count += 1
		outputforward.append(x)
		outputreverse.append(y)

for x,y in zip(choppedforward,choppedreverse):
	w = x[1]
	query = w[0:60]
	querymatch = regex.search('('+reverse_primer+'){e<=3}', query)
	if querymatch:
		count += 1
		outputforward.append(y)
		outputreverse.append(x)

print("Sum of reads passing primer checking and reorientation:" + '\t' + str(count))		

#-------------------Output the corrected reads-----------------------#		

with gzip.open(os.path.join(mycwd, "forward.fastq.gz"), "wt") as forwardoutput:
	for x in outputforward:
		for y in x:
			forwardoutput.write(y)
		
forwardoutput.close()

with gzip.open(os.path.join(mycwd, "reverse.fastq.gz"), "wt") as reverseoutput:
	for x in outputreverse:
		for y in x:
			reverseoutput.write(y)
		
reverseoutput.close()
