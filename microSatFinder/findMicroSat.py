#!/usr/bin/env python3
# Name: Immaad Mir (imir)
# Group Members: None

"""
Description:
	Program to identify and locate Microsatilitte DNA.
	input: Fasta file.
	output: User given output file or stdout of microsatilitte.
	-TO RUN WITH DEFAULT PARAMS: python3 findMicroSat.py -infile <infile> 
	-If needed to run on hummingbird cluster with SLURM workload manager, 
	 edit test.sh as needed and run with command: sbatch test.py.
"""

import sys
import argparse
import itertools
from Bio import SeqIO
import gzip
import pandas as pd

def Arguments():
	"""
	Retirive Command Line arguments:
		-infile or -i: input file
		-outfile or -o: output file
		-minMonomerLen or -min: Minimum length of repeat monomer
		-maxMonomerLen or -max: Maximum length of repeat monomer
		-minTotLen or minTot: Minimum amount of repeats to be considered a microsatelite
	"""
	parser = argparse.ArgumentParser(
			description = 'Program prolog - Program to find Microsatilitte repeats.', 
			epilog = 'Program epilog- please enter custom arguments to get wanted results',
			add_help = True, #default is True
			prefix_chars = '-',
			usage = 'python3 findMicroSat.py -infile <infile> -option1 <...> -option2 <...>')#how to run
	parser.add_argument('-infile', '-i', required=True, help='Onput file name')
	parser.add_argument('-outfile', '-o', required=False, type=argparse.FileType('w'), default=sys.stdout, help='Output file name. Default = stdout.')
	parser.add_argument('-minMonomerLen', '-min', required=False, type=int, default=1, help='Minimum length of repeate monomer')
	parser.add_argument('-maxMonomerLen', '-max', required=False, type=int, default=6, help='Maximum length of repeate monomer')
	parser.add_argument('-minTotLen', '-minTot', required=False, type=int, default=40, help='Minimum length to be considered a microsatelite' )
	return parser.parse_args()

def reverseComplement(seq):
	"""
	Returns the reverse complement of the single DNA strand given.
	Useful for when identifying microsatilittes that may be on the
	input sequence complementry strand.
	"""
	comp = ''
	for base in seq[::-1]:#revese string
		#add compliment bases
		if base is 'A':
			comp += 'T'
		if base is 'T':
			comp += 'A'
		if base is 'G':
			comp += 'C'
		if base is 'C':
			comp += 'G'
	return comp

def expandRepeatMonomers(repeat, length):
	"""
	Helper function for simulating repeat monomers.
	Expand the repeat monomer.
	"""
	expandedString = ''
	iterator = 0
	while len(expandedString) < length: #add bases
		expandedString += repeat[iterator]
		iterator += 1
		if iterator >= len(repeat): #end of repeat. Start increasing from begining of repeat
			iterator = 0
	return expandedString

def overlap(repeat):
	"""
	Helper function for simulating repeat monomers.
	Get the overlapping regions.
	"""
	cycles = []
	for iterator in range(len(repeat)):
		cycles.append(repeat[iterator:] + repeat[:iterator])
	return cycles


def simulateRepeatMonomers(minLen, maxLen):
	"""
	Create monomer List. Using itertools.product() we are
	able to find all possible monomers of A,C,G,T also known as
	the Cartisan Product. Helper function expandRepeatMonomers()
	helps us make sure there are no possible duplications of 
	monomers.
	"""
	   
	nucleotides = ['A', 'C', 'G', 'T']
	simulatedRep = []
	expanded = set()
	repeats = set()
	for length in range(minLen, maxLen+1):
		for cartesian_product in itertools.product(nucleotides, repeat=length):
			repeat = ''.join(cartesian_product)
			repeatRevComp = reverseComplement(repeat)
			expandedString = expandRepeatMonomers(repeat, maxLen)
			if expandedString not in expanded:
				RepCycles = overlap(repeat)
				for cycle in RepCycles:
					strand = "+"
					string = expandRepeatMonomers(cycle, maxLen)
					expanded.add(string)
					if cycle not in repeats:
						repeats.add(cycle)
						simulatedRep.append('\t'.join([cycle, repeat, str(len(cycle)), strand]))
				if repeatRevComp != repeat:
					RepCycles = overlap(repeatRevComp)
					for cycle in RepCycles:
						strand = '-'
						string = expandRepeatMonomers(cycle, maxLen)
						expanded.add(string)
						if cycle not in repeats:
							repeats.add(cycle)
							simulatedRep.append('\t'.join([cycle, repeat, str(len(cycle)), strand]))
	return simulatedRep


def buildDict(monomers, **kwargs):
	"""
	Creats dictionary from monomer list containing helper information for output:
	Monomer, length of monomer, which strand.
	"""
	repeats = {}
	cutoff = kwargs.get('minTotLen', None)
	if cutoff is not None:
		lenMonomers = []
		for line in monomers:
			infoDict = {}
			info = line.strip()
			info = line.split('\t')
			monomer = info[0]
			lenMonomer = int(info[2])
			if lenMonomer not in lenMonomers:
				lenMonomers.append(lenMonomer)
			iterator = 0
			lenHolder = lenMonomer
			while iterator < lenMonomer and lenHolder < cutoff:
				monomer += monomer[iterator]
				iterator += 1
				if iterator >= lenMonomer:
					iterator = 0
				lenHolder = len(monomer)
			infoDict['LenMonomer'] = lenMonomer
			infoDict['strand'] = info[3]
			repeats[monomer] = infoDict
		repeats['rep_lengths'] = [cutoff]
	return repeats


def writeToFile(headSeq, repeatDict, monomers, outFILE):
	"""
	Find microsatilites in sequences by looking for 
	matches in the simulated monomer list.
	Write all found microsatellite dna to output file.
	Each Microsatellite is outputted in the format:
		sequenceID : repeat monomer : start Index : end index : length : number of monomers : strand
	"""
	simulatedLengths = repeatDict['rep_lengths']
	inputSeq = str(headSeq.seq).upper()
	inputSeqLen = len(inputSeq)
	for cutoff in simulatedLengths:
		sub = cutoff - 1
		startCnt = 0
		stopCnt = startCnt + simulatedLengths[-1]
		while stopCnt <= inputSeqLen:
			stopCnt = startCnt + cutoff
			seq = inputSeq[startCnt:stopCnt]
			if seq in monomers: #match found
				match = True
				monomerLength = repeatDict[seq]['LenMonomer']
				offset = cutoff % monomerLength
				repeatSeq = inputSeq[startCnt+offset:startCnt+offset+monomerLength]
				iterator = 0
				while match:
					j = stopCnt
					if inputSeq[j] == repeatSeq[iterator]:
						stopCnt += 1
						iterator += 1
						if iterator >= monomerLength:
							iterator = 0
					else: #print to output file using Pandas
						match = False
						lenMatch = stopCnt - startCnt
						units = int(lenMatch/monomerLength)
						df = pd.DataFrame({'Length':[lenMatch], 'Stop':[stopCnt], 'Start': [startCnt], 'Sequence ID':[headSeq.id],\
						'#OfMonomer':[units], 'Monomer':[seq[:monomerLength]], 'Strand':[repeatDict[seq]['strand']]}, index = ['*'])
						df = df[['Sequence ID', 'Monomer', 'Start', 'Stop', 'Length', '#OfMonomer', 'Strand']]
						print(df, file=outFILE)
						print(file=outFILE)
						startCnt = stopCnt - sub
			else:
				startCnt += 1

def main():
	"""
	Open input file for reading. If file is zipped w/ gunzip, unzip file.
	Simulate the repeate monomers and build dictionary. For each sequence 
	in fasta file write found microsatillites to output file.
	"""
	#set custom options for Pandas output
	pd.set_option('display.max_rows', 500)
	pd.set_option('display.max_columns', 500)
	pd.set_option('display.width', 150)
	args = Arguments()
	#simulate repeats
	monomers = simulateRepeatMonomers(args.minMonomerLen, args.maxMonomerLen)
	#build dictionary
	infoDict = buildDict(monomers, minTotLen = args.minTotLen)
	infoSet = set(infoDict.keys())
	#open Input file
	if args.infile.endswith('.gz'):
		infile = gzip.open(args.infile, 'rt')
	else:
		infile = open(args.infile, 'r')
	#for eash sequence in fasta search for microsatilites
	records = SeqIO.parse(infile, 'fasta') 
	for record in records:
		if 0 <= len(record.seq):
			writeToFile(record, infoDict, infoSet, args.outfile)
	args.outfile.close()

if __name__ == '__main__':
	main()

