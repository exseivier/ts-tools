#!/usr/bin/env python

"""
OPEN infile FOR READ
OPEN outfile FOR WRITE
DEFINE dictionary SET None
DEFINE header SET None
DEFINE speciesID SET 0001
FOR line in file:
				{	IF line m= header:
									{	SET header line
										SET dictionary[line] ""	}
					ELSE:
									{	CONCAT dictionary[header] line	}	}

DEFINE key, values;
SET key, value FOR dictionary.items:
	WRITE key\tspeciesID\tvalue\n outfile

"""

def load_Ome(filename):
	"""(STR) -> HASH[STR:STR]
	Requires a fasta file name of the reference genome sequence
	It takes the chromosomes of a genome and their sequences
	and stores them into a HASH variable "HASH[header:sequence]"
	"""
	# Counting total headers
	IN = open(filename, "r")
	total = 0
	chrs = {}
	for line in IN:
		line = line.strip()
		try:
			if line[0] == ">":
				pass
		except Exception:
			continue
		if line[0] == ">":
			header = line
			total += 1
		else:
			try:
				chrs[header] += 1
			except Exception:
				chrs[header] = 1
	IN.close()
	# Declaring and Assigning GENOME HASH variable
	OME = {}
	counter = 0
	IN = open(filename, "r")
	for line in IN:
		line = line.strip()
		try:
			if line[0] == ">":
				pass
		except Exception:
			continue
		if line[0] == ">":
			counter += 1
			print "\rParsing %s [%f" % (line, (float(counter)/float(total))*100) + " %]"
			header = line
		else:
			try:
				OME[header] += line
			except Exception:
				OME[header] = line
	IN.close()
	# Dumping HASH variable to a file
	if len(OME) != 0:
		# Returning HASH
		return OME
	else:
		print "'Ome' is empty!"
		exit(0)

from sys import argv, exit, path
path.append('/Users/javier/scripts/ban-4.0/utils/')
from Utils import load_Ome, dna_2_rna
infile = argv[1]
speciesID = "0001"
sequences = load_Ome(infile)
outfile = "%s.ints" % ".".join(infile.split(".")[:-1])
OUT = open(outfile, "w+")
for key, value in sequences.iteritems():
	OUT.write("%s\t%s\t%s\n" % (key.split("|")[0][1:], speciesID, dna_2_rna(value)))
print "OK..."
