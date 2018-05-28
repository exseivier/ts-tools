#!/usr/bin/env python

from sys import argv, exit, path
#path.append('/Users/javier/scripts/ban-4.0/utils') Commented because Utils function was copied and pasted in this script
#from Utils import load_Ome

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

def miR_2_TS(infile):
	"""(STR) -> VOID (Writes to file)

	"""
	mirnas = load_Ome(infile)
	outfile = "%s.ints" % ".".join(infile.split(".")[:-1])
	OUT = open(outfile, "w+")
	dictionary = {}
	for mirna, seq in mirnas.iteritems():
		if seq[1:8] in dictionary:
			continue
		else:
			dictionary[seq[1:8]] = seq[1:8]
	
	for key, seed in dictionary.iteritems():
		OUT.write("%s\t%s\t0001\n" % (key, seed))
	print "OK..."

def miR_2_TSc(infile):
	"""(STR) -> VOID (Writes to file)

	"""
	mirnas = load_Ome(infile)
	outfile = "%s.intsc++" % ".".join(infile.split(".")[:-1])
	OUT = open(outfile, "w+")
	dictionary = {}
	for mirna, seq in mirnas.iteritems():
		if seq[1:8] in dictionary:
			dictionary[seq[1:8]].update({mirna : seq})
		else:
			dictionary[seq[1:8]] = {mirna : seq}
	
	for seed, miR in dictionary.iteritems():
		for mid, mature in miR.iteritems():
			OUT.write("%s\t0001\t%s\t%s\n" % (seed, mid.split(" ")[0], mature))
	
	print "OK..."

def main():
	"""
	
	"""
	infile = argv[1]
	miR_2_TS(infile)
	miR_2_TSc(infile)

if __name__ == "__main__": main()
	
