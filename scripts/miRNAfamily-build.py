#!/usr/bin/env python

from sys import argv, exit, path
path.append('/Users/javier/scripts/ban-4.0/utils')
from Utils import load_Ome

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
	
