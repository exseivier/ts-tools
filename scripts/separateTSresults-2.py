#!/usr/bin/env python

from sys import argv

infile = argv[1]

IN = open(infile, "r")
IN.readline()	# Wastes the first line
dictionary = {}
for line in IN:
	line = line.strip()
	fields = line.split("\t")
	gene_id		=	fields[0]	# Gene name
	mirbase_id	=	fields[2]	# microRNA ID, (i.e. miR-124)
	site_type	=	fields[3]	# site type (i.e. 6mer, 7mer or 8mer)
	cspp		=	float(fields[27])	# Context score
	mirna_seed	=	fields[35]	# microRNA seed (i.e. "UUAGGGC")
	if mirna_seed in dictionary:
		if gene_id in dictionary[mirna_seed]:
			dictionary[mirna_seed][gene_id] += cspp
		else:
			dictionary[mirna_seed][gene_id] = cspp
	else:
		dictionary[mirna_seed] = {gene_id : cspp}

for seed, value in dictionary.iteritems():
	OUT = open("%s.data" % seed, "w+")
	for gene_name, cspp in value.iteritems():
		OUT.write("%s\t%f\n" % (gene_name, cspp))
	OUT.close()

print "OK..."
