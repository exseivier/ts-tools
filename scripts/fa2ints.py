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
