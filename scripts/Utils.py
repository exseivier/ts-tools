#!/usr/bin/env python

def formatseq(filename, formhead):
	"""

	"""
	if formhead == "style1":
		IN = open(filename, "r")
		OUT = open("%s.out" % filename, "w+")
		init_line = IN.readline().strip()
		items = init_line.split("|")
		OUT.write(">%s|%s\n" % (items[2][1:], items[1][:-2]))
		for line in IN:
			line = line.strip()
			if line[0] == ">":
				items = line.split("|")
				new_line = "%s|%s" % (items[2][1:], items[1][:-2])
				OUT.write("\n>%s\n" % new_line)
			else:
				OUT.write("%s" % line)
		IN.close()
		OUT.close()
	elif formhead == "style2":
		IN = open(filename, "r")
		OUT = open("%s.out" % filename, "w+")
		init_line = IN.readline().strip()
		items = init_line.split("|")
		OUT.write(">%s|%s\n" % (items[3][1:], items[1]))
		for line in IN:
			line = line.strip()
			if line[0] == ">":
				items = line.split("|")
				new_line = "%s|%s" % (items[3][1:], items[1])
				OUT.write("\n>%s\n" % new_line)
			else:
				OUT.write("%s" % line)
		IN.close()
		OUT.close()

	else:
		IN = open(filename, "r")
		OUT = open("%s.out" % filename, "w+")
		OUT.write("%s\n" % IN.readline().strip())
		for line in IN:
			line = line.strip()
			if line[0] == ">":
				OUT.write("\n%s\n" % line)
			else:
				OUT.write("%s" % line)
		IN.close()
		OUT.close()
	return 1

def headsArray(filename):
	"""

	"""
	heads = []
	IN = open(filename, "r")
	for line in IN:
		line = line.strip()
		if line[0] == ">":
			heads.append(line[1:])
	IN.close()
	return heads

def geneModel(gff_file):
	"""(STR) -> HASH[(STR | INT) & HASH[(STR | INT)]]

	"""
	# Open gff file
	GFF_IN = open(gff_file, "r")
	
	# Building the gene_models data structure v1.0
	gene_models = {}

	for entry in GFF_IN:
		entry = entry.strip()
		if entry[0] == "#":
			continue
		DATA = entry.split("\t")
		if DATA[2] == "gene":
			geneID = DATA[8].split(";")[0].split("=")[1]
			geneName = DATA[8].split(";")[1].split("=")[1]
			# Building the first level
			gene_models[geneID] = 	{"name"	:	geneName,
									"start"	:	int(DATA[3]),
									"end"	:	int(DATA[4]),
									"ori"	:	DATA[6],
									"locus"	:	DATA[0]
									}

		else:
			try:
				# Building the second level
				gene_models[geneID][DATA[2]].append({"start":	int(DATA[3]),
												"end"		:	int(DATA[4]),
												"ID"		:	DATA[8].split(";")[0].split("=")[1],
												"parent"	:	DATA[8].split(";")[1].split("=")[1]	# BUG: Parent is not in index 1
												})

			except Exception:
				gene_models[geneID][DATA[2]] = [{"start"	:	int(DATA[3]),
												"end"		:	int(DATA[4]),
												"ID"		:	DATA[8].split(";")[0].split("=")[1],
												"parent"	:	DATA[8].split(";")[1].split("=")[1] # BUG: Parent is not in index 1
												}]
	GFF_IN.close()
	return gene_models

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

def rev_comp(adn):
	"""(STR) -> (STR)
	Returns the reverse complement of a DNA string
	"""
	rc_adn = ""
	nt_dict = {
		"A":"T",
		"C":"G",
		"G":"C",
		"T":"A",
		"N":"N"
		}
	for nt in adn:
		rc_adn = nt_dict[nt.upper()] + rc_adn
	return rc_adn

def rev_trans(rna):
	"""(STR) -> (STR)
	Returns the DNA string from a RNA strings: Replaces the (U)racil for (T)imine.
	"""
	dna = ""
	for nt in rna:
		if nt.upper() == "U":
			dna += "T"
		else:
			dna += nt.upper()
	return dna

def count_matches(pattern, sequence):
	"""(STR, STR) -> INT
	Counts the number of times 'pattern' matches inside 'sequence'
	"""
	counts = 0
	for i in xrange(0, len(sequence)-len(pattern)):
		if pattern.upper() == sequence[i : i + len(pattern)].upper():
			counts += 1
		else:
			pass
	return counts

def PDStruct(parentDB):
	"""(STR) -> HASH[ ID : ARRAY[STR] ]
	Creates a HASh data structure with data from parentDB.txt.
	The first column is the ID or feature name, and the data
	of the next columns will be stored in an ARRAY structure
	inside the HASH. And those data stored in the ARRAY will be
	used to replace the current header. The ID must be the same
	of the current header
	"""
	parentDic = {}
	IN = open(parentDB, "r")
	for line in IN:
		line = line.strip()
		fields = line.split("\t")
		for i in xrange(1, len(fields)):
			try:
				parentDic[fields[0]].append(fields[i])
			except Exception:
				parentDic[fields[0]] = [fields[i]]
	
	return parentDic


def replace_header(filename, output, PDStruct):
	"""(STR, STR, PDStruct) -> INT
	Replaces the current header by another of a set of fasta sequences.
	Requires the "PDStruct" object, the filename of the fasta sequences
	and the output filename. The resulting sequences will be written
	to the output file.
	"""
	IN = open(filename, "r")
	OUT = open(output, "w+")
	header = ""
	for line in IN:
		line = line.strip()
		if line[0] == ">":
			header = "|".join(PDStruct[line[1:]])
			OUT.write(">%s\n" % header)
		else:
			OUT.write("%s\n" % line)


	IN.close()
	OUT.close()


	return 1

def batchSeq_len(fastaFilename):
	"""

	"""
	formatseq(fastaFilename, None)
	IN = open("%s.out" % fastaFilename, "r")
	OUT = open("%s.len" % fastaFilename, "w+")
	for line in IN:
		line = line.strip()
		if line[0] == ">":
			header = line
			subgType = line.split("|")[0][-1]
		else:
			OUT.write("%s\t%d\t%s\n" % (header, len(line), subgType))
	IN.close()
	OUT.close()

def tab2fasta(filename, mode, delimiter):
	"""

	"""
	if mode == 1:
		# first sequence then header.
		IN = open(filename, "r")
		OUT = open("%s.fa" % filename, "w+")
		HASH = {}
		for line in IN:
			line = line.strip()
			OUT.write(">%s\n%s\n" % (line.split(delimiter)[1], line.split(delimiter)[0]))
	
	IN.close()
	OUT.close()

def extract_homologs(file1, file2):
	"""

	"""
	FH1 = open(file1, "r")
	FH2 = open(file2, "r")
	homologs = {}

	# Processing file1
	for line in FH1:
		line = line.strip()
		print line[0]
		if line[0] == ">":
			header = line
			ID = line.split("|")[0]
			print "this is the ID", ID
			homologs[ID] = []
		else:
			homologs[ID].append({"header":header, "sequence":line})
	
	# Processing file2
	for line in FH2:
		line = line.strip()
		if line[0] == ">":
			header = line
			ID = line.split("|")[0]
			ID = ".".join(ID.split(".")[:-1])
		else:
			try:
				homologs[ID].append({"header":header, "sequence":line})
			except Exception:
				pass
	
	FH1.close()
	FH2.close()

	OUT_FH1 = open("%s.hg.fa" % ".".join(file1.split(".")[:-1]), "w+")
	OUT_FH2 = open("%s.hg.fa" % ".".join(file2.split(".")[:-1]), "w+")
	CHK = open("checking_data.txt", "w+")
	CHK_ALL = open("checking_alldata.txt", "w+")

	for key, value in homologs.iteritems():
		for item in value:
			CHK_ALL.write("%s\t" % item["header"])
		CHK_ALL.write("\n")
		if len(value) == 3:
			OUT_FH1.write("%s\n%s\n" % (value[0]["header"], value[0]["sequence"]))
			OUT_FH2.write("%s\n%s\n" % (value[1]["header"], value[1]["sequence"]))
			OUT_FH2.write("%s\n%s\n" % (value[2]["header"], value[2]["sequence"]))
			CHK.write("%s\t%s\t%s\n" % (value[0]["header"], value[1]["header"], value[2]["header"]))
	
	OUT_FH1.close()
	OUT_FH2.close()

#//
#	DNA_2_RNA	#
#//
def dna_2_rna(dna):
	"""

	"""
	rna=""
	for nt in dna:
		if nt.upper() == "T":
			rna = rna + "U"
		else:
			rna = rna + nt.upper()
	return rna

#//
#	PICKUP_SEQS	#
#//
def pickup_seqs(headers_file, fasta_file):
	"""(STR, STR) -> INT & writes to file
		Requires the name of the headers file and the fasta file.
		
		Example:
		pickup_seqs("headers.txt", "three_utrs_xlae.fa").
		#>>> Returns a fasta file "headers.pckup.fa" with the sequences
		#>>> requested by each header from headers.txt file.
		<pseudo_code>
			1 open headers file for read
			2 build a dictionary with the headers as keys and empty values
			3 close headers file
			4 open sequences file for read
			5 assign the sequence to dictionary if sequence header match the dictionary key
			6 close sequences file
			7 open output file for write+
			8 write headers and sequences to output file
			9 return 1; __EOF__
		</pseudo_code>
	"""
	# Steps 1-3
	dictionary = {}
	FH_HEAD = open(headers_file, "r")
	for _IDENTIFIER_ in FH_HEAD:
		_IDENTIFIER_ = _IDENTIFIER_.strip()
		dictionary[_IDENTIFIER_] = ""
	FH_HEAD.close()
	# Steps 4-6
	FH_SEQS = open(fasta_file, "r")
	for SEQ in FH_SEQS:
		SEQ = SEQ.strip()
		if SEQ[0] == ">":
			_head_ = SEQ
		else:
			try:
				dictionary[_head_] = dictionary[_head_] + SEQ
			except Exception:
				print "ERROR: identifier not found in dictionary!"
	FH_SEQS.close()
	#Steps 7-9
	FH_OUT = open("%s.pckup.fa" % ".".join(headers_file.split(".")[:-1]), "w+")
	for IDENTIFIER, SEQUENCE in dictionary.iteritems():
		if SEQUENCE != "":
			FH_OUT.write("%s\n%s\n" % (IDENTIFIER, SEQUENCE))
		else:
			continue
	
	return 1


def help():
	"""
	HELP FUNCTION
	"""

	print """
	Library Utils.
	- formatseq(filename, formhead)
		Deletes the return carriage of fasta formatted sequences.
		Changing header style:
		if formhead == "style1": splits the header (sep="|") and
		change the format from >item1|item2|item3 to >item3[1:]|item[2][:-2]
	- headsArray(filename)
		Creates an array of sequence headers from data of fasta
		file.
	- geneModel(gff_file)
		Creates a data structure that stores genetic features data
		from a transcriptome. For more information see README.txt.
	- load_Ome(filename)
		Loads sequences of an "Ome" (genome, transcriptome, proteome)
		and stores them in a HASH variable HASH[header:sequence]
	- rev_comp(adn)
		Returns the reverse complement of a DNA string
	- rev_trans(rna)
		Returns the reverse transcription of a RNA string
	- count_matches(pattern, sequence)
		Returns the number of occurences of pattern in sequence.
	- replace_header(filename, output, PDStruct)
		Replaces headers of fasta formated sequences.
	- PDStruct(filename)
		Creates a HASH data structure with the info from filename
		file. This file consists of two columns; the first is the
		ID of the header replacement. (parentDB.txt)
	- batchSeq_len(filename)
		returns a file with 2 columns of data; the first column
		is the header ID of each sequence and the second column
		is the length of that sequence.
	- extract_homologs(file1, file2)
		Extract homologos gene sequences from two files acording to
		the header of the sequences.
	- pickup_seqs(headers_file, fasta_file)
		Pickup sequences based on the identifier. The identifier header
		of the sequence must match the header of the headers file.
	"""

if __name__ == "__main__":
	help()
