# ts-tools
Tools to automate "target scan 7" program.

----

		*******************************
		****      Target-Scan      ****
		****         PIPE          ****
		****       TUTORIAL        ****
		*******************************
	Usage:
		ts-tools -s all -i < utr.fasta > -m < mirnaome.fasta > --tscs /path/to/targetscan_context_score/parameters/ --type1 'L' --type2 'S' --window 100 --percent 25

	Options:

	-s		Stage:
			all: runs the entire pipeline.
			format: fomrats input files (UTRs and miRNA database).
			tscan: runs targetScan_70.pl.
			emulfile: emulate output file from PCT score and BCL.
			cspp: runs targetScan_context_score.pl.
			tsplot: runs Rscript for ploting target scan results.
	-i		UTR sequences in fasta format.
	-m		mature miRNA sequences in fasta format.
	--type1		Paralog gene 1 label.
	--type2		Paralog gene 2 label.
	--window	Window size in hypergeometric test of target-scan results.
	--percent	Percentage of top gene targets used for target-scan results analysis.
	--tscs		Path to mathematical model parameters files.
----

## Instalation instructions.
