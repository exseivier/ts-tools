#!/usr/bin/env bash

usage(){
	
	echo -e "\t\t*******************************"
	echo -e "\t\t****      Target-Scan      ****"
	echo -e "\t\t****         PIPE          ****"
	echo -e "\t\t****       TUTORIAL        ****"
	echo -e "\t\t*******************************"
	echo -e "\n\n"
	echo -e "\tUsage:"
	echo -e "\t\tts-tools -s all -i < utr.fasta > -m < mirnaome.fasta > --tscs /path/to/targetscan_context_score/parameters/ --type1 'L' --type2 'S' --window 100 --percent 25\n"
	echo -e "\tOptions:\n"
	echo -e "\t-s\t\tStage:"
	echo -e "\t\t\tall: runs the entire pipeline."
	echo -e "\t\t\tformat: fomrats input files (UTRs and miRNA database)."
	echo -e "\t\t\ttscan: runs targetScan_70.pl."
	echo -e "\t\t\temulfile: emulate output file from PCT score and BCL."
	echo -e "\t\t\tcspp: runs targetScan_context_score.pl."
	echo -e "\t\t\ttsplot: runs Rscript for ploting target scan results."
	echo -e "\t-i\t\tUTR sequences in fasta format."
	echo -e "\t-m\t\tmature miRNA sequences in fasta format."
	echo -e "\t--type1\t\tParalog gene 1 label."
	echo -e "\t--type2\t\tParalog gene 2 label."
	echo -e "\t--window\tWindow size in hypergeometric test of target-scan results."
	echo -e "\t--percent\tPercentage of top gene targets used for target-scan results analysis."
	echo -e "\t--tscs\t\tPath to mathematical model parameters files."
}


while getopts "hs:i:m:-:" opt
do
	if [ -z $opt ]; then
		echo "arguments were not specified"
		clear
		usage
		exit 0
	fi
	case $opt in
		h)
			clear
			usage
			exit 0
		;;
		s)
			STAGE=$OPTARG
		;;
		i)
			INPUT=$OPTARG
		;;
		m)
			MIRNA=$OPTARG
		;;
		-)
			case $OPTARG in
				type1)
					TYPE1=${!OPTIND}
					OPTIND=$(( $OPTIND + 1 ))
				;;
				type2)
					TYPE2=${!OPTIND}
					OPTIND=$(( $OPTIND + 1 ))
				;;
				window)
					WINDOW=${!OPTIND}
					OPTIND=$(( $OPTIND + 1 ))
				;;
				percent)
					PERCENT=${!OPTIND}
					OPTIND=$(( $OPTIND + 1 ))
				;;
				tscs)
					TSCS_HOME=${!OPTIND}
					OPTIND=$(( $OPTIND + 1 ))
				;;
				help)
					echo "Bad double hyphen argument"
					clear
					usage
					exit 0
				;;
			esac
		;;
		*)
			echo "Bad single hyphen argument"
			clear
			usage
			exit 1
		;;
	esac
done

INPUTFILE=$INPUT
MIRNAFILE=$MIRNA

#INPUTFILE=$(basename $INPUT '/')
#MIRNAFILE=$(basename $MIRNA '/')


case $STAGE in

	all)
		rand=$RANDOM
		# Hacer que agregue un random numero al archivo de salida para correr varias instancias de este programa y no confundir los archivos de salida si el de entrada es el mismo
		# TODO
		echo "Preprocessing input files..."
		fa2ints $INPUTFILE
		miRNAfamily-build $MIRNAFILE
		echo "Finding miRNA targets..."
		targetScan ${MIRNAFILE%.*}.ints ${INPUTFILE%.*}.ints ${INPUTFILE%.*}.outs
		echo "Emulating file..."
		target2context ${INPUTFILE%.*}.outs > ${INPUTFILE%.*}.outsc++
		echo "Linking parameters..."
		ln -s ${TSCS_HOME}/Agarwal_2015_parameters.txt Agarwal_2015_parameters.txt
		ln -s ${TSCS_HOME}/All_cell_lines.AIRs.txt All_cell_lines.AIRs.txt
		ln -s ${TSCS_HOME}/TA_SPS_by_seed_region.txt TA_SPS_by_seed_region.txt
		echo "Calculating context score..."
		targetScan-ctx++ ${MIRNAFILE%.*}.intsc++ ${INPUTFILE%.*}.ints ${INPUTFILE%.*}.outsc++ \
				 ${TSCS_HOME}/ORF_8mer_counts_sample.txt \
				 ${TSCS_HOME}/ORF_Sequences_sample.lengths.txt \
				 ${INPUTFILE%.*}.ts
		echo "Separating results..."
		targetScan-final ${INPUTFILE%.*}.ts
		echo "Preprocessing data..."
		mkdir ${INPUTFILE%.*}_out
		mv *.data ${INPUTFILE%.*}_out
		cd ${INPUTFILE%.*}_out
		targetScan-pp "$TYPE1" "$TYPE2"
		cd ..
		echo "Sattistic analysis & ploting..."
		targetScan-stats ${INPUTFILE%.*}_out $WINDOW $PERCENT
	;;
	format)
		echo "Preprocessing input files..."
		fa2ints $INPUTFILE
		miRNAfamily-build $MIRNAFILE
	;;
	tscan)
		echo "Finding miRNA targets..."
		targetScan ${MIRNAFILE%.*}.ints ${INPUTFILE%.*}.ints ${INPUTFILE%.*}.outs
	;;
	emulfile)
		echo "Emulating file..."
		target2context ${INPUTFILE%.*}.outs > ${INPUTFILE%.*}.outsc++
	;;
	cspp)
		echo "Linking parameters..."
		ln -s ${TSCS_HOME}/Agarwal_2015_parameters.txt Agarwal_2015_parameters.txt
		ln -s ${TSCS_HOME}/All_cell_lines.AIRs.txt All_cell_lines.AIRs.txt
		ln -s ${TSCS_HOME}/TA_SPS_by_seed_region.txt TA_SPS_by_seed_region.txt
		echo "Calculating context score..."
		targetScan-ctx++ ${MIRNAFILE%.*}.intsc++ ${INPUTFILE%.*}.ints ${INPUTFILE%.*}.outsc++ \
				 ${TSCS_HOME}/ORF_8mer_counts_sample.txt \
				 ${TSCS_HOME}/ORF_Sequences_sample.lengths.txt \
				 ${INPUTFILE%.*}.ts
		echo "Separating results..."
		targetScan-final ${INPUTFILE%.*}.ts
		echo "Preprocessing data..."
		mkdir ${INPUTFILE%.*}_out
		mv *.data ${INPUTFILE%.*}_out
		cd ${INPUTFILE%.*}_out
		targetScan-pp "$TYPE1" "$TYPE2"
		cd ..
	;;
	tsplot)
	echo "Sattistic analysis & ploting..."
	targetScan-stats ${INPUTFILE%.*}_out $WINDOW $PERCENT
	;;
esac
echo "Unliking parameters..."
unlink Agarwal_2015_parameters.txt
unlink All_cell_lines.AIRs.txt
unlink TA_SPS_by_seed_region.txt

echo "Everything is OK... ;)"

