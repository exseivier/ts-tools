#!/usr/bin/env bash


echo -e "Postprocessing TargetScan Results..."
echo -e "\n\n"
echo -e "Sorting homologues genes & extracting 25% of top and bottom data..."
for file in *.data;
do
	sort -nk 2 $file > tmp
	mv tmp $file
	ehecatl -p join-homologs -i $file --bool-single --delim "tab"
	head -n $(echo "25/100 * $(wc -l $file | tr -cd [:digit:])" | bc -l | cut -d "." -f1) $file > ${file}.fu
	tail -n $(echo "25/100 * $(wc -l $file | tr -cd [:digit:])" | bc -l | cut -d "." -f1) $file > ${file}.fd
done


echo -e "Calculating the number of total, $1 ans $2 genes in dataset..."
echo -e "UP.Name\tUP.$1_genes\tUP.$2_genes\tUP.Total" > filtered-up.txt
echo -e "DOWN.Name\tDOWN.$1_genes\tDOWN.$2_genes\tDOWN.Total" > filtered-down.txt
echo -e "miR-Family\tUP.Intersection\tUP.$1_genes\tUP.$2_genes" > up.inter.txt
echo -e "miR-Family\tDOWN.Intersection\tDOWN.$1_genes\tDOWN.$2_genes" > down.inter.txt
for file in *.{fu,fd};
do
	file_ext=$(echo $(echo $file | cut -d"." -f3))
	if [ "$file_ext" == "fu" ]; then
		output="filtered-up.txt"
		output2="up.inter.txt"
	elif [ "$file_ext" == "fd" ]; then
		output="filtered-down.txt"
		output2="down.inter.txt"
	else
		echo "Something was wrong with selecting outputfile..."
		exit 0
	fi
	L=$(grep "\.$1" $file | wc -l)
	if [ "$2" == "" ]; then
		S=$(grep -v "\.$1" $file | wc -l)
	else
		S=$(grep "\.$2" $file | wc -l)
	fi
	INTER=$(cut -f1 $file | cut -d"." -f1 | sort | uniq -c | grep "2 " | wc -l)
	Lm=$(echo "$L - $INTER" | bc -l)
	Sm=$(echo "$S - $INTER" | bc -l)
	TOTAL=$(echo "$S + $L" | bc -l)
	# WRITE DATA TO FILE.. TODO
	echo -e "$file\t$L\t$S\t$TOTAL" >> $output
	echo -e "$file\t$INTER\t$Lm\t$Sm" >> $output2
done

echo -e "Pasting filtered files..."

paste filtered-up.txt filtered-down.txt > filtered.txt
