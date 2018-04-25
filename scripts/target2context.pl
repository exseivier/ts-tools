#!/usr/bin/env perl

my $file=$ARGV[0];
my $counter=0;
open (SEEDS,"$file");
while (my $line=<SEEDS>){
	my @seeds=split(/\t/,$line);
	#print "$seeds[0]\t$seeds[1]\t$seeds[2]\t$seeds[3]\t$seeds[4]\t$seeds[5]\t$seeds[6]\t$seeds[7]\t$seeds[8]\t$seeds[9]\t$seeds[10]\t$seeds[11]\t$seeds[12]\t$seeds[13]" if ($counter<=10);
	print "$seeds[0]\t$seeds[1]\t$seeds[2]\t$seeds[3]\t$seeds[4]\t$seeds[5]\t$seeds[6]\t$seeds[7]\t$seeds[8]\t$seeds[9]\t$seeds[10]\tBranchLS\tPct\tconserved\n" if ($counter<=0);
	print "$seeds[0]\t$seeds[1]\t$seeds[2]\t$seeds[3]\t$seeds[4]\t$seeds[5]\t$seeds[6]\t$seeds[7]\t$seeds[8]\t$seeds[9]\t$seeds[10]\t0\tNA\t0\n" if ($counter>=1);
	$counter++;
}
