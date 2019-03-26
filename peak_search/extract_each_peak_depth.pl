#!/usr/bin/perl
use warnings;
use strict;

my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}

my $bed="100_peak.bed";
open(IN,"100_peak_list.tsv");
open(BED,">$bed");
<IN>;
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		if($line[1] >0){$line[1]--;}
		print BED "$line[0]\t$line[1]\t$line[2]\n";
}
close IN;
close BED;

my @sample_num=qw(01 02 03 04 05 06 07 08 09 10 11 12);
open(OUT,">all_peak_table_tidy.tsv");
#GFP
foreach my $num(@sample_num){
		my $sample_id ="GFP-$num";
		&main($sample_id);
		print "done $sample_id\n";
}
#mCherry
foreach my $num(@sample_num){
		my $sample_id ="mCherry-$num";
		&main($sample_id);
		print "done $sample_id\n";
}
close OUT;
exit;

sub main($){
		my $sample =$_[0];
		my $bam="$sample/$sample.bam";
		open(DP,"samtools depth -a -b $bed $bam|");
		my($chr,$start,$bef_position,$width,$sum_depth)=("",0,0,0,0);
		while(<DP>){
				chomp;
				my @line = split(/\t/,);
				if(($line[0] eq $chr)&&($bef_position +1 ==$line[1])){
						$width++;
						$sum_depth += $line[2];
						$bef_position++;
				}else{
						if($width ==0){
								($chr,$start,$bef_position,$width,$sum_depth)=($line[0],$line[1],$line[1],1,$line[2]);
								next;
						}
						my $average_depth = $sum_depth/$width;
						print OUT "$sample\t$chr\t$start\t$bef_position\t$average_depth\n";
						($chr,$start,$bef_position,$width,$sum_depth)=($line[0],$line[1],$line[1],1,$line[2]);
				}
		}
		my $average_depth = $sum_depth/$width;
		print OUT "$sample\t$chr\t$start\t$bef_position\t$average_depth\n";
}
