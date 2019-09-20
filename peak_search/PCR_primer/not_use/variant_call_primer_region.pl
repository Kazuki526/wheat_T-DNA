#!/usr/bin/perl
use warnings;
use strict;


my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}


open(PRIMER,"$ENV{HOME}/Dropbox/cooperative/plant_gen/T-DNA/count_all_read/primer_comp.txt");
my $header = <PRIMER>;chomp $header;
my %col = &header2hash($header);
while(<PRIMER>){
		chomp;
		my @line = split(/\t/,);
		my($sample,$chr,$start,$end)=($line[$col{max_sample}],$line[$col{chr}],$line[$col{primer_start}],$line[$col{primer_end}]);
		my $region = "$chr:$start-$end";
		my $bam = "$sample/$sample"."_grouped.bam";
		my $vcf = "PCR_primer/vcf/$sample"."_$region.vcf";
		`$ENV{HOME}//gatk-4.1.2.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R /Volumes/kazcancer/T-DNA_search/genome/all.fa -I $bam -O $vcf --dont-use-soft-clipped-bases true -L $region -DF NotSecondaryAlignmentReadFilter -DF MappingQualityReadFilter`;
}
close PRIMER;




sub header2hash ( $ ){
		my $header = $_[0];
		my @colm = split(/\t/,$header);
		my %out = ();
		for(my $i=0; $i < @colm; $i++){
				$out{$colm[$i]}=$i;
		}
		return(%out);
}
