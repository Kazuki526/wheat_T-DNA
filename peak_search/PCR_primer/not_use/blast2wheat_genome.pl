#!/usr/bin/perl
use strict;
use warnings;

my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}

my $dir = "PCR_primer";
my $peak_file = "all_peak_table_with_p.tsv";
my $all_peak_fasta = "$dir/all_peak.fasta";
my $genome = "../genome/all.fa";
(-e $genome) or die "ERROR::$genome is not exist\n";

#exrract around (+-25bp) peak sequence
open(PEAK,"$peak_file") or die "ERROR::cannot open $peak_file\n";
my $header = <PEAK>;
my %col = &header2hash($header);
while(<PEAK>){
		chomp;
		my @line = split(/\t/,);
		if($line[$col{max_p}] > 0.00001){last;}
		my ($start,$end)=($line[1]-25,$line[2]+25);
		my $region = "$line[0]:$start-$end";
		`samtools faidx $genome $region >>$all_peak_fasta`;
}
close PEAK;

#BLAST
#`makeblastdb -in $genome -out $genome -dbtype nucl -parse_seqids`;
system("blastn -query $all_peak_fasta -db $genome -out $dir/signif_peak_blast.tsv -outfmt 6");






sub header2hash ( $ ){
		my $header = $_[0];
		my @colm = split(/\t/,$header);
		my %out = ();
		for(my $i=0; $i < @colm; $i++){
				$out{$colm[$i]}=$i;
		}
		return(%out);
}
