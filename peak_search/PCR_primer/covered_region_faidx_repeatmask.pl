#!/usr/bin/perl
use strict;
use warnings;

my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}
my $genome_file = "../genome/all.fa";
-e $genome_file or die "ERROR::not exist $genome_file\n";

open(IN,"PCR_primer/primer_comp_homologous_removed.txt") or die "ERROR::cannot open primer_comp.txt\n";
<IN>; #header
my %primer=();
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		my $region = "$line[1]:$line[8]-$line[9]";
		my $id = "$line[4]_$region";
		$primer{$id}{sample}=$line[4];
		$primer{$id}{region}=$region;
		$primer{$id}{chr}=$line[1];
		$primer{$id}{p_start}=$line[8];
		$primer{$id}{p_end}=$line[9];
		$primer{$id}{peak}="$line[1]:$line[2]-$line[3]";
		$primer{$id}{side}=$line[5];
		$primer{$id}{seq}=&faidx_seq($region);
}
close IN;
$|=1;
my %variant=();
open(VCF,"PCR_primer/manual_all_mutation_uchino_ym-check.txt") or die "ERROR::cannot open manual_vcf file\n";
<VCF>;#header
while(<VCF>){
		chomp;
		my @line = split(/\t/,);
		if(!defined $variant{$line[0]}){
				$variant{$line[0]}="$line[2]:$line[3]:$line[4]";
		}else{
				$variant{$line[0]}="$line[2]:$line[3]:$line[4],$variant{$line[0]}";
		}
}
close VCF;


foreach my $id (keys %variant){
		if(!defined $primer{$id}{side}){next;}
		my @vars=split(/,/,$variant{$id});
		foreach my $var (@vars){
				my ($posi,$ref,$alt) = split(/:/,$var);
				if(($posi < $primer{$id}{p_start})||($posi > $primer{$id}{p_end})){next;}
				my $ref_=substr($primer{$id}{seq},$posi-$primer{$id}{p_start},length($ref),$alt);
				if($ref_ ne $ref){die "ERROR::$id $var substr $ref_\n";}
		}
}

print "done faidx and substring variant\n";

my $fasta = "PCR_primer/all_primer_homologous_removed/all_primer_homologous_removed_region.fasta";
open(OUT,">$fasta");
foreach my $id(sort keys %primer){
		if(!defined $primer{$id}{seq}){next;}
		print OUT ">$id\n$primer{$id}{seq}\n";
}
close OUT;

print "repeat masker\n";
system("repeatmasker -lib /../../mipsREdat_9.3p_ALL.fasta $fasta");

print "blast\n";
system("blastn -db ../genome/all.fa -query $fasta.masked -out PCR_primer/all_primer_homologous_removed/all_primer_region_blast.tsv -outfmt 7");





sub header2hash ( $ ){
		my $header = $_[0];
		my @colm = split(/\t/,$header);
		my %out = ();
		for(my $i=0; $i < @colm; $i++){
				$out{$colm[$i]}=$i;
		}
		return(%out);
}
sub faidx_seq( $ ){
		my $region = $_[0];
		open(FA,"samtools faidx $genome_file $region|");
		my $name = <FA>;
		if($name !~/^>/){die "ERROR::primer region $region is not exist\n";}
		my $seq = "";
		while(<FA>){
				chomp;
				$seq.=uc $_;
		}
		close FA;
		return($seq);
}
