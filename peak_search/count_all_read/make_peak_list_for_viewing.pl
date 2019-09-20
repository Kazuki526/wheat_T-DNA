#!/usr/bin/perl
use warnings;
use strict;

my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}

my $fasta = "../genome/all.fa";

# read all_peak_table
my %peak_tbl =();
my %for_peak_arrange=();
open(IN,"count_all_read/all_peak_table_with_p.tsv") or die "ERROR::cannot open all_peak_table_with_p.tsv\n";
my $header=<IN>;chomp $header;
my %col =&header2hash($header);
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		if($line[$col{max_p}] <10**-5){$for_peak_arrange{$line[0]}{$line[1]}=$line[2];}
		my $region = "$line[0]:$line[1]-$line[2]";
		$peak_tbl{$region}{p} = $line[$col{max_p}];
		$peak_tbl{$region}{sample} = $line[$col{max_sample}];
		$peak_tbl{$region}{all_depth} = $line[$col{$line[$col{max_sample}]}];
		$line[1] -=100; $line[2] +=100;
		$peak_tbl{$region}{around} = "$line[0]:$line[1]-$line[2]";
}
close IN;

#only new peak
my %for_peak_arrange_new=();
open(IN,"count_all_read/new_peak/new_peak_with_p.tsv") or die "ERROR::cannot open new_peak_with_p.tsv\n";
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		$for_peak_arrange_new{$line[0]}{$line[1]}=$line[2];
}
close IN;


open(IN,"count_all_read/all_peak_table_paired_end_plasmid.tsv") or die "ERROR::cannot open all_peak_table_paired_end_plasmid.tsv\n";
$header=<IN>;chomp $header;
%col =&header2hash($header);
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		my $region = "$line[0]:$line[1]-$line[2]";
		$peak_tbl{$region}{paired_depth} = $line[$col{max_depth_pairedend}];
		$peak_tbl{$region}{paired_second_d}="";
		$peak_tbl{$region}{paired_second_s}="";
		if($line[$col{second_depth_pairedend}] eq "NA"){next;}
		if($line[$col{second_depth_pairedend}] >10){
				$peak_tbl{$region}{paired_second_d}=$line[$col{second_depth_pairedend}];
				$peak_tbl{$region}{paired_second_s}=$line[$col{second_sample_pairedend}];
		}
}
close IN;

open(IN,"count_all_read/all_peak_table_chimeric_align_plasmid.tsv") or die "ERROR::cannot open all_peak_table_chimeric_align_plasmid.tsv\n";
$header=<IN>;chomp $header;
%col =&header2hash($header);
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		my $region = "$line[0]:$line[1]-$line[2]";
		$peak_tbl{$region}{chim_depth} = $line[$col{max_depth_chimeric}];
		$peak_tbl{$region}{chim_second_d}="";
		$peak_tbl{$region}{chim_second_s}="";
		if($line[$col{second_depth_chimeric}] eq "NA"){next;}
		if($line[$col{second_depth_chimeric}] >10){
				$peak_tbl{$region}{chim_second_d}=$line[$col{second_depth_chimeric}];
				$peak_tbl{$region}{chim_second_s}=$line[$col{second_sample_chimeric}];
		}
}
close IN;

open(IN,"count_all_read/homologous_peak_list.tsv") or die "ERROR::cannot open homologous_peak_list.tsv\n";
$header=<IN>;chomp $header;
%col =&header2hash($header);
while(<IN>){
		chomp;
		my @line =split(/\t/,);
		$peak_tbl{$line[$col{peak_region}]}{homo}=$line[$col{homologous_peak_list}];
}
close IN;

my $all_seq = "count_all_read/fasta_dir_signif/all_peak_region_for_blast.fasta";
if(-e $all_seq){unlink $all_seq;}
open(OUT,">count_all_read/all_peak_for_viewing.tsv");
print OUT "chr\tstart\tend\tmax_sample\tp_value\tPrimary_or_Secondary\tall_depth\tpaired_depth\tchimeric_depth\tpaired_second\tchimeric_second\thomollgous_list\n";
foreach my $chr(sort keys%for_peak_arrange){
		foreach my $start (sort {$a <=> $b} keys %{$for_peak_arrange{$chr}}){
				my $end = $for_peak_arrange{$chr}{$start};
				my $region = "$chr:$start-$end";
				my $primary_secondary = "Primary";
				if(defined $for_peak_arrange_new{$chr}{$start}){$primary_secondary="Secondary";}
				print OUT "$chr\t$start\t$end\t";
				print OUT "$peak_tbl{$region}{sample}\t$peak_tbl{$region}{p}\t$primary_secondary\t";
				print OUT "$peak_tbl{$region}{all_depth}\t$peak_tbl{$region}{paired_depth}\t$peak_tbl{$region}{chim_depth}\t";
				if($peak_tbl{$region}{paired_second_s} ne ""){
						print OUT "$peak_tbl{$region}{paired_second_s}:$peak_tbl{$region}{paired_second_d}\t";
				}else{print OUT "NA\t";}
				if($peak_tbl{$region}{chim_second_s} ne ""){
						print OUT "$peak_tbl{$region}{chim_second_s}:$peak_tbl{$region}{chim_second_d}\t";
				}else{print OUT "NA\t";}
				if(defined $peak_tbl{$region}{homo}){
						print OUT "$peak_tbl{$region}{homo}\n";
				}else{print OUT "NA\n";}
				$start -=500; $end +=500;
				my $around = "$chr:$start-$end";
				`samtools faidx $fasta $around >>$all_seq`;
		}
}
close OUT;
#only new peak list
open(OUT,">count_all_read/new_peak/new_peak_for_viewing.tsv");
print OUT "chr\tstart\tend\tmax_sample\tp_value\tall_depth\tpaired_depth\tchimeric_depth\tpaired_second\tchimeric_second\thomollgous_list\n";
foreach my $chr(sort keys%for_peak_arrange_new){
		foreach my $start (sort {$a <=> $b} keys %{$for_peak_arrange_new{$chr}}){
				my $end = $for_peak_arrange_new{$chr}{$start};
				my $region = "$chr:$start-$end";
				print OUT "$chr\t$start\t$end\t";
				print OUT "$peak_tbl{$region}{sample}\t$peak_tbl{$region}{p}\t$peak_tbl{$region}{all_depth}\t";
				print OUT "$peak_tbl{$region}{paired_depth}\t$peak_tbl{$region}{chim_depth}\t";
				if($peak_tbl{$region}{paired_second_s} ne ""){
						print OUT "$peak_tbl{$region}{paired_second_s}:$peak_tbl{$region}{paired_second_d}\t";
				}else{print OUT "NA\t";}
				if($peak_tbl{$region}{chim_second_s} ne ""){
						print OUT "$peak_tbl{$region}{chim_second_s}:$peak_tbl{$region}{chim_second_d}\t";
				}else{print OUT "NA\t";}
				if(defined $peak_tbl{$region}{homo}){
						print OUT "$peak_tbl{$region}{homo}\n";
				}else{print OUT "NA\n";}
				$start -=500; $end +=500;
				my $around = "$chr:$start-$end";
		}
}
close OUT;


sub header2hash ( $ ){
		my $header = $_[0];
		my @colm = split(/\t/,$header);
		my %out = ();
		for(my $i=0; $i < @colm; $i++){
				$out{$colm[$i]}=$i;
		}
		return(%out);
}

