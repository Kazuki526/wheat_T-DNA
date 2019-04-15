#!/usr/bin/perl
use warnings;
use strict;

my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}

my $fasta = "../genome/all.fa";

# read all_peak_table
my %peak_tbl =();
open(IN,"all_peak_table_with_p.tsv") or die "ERROR::cannot open all_peak_table_with_p.tsv\n";
my $header=<IN>;chomp $header;
my %col =&header2hash($header);
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		if($line[$col{max_p}] >10**-5){last;}
		my $region = "$line[0]:$line[1]-$line[2]";
		$peak_tbl{$region}{p} = $line[$col{max_p}];
		$peak_tbl{$region}{all_depth} = $line[$col{$line[$col{max_sample}]}];
		$line[1] -=500; $line[2] +=500;
		$peak_tbl{$region}{around} = "$line[0]:$line[1]-$line[2]";
}
close IN;

open(IN,"all_peak_table_paired_end_plasmid.tsv") or die "ERROR::cannot open all_peak_table_paired_end_plasmid.tsv\n";
$header=<IN>;chomp $header;
%col =&header2hash($header);
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		my $region = "$line[0]:$line[1]-$line[2]";
		if(!defined $peak_tbl{$region}){last;}
		$peak_tbl{$region}{paired_depth} = $line[$col{max_depth}];
		$peak_tbl{$region}{paired_second_d}="";
		$peak_tbl{$region}{paired_second_s}="";
		if($line[$col{second_depth}] >10){
				$peak_tbl{$region}{paired_second_d}=$line[$col{second_depth}];
				$peak_tbl{$region}{paired_second_s}=$line[$col{second_sample_pairedend}];
		}
}
close IN;

open(IN,"all_peak_table_chimeric_align_plasmid.tsv") or die "ERROR::cannot open all_peak_table_chimeric_align_plasmid.tsv\n";
$header=<IN>;chomp $header;
%col =&header2hash($header);
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		my $region = "$line[0]:$line[1]-$line[2]";
		if(!defined $peak_tbl{$region}){last;}
		$peak_tbl{$region}{chim_depth} = $line[$col{max_depth}];
		$peak_tbl{$region}{chim_second_d}="";
		$peak_tbl{$region}{chim_second_s}="";
		if($line[$col{second_depth}] >10){
				$peak_tbl{$region}{chim_second_d}=$line[$col{second_depth}];
				$peak_tbl{$region}{chim_second_s}=$line[$col{second_sample_chimeric}];
		}
}
close IN;

if(!-e "fasta_dir_signif/group"){mkdir "fasta_dir_signif/group";}
my @sample_num=qw(01 02 03 04 05 06 07 08 09 10 11 12);
open(OUT,">homologous_peak_list.tsv");
print OUT "sample\tgroup\tpeak_region\tp_value\tall_depth\t/paired_depth\tchimeric_depth\tpaired_second\tchimeric_second\n";
my $group=1;
#GFP
foreach my $num(@sample_num){
		my $sample_id ="GFP-$num";
		my ($out,$group_) = &main($sample_id,$group);
		$group = $group_;
		print OUT "$out";
}
#mCherry
foreach my $num(@sample_num){
		my $sample_id ="mCherry-$num";
		my ($out,$group_) = &main($sample_id,$group);
		$group = $group_;
		print OUT "$out";
}
exit;

sub header2hash ( $ ){
		my $header = $_[0];
		my @colm = split(/\t/,$header);
		my %out = ();
		for(my $i=0; $i < @colm; $i++){
				$out{$colm[$i]}=$i;
		}
		return(%out);
}


sub main( $ $){
		my ($sample, $group)=@_;
		my $hom_file = "homologous_peak_list_signif/$sample.tsv";
		my %hom_list=();
		my %group=();
		my $out="";
		open(HOM,"$hom_file")or die "ERROR::cannot open $hom_file\n";
		<HOM>;
		my $region="";
		my @list=();
		while(<HOM>){
				chomp;
				if($_ =~ /^##/){
						if($_ =~ /from\s(chr\S+:\d+-\d+)\s/){$region =$1;}else{die "ERROR:cannot extract region $hom_file:$_\n";}
						next;
				}elsif($_ eq ""){
						if(scalar(@list) ==0){next;}
						$out .="$sample\t$group\t$region\t$peak_tbl{$region}{p}\t$peak_tbl{$region}{all_depth}\t$peak_tbl{$region}{paired_depth}\t$peak_tbl{$region}{chim_depth}\t";
						if($peak_tbl{$region}{paired_second_s} ne ""){
								$out .="$peak_tbl{$region}{paired_second_s}:$peak_tbl{$region}{paired_second_d}\t";
						}else{$out .="NA\t";}
						if($peak_tbl{$region}{chim_second_s} ne ""){
								$out .= "$peak_tbl{$region}{chim_second_s}:$peak_tbl{$region}{chim_second_d}\n";
						}else{$out .= "NA\n";}
						`samtools faidx $fasta $peak_tbl{$region}{around} >fasta_dir_signif/group/group$group.fasta`;
						foreach $region(@list){
								$out .="$sample\t$group\t$region\t$peak_tbl{$region}{p}\t$peak_tbl{$region}{all_depth}\t$peak_tbl{$region}{paired_depth}\t$peak_tbl{$region}{chim_depth}\t";
								if($peak_tbl{$region}{paired_second_s} ne ""){
										$out .="$peak_tbl{$region}{paired_second_s}:$peak_tbl{$region}{paired_second_d}\t";
								}else{$out .="NA\t";}
								if($peak_tbl{$region}{chim_second_s} ne ""){
										$out .= "$peak_tbl{$region}{chim_second_s}:$peak_tbl{$region}{chim_second_d}\n";
								}else{$out .= "NA\n";}
								`samtools faidx $fasta $peak_tbl{$region}{around} >>fasta_dir_signif/group/group$group.fasta`;
						}
						@list=();
						$group++;
				}else{
						my @line = split(/\t/,);
						if(defined $hom_list{$line[0]}){next;}
						$hom_list{$line[1]}="exist";
						push(@list,$line[1]);
				}
		}
		close HOM;
		return($out,$group);
}


