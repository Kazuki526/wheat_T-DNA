#!/usr/bin/perl
use warnings;
use strict;

my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}

open(PEAK,"count_all_read/all_peak_table_with_p.tsv") or die "ERROR::cannot open PEAK\n";
my $header = <PEAK>;
my %col = &header2hash($header);
my %peak_info =();
my @peak_list =();
my %sample_peak =();
while(<PEAK>){
		chomp;
		my @line = split (/\t/,);
		if($line[$col{max_p}] > 0.00001){last;}
		my $region = "$line[0]:$line[1]-$line[2]";
		my $sample = $line[$col{max_sample}];
		push(@peak_list,$region);
		$peak_info{$region}=$sample;
		my ($start,$end) = ($line[1]-1000,$line[2]+1000);
		$sample_peak{$sample}{$line[0]}="$line[1]-$line[2]";
}
close PEAK;

open(OUT,">count_all_read/signif_peak_mate&chimeric_grouping.tsv");
print OUT "#group1 same chromosome or not chimeric aligned\n";
print OUT "#group2 palsmid\n";
print OUT "#group3 same sample other signif peak\n";
print OUT "#group4 other\n";
print OUT "chr\tstart\tend\tsample\tall_read\tmate_group1\tmate_group2\tmate_group3\tmate_group4\tchimeric_group1\tchimeric_group2\tchimeric_group3\tchimeric_group4\n";
foreach my $region (@peak_list){
		my ($chr,$start,$end);
		if($region =~ /^(chr.+):(\d+)-(\d+)$/){($chr,$start,$end)=($1,$2,$3);
		}else{die "ERROR::what region $region\n";}
		my ($start_,$end_)=($start-1000,$end+1000);
		my $sample = $peak_info{$region};
		my $all_read=0;
		my @mate_info =(0,0,0,0);
		my @chimeric_info=(0,0,0,0);
		open(BAM,"samtools view $sample/$sample.bam $region|");
		while(<BAM>){
				chomp;
				my @line = split(/\t/,);
				$all_read++;
#mate grouping
				if(($line[6] eq "=")&&($line[7] > $start_)&&($line[7] < $end_)){
						$mate_info[0]++;
				}elsif($line[6] =~ /plasmid/){
						$mate_info[1]++;
				}elsif(&other_peak_search($sample,$line[6],$line[7])){
						$mate_info[2]++;
				}else{
						$mate_info[3]++;
				}
#chimeric info
				my ($cchr,$cposi);
				if($_ =~ /SA:Z:([^,]+),(\d+),/){($cchr,$cposi)=($1,$2);
				}else{
						$chimeric_info[0]++;next;
				}
				if($cchr =~ /plasmid/){
						$chimeric_info[1]++;
				}elsif(&other_peak_search($sample,$cchr,$cposi)){
						$chimeric_info[2]++;
				}else{
						$chimeric_info[3]++;
				}
		}
		close BAM;
		my $out = "$chr\t$start\t$end\t$sample\t$all_read\t" . join("\t",@mate_info) . "\t". join("\t",@chimeric_info) . "\n";
		print OUT $out;
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

sub other_peak_search( $ $ $ ){
		my($sample, $chr, $posi) = @_;
		my $other_peak=0;
		if(defined $sample_peak{$sample}{$chr}){
				foreach my $sted (split(/;/,$sample_peak{$sample}{$chr})){
						my($mst,$med) = split(/-/,$sted);
						if(($mst < $posi)&&($med >$posi)){$other_peak++;}
				}
		}
		return($other_peak);
}

