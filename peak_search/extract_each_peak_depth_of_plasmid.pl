#!/usr/bin/perl
use warnings;
use strict;

my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}

my $bed="100_peak.bed";





my @sample_num=qw(01 02 03 04 05 06 07 08 09 10 11 12);
open(OUTP,">paired_plasmid_table_tidy.tsv");
open(OUTC,">chimeric_align_talbe_tidy.tsv");

#GFP
foreach my $num(@sample_num){
		my $sample_id ="GFP-$num";
		my $bam = "$sample_id/$sample_id.bam";
		my ($pbam,$cbam) = &extract_plasmid_related_read($bam);
		open(COV,"samtools bedcov $bed $pbam|");
		while(<COV>){
				print OUTP "$sample_id\t$_";
		}
		close COV;
		open(COV,"samtools bedcov $bed $cbam|");
		while(<COV>){
				print OUTC "$sample_id\t$_";
		}
		close COV;
		print "done $sample_id\n";
}
#mCherry
foreach my $num(@sample_num){
		my $sample_id ="mCherry-$num";
		my $bam = "$sample_id/$sample_id.bam";
		my ($pbam,$cbam) = &extract_plasmid_related_read($bam);
		open(COV,"samtools bedcov $bed $pbam|");
		while(<COV>){
				print OUTP "$sample_id\t$_";
		}
		close COV;
		open(COV,"samtools bedcov $bed $cbam|");
		while(<COV>){
				print OUTC "$sample_id\t$_";
		}
		close COV;
		print "done $sample_id\n";
		unlink $pbam;
		unlink $cbam;
		unlink "$pbam.bai";
		unlink "$cbam.bai";
}
close OUTP;
close OUTC;
exit;


sub extract_plasmid_related_read( $ ){
		my ($bam)=$_[0];
		my ($psam,$csam)=("paired_plasmid.sam","chimeric_align.sam");
		my ($pbam,$cbam)=("paired_plasmid.bam","chimeric_align.bam");
		open(BAM,"samtools view -h $bam -L $bed|");
		open(PS,">$psam");
		open(CS,">$csam");
		while(<BAM>){
				if($_ =~ /^\@/){print PS "$_";print CS "$_";next;}
				my @line = split(/\t/,);
				if(($line[6] eq "GFP_plasmid")||($line[6] eq "mCherry_plasmid")){
						print PS "$_";
				}
				if($_ =~ /SA:Z:(.+);\s/){
						my@chimeric = split(/,/,$1);
						if(($chimeric[0] eq "GFP_plasmid")||($chimeric[0] eq "mCherry_plasmid")){
								print CS "$_";
						}
				}
		}
		close BAM;
		close PS;
		close CS;
		`samtools view -b $psam -o $pbam`;
		`samtools index $pbam`;
		`samtools view -b $csam -o $cbam`;
		`samtools index $cbam`;
		unlink $csam;
		unlink $psam;
		return($pbam,$cbam);
}


