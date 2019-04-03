#!/usr/bin/perl
use strict;
use warnings;

my $kmer = 25;

my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}

my $fasta = "../genome/all.fa";

my $peak_file="all_peak_table_with_p.tsv";
my %peak_region=();
my %pvalue=();
open(IN,"$peak_file") or die "ERROR::cannot open $peak_file\n";
<IN>; #header line
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		if($line[29] > 0.00001){last;}
		$peak_region{$line[3]} .= "$line[0]:$line[1]-$line[2]\t";
		$pvalue{"$line[0]:$line[1]-$line[2]"}=$line[29];
}
close IN;

my $fadir="fasta_dir_signif";
if(!-e $fadir){mkdir $fadir;}
my @sample=();
foreach my $sample (keys %peak_region){
		my @region = split(/\t/,$peak_region{$sample});
		my $out_fa = "$fadir/$sample"."_peak.fa";
		foreach my $region (@region){
				`samtools faidx $fasta $region >>$out_fa`;
		}
}

# pick homologous peak region
my $homologdir= "homologous_peak_list_signif";
if(!-e $homologdir){mkdir $homologdir;}
foreach my $sample (keys %peak_region){
		my $fa = "$fadir/$sample"."_peak.fa";
		my %peak_fa = &fasta2hash($fa);
		my $out = "$homologdir/$sample.tsv";
		open(OUT,">$out") or die "ERROR::cannot open $out\n";
		print OUT "k-mer k=$kmer\n";
		foreach my $db ( keys %peak_fa){
				print OUT "##### from $db\tp_value = $pvalue{$db}\n";
				my %homolog=();
				for(my $dbposi=0; $dbposi <= length($peak_fa{$db}{sequence})-$kmer+1; $dbposi++){
						my $seq = substr($peak_fa{$db}{sequence},$dbposi,$kmer);
						my $reverse = $seq;
						$reverse =~ tr/AaTtGgCc/TtAaCcGg/;
						$reverse = reverse($reverse);
						foreach my $query (keys %peak_fa){
								if($db eq $query){next;}
								if($peak_fa{$query}{sequence} =~ /$seq/){$homolog{$query}.="$dbposi,";}
								if($peak_fa{$query}{sequence} =~ /$reverse/){$homolog{$query}.="-$dbposi,";}
						}
				}
				foreach my $query( sort keys %homolog){
						my $count=0;
						$count++ while($homolog{$query} =~ m/,/g);
						print OUT "$db\t$query\t$pvalue{$query}\t$homolog{$query}\t$count\n";
				}
				print OUT "\n";
		}
		close OUT;
}
exit;
#------------------------------------------------------------------------------
sub fasta2hash ( $ )
{
		my ($file,$key,$value,$comment);
		my (%fasta_hash);
		$file=$_[0];
		if ($file =~ /\.gz$/) {
				open (IN,"gunzip -c $file |") || die "problem with $file\n";
		} else {
				open (IN,$file) || die "problem with $file\n";
		}
		while (<IN>)
		{
				chomp;
				if (/^>(\S+)\s*(.*)/)
				{
						$key = $1;
						$comment = $2;
						if ($key =~ /^.+\|.+\|.+\|(.+)\|$/) {
								$comment = "$key $comment";
								$key = $1;
						} 
						$fasta_hash{$key}{"comment"}=$comment;
		$fasta_hash{$key}{"sequence"}='';
				} #if (/^>(\w)$/)
				else 
				{
						$key || die "File $file is not a fasta file!\n$key\n$_\n";
						s/\s+//g;
						$fasta_hash{$key}{"sequence"}.=$_;
				} #else 
		} #while (<IN>)
		close IN;
		return (%fasta_hash); 
} #fasta2hash ( $ )

