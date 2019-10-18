#!/usr/bin/perl
use strict;
use warnings;

my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1/PCR_primer"){die "ERROR:working on wrong dir\n";}

my $fasta_file = "all_primer_homologous_removed/all_primer_homologous_removed_region.fasta";
my %fasta = &fasta2hash($fasta_file);

my $vcf_file = "manual_all_mutation_uchino_ym-check.txt";
open(VCF,"$vcf_file") or die "ERROR::cannot open $vcf_file";
open(OUT,">read_mutation.vcf");
my($bef_peak,$posi_adj) = ("",0);
<VCF>;
while(<VCF>){
		chomp;
		my @line =split(/\t/,);
		if(!defined$fasta{$line[0]}){next;}
		if($line[0] ne $bef_peak){
				$posi_adj = 0;
				$bef_peak = $line[0];
		}
		my ($sample,$peak_start,$peak_end);
		if($line[0] =~ /(\w+-\d\d)_chr[^:]+:(\d+)-(\d+)$/){
				($sample,$peak_start,$peak_end) = ($1, $2, $3);
		}else{die "what line $_\n";}
		if(($line[2] < $peak_start)||($line[2] > $peak_end)){next;}
		my $posi = $line[2] - $peak_start +1 + $posi_adj;
		my $seq = substr($fasta{$line[0]}{seq},$posi-1,length($line[4]));
		if($seq ne $line[4]){die "pick error $_ : $seq : $posi\n";}
		print OUT "$line[0]\t$posi\t$line[4]\t$line[3]\tORIGIN=$line[1]:$line[2]:$line[3]:$line[4]\n";
		if(length($line[3]) != length($line[4])){
				$posi_adj += length($line[4])-length($line[3]);
		}
}
close VCF;
close OUT;


#------------------------------------------------------------------------------
sub fasta2hash ( $ ){
  my ($file,$key,$value,$comment);
  my (%fasta_hash);
  $file=$_[0];
  if ($file =~ /\.gz$/) {
   open (IN,"gunzip -c $file |") || die "problem with $file\n";
  } else {
   open (IN,$file) || die "problem with $file\n";
  }
  while (<IN>){
    chomp;
    if (/^>(\S+)\s*(.*)/)
     {
	 #my @rec = split(m/\|/,$1);
	# while(scalar(@rec)>0)
	# {
      	#$key = pop(@rec);
	#	last if(length($key)>3);
	# }
      $key = $1;
      $comment = $2;
      if ($key =~ /^.+\|.+\|.+\|(.+)\|$/) {
	$comment = "$key $comment";
      	$key = $1;
      }
      $fasta_hash{$key}{"comment"}=$comment;
		#if(defined($fasta_hash{$key}{"sequence"}) || $fasta_hash{$key}{"sequence"} ne '')
		#{
		#	print STDERR "sequence $key already exists\n";
			$fasta_hash{$key}{"seq"}='';
		#}
     } #if (/^>(\w)$/)
    else
     {
      $key || die "File $file is not a fasta file!\n$key\n$_\n";
      s/\s+//g;
	  $_ =uc $_;
      $fasta_hash{$key}{"seq"}.=$_;
     } #else
   } #while (<IN>)
  close IN;
  return (%fasta_hash);
} #fasta2hash ( $ )
