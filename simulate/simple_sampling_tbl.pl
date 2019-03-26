#!/usr/bin/perl
use warnings;
use strict;

my $t = time();
chdir "/Volumes/kazcancer/T-DNA_search";

my @nonN=qw(only_original perfect_hit multi_hit);
my @N   =qw(Nonly_original Nperfect_hit Nmulti_hit N_non_hit);
open(OUT,">simple_sampling/result_tbl.tsv");
print OUT "read_length\t" . join("\t",(@nonN,@N)) . "\n";
for(my $read_length=100;300>=$read_length;$read_length+=10){
		print "now doin $read_length !\n";
		my $sample_dir="simple_sampling/read$read_length/";
		my %result = ();
		foreach my$i(@nonN,@N){$result{$i}=0;}
# not containing N
		open(BAM ,"samtools view $sample_dir/out.bam|");
		while(<BAM>){
				chomp;
				my @line =split(/\t/,);
				if($line[$#line] !~ /^XA:/){
						$result{only_original}++;
				}else{
						my $xa_ =$line[$#line];
						$xa_ =~ s/XA:Z://;
						my @xa = split(/;/,$xa_);
						my $n_match=0;
						foreach my $xa (@xa){
								my @ali = split(/,/,$xa);
								if($ali[2] eq ($read_length ."M")){$n_match++;}
						}
						if($n_match >=2){$result{multi_hit}++;
						}else{$result{perfect_hit}++;}
				}
		}
		close BAM;
# containing N
		open(BAMN,"samtools view $sample_dir/outN.bam|");
		while(<BAMN>){
				chomp;
				my @line = split(/\t/,);
				if($line[1] == 4){$result{N_non_hit}++;next;}
				my ($chr,$start,$cigar)=($line[2],$line[3],$line[5]);
				if($cigar =~ /H/){next;} # NをHard clipとして扱わずalignしているパターン
				if($cigar =~ /^(\d+)S/){$start-=$1;}
				my$end=$start + $read_length - 1;
				if($line[0] ne "$chr:$start:$end"){
						if($line[$#line] =~/^XA:Z/){$result{Nmulti_hit}++;
						}else{$result{N_non_hit}++;}
				}else{
						if($line[$#line] !~ /^XA:Z:/){
								$result{Nonly_original}++;
						}else{
								my $xa_ =$line[$#line]; $xa_=~s/^XA:Z://;
								my @xa = split(/;/,$xa_);
								my $n_match=0;
								foreach my $xa (@xa){
										my @ali =split(/,/,$xa);
										my ($xa_strand,$xa_cigar)=($ali[2],"");
										if($ali[1] =~ /^([-+])/){$xa_strand=$1;}
										if($xa_strand eq "+"){
												if($xa_cigar eq $cigar){$n_match++}
										}else{
												if($xa_cigar =~/(\d+[MS]).*(\d+[MS])/){
														if("$2$1" eq $cigar){$n_match++;}
												}
										}
								}
								if($n_match >=2){$result{Nmulti_hit}++;
								}else{$result{Nperfect_hit}++;}
						}
				}
		}
		close BAMN;
		print OUT $read_length;
		foreach my $t(@nonN,@N){
				print OUT "\t$result{$t}";
		}
		print OUT "\n";
}
close OUT;
exit;

