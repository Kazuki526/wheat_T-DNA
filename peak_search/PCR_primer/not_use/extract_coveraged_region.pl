#!/usr/bin/perl
use warnings;
use strict;


my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}


open(PRIMER,"$ENV{HOME}/Dropbox/cooperative/plant_gen/T-DNA/count_all_read/primer.txt") or die "ERROR::cannot open primer.txt\n";
open(OUT,">$ENV{HOME}/Dropbox/cooperative/plant_gen/T-DNA/count_all_read/primer_comp.txt");
my $header = <PRIMER>;chomp $header;
print OUT "$header\n";
my %col = &header2hash($header);
while(<PRIMER>){
		chomp;
		my @line = split(/\t/,);
		my($start,$end);
		if(defined($line[$col{primer_end}])){($start,$end)=($line[$col{primer_start}],$line[$col{primer_end}]);
		}elsif($line[$col{side}] eq "left"){($start,$end)=($line[$col{junction}]-800,$line[$col{junction}]);
		}elsif($line[$col{side}] eq "right"){($start,$end)=($line[$col{junction}],$line[$col{junction}]+800);
		}else{die "ERROR::what line \'side\' colum is error??\n";}
		my %depth = &count_all_dpeth($line[$col{max_sample}],$line[$col{chr}],$start,$end);
		if($line[$col{side}] eq "left"){
				for(my $posi=$end;$posi>=$start;$posi--){
						if(!defined($depth{$posi})){$start=$posi+1;last;
						}elsif($depth{$posi}<10){$start=$posi+1;last;}
				}
				print OUT join("\t",@line[0..$col{junction}])."\t$start\t$end\n";
		}else{
				for(my $posi=$start;$posi<=$end;$posi++){
						if(!defined($depth{$posi})){$end=$posi-1;last;
						}elsif($depth{$posi}<10){$end=$posi-1;last;}
				}
				print OUT join("\t",@line[0..$col{junction}])."\t$start\t$end\n";
		}
}
close PRIMER;
close OUT;
				


sub count_all_dpeth( $ $ $ $){
		my($sample,$chr,$start,$end) = @_;
		my $bam_file = "$sample/$sample.bam";
		open(SAM,"samtools view $bam_file $chr:$start-$end|") or die "ERROR:: cannot open $bam_file";
		my %all_depth=();
		while(<SAM>){
				chomp;
				my @line = split(/\t/,);
				if($line[2] eq "*"){next;}
				my ($cigar,$posi) = ($line[5],$line[3]);
				while($cigar=~/./){
						$cigar =~ s/^(\d+)([HSIDM])//;
						my($length,$operator)=($1,$2,);
						if($operator =~ /[IHS]/){next;}
						while($length>0){
								$all_depth{$posi}++;
								$posi++;$length--;
						}
				}
		}
		close SAM;
		return(%all_depth);
}

sub header2hash ( $ ){
		my $header = $_[0];
		my @colm = split(/\t/,$header);
		my %out = ();
		for(my $i=0; $i < @colm; $i++){
				$out{$colm[$i]}=$i;
		}
		return(%out);
}
