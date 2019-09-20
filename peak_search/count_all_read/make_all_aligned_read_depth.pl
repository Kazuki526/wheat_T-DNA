#!/usr/bin/perl
use warnings;
use strict;

my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}

my @id_num = qw(01 02 03 04 05 06 07 08 09 10 11 12);
foreach my $num (@id_num){
		my $sample = "GFP-$num";
		print "start $sample =>";
		my $bam = "$sample/$sample.bam";
		if(!-e "count_all_read/$sample"){mkdir "count_all_read/$sample";}
		my $out = "count_all_read/$sample/$sample"."_depth.tsv.gz";
		&print_all_depth($bam,$out);
		print " end\n";
}
=pod
foreach my $num (@id_num){
		my $sample = "mCherry-$num";
		print "start $sample =>";
		my $bam = "$sample/$sample.bam";
		if(!-e "count_all_read/$sample"){mkdir "count_all_read/$sample";}
		my $out = "count_all_read/$sample/$sample"."_depth.tsv.gz";
		&print_all_depth($bam,$out);
		print " end\n";
}
=cut

#====================================================================================================================-
sub print_all_depth( $ $ ){
		my ($bam_file, $out_file) = @_;
		open(SAM,"samtools view $bam_file|") or die "ERROR:: cannot open $bam_file\n";
		open(OUT,"|gzip -c >$out_file") or die "ERROR:: cannot open $out_file\n";
		my ($chr,$window) = ("",0);
		my %all_read = ();
		my %chimeric = ();
		my %paired = ();
		while(<SAM>){
				chomp;
				my @line = split(/\t/,);
#				if($chr ne $line[2]){print "start $line[2]\n";}
				if($window != int($line[3]/25)){
						while(($window != int($line[3]/25))&&(scalar(keys%all_read)!=0)){
								my ($all,$chim,$pair,$start,$end)=(0,0,0,$window*25,$window*25+24);
								for(my$posi = $start;$posi <= $end;$posi++){
										if(defined $all_read{$posi}){$all +=$all_read{$posi};delete$all_read{$posi};}
										if(defined $chimeric{$posi}){$chim+=$chimeric{$posi};delete$chimeric{$posi};}
										if(defined $paired{$posi}  ){$pair+=$paired{$posi};  delete$paired{$posi};  }
								}
								if($all>0){
										if($chr eq ""){$chr=$line[2];}
										print OUT "$chr\t$start\t$end\t$all\t$chim\t$pair\n";
								}
								$window++;
						}
						($chr,$window)=($line[2],int($line[3]/25));
				}
				if($chr eq "*"){last;}
				if($line[5] eq "*"){next;}
				my ($cigar,$posi) = ($line[5],$line[3]);
				my ($chimeric_focal,$paired_focal) = (0,0);
				if($line[6] =~ /plasmid/){$paired_focal=1;}
				if($_ =~ /SA:Z:[GFPmChery]+_plasmid,/){$chimeric_focal=1;}
				while($cigar =~ /./){
						$cigar =~ s/^(\d+)([HSIDM])//;
						my($length,$operator)=($1,$2);
						if($operator =~ /[IHS]/){next;}
						while($length>0){
								$all_read{$posi}++;
								if($chimeric_focal){$chimeric{$posi}++;}
								if($paired_focal){$paired{$posi}++;}
								$posi++;$length--;
						}
				}
		}
		close SAM;
		close OUT;
}



