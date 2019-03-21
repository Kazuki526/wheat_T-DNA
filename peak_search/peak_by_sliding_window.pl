#!/usr/bin/perl
use warnings;
use strict;

my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}

my @sample_num=qw(01 02 03 04 05 06 07 08 09 10 11 12);
#GFP
foreach my $num(@sample_num){
		my $sample_id ="GFP-$num";
		&main($sample_id);
		print "done $sample_id\n";
}
#mCherry
foreach my $num(@sample_num){
		my $sample_id ="mCherry-$num";
		&main($sample_id);
		print "done $sample_id\n";
}
exit;



sub main( $ ){
		my $sample = $_[0];
		my ($in_file,$out_file) = ("$sample/$sample.bam","$sample/$sample"."_allpeak.tsv");
		open(IN,"samtools depth $in_file|");
		open(OUT,"|gzip -c >$out_file");
		print OUT "chr\tstart\tend\tdepth\taverage_depth\n";
		my ($chr,$window_num,$posi_num) =("",0,0);
		my @window_depth =();
		while(<IN>){
				chomp;
				my @line = split(/\t/,);
				my($now_wn,$now_posi)=(int($line[1]/25),$line[1]%25);
				if(($chr eq $line[0])&&($window_num==$now_wn)){
						for(my$interval=$posi_num;$interval<$now_posi-1;$interval++){push(@window_depth,0);}
						push(@window_depth,$line[2]);
						$posi_num = $now_posi;
				}else{
				# calculate average depth & print out if peak of before window
						if(scalar(@window_depth) !=0){ 
								# push depth 0 to interval
								for(my$interval=$posi_num;$interval<24;$interval++){push(@window_depth,0);} 
								if(scalar(@window_depth) != 25){die "ERROR::window size is not 25 at $_\n@window_depth\n";}
								my $average_depth=0;
								foreach my $d(@window_depth){$average_depth += $d/25;}
								# print out if peak (average_depth>100)
								if($average_depth > 100){  
										my($start,$end)=($window_num*25,$window_num*25+24);
										print OUT "$chr\t$start\t$end\t".join(",",@window_depth)."\t$average_depth\n";
								}
								@window_depth=();
						}
						# next window
						for(my$interval=0;$interval<$now_posi;$interval++){push(@window_depth,0);}
						push(@window_depth,$line[2]);
						($chr,$window_num,$posi_num)=($line[0],$now_wn,$now_posi);
				}
		}
		# caldulate and print out last window
		for(my$interval=scalar(@window_depth); $interval<25; $interval++){push(@window_depth,0);} 
		if(scalar(@window_depth) != 25){die "ERROR::window size is not 25 at $_\n@window_depth\n";}
		my $average_depth=0;
		foreach my $d(@window_depth){$average_depth += $d/25;}
		# print out if peak (average_depth>100)
		if($average_depth > 100){  
				my($start,$end)=($window_num*25,$window_num*25+24);
				print OUT "$chr\t$start\t$end\t".join(",",@window_depth)."\t$average_depth\n";
		}

		close IN;
		close OUT;
}

