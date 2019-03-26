#!/usr/bin/perl
use warnings;
use strict;

my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}

my $peak_depth = 100;
my @sample_num=qw(01 02 03 04 05 06 07 08 09 10 11 12);
my %peak_list=();
#GFP
foreach my $num(@sample_num){
		my $sample_id ="GFP-$num";
		my @peak_ls = &peak_connect($sample_id);
		foreach my $peak(@peak_ls){
				my($chr,$start,$end);
				if($peak =~ /^(chr.+):(\d+)-(\d+)$/){($chr,$start,$end)=($1,$2,$3);}else{die "ERROR::what peak ?? $peak \n";}
				my $existed =0;
				if(!defined $peak_list{$chr}){
						$peak_list{$chr}{$start}{end}=$end;
						$peak_list{$chr}{$start}{ids}=$sample_id;
						next;
				}
				foreach my $other_start(sort{$a<=>$b} keys%{$peak_list{$chr}}){
						my ($other_end,$others) = ($peak_list{$chr}{$other_start}{end},$peak_list{$chr}{$other_start}{ids});
						if(!defined $other_end){die"$other_start:at $start\n";}
						if((($other_start<=$start)&&($other_end >=$start))||
						   (($other_start<=$end)&&($other_end >=$end))){
								$existed++;
								$others.="\t$sample_id";
								if($start < $other_start){
										delete $peak_list{$chr}{$other_start};
										$peak_list{$chr}{$start}{end}=$other_end;
										$other_start=$start;
								}
								if($end > $other_end){$peak_list{$chr}{$other_start}{end}=$end;}
								$peak_list{$chr}{$other_start}{ids}=$others;
						}
				}
				if($existed==0){
						$peak_list{$chr}{$start}{end}=$end;
						$peak_list{$chr}{$start}{ids}=$sample_id;
				}
		}
		foreach my $chr(keys%peak_list){
				foreach my $start(sort{$a<=>$b}keys%{$peak_list{$chr}}){
						if(!defined $peak_list{$chr}{$start}){next;}
						foreach my $other_start(sort{$a<=>$b}keys%{$peak_list{$chr}}){
								if(($start >= $other_start)||($peak_list{$chr}{$start}{end} <$other_start)){next;}
								if($peak_list{$chr}{$start}{end} < $peak_list{$chr}{$other_start}{end}){
										$peak_list{$chr}{$start}{end} = $peak_list{$chr}{$other_start}{end};
								}
								my($ids,$other_ids)=($peak_list{$chr}{$other_start}{ids},$peak_list{$chr}{$other_start}{ids});
								my @other_ids=split(/\t/,$other_ids);
								foreach my $oid(@other_ids){
										if($ids !~ /$oid/){$ids.="\t$oid";}
								}
								$peak_list{$chr}{$start}{ids}=$ids;
								delete $peak_list{$chr}{$other_start};
						}
				}
		}
		print "done $sample_id\n";
}
#mCherry
foreach my $num(@sample_num){
		my $sample_id ="mCherry-$num";
		my @peak_ls = &peak_connect($sample_id);
		foreach my $peak(@peak_ls){
				my($chr,$start,$end);
				if($peak =~ /^(chr.+):(\d+)-(\d+)$/){($chr,$start,$end)=($1,$2,$3);}else{die "ERROR::what peak ?? $peak \n";}
				my $existed =0;
				if(!defined $peak_list{$chr}){
						$peak_list{$chr}{$start}{end}=$end;
						$peak_list{$chr}{$start}{ids}=$sample_id;
						next;
				}
				foreach my $other_start(sort{$a<=>$b} keys%{$peak_list{$chr}}){
						my ($other_end,$others) = ($peak_list{$chr}{$other_start}{end},$peak_list{$chr}{$other_start}{ids});
						if((($other_start<=$start)&&($other_end >=$start))||
						   (($other_start<=$end)&&($other_end >=$end))){
								$existed++;
								$others.="\t$sample_id";
								if($start < $other_start){
										delete $peak_list{$chr}{$other_start};
										$peak_list{$chr}{$start}{end}=$other_end;
										$other_start=$start;
								}
								if($end > $other_end){$peak_list{$chr}{$other_start}{end}=$end;}
								$peak_list{$chr}{$other_start}{ids}=$others;
						}
				}
				if($existed==0){
						$peak_list{$chr}{$start}{end}=$end;
						$peak_list{$chr}{$start}{ids}=$sample_id;
				}
		}
		foreach my $chr(keys%peak_list){
				foreach my $start(sort{$a<=>$b}keys%{$peak_list{$chr}}){
						if(!defined $peak_list{$chr}{$start}){next;}
						foreach my $other_start(sort{$a<=>$b}keys%{$peak_list{$chr}}){
								if(($start >= $other_start)||($peak_list{$chr}{$start}{end} <$other_start)){next;}
								if($peak_list{$chr}{$start}{end} < $peak_list{$chr}{$other_start}{end}){
										$peak_list{$chr}{$start}{end} = $peak_list{$chr}{$other_start}{end};
								}
								my($ids,$other_ids)=($peak_list{$chr}{$other_start}{ids},$peak_list{$chr}{$other_start}{ids});
								my @other_ids=split(/\t/,$other_ids);
								foreach my $oid(@other_ids){
										if($ids !~ /$oid/){$ids.="\t$oid";}
								}
								$peak_list{$chr}{$start}{ids}=$ids;
								delete $peak_list{$chr}{$other_start};
						}
				}
		}
		print "done $sample_id\n";
}

open(OUT,">$peak_depth"."_peak_list.tsv");
print OUT "chr\tstart\tend\tsamples\tpeak_sample_num\n";

#間が100bp以内のpeakは連結
my ($bef_start,$bef_end,$bef_ids)=(0,0,"");
foreach my $chr(sort keys %peak_list){
		($bef_start,$bef_end,$bef_ids)=(0,0,"");
		foreach my $start(sort{$a<=>$b} keys %{$peak_list{$chr}}){
				my($end,$ids)=($peak_list{$chr}{$start}{end},$peak_list{$chr}{$start}{ids});
				if($start > $bef_end+101){
						($bef_start,$bef_end,$bef_ids)=($start,$end,$ids);
						next;
				}else{
						my @ids =split(/\t/,$ids);
						foreach my $oids (@ids){
								if($bef_ids !~/$oids/){$bef_ids.="\t$oids";}
						}
						$peak_list{$chr}{$bef_start}{end}=$end;
						$peak_list{$chr}{$bef_start}{ids}=$bef_ids;
						delete $peak_list{$chr}{$start};
						$bef_end=$end;
				}
		}
}
#print out
foreach my $chr(sort keys %peak_list){
		foreach my $start(sort{$a<=>$b} keys %{$peak_list{$chr}}){
				my($end,$ids)=($peak_list{$chr}{$start}{end},$peak_list{$chr}{$start}{ids});
				my @ids = sort(split(/\t/,$ids));
				my $num=scalar(@ids);
				print OUT "$chr\t$start\t$end\t".join(",",@ids)."\t$num\n";
		}
}
close OUT;

exit;



sub peak_connect( $ ){
		my $sample=$_[0];
		my $peak_file="$sample/$sample"."_allpeak.tsv.gz";
		my @peak_list =();
		open(IN,"gunzip -c $peak_file|");
#		open(OUT,">$sample/$sample"."_allconnect.tsv");
		<IN>;#print OUT "chr\tstart\tend\taverage_depth\n";
		my($bef_end,$chr,$start,$average_depth)=(0,"",0,0);
		while(<IN>){
				chomp;
				my @line = split(/\t/,);
				if(($line[0] eq "GFP_plasmid")||($line[0] eq "mCherry_plasmid")){next;}
				if(($chr ne $line[0]) || ($bef_end+1 != $line[1])){
						my $window_num = ($bef_end - $start +1)/25;
						$average_depth/=$window_num;
						if(($chr ne "")&&($average_depth > $peak_depth)){
#								print OUT "$chr\t$start\t$bef_end\t$average_depth\n";
								push(@peak_list,"$chr:$start-$bef_end");
						}

						($chr,$start,$bef_end,$average_depth) = ($line[0],$line[1],$line[2],$line[4]);
				}else{
						$bef_end=$line[2];
						$average_depth += $line[4];
				}
		}
		my $window_num = ($bef_end - $start +1)/25;
		$average_depth/=$window_num;
#		print OUT "$chr\t$start\t$bef_end\t$average_depth\n";
		push(@peak_list,"$chr:$start-$bef_end");
		return(@peak_list);
		close IN;
#		close OUT;
}
		


