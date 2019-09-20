#!/usr/bin/perl
use warnings;
use strict;
my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}

#make peak bed 
my @id_num = qw(01 02 03 04 05 06 07 08 09 10 11 12);
my $out = "count_all_read/all_peak.bed";
if(-e $out){unlink $out;}
foreach my $num (@id_num){
		my $sample = "GFP-$num";
		my $depth = "count_all_read/$sample/$sample"."_depth.tsv.gz";
		&extract_peak($depth,$out);
}
foreach my $num (@id_num){
		my $sample = "mCherry-$num";
		my $depth = "count_all_read/$sample/$sample"."_depth.tsv.gz";
		&extract_peak($depth,$out);
}

sub extract_peak( $ $ ){
		my ($in,$bed)=@_;
		open(IN,"gunzip -c $in|") or die "EROOR::cannot open $in\n";
		open(OUT,">>$bed");
		my ($chr,$start);
		while(<IN>){
				chomp;
				my @line = split(/\t/,);
				if($line[0] =~ /plasmid/){last;}
				if($line[3]>2500){
						$line[2]++;
						print OUT "$line[0]\t$line[1]\t$line[2]\n";
				}
		}
		close IN;
		close OUT;
}
print "printed out all peak bed\n";
# merge peak by bedtools 
my $sorted_bed = "count_all_read/all_peak_sorted.bed";
my $merged_bed_ = "count_all_read/all_peak_merged_.bed";
my $merged_bed = "count_all_read/all_peak_merged.bed";
`bedtools sort -i $out >$sorted_bed`;
system("bedtools merge -i $sorted_bed >$merged_bed_");

#merge peak interval is 100bp
open(BED,"$merged_bed_")or die "ERROR::cannot open $merged_bed_\n";
open(OUT,">$merged_bed") or die "ERROR::cannot open $merged_bed\n";
my ($chr,$start,$end)=("",0,0);
while(<BED>){
		chomp;
		my @line = split(/\t/,);
		if(($chr eq $line[0])&&($end+100 >= $line[1])){
				$end =$line[2];
		}else{
				if($chr ne ""){print OUT "$chr\t$start\t$end\n";}
				($chr,$start,$end)=($line[0],$line[1],$line[2]);
		}
}
close BED;
close OUT;

#make peak table tidy
my %peak =();
open(BED,"$merged_bed") or die "ERROR:cannot open $merged_bed\n";
while(<BED>){
		chomp;
		my @line = split(/\t/,);
		$peak{$line[0]}.="$line[1]-$line[2];";
}
close BED;

my $out_a = "count_all_read/all_peak_table_tidy.tsv";
my $out_c = "count_all_read/all_peak_table_chimeric_tidy.tsv";
my $out_p = "count_all_read/all_peak_table_paired_tidy.tsv";
if(-e $out_a){unlink $out_a;}
if(-e $out_c){unlink $out_c;}
if(-e $out_p){unlink $out_p;}

foreach my $num (@id_num){
		my $sample = "GFP-$num";
		print "start $sample =>";
		&print_peak_tidy($sample,$out_a,$out_c,$out_p);
		print " end\n";
}
@id_num = qw(01 02 03 04 05 06 07 08 09 10 11 12);
foreach my $num (@id_num){
		my $sample = "mCherry-$num";
		print "start $sample =>";
		&print_peak_tidy($sample,$out_a,$out_c,$out_p);
		print " end\n";
}

sub print_peak_tidy( $ $ $ $ ){
		my ($sample,$out_all,$out_chim,$out_pair) = @_;
		my $depth = "count_all_read/$sample/$sample"."_depth.tsv.gz";
		open(DEP,"gunzip -c $depth|");
		open(OUTA,">>$out_all");
		open(OUTC,">>$out_chim");
		open(OUTP,">>$out_pair");
		my %peak_depth =();
		while(<DEP>){
				chomp;
				my @line =split(/\t/,);
				my $focal = &peak_overlap($line[0],$line[1],$line[2]);
				if($focal){
						$peak_depth{$focal}{all}  +=$line[3];
						$peak_depth{$focal}{chim} +=$line[4];
						$peak_depth{$focal}{pair} +=$line[5];
				}
		}
		close DEP;
		foreach my $region (keys %peak_depth){
				my ($chr,$start,$end) = split(/:/,$region);
				print OUTA "$sample\t$chr\t$start\t$end\t$peak_depth{$region}{all}\n";
				print OUTC "$sample\t$chr\t$start\t$end\t$peak_depth{$region}{chim}\n";
				print OUTP "$sample\t$chr\t$start\t$end\t$peak_depth{$region}{pair}\n";
		}
		close OUTA;
		close OUTC;
		close OUTP;
}



sub peak_overlap( $ $ $ ){
		my ($chr,$start,$end)=@_;
		my $focal =0;
		if($chr =~ /plasmid/){
				return($focal);
		}else{
				my @region = split(/;/,$peak{$chr});
				foreach my $region (@region){
						my($peak_start,$peak_end)=split(/-/,$region);
						if(($peak_start <=$start)&&($peak_end >= $end)){
								$peak_end--;
								$focal = "$chr:$peak_start:$peak_end";
						}
				}
				return($focal);
		}
}
