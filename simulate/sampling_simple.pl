#!/usr/bin/perl
use warnings;
use strict;

my $t =time();
chdir "/Volumes/kazcancer/T-DNA_search";

# CS fasta file
my $cs_fa_file = "Chinese_Spring/all_chr_v1.0.fasta";
# each chr length file
my $chr_length_file = "Chinese_Spring/CSrefseq_detail.txt";
(-f $cs_fa_file and -f $chr_length_file) or die "ERROR::there are not exist Chinese spring fasta or chr length file\n";

my %num2chr = (
	1  => "chr1A", 2  => "chr2A", 3  => "chr3A", 4  => "chr4A", 5  => "chr5A", 6  => "chr6A", 7  => "chr7A",
	8  => "chr1B", 9  => "chr2B", 10 => "chr3B", 11 => "chr4B", 12 => "chr5B", 13 => "chr6B", 14 => "chr7B",
	15 => "chr1D", 16 => "chr2D", 17 => "chr3D", 18 => "chr4D", 19 => "chr5D", 20 => "chr6D", 21 => "chr7D");

#read chr length file
my %chr_length = ();
open(IN,"$chr_length_file");
while(<IN>){
		chomp;
		my @line = split(/\t/,);
		$chr_length{$line[0]} = $line[1];
}
close IN;

#check single_end dir exist
my $sample_dir = "simple_sampling";
if(! -d $sample_dir){mkdir $sample_dir;}

my @read_length = (); #100 ~ 300 by 10
for(my $i =10;30 >=$i;$i++){
		push(@read_length,$i*10);
}
foreach my $read_length(@read_length){
		my $out_dir = "$sample_dir/read$read_length";
		mkdir $out_dir;
		open(OUTFA, ">$out_dir/sample.fa");  # sequence not including N
		open(OUTFAN,">$out_dir/sampleN.fa"); # sequence containing N
		for(my $times = 1;1000000 >=$times;$times++){
				my $n=0;
				my $chr = $num2chr{int(rand(21))+1};
				my $start = int(rand($chr_length{$chr} - $read_length)) +1;
				my $end = $start + $read_length -1;
				my $samtools_out = `samtools faidx $cs_fa_file $chr:$start-$end`;
				$samtools_out =~ s/^>$chr:$start-$end//;
				$samtools_out =~ s/\n//g;
				my $seq = uc $samtools_out;
				if($seq =~ /N/){
						print OUTFAN ">$chr:$start:$end\n$seq\n";
				}else{
						print OUTFA ">$chr:$start:$end\n$seq\n";
				}
		}
		close OUTFA;
		close OUTFAN;
		# mapping 
		system("bwa mem -M $cs_fa_file $out_dir/sample.fa >$out_dir/out.bam");
		system("bwa mem -M $cs_fa_file $out_dir/sampleN.fa >$out_dir/outN.bam");
		print "finish read length = $read_length\t:time spent " . (time() - $t) . "sec\n";
}
