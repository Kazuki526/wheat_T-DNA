#!/usr/bin/perl
use strict;
use warnings;

my $pwd = `pwd`;chomp $pwd;
if($pwd !~ /T-DNA_search$/){die "ERROR:please mapping on T-DNA_search dir\n";}

# refference is merge T-DNA sequence (GFP and mCherry) and 161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta"
my $ref="genome/all.fa";
if(!-e $ref){die "ERROR::$ref isnot exist";}
my @gfp_list=();
my @mcherry_list=();
for(my $i=7;12>=$i;$i++){
		if(length($i)>1){push(@gfp_list,"GFP-$i");
		}else{push(@gfp_list,"GFP-0$i");}
		if(length($i)>1){push(@mcherry_list,"mCherry-$i");
		}else{push(@mcherry_list,"mCherry-0$i");}
}
foreach my $sample(@gfp_list){
		my $adapter = "HN00100341_hdd1/Adapter/Adapter_$sample.fasta";
		my $fastq1 = "HN00100341_hdd1/$sample/$sample"."_R1.fastq.gz";
		my $fastq2 = "HN00100341_hdd1/$sample/$sample"."_R2.fastq.gz";
		my ($trimed1,$trimed2) = ("HN00100341_hdd1/$sample/trimed1.fastq","HN00100341_hdd1/$sample/trimed2.fastq");
		if((!-e $adapter)||(!-e $fastq1)||(!-e $fastq2)){die "ERROR::GFP adapter or fastq.gz was not exit!\n";}
		system("trimmomatic PE -phred33 $fastq1 $fastq2 $trimed1 /dev/null $trimed2 /dev/null ILLUMINACLIP:$adapter:2:20:10 SLIDINGWINDOW:4:25 LEADING:20 TRAILING:20 MINLEN:80");
		my $out_bam = "HN00100341_hdd1/$sample/$sample.bam";
		system("bwa mem -M $ref $trimed1 $trimed2|samtools sort >$out_bam");
		system("samtools index $out_bam");
}
=pod
foreach my $sample(@mcherry_list){
		my $adapter = "HN00100341_hdd1/Adapter/Adapter_$sample.fasta";
		my $fastq1 = "HN00100341_hdd1/$sample/$sample"."_R1.fastq.gz";
		my $fastq2 = "HN00100341_hdd1/$sample/$sample"."_R2.fastq.gz";
		my ($trimed1,$trimed2) = ("HN00100341_hdd1/$sample/trimed1.fastq","HN00100341_hdd1/$sample/trimed2.fastq");
		if((!-e $adapter)||(!-e $fastq1)||(!-e $fastq2)){die "ERROR::mCherry adapter or fastq.gz was not exit!\n";}
		system("trimmomatic PE -phred33 $fastq1 $fastq2 $trimed1 /dev/null $trimed2 /dev/null ILLUMINACLIP:$adapter:2:20:10 SLIDINGWINDOW:4:25 LEADING:20 TRAILING:20 MINLEN:80");
		my $out_bam = "HN00100341_hdd1/$sample/$sample.bam";
		system("bwa mem -M $ref $trimed1 $trimed2|samtools sort >$out_bam");
		system("samtools index $out_bam");
}
=cut
