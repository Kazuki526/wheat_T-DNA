#!/usr/bin/perl
use warnings;
use strict;

my $pwd =`pwd`;chomp $pwd;
if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}
my $genome_file = "../genome/all.fa";
-e $genome_file or die "ERROR::not exist $genome_file\n";

my $fasta = "PCR_primer/all_primer_homologous_removed/all_primer_homologous_removed_region.fasta";
my $masked = "PCR_primer/all_primer_homologous_removed/all_primer_homologous_removed_region.fasta.masked";
my $blast = "PCR_primer/all_primer_homologous_removed/all_primer_region_blast.tsv";

my %masked = &fasta2hash($masked);

my %blast =();
open(BLAST,"$blast") or die "ERROR::$blast cannot oepnt\n";
while(<BLAST>){
		if($_=~/^#/){next;}
		chomp;
		my @line =split(/\t/,);
		my ($sample,$chr,$start,$end);
		if($line[0] =~ /^(\w+-\d\d)_(chr[^:]+):(\d+)-(\d+)$/){
				($sample,$chr,$start,$end) = ($1,$2,$3,$4);
		}else{die "ERROR::$line[0] not match id regex\n";}
		my($blast_st,$blast_ed,$strand)=($line[8],$line[9],"+");
		if($blast_st>$blast_ed){($blast_st,$blast_ed,$strand)=($line[9],$line[8],"-");}
		if(($line[1] eq $chr)&&($blast_st >=$start)&&($blast_ed <=$end)){ #query region
				$blast{$line[0]} .= "self;";
		}else{$blast{$line[0]}.="$line[1]:$blast_st-$blast_ed$strand;";
		}
}
close BLAST;

my $align_dir = "PCR_primer/covered_align";
mkdir $align_dir;
my %out = ();
my @align_id=();
foreach my $id  (sort keys %masked){
		my $seq = $masked{$id}{seq};
		my $noN = $seq;
		$noN =~ s/N//g;
		$out{$id}=1- length($noN)/length($seq);
		if(!defined $blast{$id}){
				$out{$id} .= "\t0\tNA\tNA\n";
		}else{
				my @homolog = split(/;/,$blast{$id});
				if(scalar(@homolog) == ($blast{$id} =~ tr/self/self/)/4){$out{$id} .= "\t1\t".length($noN)."\t0\n";
				}else{
						push(@align_id, $id);
						open(OUT,">$align_dir/$id.fasta");
						print OUT ">$id\n$seq\n";
						close OUT;
						my $homolog_n=1;
						foreach my $homolog(@homolog){
								if($homolog eq "self"){next;}
								my $strand;
								if($homolog =~ /([-+])$/){$strand=$1;}
								$homolog =~ s/[+-]$//;
								if($strand eq "+"){
										`samtools faidx $genome_file $homolog >>$align_dir/$id.fasta`;
								}else{
										`samtools faidx $genome_file -i $homolog >>$align_dir/$id.fasta`;
								}
								$homolog_n++;
						}
						$out{$id}.="\t$homolog_n";
				}
		}
}

#align and count specific nucl
print "align and count specific nucl\n";
foreach my $id (@align_id){
		my $align_file = "$align_dir/$id.align.fasta";
		`mafft --auto --quiet $align_dir/$id.fasta > $align_file`;
		if(-e "$align_dir/$id.manual_align.fasta"){$align_file = "$align_dir/$id.manual_align.fasta";}
		my %align = &fasta2hash($align_file);
		my $align_n =scalar(keys %align);
		if(!defined$align{$id}){die "ERROR::there are not $id sequence at line 64\n";}
		my($no_homolog_sum,$specific_sum)=(0,0);
		for(my $posi=0; $posi<length($align{$id}{seq}); $posi++){
				my ($focal_nucl,$no_homolog,$specific) = (substr($align{$id}{seq},$posi,1),1,1);
				if($focal_nucl eq "N"){next;}
				foreach my $region(keys %align){
						if($region eq $id){next;}
						my $nucl = substr($align{$region}{seq},$posi,1);
						if($nucl ne $focal_nucl){
								if($nucl eq "-"){$no_homolog++;}
								$specific++;
						}
				}
				if($no_homolog == $align_n){$no_homolog_sum++;
				}elsif($specific == $align_n){$specific_sum++;}
		}
		$out{$id}.="\t$no_homolog_sum\t$specific_sum\n";
		print "done $id\n";
}

#print out
open(OUT,">PCR_primer/covered_region_info.tsv");
print OUT "primer_region_id\trepeat_ratio\tblast_hit\tno_homologous_nucl\tspecific_nucl\n";
foreach my $id (sort keys %masked){
		if(!defined $masked{$id}){die "ERROR::not defined $id at print out block\n";}
		print OUT "$id\t$out{$id}";
}
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
