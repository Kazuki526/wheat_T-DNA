#!/usr/bin/perl
use warnings;
use strict;

my $pwd =`pwd`;chomp $pwd;
#if($pwd ne "/Volumes/kazcancer/T-DNA_search/HN00100341_hdd1"){die "ERROR:working on wrong dir\n";}
if($pwd ne "/Users/kaz/Dropbox/cooperative/plant_gen/T-DNA/forPCR"){die "ERROR:working on wrong dir\n";}

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
				$blast{$line[0]}{self} .= "$line[6]-$line[7];";
		}else{$blast{$line[0]}{homolog}.="$line[1]:$blast_st-$blast_ed$strand;";
		}
}
close BLAST;

my $primer_file = "PCR_primer/primer_comp.txt";-e $primer_file or die "ERROR::$primer_file is not exist\n";
my %primer_comp = ();
open(PR,"$primer_file");
<PR>;
while(<PR>){
		chomp;
		my @line = split(/\t/,);
		my $primer_id = "$line[4]_$line[1]:$line[8]-$line[9]";
		my $peak_id = "$line[4]_$line[1]:$line[2]-$line[3]";
		$primer_comp{$primer_id}{side} = $line[5];
		$primer_comp{$primer_id}{peak_type} = $line[0];
		$primer_comp{$primer_id}{peak_id} = $peak_id;
}
close PR;

my $read_mutation_file = "PCR_primer/read_mutation.vcf"; -e $read_mutation_file or die "ERROR::$read_mutation_file is not exist\n";
my %read_mutation =();
open(VCF,"$read_mutation_file");
while(<VCF>){
		chomp;
		my @line = split(/\t/,);
		my $origin;
		if($line[4] =~ /^ORIGIN=(.+)$/){$origin = $1;}
		my $end=$line[1]+length($line[2])-1;
		$read_mutation{$line[0]}{mut}.="$line[1]-$end:$line[2]>$line[3];";
		$read_mutation{$line[0]}{origin}.="$origin;";
}
close VCF;

my $primer_dir = "PCR_primer/primer_tbl";
mkdir $primer_dir;
open(TBL,">$primer_dir/all_peak_primer_info.tsv");
print TBL "peak_id\tprimer_id\tcomment\tN(repeat)_freq\tblast_hit_n\tprimer_number\n";
open(TWO,">$primer_dir/two-sided_peak_primer.tsv");
print TWO "peak_id\tprimer_id\tchr\tstart_posi\tend_posi\tlength\tstrand\tseq\tread_mutation\tread_mutation_from_refference\n";
open(ONE,">$primer_dir/one-sided_peak_primer.tsv");
print ONE "peak_id\tprimer_id\tchr\tstart_posi\tend_posi\tlength\tstrand\tseq\tread_mutation\tread_mutation_from_refference\n";

open(INDEL,">$primer_dir/indel_info.tsv");

foreach my $primer_id(sort keys %primer_comp){
		if(!defined$masked{$primer_id}){
				print TBL "$primer_comp{$primer_id}{peak_id}\t$primer_id\tpart_of_other_peak\n";
		}else{
				my $seq = $masked{$primer_id}{seq};
				my $noN = $seq; $noN =~ s/N//g;
				my $Nfreq = length($noN)/length($seq);
				if(!defined $blast{$primer_id}){
						print TBL "$primer_comp{$primer_id}{peak_id}\t$primer_id\tno_blast_hit\t$Nfreq\t0\n";
				}else{
						if(!defined $blast{$primer_id}{homolog}){ #no homolog
								my @blast_covered = split(/;/,$blast{$primer_id}{self});
								my $out = "";
								foreach my $cover_region(@blast_covered){
										$out .= &get_primer_tbl($primer_id,$cover_region);
								}
								if($primer_comp{$primer_id}{peak_type} eq "both"){
										print TWO "$out";
								}else{
										print ONE "$out";
								}
								my $count = $out =~ tr/\n/\n/;
								print TBL "$primer_comp{$primer_id}{peak_id}\t$primer_id\t\t$Nfreq\t$count\n"
						}else{ #have some homolog
								my @blast_covered =();
								if(defined $blast{$primer_id}{self}){@blast_covered=split(/;/,$blast{$primer_id}{self});}
								my $align_file="PCR_primer/covered_align/$primer_id.manual_align.fasta";
								if(!-e $align_file){$align_file="PCR_primer/covered_align/$primer_id.align.fasta";}
								my %specific_mutation = &pick_specific_mutation($primer_id,$align_file);
								my $out="";
								foreach my $cover_region(@blast_covered){
										$out .= &get_specific_primer_tbl($primer_id,$cover_region,\%specific_mutation);
								}
								if($primer_comp{$primer_id}{peak_type} eq "both"){
										print TWO "$out";
								}else{
										print ONE "$out";
								}
								my $count = $out =~ tr/\n/\n/;
								print TBL "$primer_comp{$primer_id}{peak_id}\t$primer_id\t\t$Nfreq\t$count\n"
						}
				}
		}
}
close TBL;
close TWO;
close ONE;
close INDEL;



#------------------------------------------------------------------------------
sub get_specific_primer_tbl( $ $ $ ){
		my ($primer_id,$region)=($_[0],$_[1]);
		my %specific_mutation = %{$_[2]};
		my ($sample,$chr,$primer_start,$primer_end,$region_start,$region_end);
		if($primer_id =~ /^(\w+-\d\d)_(chr[^:]+):(\d+)-(\d+)$/){
				($sample,$chr,$primer_start,$primer_end) = ($1, $2, $3, $4);
		}
		if($region =~ /^(\d+)-(\d+)$/){
				($region_start,$region_end) = ($1, $2);
		}
		my $read_mutation = "";
		my @read_mutation=();
		my @read_mutation_orig=();
		if(defined$read_mutation{$primer_id}{mut}){$read_mutation = $read_mutation{$primer_id}{mut};
				@read_mutation = split(/;/,$read_mutation{$primer_id}{mut});
				@read_mutation_orig = split(/;/,$read_mutation{$primer_id}{origin});
		}
		my %ref_posi=&position_translation_ref($primer_start,$masked{$primer_id}{seq},$read_mutation);
		my $seq =$masked{$primer_id}{seq};
		my $out="";
		for(my $length=19;$length<=25;$length++){
				for(my$start=$region_start;$start<=$region_end-$length+1;$start++){
						if($primer_comp{$primer_id}{side} eq "left"){
								if(!defined$specific_mutation{$start}){next;}
						}else{
								if(!defined$specific_mutation{$start+$length-1}){next;}
						}
						my($specific,$specific_ref)=("","");
						my $out_seq = substr($seq,$start-1,$length);
						if($out_seq =~ /N/){next;}
						my ($start_posi,$end_posi,$strand)=($primer_start+$start-1,$primer_start+$start+$length-2,"+");
						if($primer_comp{$primer_id}{side} eq "right"){
								($start_posi,$end_posi)=($end_posi,$start_posi);
								$out_seq =~ tr/atgcATGC/tacgTACG/;
								$out_seq = reverse($out_seq);
								$strand = "-";
						}
						my ($out_read_mut,$out_read_mut_orig) = ("","");
						for(my $i=0;$i<scalar(@read_mutation);$i++){
								if($read_mutation[$i] =~ /^(\d+)-(\d+):(.+)$/){
										if( (($1>=$start)&&($1<$start+$length)) || (($2>=$start)&&($2<$start+$length)) ){
												$out_read_mut.="$1:$3;";
												$out_read_mut_orig.="$read_mutation_orig[$i];";
										}
								}
						}
						$out .= "$primer_comp{$primer_id}{peak_id}\t$primer_id\t$chr\t$start_posi\t$end_posi\t$length\t$strand\t";
						$out .= "$out_seq\t$specific\t$specific_ref\t$out_read_mut\t$out_read_mut_orig\n";
				}
		}
		return($out);
}

#------------------------------------------------------------------------------
sub pick_specific_mutation ( $ ){
		my($primer_id,$file) = @_;
		my ($sample,$chr,$primer_start,$primer_end,$region_start,$region_end);
		my %align_fasta = &fasta2hash($file);
		my $homolog_num = scalar(keys %align_fasta)-1;
		my ($nucl_posi,$posi)=(1,1);
		my %specific_mutation =();
		while($posi<=length($align_fasta{$primer_id}{seq})){
				if(substr($align_fasta{$primer_id}{seq},$posi-1,1) eq "-"){
						$posi++;
				}else{
						my($dif_n,$indel_n) = &check_dif($primer_id,$posi-1,\%align_fasta);
						if($indel_n ==$homolog_num){
								my $indel_length =0;
								while($indel_n ==$homolog_num){
										$indel_length++;
										($dif_n,$indel_n)=&check_dif($primer_id,$posi-1+$indel_length,\%align_fasta);
								}
								print INDEL "$primer_id\t$posi\t$nucl_posi\t$indel_length\n";
								$specific_mutation{$nucl_posi}{seq}=substr($align_fasta{$primer_id}{seq},$posi,$indel_length);
								$posi += $indel_length;
								$nucl_posi += $indel_length;
#####################
						}elsif($dif_n == $homolog_num){
								$specific_mutation{$nucl_posi}{seq}=substr($align_fasta{$primer_id}{seq},$posi,1);
								$posi++;
								$nucl_posi++;
						}else{
								$posi++;
								$nucl_posi++;
						}
				}
		}
		return(%specific_mutation);
}




#------------------------------------------------------------------------------
sub check_dif( $ $ $){
		my ($primer_id,$posi) = ($_[0],$_[1]);
		my %fasta = %{$_[2]};
		my ($primer_nuc,$dif,$indel) = (substr($fasta{$primer_id}{seq},$posi,1),0,0);
		foreach my $id(keys %fasta){
				if($id eq $primer_id){next;}
				my $nuc = substr($fasta{$id}{seq},$posi,1);
				if($primer_nuc ne $nuc){$dif++;}
				if($nuc eq "-"){$indel++;}
		}
		return($dif,$indel);
}

#------------------------------------------------------------------------------
sub get_primer_tbl ( $ $ ){
		my ($primer_id,$region)=@_;
		my ($sample,$chr,$primer_start,$primer_end,$region_start,$region_end);
		if($primer_id =~ /^(\w+-\d\d)_(chr[^:]+):(\d+)-(\d+)$/){
				($sample,$chr,$primer_start,$primer_end) = ($1, $2, $3, $4);
		}
		if($region =~ /^(\d+)-(\d+)$/){
				($region_start,$region_end) = ($1, $2);
		}
		my $read_mutation = "";
		my @read_mutation=();
		my @read_mutation_orig=();
		if(defined$read_mutation{$primer_id}{mut}){$read_mutation = $read_mutation{$primer_id}{mut};
				@read_mutation = split(/;/,$read_mutation{$primer_id}{mut});
				@read_mutation_orig = split(/;/,$read_mutation{$primer_id}{origin});
		}
		my %ref_posi=&position_translation_ref($primer_start,$masked{$primer_id}{seq},$read_mutation);
		my $seq =$masked{$primer_id}{seq};
		my $out="";
		for(my $length=19;$length<=25;$length++){
				for(my$start=$region_start;$start<=$region_end-$length+1;$start++){
						my $out_seq = substr($seq,$start-1,$length);
						if($out_seq =~ /N/){next;}
						my ($start_posi,$end_posi,$strand)=($primer_start+$start-1,$primer_start+$start+$length-2,"+");
						if($primer_comp{$primer_id}{side} eq "right"){
								($start_posi,$end_posi)=($end_posi,$start_posi);
								$out_seq =~ tr/atgcATGC/tacgTACG/;
								$out_seq = reverse($out_seq);
								$strand = "-";
						}
						my ($out_read_mut,$out_read_mut_orig) = ("","");
						for(my $i=0;$i<scalar(@read_mutation);$i++){
								if($read_mutation[$i] =~ /^(\d+)-(\d+):(.+)$/){
										if( (($1>=$start)&&($1<$start+$length)) || (($2>=$start)&&($2<$start+$length)) ){
												$out_read_mut.="$1:$3;";
												$out_read_mut_orig.="$read_mutation_orig[$i];";
										}
								}
						}
						$out .= "$primer_comp{$primer_id}{peak_id}\t$primer_id\t$chr\t$start_posi\t$end_posi\t$length\t$strand\t";
						$out .= "$out_seq\tno_homolog\tno_homolog\t$out_read_mut\t$out_read_mut_orig\n";

				}
		}
		return($out);
}

#------------------------------------------------------------------------------
sub position_translation_ref( $ $ $ ){
		my($ref_start,$seq,$read_mut)=@_;
		my @read_mutation=split(/;/,$read_mut);
		my %indel = ();
		foreach my $read_mutation (@read_mutation){
				if($read_mutation =~ /(\d+)-(\d+):(\w+)>(\w+)$/){
						if($1!=$2){$indel{$1}=length($3)-length($4);}
				}else{die "ERROR::waht indel $read_mutation in start:$ref_start peak\n";}
		}
		my $ref_posi=$ref_start;
		my %translation=();
		for(my $posi=1;$posi<length($seq);$posi++){
				if(defined $indel{$posi}){
						if($indel{$posi}<0){
								$translation{$posi}=$ref_posi;
								$ref_posi+=-$indel{$posi}+1;
						}else{#$indel{$posi}>0
								$translation{$posi}=$ref_posi;
								while($indel{$posi}>0){
										$translation{$posi}="$ref_posi-";
										$posi++;$indel{$posi}--;
								}
								$ref_posi++;
						}
				}else{
						$translation{$posi}=$ref_posi;
						$ref_posi++;
				}
		}
		return(%translation);
}

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
