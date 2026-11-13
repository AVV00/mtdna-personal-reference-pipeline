#!/usr/local/bin/perl
use strict;
use Getopt::Std;
use Statistics::Descriptive;


  

use vars qw($ref $opt_i $opt_m $opt_q $opt_d $q_t $opt_r $opt_s $NO_read_bin $opt_a);

&init();

if ($opt_i=~/\.gz$/){
	open IN,"gzip -dc $opt_i|"||die $!;
}else {
	open IN,$opt_i ||die $!;
}	

my $content;
my $quality;
my $position;
my $mapping;
my $difference;
my ($pos, $ref, $cov, $chr);
my $indel;

$NO_read_bin=int($opt_r/$opt_s-0.0000001)+1;

my $range_start;
my $range_end;
my $range_chr;

if($opt_a){
	if($opt_a=~/(\w+):(\d+)-(\d+)/){
		$range_chr=$1;
		$range_start=$2;
		$range_end=$3;
	}
	else{
		die "the target region should be in the format of chrz:x-y\n";
	}
}

while(chomp(my $l=<IN>)){
	my @line=split /\s/,$l;
	$chr=$line[0];
	$pos=$line[1];
	$ref=$line[2];
	next if ($chr ne $range_chr and $opt_a);
	next if ($pos<$range_start and $opt_a);
	next if ($pos>$range_end and $opt_a);
	
	$cov=0;#$line[3];
	$content=&_content($line[4]);  #content is not in the pileup format! actually .<>, . reversed ,forward
	$quality=&_quality($line[5]);
	$mapping=&_quality($line[6]); # mapping QS is in the same format with QS
	$position=&_position($line[7]);
	$difference=&_mismatch($line[8]);
#	print "@$mapping\n";
#	if($opt_m){
#		($content,$quality,$mapping,$position,$difference)=&_filter_great($content,$quality,$mapping,$position,$difference,$mapping,$opt_m);
#	}
#	if($opt_q){
		#my @quality_tmp=@$quality;
#		($content,$quality,$mapping,$position,$difference)=&_filter_great($content,$quality,$mapping,$position,$difference,$quality,$opt_q);
#	}
#	if($opt_d){
#		($content,$quality,$mapping,$position,$difference)=&_filter_small($content,$quality,$mapping,$position,$difference,$difference,$opt_d);
#	}	
	($content,$quality,$mapping,$position,$difference)=&_filter_small($content,$quality,$mapping,$position,$difference);
	
	my $basic=&_basic($content, $quality, $position);
	print $pos,"\t",$ref,"\t",$cov,"\t",$basic->{'first_allele'},"\t",$basic->{'second_allele'},"\t",$basic->{'other_allele'},"\t",$basic->{'first_count'},"\t",$basic->{'second_count'},"\t",$basic->{'other_count'},"\t",$basic->{'first_double'},"\t",$basic->{'second_double'},"\t",$basic->{'other_double'},"\t",$basic->{'first_deviation'},"\t",$basic->{'second_deviation'},"\t",$basic->{'other_deviation'},"\t",$basic->{'first_10'},"\t",$basic->{'second_10'},"\t",$basic->{'other_10'},"\t",$basic->{'flag'},"\n";
}

sub _filter_small(){
	my $content=shift;
	my $quality=shift;
	my $mapping=shift;
	my $position=shift;
	my $difference=shift;
	
#	my $flag=shift;
#	my $value=shift;
	
	my @content_new;
	my @quality_new;
	my @mapping_new;
	my @position_new;
	my @difference_new;
=head 	
	for(my $i=0;$i<@$flag;$i++){	
		if ($flag->[$i] <= $value){
			push @content_new,$content->[$i];
			push @quality_new,$quality->[$i];
			push @mapping_new,$mapping->[$i];
			push @position_new,$position->[$i];
			push @difference_new,$difference->[$i];
		}
		
	}
=cut 
	my ($m_check,$q_check,$d_check)= (0)x3;
	$m_check = 1 if $opt_m;
	$q_check = 2 if $opt_q;
	$d_check = 4 if $opt_d;
	my $sum = $m_check + $q_check + $d_check;
	my %check;
	if ($sum==7){			
		for (my $i=0;$i<@$content;$i++){
			next if ($mapping->[$i]<$opt_m);
			next if ($quality->[$i]<$opt_q);
			next if ($difference->[$i]>$opt_d);
			push @content_new,$content->[$i];
			push @quality_new,$quality->[$i];
			push @mapping_new,$mapping->[$i];
			push @position_new,$position->[$i];
			push @difference_new,$difference->[$i];
		}
	}elsif ($sum==6){	
		for (my $i=0;$i<@$content;$i++){
#			next if ($mapping->[$i]<$opt_m);
			next if ($quality->[$i]<$opt_q);
			next if ($difference->[$i]>$opt_d);
			push @content_new,$content->[$i];
			push @quality_new,$quality->[$i];
			push @mapping_new,$mapping->[$i];
			push @position_new,$position->[$i];
			push @difference_new,$difference->[$i];
		}
	}elsif($sum==5){
		for (my $i=0;$i<@$content;$i++){
			next if ($mapping->[$i]<$opt_m);
#			next if ($quality->[$i]<$opt_q);
			next if ($difference->[$i]>$opt_d);
			push @content_new,$content->[$i];
			push @quality_new,$quality->[$i];
			push @mapping_new,$mapping->[$i];
			push @position_new,$position->[$i];
			push @difference_new,$difference->[$i];
		}
	}elsif ($sum==4){
		for (my $i=0;$i<@$content;$i++){
#			next if ($mapping->[$i]<$opt_m);
#			next if ($quality->[$i]<$opt_q);
			next if ($difference->[$i]>$opt_d);
			push @content_new,$content->[$i];
			push @quality_new,$quality->[$i];
			push @mapping_new,$mapping->[$i];
			push @position_new,$position->[$i];
			push @difference_new,$difference->[$i];
		}
	}elsif ($sum==3){
		for (my $i=0;$i<@$content;$i++){
			next if ($mapping->[$i]<$opt_m);
			next if ($quality->[$i]<$opt_q);
#			next if ($difference->[$i]>$opt_d);
			push @content_new,$content->[$i];
			push @quality_new,$quality->[$i];
			push @mapping_new,$mapping->[$i];
			push @position_new,$position->[$i];
			push @difference_new,$difference->[$i];
		}
	}elsif($sum==2){
		for (my $i=0;$i<@$content;$i++){
#			next if ($mapping->[$i]<$opt_m);
			next if ($quality->[$i]<$opt_q);
#			next if ($difference->[$i]>$opt_d);
			push @content_new,$content->[$i];
			push @quality_new,$quality->[$i];
			push @mapping_new,$mapping->[$i];
			push @position_new,$position->[$i];
			push @difference_new,$difference->[$i];
		}
	}elsif($sum==1){
		for (my $i=0;$i<@$content;$i++){
			next if ($mapping->[$i]<$opt_m);
#			next if ($quality->[$i]<$opt_q);
#			next if ($difference->[$i]>$opt_d);
			push @content_new,$content->[$i];
			push @quality_new,$quality->[$i];
			push @mapping_new,$mapping->[$i];
			push @position_new,$position->[$i];
			push @difference_new,$difference->[$i];
		}
	}elsif ($sum==0){
		@content_new=@$content;
		@quality_new=@$quality;
		@mapping_new=@$mapping;
		@position_new=@$position;
		@difference_new=@$difference;
	}	
#		for (my $i=0;$i<@$content;$i++){
#			next if ($mapping->[$i]<$opt_m);
#			next if ($quality->[$i]<$opt_q);
#			next if ($difference->[$i]>$opt_d);
#			push @content_new,$content->[$i];
#			push @quality_new,$quality->[$i];
#			push @mapping_new,$mapping->[$i];
#			push @position_new,$position->[$i];
#			push @difference_new,$difference->[$i];
#		}
	return (\@content_new,\@quality_new,\@mapping_new,\@position_new,\@difference_new);
}

sub _filter_great(){
	my $content=shift;
	my $quality=shift;
	my $mapping=shift;
	my $position=shift;
	my $difference=shift;
	
	my $flag=shift;
	my $value=shift;
	
	my @content_new;
	my @quality_new;
	my @mapping_new;
	my @position_new;
	my @difference_new;
	
	for(my $i=0;$i<@$flag;$i++){	
		if ($flag->[$i] >= $value){
			push @content_new,$content->[$i];
			push @quality_new,$quality->[$i];
			push @mapping_new,$mapping->[$i];
			push @position_new,$position->[$i];
			push @difference_new,$difference->[$i];
		}
		
	}
	return (\@content_new,\@quality_new,\@mapping_new,\@position_new,\@difference_new);
}


sub _basic(){
	my $content=shift;
	my $quality=shift;
	my $position=shift;
	my %hash;
	my %hash_quality;
	my %hash_strand;
	my %hash_position;
	my $flag;
	my %result;
	
	for (my $i=0;$i<@$content;$i++){
		$hash{$content->[$i]}++;
		push @{$hash_quality{$content->[$i]}},$quality->[$i];
		push @{$hash_position{$content->[$i]}},$position->[$i];  #postion added
		my $tmp=uc $content->[$i];
		$hash_strand{$tmp}++;
	}
	
	my $count;
	my @for_10=(0)x$NO_read_bin;
	my @rev_10=(0)x$NO_read_bin;
	foreach my $base ( sort {$hash_strand{$b}<=>$hash_strand{$a}}  keys %hash_strand){
		$count++;
		my ($for,$rev,$qua_f,$qua_r,$for_u,$rev_u,$for_d,$rev_d,$for_m,$rev_m);
		$for=0;$rev=0;$for_u=0;$rev_u=0;
		my @first_for_10=(0)x$NO_read_bin;
		my @first_rev_10=(0)x$NO_read_bin;
		
		my $tmp=uc $base;
		$for=$hash{$tmp}||0;
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@{$hash_quality{$tmp}});
		$qua_f=$stat->median() if (@{$hash_quality{$tmp}}>0);
		$qua_f=$qua_f||0;
		
		$for_u=&_unique_read([@{$hash_position{$tmp}}]) if (defined(@{$hash_position{$tmp}}));
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@{$hash_position{$tmp}});
		$for_d=$stat->standard_deviation();
		$for_d=sprintf("%.".(0)."f", $for_d);
		$for_m=$stat->mean();
		$for_m=sprintf("%.".(0)."f", $for_m);
		
		if($count>1){
			@for_10=@{&_10_divi([@{$hash_position{$tmp}}],[@for_10])};
		}
		else{
			@first_for_10=@{&_10_divi([@{$hash_position{$tmp}}],[@first_for_10])};
		}
				
		$tmp=lc $base;
		$rev=$hash{$tmp}||0;
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@{$hash_quality{$tmp}});
		$qua_r=$stat->median() if (@{$hash_quality{$tmp}}>0);
		$qua_r=$qua_r||0;
		
		$rev_u=&_unique_read([@{$hash_position{$tmp}}]) if (defined (@{$hash_position{$tmp}}));
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@{$hash_position{$tmp}});
		$rev_d=$stat->standard_deviation();
		$rev_d=sprintf("%.".(0)."f", $rev_d);
		$rev_m=$stat->mean();
		$rev_m=sprintf("%.".(0)."f", $rev_m);
		
		if($count>1){
			@rev_10=@{&_10_divi([@{$hash_position{$tmp}}],[@rev_10])};
		}
		else{
			@first_rev_10=@{&_10_divi([@{$hash_position{$tmp}}],[@first_rev_10])};
		}
		
		$flag.=$base.':'.$for.':'.$rev.':'.$qua_f.':'.$qua_r.':'.$for_u.':'.$rev_u."|";
		$hash{'flag'}=$flag;
		
		$cov+=$rev+$for;
		
		if ($count==1){
			$hash{'first_allele'}=$base;
			$hash{'first_count'}=$for+$rev;
			if($for && $rev){
				$hash{'first_double'}="FC:2".":($for-$rev)"."($for_u-$rev_u)";
			}
			else{
				$hash{'first_double'}="FC:1".":($for-$rev)"."($for_u-$rev_u)";
			}
			$hash{'first_deviation'}="FP:"."($for_m±$for_d:$rev_m±$rev_d)";
			$hash{'first_10'}="FD:".&_string_10([@first_for_10],[@first_rev_10]);
		}
		else{
			if($count==2){
				$hash{'second_allele'}=$base;
				$hash{'second_count'}=$for+$rev;
				if($for && $rev){
					$hash{'second_double'}="SC:2".":($for-$rev)"."($for_u-$rev_u)";
				}
				else{
					$hash{'second_double'}="SC:1".":($for-$rev)"."($for_u-$rev_u)";
				}
				$hash{'second_deviation'}="SP:"."($for_m±$for_d:$rev_m±$rev_d)";
				my $string_10=&_string_10([@for_10],[@rev_10]);
				$hash{'second_10'}="SD:".$string_10;
			}
			if($count >2) {#revised  2020.08.03  newly add 
			  $hash{'other_allele'}.=$base;
			  $hash{'other_count'}+=$for+$rev;
		  	if($for && $rev){
			  	$hash{'other_double'}.="2".":($for-$rev)"."($for_u-$rev_u)";
			  }
			  else{
				  $hash{'other_double'}.="1".":($for-$rev)"."($for_u-$rev_u)";
			  }
			  $hash{'other_deviation'}.="($for_m±$for_d:$rev_m±$rev_d)";
			}   #revised  2020.08.03 newly add 
		}
		
	}
	$hash{'other_10'}="OD:".&_string_10([@for_10],[@rev_10]) if ($hash{'other_allele'});
	$hash{'other_double'}="OC:".$hash{'other_double'} if($hash{'other_double'});
	$hash{'other_deviation'}="OD:".$hash{'other_deviation'} if($hash{'other_deviation'});
	
	return \%hash;
}

sub _string_10(){
	my $for=shift;
	my $rev=shift;
	my $string_10;
	for(my $i=0;$i<$NO_read_bin;$i++){
		my $prefix=1+$opt_s*$i;
		$string_10.=$prefix.":"."$for->[$i]"."-"."$rev->[$i]"."|";
	}
	return $string_10;
}

sub _10_divi(){
	my $pos=shift;
	my $array=shift;
#	foreach my $int (@$pos){
#		my $tmp=int(($int-1)/$opt_s);
#		$array->[$tmp]++;
#	}
	map {my $tmp=int(($_-1)/$opt_s);$array->[$tmp]++;} @$pos;
	return $array;
}

sub _unique_read(){
	my $position=shift;
	my %hash;
#	foreach my $pos (@$position){
#		$hash{$pos}++;
#	}
	map {$hash{$_}++} @$position;
	my $result=keys %hash;
	return $result;
}





############debugged;
sub _mismatch(){
	my $mismatch=shift;
	my @result;
#	my @line=split //,$mismatch;
#	foreach my $base (@line){
#		push @result,&_QS_convert_2($base);
#	}
	map {push @result,&_QS_convert_2($_)} (split //,$mismatch);
	return (\@result);
}

sub _position(){
	my $position=shift;
	my @result;
#	my @line=split /\,/,$position;
#	foreach my $base (@line){
#		push @result,$base;
#	}
	map {push @result,$_} (split /\,/,$position);
	return (\@result);
}

sub _quality(){
	my $quality=shift;
	
	my @result;
	my @line=split //,$quality;
	shift @line;
#	foreach my $base (@line){
#		push @result,&_QS_convert($base);
#	}
	map {push @result,&_QS_convert($_)} @line;
	return (\@result);
}

sub _content(){
	 my $content=shift;
	 
	 my @result;
	 my @line=split //,$content;
	 shift @line;
	 foreach my $base (@line){
		if ($base=~/[A-Z]/i){
			push @result, $base;
		}
		elsif($base=~/\,/){
			push @result, uc $ref;
		}
		elsif($base=~/\./){
			push @result, lc $ref;
		}
		else{
			print "warning: $base\n";
		}
	 }
	 return (\@result);
}

sub _QS_convert(){
	my $int=shift;
	my $score;
	$score=ord($int)-33 if($q_t eq 'phred');
	$score= 10 * log(1 + 10 ** ((ord($int) - 64) / 10.0)) / log(10) if($q_t eq 'solexa');
	return int($score);
}

sub _QS_convert_2(){
	my $int=shift;
	my $score;
	$score=ord($int)-48;
	return int($score);
}

sub init(){
	getopts( 'm:i:q:d:r:s:a:' );
	unless(defined ($opt_i)){
		die "perl filter_and_summary.pl\t-i [pileup file]\n\t\t\t\t-d [number of difference allowed]\n\t\t\t\t-q [minimum quality score]\n\t\t\t\t-m [minimum mapping quality score]\n\t\t\t\t-r [length of the read]\n\t\t\t\t-s [read bin length]\n\t\t\t\t-a [your target region chr:a-b on your reference]\n";
	}
	$opt_r=$opt_r||76;
	$opt_s=$opt_s||10;
	$q_t="phred";
}
