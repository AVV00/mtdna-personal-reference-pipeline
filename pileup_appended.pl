#!/usr/lcoal/bin/perl
use strict;

my $file=shift||die "perl pileup_appended.pl [file]\n";


my %hash;
my $num = `wc -l $file` ; chomp $num ; 
my $size ; 
if ($num <= (16569 + 550) and $num > (16299 + 550)) { 
	$size = 16569 ; ## human 
}else { 
	$size = 16299; ## mouse 
} 
open IN,$file||$!;
while(<IN>){
	my @array=split /\s+/,$_;
	next if ($array[1]<=$size);
	chomp;
	my $pos=$array[1]-$size;
	$hash{$pos}=$_;
}
close IN;

open IN,$file||$!;
while(<IN>){
	chomp;
	my @array=split /\s+/,$_;
	if ($array[1]<=550){
		my $quality=&_quality($array[6]);
		#die "can not find $array[1]+16569" unless (defined($hash{$array[1]}));   # it is not necessary that all positions have reads in appended region?
		my @array_appended=split /\s+/,$hash{$array[1]};
		$array[4]=~s/^\@//;
		$array[5]=~s/^\@//;
		$array_appended[4]=~s/^\@//;
		$array_appended[5]=~s/^\@//;			
		$array_appended[6]=&_quality($array_appended[6]);
		my $new_3=$array[3]+$array_appended[3];
		my $new_4="@".$array[4].$array_appended[4];
		my $new_5="@".$array[5].$array_appended[5];
		my $new_6="@".$quality.$array_appended[6];
		if($array[7]){
			$array[7]=$array[7].",";
		}
		my $new_7=$array[7].$array_appended[7];
		$new_4=~s/^@$//;
		$new_5=~s/^@$//;
		$new_6=~s/^@$//;
		$new_7=~s/\,$//;
		my $new_8=$array[8].$array_appended[8];
		my $count_ID;
		while($new_4=~/I/ig){
			$count_ID++;
		}
		if(length($new_4)==$new_3+1+$count_ID and length($new_5)==$new_3+1+$count_ID and length($new_6)==$new_3+1+$count_ID and length($new_8)==$new_3+$count_ID){
			print $array[0]," ",$array[1]," ",$array[2]," ",$new_3," ",$new_4," ",$new_5," ",$new_6," ",$new_7," ",$new_8,"\n";
		}
		else{
			die "length of different elements are different";
		}	
	}
	elsif($array[1]>550 and $array[1]<=$size){
		print $_,"\n";
	}
	else{
		exit;
	}
}
close IN;

sub _quality(){
	my $string=shift;
	my $new;
	my $count;
	my @array=split //,$string;
	foreach my $i (@array){
		$count++;
		next if ($count==1);
		if(&_QS_convert($i)<=20){
			$new.="]";
		}
		else{
			$new.=$i;
		}
	}
	return $new;
}

sub _QS_convert(){
        my $int=shift;
        my $score;
        $score=ord($int)-33;
        return int($score);
}



