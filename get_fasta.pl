#!/usr/bin/perl -w 
use strict;
use File::Basename; 

my @ssp = glob "./*/*ssp" ; 
foreach my $ssp (@ssp) { 
	my $name = basename ($ssp,".ssp") ; 
#	print join "\t",$ssp,"\n" ; 
	my $path = "fasta"; 
	mkdir $path unless -e $path; 
	my $seq; 
	open IN,$ssp ||die $!; 
	while (<IN>) { 
		chomp;
		my @aa = split /\t/,$_; 
		my $major = $aa[1] ; 
		if (defined $aa[3] and $aa[3]) { 
			$major = $aa[3] ; 
		}
		$seq .= $major; 
	} 
	close IN; 

	open OT,">$name/$name.fa" ||die $!;
	print OT ">$name\n" ;
	for (my $i = 0 ; $i <= length $seq; $i += 70) { 
		my $temp= substr($seq,$i,70) ; 
		print OT "$temp\n" ; 
	} 
	close OT; 
} 
