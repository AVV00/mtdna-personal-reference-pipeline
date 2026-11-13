#! /usr/bin/perl -w 
use strict;
use Excel::Writer::XLSX;
use File::Basename;
use Getopt::Std;
use vars qw($opt_a $opt_b $opt_c);
getopts('a:b:c:');

unless ($opt_a or $opt_b) {
	print
qq(Usage::perl $0  transform txt file into excel 2007
	-a <string> txt file 
	-b <string> file name
	-c <string> sheet prefix, defaul .txt .xls
);
	exit;
}
my @prefix ;

if ($opt_c) {
	push @prefix,$opt_c;
}else{
	push @prefix , ".txt",".xls";
}	
my @file; 
if (-d $opt_a) { 
	my $file = `ls $opt_a/*$opt_c` ; 
	chomp $file; 
	@file = split /\n/,$file; 
}else { 
	@file = split /[\,|\s]+/,$opt_a;
}
open my $fh,">$opt_b.xlsx" or die $!;
my $workbook = Excel::Writer::XLSX->new("$opt_b.xlsx" );

foreach my $sample (@file){
	my $name = basename ($sample,@prefix); 
	if ((length $name) > 30) {
		$name = substr ($name,0,30) ; 
	}
	my $worksheet= $workbook->add_worksheet("$name");
	my $format = $workbook->add_format();
	$format->set_align( 'left' );
	my $row = 0;
	open my $in,$sample ||die $!;
	while (<$in>){
		chomp;
		my @aa = split /\t/,$_;
		foreach my $i (0..@aa-1){
			$worksheet->write($row,$i,$aa[$i],$format);
			binmode STDOUT;
#	print $fh 
		}
		$row++;
	}
	close $in;
}	
$workbook->close();
