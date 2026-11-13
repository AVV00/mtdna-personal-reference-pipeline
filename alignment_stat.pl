#!/usr/bin/perl -w 
use strict;
use Data::Dumper;
use File::Basename;
use threads; 
use Getopt::Std;
use vars qw($opt_a $opt_b $opt_c $opt_o $opt_e);
getopts('a:b:c:o:e');
use Text::CSV_XS qw( csv );

$opt_a ||= "./" ; 
unless ($opt_a && $opt_b ) { 
	print "Usage :: perl $0 basic alignment statistics for the whole lane
	-a	<string>	the analysis file,default is ./data/fiter/analysis
	-b 	<string>	haplogroup file
	-c	<num>		threads number,default is 10
	-o 	<string>	the outfile name,default is data_summaray\$data.txt
	-e			whether transform the output into excel format\n\n"  ; 
	exit() ; 
}  
my $file = $opt_a ; 
#$file ||= "data/filter/analysis/" ; 
my $thread_num = $opt_c || 10;
my $haplogroup ; 
if (-e $opt_b) { 
	$haplogroup =  &_haplo($opt_b) ; 
}  
my $time = `date +%Y%m%d` ; chomp $time; 
$opt_o ||= "data_summary_$time.txt" ; 
my @ssp = glob "$file/*/*ssp" ; 
#push @ssp,glob "../data/filter/analysis/*plasma*/*ssp" ; 
my $i = 0 ;
my @data; 
foreach my $ssp (@ssp) { 
	my $sample = basename ($ssp,".ssp") ; 
#	print $ssp,"\n" ; 
	my $j =int ($i / $thread_num) ;
#	print join "\t",$i,$j,"\n" ; 
	$i ++; 
	push @{$data[$j]},$ssp; 
} 

open OT,">$opt_o" ||die $! ; 
print OT join "\t","#name\ttotal_reads\tQ20\tQ30\trmdup_reads\tduplication_rate\talign_hg19_reads\thg19_mapping_rate\talign_chrM_reads\tchrM_mapping_rate\tdepth\thetroplasmy_num\thomoplasmy_num\thaplogroup" ;
print OT "\n" ; 
foreach my $data (@data) { 
#	print join "\t",@$data,"\n" ; 	
	my @threads;
	my $j =0;
	foreach my $shell(@$data){
		$threads[$j]=threads->create(\&_analysis,$shell);
		$j++;
	}
	foreach my $thread (@threads){
		$thread->join();
	}

} 
my $basename = basename ($opt_o,".txt") ; 
if ($opt_e ) { 
	system ("perl /mnt/hdd/softs/mt_analysis/transform.pl -a $opt_o -b $basename")  ; 
} 

sub _analysis { 
	my $ssp = shift; 
	my $name = basename ($ssp,".ssp") ; 
	my $dir = dirname $ssp ;
	unless (-e "$dir/$name.depth.sample_summary") { 
		&_depth ("$dir/$name.depth") ; 
	} 
	my $depth = `cat $dir/$name.depth*ry |sed -n '2p' |cut -f 3` ; chomp $depth; 
	my $total_bam = &_count_bam ("$dir/$name.sort.picard.bam") ; 	
	my $rmdup = `samtools view -c $dir/$name.rmdup.bam`; chomp $rmdup; 
	my $hg19 = `samtools view -c -f 0x2 $dir/$name.rmdup.bam` ; chomp $hg19 ; 
	my $mt = `samtools view -c -f 0x2 $dir/$name.rmdup.bam chrM ` ; chomp $mt ; 
	my $hetro = `cat $dir/$name\.hetro_2.txt |wc -l ` ; chomp $hetro;  
	my $homo = `cat $dir/$name\.homo_2.txt |wc -l ` ; chomp $homo;  
	my $hap = (keys %{$haplogroup} && (exists $haplogroup->{$name}) ) ? $haplogroup->{$name}: "" ; 
	print OT join "\t",$name,@$total_bam,$rmdup,(sprintf "%.4f",(1- $rmdup/$total_bam->[0])),$hg19,(sprintf "%.4f",$hg19 / $rmdup) , $mt, (sprintf "%.4f",$mt/$rmdup), $depth,$hetro,$homo,$hap; 
	print OT "\n" ; 
} 

sub _count_bam { 
	my $bam = shift; 

	my ($total,$base,$Q20,$Q30) = (0) x 4 ; 

	open IN,"samtools view $bam|" ||die $! ; 
	while (<IN>) { 
		chomp;
		my @aa = split /\t/,$_; 
		$total ++; 
		$base+= (length $aa[10]) ; 

		my @quality = map{ord ($_) - 33} (split //,$aa[10]) ; 
		my $above_20 = grep {$_ >= 20} @quality ; 
		my $above_30 = grep {$_>=30} @quality ; 
		$Q20 += $above_20; 
		$Q30 += $above_30; 

		#print join "\t",@quality,$aa[10],"\n" ; 
	} 
	close IN; 

	my $Q20_ratio = sprintf "%.4f",$Q20 / $base; 
	my $Q30_ratio = sprintf "%.4f",$Q30 / $base; 
	
	return ([$total,$Q20_ratio,$Q30_ratio]) ; 
} 
sub _depth { 
	my $depth = shift ; 
	my $name = basename ($depth,".depth") ; 
	my $sum =0  ; 
	open IN,$depth ||die $!; 
	while (<IN>) { 
		chomp;
		my @aa = split; 
#	print join "\t",@aa,"\n" ;
		next if $aa[0]=~/Locus/;
		$sum += $aa[-1] ; 
	} 
	close IN; 
	my $average = sprintf "%.2f",$sum / 16569 ; 
	my $dir = dirname $depth ;
	open OUT,">$dir/$name.depth.sample_summary" ||die $!; 
	print OUT join "\t",'sample_id','total','mean','granular_third_quartile','granular_median','granular_first_quartile','%_bases_above_15' ;
	print OUT "\n"; 
	print OUT join "\t",$name,$sum,$average,500,500,500,99.9; 
	print OUT "\n";
	print OUT join "\t","Total",$sum,$average,"N/A","N/A","N/A" ; 
	print OUT "\n";
}

sub _haplo{ 
	my $data = shift ; 
	my %hash ; 
	my $csv = Text::CSV_XS->new ({ binary => 1, auto_diag => 1 }); 
	open my $fh, $data  or die "test.csv: $!";
	while (my $row = $csv->getline ($fh)) {
		next if $row->[0]=~/name/i; 
		$hash{$row->[0]} = $row->[1] ; 
	}
	close $fh ; 
	return (\%hash) ; 
} 	
