#! /usr/bin/perl -w 
use strict;
use File::Basename;
use Cwd 'abs_path' ;
use Data::Dumper;
use threads; 
use Getopt::Std;
use vars qw($opt_f $opt_t); 
getopts('f:t:'); 

$opt_f ||= "data/filter/analysis" ; 
$opt_t ||=2;
my @ssp =  glob "$opt_f/*/*ssp" ; 
unless (@ssp) { 
  print 
  qq(usage::perl $0 automatic run the shell
     \t-f 	<string>	the analysis file,defaut is data/filter/analysis,the ssp files are in <file>/*/*ssp;
     \t-t 	<string>	the heteroplasmy level,1 mean 1%,2 mean 2%,and default theashold is 2%  
     Please check the parameter for -f 
     
  );	
  exit();
} 
#push @ssp, glob "14_0311/*/analysis/*/*ssp" ; 
my $thresh = $opt_t; 
#push @ssp, glob "14_0311/*/analysis/*/*ssp" ; 
#my %hash; 
my @group; 
foreach my $ssp (@ssp) {
  $ssp = abs_path $ssp;
  #	next unless $ssp=~/A927P1/; 
  my $num = `wc -l $ssp |awk '{print \$1}'` ; chomp $num;
  #	next unless $num == 16569;  
  #	next unless $ssp =~/A225_B/ ; 
  #	next unless $ssp=~/unfilter/ and $ssp=~/CC667/ ;  
  push @group,$ssp;
}

my @array ; 
foreach my $i (0 .. @group-1) { 
  my $j =int ($i / 10) ;
  push @{$array[$j]}, $group[$i] ; 
} 

foreach my $array (@array) { 
  #	my @info = @{$array}; 
  my @threads;
  my $i =0;
  #	print Dumper ($array) ; 
  foreach my $shell(@$array){
    #		print join "\t",$shell,"\n" ; 
    $threads[$i]=threads->create(\&_analysis,$shell);
    $i++;
  }   
  foreach my $thread (@threads){
    $thread->join();
  }
} 
sub _analysis {
  my $ssp = shift;
  my $indel = &_major($ssp) ; 
  my $name = basename ($ssp,".ssp"); 
  my $path = dirname $ssp;
  #	my $depth =  0 ; 
  #`cat $path/*.depth.sample_summary|sed -n '2p' |cut -f 3` ; chomp $depth; 
  open IN,$ssp ||die $!;
  open OT, ">$path/$name.hetro_$thresh.txt" ||die $!;
  open OT1,">$path/$name.homo_$thresh.txt" ||die $!; 
#  print "$name\n";
   my $i=0;
  while (<IN>){
    chomp;
    next if /^#/;
      next if /^\s+$/;
    #my @aa =split;            #raw  
    my @aa =split /\t/;       #revised  2020.08.19  split不给分隔符默认会将每一行中的空格和\t当成分隔符  影响aa[16],出现Use of uninitialized value $data in substitution (s///) . Use of uninitialized value $data in split.
    next if $aa[2]==0;
    next if &_region($aa[0]);
    ##		C:455:233:65:66:19:6|G:1:0:66:0:1:0|
    my ($pos,$ref,$total) = @aa[0,1,2];
    my ($major_fw,$major_rv,$minor_fw,$minor_rv,$other_fw,$other_rv,$other_count)=(0)x7;
    my ($major_base,$minor_base,$other_base) ; 
    $aa[-1]=~s/\|$//;
    my @data = split /\|/,$aa[-1];
    my %hash;
    
    #new 2020.08.10
    #创建了%base_count,用ssp中最后一列的信息，来计算每个位点中所有类型的碱基在正链和负链的个数，大写碱基为正链，小写碱基为负链
    my %base_count;
    $base_count{"A"}=$base_count{"T"}=$base_count{"C"}=$base_count{"G"}=0;
    $base_count{"a"}=$base_count{"t"}=$base_count{"c"}=$base_count{"g"}=0;
    foreach my $i (0 .. @data-1){
      my ($base,$fw,$rv)=(split /:/,$data[$i])[0,1,2];
      my $base_fw=uc $base;
      my $base_rv=lc $base;
      $base_count{$base_fw}=$fw;
      $base_count{$base_rv}=$rv;
    }
    
    ($major_base,$major_fw,$major_rv)=(split /:/,$data[0])[0,1,2];
    ($minor_base,$minor_fw,$minor_rv)=(split /:/,$data[1])[0,1,2] if defined $data[1];
    if (scalar @data >2){
      foreach my $i(2..$#data){
                      my @temp= split /:/,$data[$i];
                    $other_base .= "$temp[0]\t";
                    $other_fw .= "$temp[1]\t";
                    $other_rv .= "$temp[2]\t";
    }
  }
  my $i = 0 ;   
  $other_count = &_sum($other_fw) + &_sum ($other_rv) ; 
  my $major_count = $major_fw + $major_rv ; 
  my $minor_count = $minor_fw + $minor_rv ; 
  my $flank_base =1  ; 
  foreach my $j (-2,-1,1,2) {
    if (defined $indel->{$pos+$j} and $indel->{$pos+$j} =~/D|I/) {
      $flank_base = 0 ; 
    }
  }  
=head 		
  if ($pos == 4135) {
    print join "\t",$total,$minor_fw,$major_fw,$minor_fw,$major_fw,$minor_fw/($minor_fw + $major_fw),$minor_rv/($minor_rv+$major_rv),$other_count/$total,$major_base,$minor_base, "\n" ; 
  }
=cut 
  if ($total >=100 and ($minor_fw + $major_fw) and $minor_fw >=3 and $minor_fw/($minor_fw + $major_fw) > $thresh/100 and ($minor_rv+$major_rv) and $minor_rv >=3 and $minor_rv/($minor_rv+$major_rv) > $thresh/100  and $other_count/$total < $thresh/100 and $minor_base!~/I|D/ and $major_base!~/I|D/ and $flank_base ) { 
    #			print "$_\n"; 
    my $j = 0 ; 
    #		if ($pos >=500 and $pos <= 16569-500) { 
    my $end = &_end ($aa[16]) ;  
    my $minor_ratio = $minor_count/$total; 
    my $bias  = &_bias([$minor_fw,$minor_rv,$major_fw,$major_rv]) ; 
    if ( $bias>1) {
      $j ++; 
    }
    $j = 0 if $minor_ratio >=0.2; 
    if($end) { 
      $j++;
    } 
    
    #new 2020.08.10
    #创建了三个变量，输出主次等位碱基、参考基因组碱基的前后一个碱基的信息。
    if((not $minor_base)){$minor_base=' '}; 
    my ($three_base_major,$three_base_minor,$three_base_ref)=(0)x3;
    if($pos==1){
      $three_base_major=" ".$major_base.&ref_base($aa[0]+1);
      $three_base_minor=" ".$minor_base.&ref_base($aa[0]+1);
      $three_base_ref=" ".$ref.&ref_base($aa[0]+1);
    }else{
      $three_base_major=&ref_base($aa[0]-1).$major_base.&ref_base($aa[0]+1);
      $three_base_minor=&ref_base($aa[0]-1).$minor_base.&ref_base($aa[0]+1);
      $three_base_ref=&ref_base($aa[0]-1).$ref.&ref_base($aa[0]+1);
    }  
    
    if ($j ==0) {
      $i ++ ; 
      if ($major_base eq $ref) {
        #print OT join "\t",$aa[0],$major_base,$minor_base,$major_fw,$major_rv,$minor_fw,$minor_rv,(sprintf "%.4f",$minor_count/$total),"\n";                                                                 #raw
        print OT join "\t",$aa[0],$ref,$minor_base,$three_base_ref,$three_base_minor,$base_count{uc $ref},$base_count{lc $ref},$minor_fw,$minor_rv,(sprintf "%.4f",$minor_count/$total),"\n";          #revised 2020.08.10 将结果中第一列碱基信息，变成参考基因组碱基的信息 新增两列碱基前后碱基信息
      }else {
        #print OT join "\t",$aa[0],$minor_base,$major_base,$minor_fw,$minor_rv,$major_fw,$major_rv,(sprintf "%.4f",$major_count/$total),"\n";                                                                 #raw
        print OT join "\t",$aa[0],$ref,$major_base,$three_base_ref,$three_base_major,$base_count{uc $ref},$base_count{lc $ref},$major_fw,$major_rv,(sprintf "%.4f",$major_count/$total),"\n";            #revised 2020.08.10 同上
      } 
    }
  }		
  if($i==0 ) {
    if ($major_base ne $ref) { 
      if($total >= 50 or ((not $minor_base) and $total>8)) {
        next if $major_base=~/I|D/ ;
        
        #new 2020.08.10
        #创建了两个变量，输出主等位碱基、参考基因组碱基的前后一个碱基的信息。
        my ($three_base_major,$three_base_ref)=(0)x3;
        if($pos==1){
          $three_base_major=" ".$major_base.&ref_base($aa[0]+1);
          $three_base_ref=" ".$ref.&ref_base($aa[0]+1);
        }else{
          $three_base_major=&ref_base($aa[0]-1).$major_base.&ref_base($aa[0]+1);
          $three_base_ref=&ref_base($aa[0]-1).$ref.&ref_base($aa[0]+1);
        }  
        
        #my $ratio = sprintf "%.2f",$major_count/($major_count + $minor_count) ;                                                                                                #raw
        my $ratio = sprintf "%.2f",$major_count/($total) ;                                                                                                                    #revised 2020.08.10  分母统一成total
        #print OT1 join "\t",$aa[0],$ref,$major_base,$minor_fw,$minor_rv,$major_fw,$major_rv,$ratio,"\n" ;                                                                      #raw
        print OT1 join "\t",$aa[0],$ref,$major_base,$three_base_ref,$three_base_major,$base_count{uc $ref},$base_count{lc $ref},$major_fw,$major_rv,$ratio,"\n" ;           #revised 2020.08.10  同hetro
      }
    }
  } 		
}
close IN;		
close OT;
close OT1;
}	
sub _major {
  my $file = shift;
  my %major; 
  
  open FILE,$file ||die$!; 
  while (<FILE>) { 
    chomp; 
    my @aa = split /\t/,$_; 
    if ($aa[3]) {
      $major{$aa[0]} = $aa[3] ; 
    }
  }  
  close FILE; 
  return (\%major); 
}

sub _sum {
  my $file = shift;
  $file=~s/\t$// if $file=~/\t$/ ; 
  my @file = split /\s+/,$file;
  my $num =0;
  map {$num+=$_} @file;
  return $num;
}
sub _region {
  my $pos=shift;
  my @region;
  push @region,302..316,513..526,566..573,16181..16194;
  my @result = grep {$pos == $_} @region;
  if (@result){
    return 1;
  }else{
    return 0;
  }
}
sub _bias {
  my $data = shift;
  my @data = @$data; 
  #	print join "\t",@data,"\n" ; 
  if (($data[0] +$data[1])==0 or ($data[2]+$data[3])== 0 ) { 
    return 0; 
  }else { 
    my $result= `strand_bias.pl @data`; chomp $result;
    return abs($result);
  } 
} 
sub _end { 
  my $data = shift;
  #	print "$data\n" ;
  ##	SD:1:0-0|6:0-0|11:0-0|16:0-0|21:0-0|26:0-0|31:0-0|36:0-0|41:0-0|46:0-0|51:0-0|56:0-0|61:0-0|66:0-0|71:0-0|76:0-0|81:0-0|86:0-0|91:0-0|96:0-0|101:0-1|106:0-0|111:0-0|116:0-0| 
  $data=~s/^SD://; 
  my @data = split /\|/,$data; 
  my ($start,$end,$total) = (0) x 3 ; 
  foreach my $i (0 .. @data-1) { 
    my @array = split /[:-]/,$data[$i] ; 
    if ($i == 0 ) { 
      $start = $array[1] + $array[2] ; 
    }elsif($i == @data-1){
      $end = $array[1] + $array[2] ; 
    } 
    $total += $array[1] + $array[2] ; 
  } 
  if ($total and ($start + $end)/$total > 0.8) {
    return 1 ; 
  }else {
    return 0 ; 
  } 
} 

#new 2020.08.10
#输入：位点的位置
#输出：返回参考基因组的位置对应的碱基
#过程：使用17119个碱基的参考基因组human_mtDNA.fasta，构建1:17119对应碱基的哈希，key为position，value为碱基
sub ref_base(){
  my $pos=shift;
  open(FILE,"<",'/mnt/hdd/ref_data/modified_genome/human_mtDNA.fasta')||die"cannot open the human_mtDNA.fasta: $!\n";
  readline FILE;
  chomp(my @data=<FILE>);
  close FILE;
  
  my @seq = split('',join('',@data));
  my %base;
  
  my $position=0;
  foreach my $i (@seq){
    $position=$position+1;
    $base{$position}=$i;
  }
  return($base{$pos});
}
