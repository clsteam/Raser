#!/usr/bin/perl -w
use strict;

#This program can print fastq format's reads quality score, and help to decide it's encode in phred33 or phred64
#Usage:
#  perl fastq_phred_decide.pl fastq_file_name [reads number,default 1000]
#

if ($#ARGV < 0){
print "faild paramater!\n";
exit;
}

my ($Q, @IN, $i, $read_num);
my ($count,$lt_zero,$base_count) = (0,0,0);

#Default Phred64 to test
$Q=64;
($#ARGV == 1 && $ARGV[1] >= 1) ? ($read_num = $ARGV[1] ) : ($read_num = 1000);

open IN,"<$ARGV[0]" or die "Can not open $ARGV[0]:$!\n";

while($count < $read_num){
$count++;
	@IN=();	
	#read four lines from IN
	for($i=0; $i<=3; $i++){
		$IN[$i]=<IN>;
	}
last if $IN[0] eq "" ;
	&print_score($IN[3],$Q);
	}


sub print_score{
	my ($read, $Q) = @_;
	my ($j,$score);
	for($j=0; $j<=length($read)-2; $j++){
		$score = ord(substr($read, $j,1))-$Q;
	$base_count++;
	if ($score < 0){ $lt_zero += 1; }
	}
}
if ($lt_zero > 0 ){
  #print "Negative value appear $lt_zero times in $base_count base quality score, it should be Phred33.\n\n";
  print "-phred33";
}else{
  #print "Negative value appear $lt_zero times in $base_count base quality score, it should be Phred64.\n\n";
  print "-phred64";
}
