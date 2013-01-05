#!/usr/bin/perl
# paired_min_qc_filter.pl filters out paired fastq entries with a given number of phred scores below the given thresshold.
# only pairs for which both entries pass the filter will be written to the output files.
# usage: paired_min_qc_filter.pl in_file out_file score_thresh num_thresh

if(@ARGV != 6) {
	print "usage: paired_min_qc_filter.pl in_file1 out_file1 in_file2 out_file2 score_thresh num_thresh\n";
} else {

$in_file1 = @ARGV[0];
$out_file1 = @ARGV[1];
$in_file2 = @ARGV[2];
$out_file2 = @ARGV[3];
$score = @ARGV[4];
$num = @ARGV[5];

open(INPUT1, "<$in_file1");
open(INPUT2, "<$in_file2");

open(OUTPUT1, ">$out_file1");
open(OUTPUT2, ">$out_file2");

$x = 0;

while(<INPUT1>) {
	
	$line = $_;
	chomp($line);
	$line2 = <INPUT2>;
	chomp($line2);
	
	if($x == 3) {
		
		push(@entry, $line);
		push(@entry2, $line2);
		
		@checkscores = split(//,@entry[3]);
		@checkscores2 = split(//,@entry2[3]);
		
		$below = 0;
		$below2 = 0;
		
		for(@checkscores) {
			
			$check = ord($_) - 33;
			
			if($check < $score) {
				$below++;
			}
			
			if($below > $num) {
				last;
			}
			
		}
		
		for(@checkscores2) {
			
			$check = ord($_) - 33;
			
			if($check < $score) {
				$below2++;
			}
			
			if($below2 > $num) {
				last;
			}
			
		}
		
		if(($below <= $num) && ($below2 <= $num)) {
			for(@entry) {
				print OUTPUT1 $_ . "\n";
			}
			for(@entry2) {
				print OUTPUT2 $_ . "\n";
			}
		}
		
		$x = 0;
		
		@entry = ();
		@entry2 = ();
		
	} else {
		push(@entry, $line);
		push(@entry2, $line2);
		$x++;
	}
	
}

close(INPUT1);
close(INPUT2);
close(OUTPUT1);
close(OUTPUT2);


}