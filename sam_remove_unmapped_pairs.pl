#!/usr/bin/perl
#sam_remove_unmapped_pairs.pl - removes SAM entries with flags 77 or 141 which indicate that
#the read was not aligned.
#usage: sam_remove_unmapped_pairs.pl in_file out_file

if(@ARGV != 2) {
	print "usage: sam_remove_unmapped_pairs.pl in_file out_file\n";
}

$in_file = @ARGV[0];
$out_file = @ARGV[1];

open(IN, "<$in_file");
open(OUTPUT, ">$out_file");


while(<IN>) {

	$line = $_;
	chomp($line);
	@columns = split(/\t/, $line);

	if ((@columns[1] != 77) && (@columns[1] != 141)) {

		print OUTPUT $line . "\n";
	
	}

}

close(IN);
close(OUTPUT);