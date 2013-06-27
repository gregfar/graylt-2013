#!/usr/bin/perl
# bed_sort.pl - sorts a BED file based on chromosome and start position.
# Because hashes are used for sorting, !!identical lines will be removed!!
# usage: bed_sort.pl in_bed out_bed

if (@ARGV != 2) {

	print "usage: bed_sort.pl in_bed out_bed\n";

} else {
	
	$in_file = @ARGV[0];
	$out_file = @ARGV[1];
	
	open(INPUT, "<$in_file");
	open(OUTPUT, ">$out_file");
	while(<INPUT>) {
		
		$line = $_;
		chomp($line);
		
		@hash_check = split(//,$line);
		if(@hash_check[0] eq "#") {
		} elsif (@hash_check[0] eq "t") {
			print OUTPUT $line . "\n";
		} else {
			@line_split = split(/\t/,$line);
			
			#hash of unique chromosomes
			$chr_list{@line_split[0]} = 0;
			
			#put the line in a hash for its chromosome
			${"@line_split[0]"}{$line} = @line_split[1];
		}	
	}
	
	close(INPUT);
	
	
	
	#sort the chromosomes
	
	@sorted_chr = sort keys %chr_list;
	
	foreach(@sorted_chr) {
		
		foreach $key (sort { ${"$_"}{$a} <=> ${"$_"}{$b} } keys %{"$_"}) {
			print OUTPUT $key . "\n";
		}
		
		undef %{"$_"};
		
	}
	
	close(OUTPUT);
	
}