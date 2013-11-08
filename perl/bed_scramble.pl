#!/usr/bin/perl
# bed_scramble.pl - randomizes the positions of BED regions given a set of valid
# regions for each chromosome. This preserves the genome coverage while randomly varying start
# position. The antigap.bed file used for GREAT is a good source of valid regions.
# usage: bed_scramble.pl in_bed in_regions chr_weights out_bed

if(@ARGV != 4) {
	print "usage: bed_scramble.pl in_bed in_regions chr_weights out_bed\n";
} else {

$in_file = @ARGV[0];
$in_regions = @ARGV[1];
$chr_weights = @ARGV[2];
$out_file = @ARGV[3];

open(REGIONS, "<$in_regions");

while(<REGIONS>) {
	
	$line = $_;
	chomp($line);
	@region_split = split(/\t/, $line);
	$region = @region_split[1] . "\t" . @region_split[2];
	push(@{"@region_split[0]"}, $region);
	
}
close(REGIONS);

open(WEIGHTS, "<$chr_weights");

while(<WEIGHTS>) {
	
	$line = $_;
	chomp($line);
	
	push(@weights, $line);
	
}
close(WEIGHTS);

sub pick_chr {
	#pick a new chromosome
	$chr_rand = rand(1);
	$weight_start = 0;
	foreach(@weights) {
		
		@weight_split = split(/\t/, $_);
		$weight_end = @weight_split[1];
		
		if(($chr_rand > $weight_start) && ($chr_rand < $weight_end)) {
			$new_chr = @weight_split[0];
			last;
		} else {
			$weight_start = @weight_split[1];
		}
		
	}
	
	#determine the end of the chromosome
	($z, $chr_end) = split(/\t/, @{"$new_chr"}[-1]);
}



open(BED, "<$in_file");
open(OUTPUT, ">$out_file");

while(<BED>) {
	
	$line = $_;
	chomp($line);
	
	@bed_split = split(/\t/, $line);
	
	#find out the length of the regulatory domain
	$bed_length = @bed_split[2] - @bed_split[1];
	
	#select a random chromosome, get the end of the chromosome as a boundary for random selection
	pick_chr();
	
	$valid = 0;
	$check_counter = 0;
	
	#pick a random position, and check to see if it fits in a mappable region
	while($valid == 0) {
		
		#if the check fails too many times, this chromosome may not have a large enough mappable
		#region (frequent for chr_random), so choose another.
		if($check_counter == 20) {
			$check_counter = 0;
			pick_chr();
		} else {
			$check_counter++;
			$random_start = int(rand($chr_end));
			$random_end = $random_start + $bed_length;
			
			foreach(@{"$new_chr"}) {
				@check_split = split(/\t/, $_);
				
				if(($random_start >= @check_split[0]) && ($random_end <= @check_split[1])) {
					
					print OUTPUT $new_chr . "\t" . $random_start . "\t" . $random_end;
					
					if(@bed_split > 3) {
						for($x = 3; $x < @bed_split; $x++) {
							print OUTPUT "\t" . @bed_split[$x];
						}
					}
					
					print OUTPUT "\n";
					
					$valid = 1;
					last;
					
				}
					
			}
			
		}
		
	}
	

	
}

close(OUTPUT);
close(BED);

}