#!/usr/bin/perl
#bed_subtract.pl - compares two bed files and returns lines of the first file that do not overlap the second.
#Both bed files should be sorted first. Should work for any files with chr start end as the
#first 3 tab-separated columns. Faster to use the shorter of the two bed files as the second listed.
# usage: bed_overlap.pl bed1_in bed2_in overlap_out

if (@ARGV != 3) {
	print "usage: bed_overlap.pl bed1_in bed2_in overlap_out\n";
} else {
	
	$peaks_in = @ARGV[0];
	$bed_in = @ARGV[1];
	$bed_out = @ARGV[2];
		
	# first, read in the peaks to arrays for searching each chromosome
	
	open(PEAKS, "<$peaks_in");
	while(<PEAKS>) {
		
		$line = $_;
		chomp($line);
		
		@line_split = split(/\t/, $line);
		
		$test_line = $line . "\t0";
		
		push(@{"@line_split[0]"}, $test_line);

	}
	close(PEAKS);
	
	#Then, look through the bed file for positions that are in the regions
	$display_counter = 0;
	$total_counter = 0;
	
	$cur_chr = "Regions";
	
	open(BED, "<$bed_in");
	open(OUTPUT, ">$bed_out");
	
	while(<BED>) {
		
		$line = $_;
		chomp($line);
		
		#check through the regions for the current chromosome
		
		@split_bed = split(/\t/, $line);
		
		#If this is a new chromosome
		if(@split_bed[0] ne $cur_chr) {
			
			#Check to see if the array of locations from the previous chromosome still has
			#locations that need to be checked.
			
			if(@{"$cur_chr"} > 0) {
			
				#If so, for each location, split and look for the hit flag
				foreach(@{"$cur_chr"}) {
					
					@region_split = split(/\t/, $_);
					$hit_check = @region_split[-1];
					
					#if the hit flag is a 1, output the region.
					if($hit_check == 0) {
						$out_line = substr(@{"$cur_chr"}[$i], 0, -2);
						print OUTPUT $out_line . "\n";
					}
					
				}
				
			}
			
			#Then, switch over to the new chromosome.
			$cur_chr = @split_bed[0];
			
		}
		
		for ($i = 0; $i < @{"$cur_chr"}; $i++) {
			
			@region_split = split(/\t/, @{"$cur_chr"}[$i]);
			
			$start = @region_split[1];
			$end = @region_split[2];
			$hit_check = @region_split[-1];
			
			#check to see if the bed location falls in the region
			
			if ( @split_bed[2] < $start ) {
				#if the bed position is below the start of the current region
				#stop looping through regions. This should improve speed.
				last;
				
			} elsif (@split_bed[1] > $end) {
				#if the bed position is above the end of the current region, check to see
				#if it overlapped a previous region. If so, remove the hit check and output
				#the line. 
				
				if($hit_check == 0) {
					$out_line = substr(@{"$cur_chr"}[$i], 0, -2);
					print OUTPUT $out_line . "\n";
				}
				
				#Then, remove the region to speed up further comparisons.
				shift(@{"$cur_chr"});
				$i--;
				
			} elsif ( @split_bed[1] <= $end ) {
				#if it does fall in the region, change the hit flag at the end of the line.
				
				@{"$cur_chr"}[$i] = substr(@{"$cur_chr"}[$i],0,-1) . "1";
				
			} 
			
		}
		
	}
	
	#Before we're all done, we need to do one last check of the final chromosome to make sure
	#everything is properly output.
	
	if(@{"$cur_chr"} > 0) {
			
		#If so, for each location, split and look for the hit flag
		foreach(@{"$cur_chr"}) {
			
			@region_split = split(/\t/, $_);
			$hit_check = @region_split[-1];
			
			#if the hit flag is a 1, output the region.
			if($hit_check == 0) {
				$out_line = substr(@{"$cur_chr"}[$i], 0, -2);
				print OUTPUT $out_line . "\n";
			}
			
		}
		
	}
	
	close(BED);
	close(OUTPUT);
	
	
}