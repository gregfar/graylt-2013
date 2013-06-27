#!/usr/bin/perl
# makemotif.pl - generates a MEME minimal motif file for the input sequence
# usage: makemotif.pl motif motif_name output_file

if(@ARGV < 3) {
	print "usage: makemotif.pl motif motif_name output_file\n";
} else {

#read in the arguments
$motif = @ARGV[0];
$motif_name = @ARGV[1];
$out_file = @ARGV[2];

@motif_array = split(//, $motif);

open OUTPUT, ">", $out_file;

print OUTPUT "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from\nA 0.250 C 0.250 G 0.250 T 0.250\n\nMOTIF " . $motif_name . "\n";

print OUTPUT "letter-probability matrix: alength= 4 w= " . @motif_array . " nsites= " . @motif_array . " E= 0\n";

foreach(@motif_array) {
	if( $_ eq "A" || $_ eq "a") {
		print OUTPUT " 1.000000  0.000000  0.000000  0.000000 \n"
	} elsif ($_ eq "C" || $_ eq "c") {
		print OUTPUT " 0.000000  1.000000  0.000000  0.000000 \n"
	} elsif ($_ eq "G" || $_ eq "g") {
		print OUTPUT " 0.000000  0.000000  1.000000  0.000000 \n"
	} elsif ($_ eq "T" || $_ eq "t") {
		print OUTPUT " 0.000000  0.000000  0.000000  1.000000 \n"
	} elsif ($_ eq "N" || $_ eq "n") {
		print OUTPUT " 0.250000  0.250000  0.250000  0.250000 \n"
	} elsif ($_ eq "W" || $_ eq "n") {
		print OUTPUT " 0.500000  0.000000  0.000000  0.500000 \n"
	} elsif ($_ eq "S" || $_ eq "n") {
		print OUTPUT " 0.000000  0.500000  0.500000  0.000000 \n"
	} elsif ($_ eq "M" || $_ eq "n") {
		print OUTPUT " 0.500000  0.500000  0.000000  0.000000 \n"
	} elsif ($_ eq "K" || $_ eq "n") {
		print OUTPUT " 0.000000  0.000000  0.500000  0.500000 \n"
	} elsif ($_ eq "R" || $_ eq "n") {
		print OUTPUT " 0.500000  0.000000  0.500000  0.000000 \n"
	} elsif ($_ eq "Y" || $_ eq "n") {
		print OUTPUT " 0.000000  0.500000  0.000000  0.500000 \n"
	} elsif ($_ eq "B" || $_ eq "n") {
		print OUTPUT " 0.000000  0.333333  0.333333  0.333333 \n"
	} elsif ($_ eq "D" || $_ eq "n") {
		print OUTPUT " 0.333333  0.000000  0.333333  0.333333 \n"
	} elsif ($_ eq "H" || $_ eq "n") {
		print OUTPUT " 0.333333  0.333333  0.000000  0.333333 \n"
	} elsif ($_ eq "V" || $_ eq "n") {
		print OUTPUT " 0.333333  0.333333  0.333333  0.000000 \n"
	}
	
}

close(OUTPUT);

}