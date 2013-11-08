#!/usr/bin/perl


if(@ARGV != 2) {
	print "usage: fasta_gc_percent.pl in.fa out_gc.txt\n"
} else {

$in_fa = @ARGV[0];
$out_gc = @ARGV[1];

open(INPUT,"<$in_fa");
open(OUTPUT,">$out_gc");

$z = 0;

while(<INPUT>) {
	$line = $_;
	chomp($line);
	$first = substr($line,0,1);

	if($first eq ">" && $z == 1) {
		$l = length($seq);
		$c = grep /[GgCc]/,split(//,$seq);
		$out_percent = $c / $l;

		$n = grep /[Nn]/, split(//,$seq);

		if($n > 0) {
			print OUTPUT $name . "\tN\n"
		} else {
			print OUTPUT $name . "\t" . $out_percent . "\n";
		}

		$name = substr($line,1);
		$seq = "";
	} else {
		$seq .= $line;
	}
	if($z == 0) {
		$name = substr($line,1);
		$z = 1;
	}

}

$l = length($seq);
$c = grep /[GgCc]/,split(//,$seq);
$out_percent = $c / $l;
print OUTPUT "$name\t$out_percent\n";

}