#!/usr/bin/perl -w

use strict;
use warnings;

# exit with the proper diagnostic message in case that program 
# has not been called with the proper arguments (filename)

if ($#ARGV != 1) {
  die "Incorrect arguments given!\n";
}

my $input_file = $ARGV[0]; 			# input file
my $output_file = $ARGV[1];			# the output file
my $def_sep = "\n";                      	# default input record separator
my @entries = ();                        	# array for keeping the input file

open(IN, "<$input_file") || die $!;
# change the input record separator so that each element of @entries is an entry in the input file
@entries = <IN>;
close(IN);

open(OUT, ">$output_file") || die $!;
foreach (@entries) {
	if (/^(\d+)[ \t]+(.*)/) {
		print OUT $1." : ".$2."\n";
	}
}
close(OUT);

