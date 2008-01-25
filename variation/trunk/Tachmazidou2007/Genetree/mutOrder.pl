#!/usr/bin/perl -w

use strict;
use warnings;

# exit with the proper diagnostic message in case that program 
# has not been called with the proper arguments (filename)


my $input_file = $ARGV[0];                      # input file
my $output_file = "SNPsOrder.txt";			# output file
my $def_sep = "\n";                      	      # default input record separator
my @entries = ();                        	      # array for keeping the input file
my $site_labels = "";                           # the section with the site labels

unlink "SNPsOrder.txt";
open(IN, "<$input_file") || die "Can't open : $!";
# change the input record separator so that each element of @entries is an entry in the input file
$/ = "\n/";
@entries = <IN>;
close(IN);
# change back to default input record separator
$/ = $def_sep;

open(OUT, ">$output_file") || die "Can't write : $!";
$site_labels = $entries[2];
while ($site_labels =~ /[\d\.]+ ([\d\.]+) moveto\n\(([^\)]+)\) show\n/gm) {
  print OUT $2."\t".$1."\n";
}
close(OUT);
