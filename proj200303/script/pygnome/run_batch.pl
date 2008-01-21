#!/usr/bin/perl -w

#start gimp as
# gimp  --no-interface --batch '(extension-perl-server 0 0 0)' &

use strict;

#read a file with bubbles definitions, ie.
#B_lectin|BNEWLINElectin|115|#9b1e00|rect|blue

open (IN, '/tmp/bubbles.def');

while (<IN>) {
  chomp;
  my @line = split (/\|/);
  my $name = $line[0];
  my $text = $line[1];
  my $width= $line[2];
  my $color = $line[3];
  my $shape= $line[4];
#  my $text_color = $line[5];
#  my $text_color2 = $text_color;
#  if ($text =~ /c$/ or $name eq "AAA" or $name eq "A1pp" or $name eq 'EXOIII') { #catalitic domains are red
#    if ($text_color eq  'blue') { #extracelular catalitic are two color
#      $text_color2 = 'red';
#    } else {
#      $text_color2= $text_color = 'red';  
#    }
#  }
  #run bubbles.pl and generate bubble in /tmp/bubbles/

 `./bubbles.pl -o $name -bgcolor "$color" -shape $shape -width $width -height 50 -text "$text" `;
}
# -text_color $text_color -text_color2 $text_color2`;
#}

