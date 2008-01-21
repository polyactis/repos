#!/usr/bin/perl -w

=head1 NAME

  bubbles.pl - Draw SMART domain bubbles

=head1 SYNOPSIS

  bubbles.pl

=head1 DESCRIPTION

 Draw SMART domain bubbles with gimp. Gimp must have a perl server installed and 
 started (via menu Xtns->Perl->Server). For batch running start gimp as:

 gimp  --no-interface --batch '(extension-perl-server 0 0 0)' &

 and run bubbles.pl with -output FILENAME
 
 In batch mode output format depends on file extension (ie. if you want PNG output
 use bubbles.pl -output /tmp/test.png)

 If you want bubble text in two lines, separate them with string NEWLINE, ie.  
 CalxNEWLINEBeta will put 'Calx' in first line and 'Beta' in second.
 
 See bubbles.pl --help for full list of command line arguments.

=cut

=head1 AUTHOR

Ivica Letunic <letunic@embl-heidelberg.de>
http://www.embl-heidelberg.de/~letunic

=cut

=head1 COPYRIGHT

Copyright (c) 2001 Ivica Letunic

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

=cut

use Gimp ":auto";
use Gimp::Fu;

sub start {

  my ($shape, $bgcolor, $text_color, $text_color2, $text, $width, $height, $font)= @_;
  my $gradient_type=0; #0=linear 2=radial
  my $fgcolor="#ffffff";
  my $bgcolor2= "White";
  my $border=5;
  my $corner_radius=0;
  my $ellipse=0;
  my ($no_points, @points_list);
  my $rotate=0; # rotate the text 90 deg? (used for narrow bubbles)
  my $offset_x=3; #for precise positioning of text
  my $offset_y=3;
  my $antialias=1;
  my ($text1, $text2);

  #define text box size
  #if the width is smaller than 30, force "rect" shape and rotate the text

  if ($width <= 30 or ($height == 60 and $width <= 50)) {
    $shape="rect";
    $tx1= 0; $ty1=0; 
    $tx2= int($height*4/5); $ty2= int($width*4/5);
    $rotate=90;
    $offset_x= -5;
    $antialias=0;
  } else {
    $tx1= int($width/6); $ty1=int($height/6); 
    $tx2= int($width*5/6); $ty2= int($height*5/6);
    $ty2=$ty1+12 if (($ty2 - $ty1) > 12);
  }
  
  #define the shapes of the bubbles
  #you can new shapes here

  if ($shape eq "dia") { #diamond
    $no_points=4; 
    @points_list=[0,$height/2, $width/2,0, $width,$height/2, $width/2, $height]; 
    $tx1*= 0.9; $ty1*=0.9; $tx2*=0.9; $ty2*=0.9; #adjust text box size
  } 
  if ($shape eq "hex_h") {#Hexagon (horiz)
    $no_points=6; 
    @points_list=[0,$height/2, $width/3,0, $width*2/3,0, $width, $height/2, $width*2/3,($height-1), $width/3,($height-1)]; 
  } 

  if ($shape eq "hex_v") {  #Hexagon (vert)
    $no_points=6; 
    @points_list=[0,$height/3, $width/2,0, $width,$height/3, $width, $height*2/3, $width/2,$height, 0,$height*2/3 ];  
  }

  if ($shape eq "oct") {  #Octagon
    $no_points=8; 
    @points_list=[0,$height*2/3, 0,$height/3, $width/3,0, $width*2/3,0, $width,$height/3, $width,$height*2/3, $width*2/3,($height-1), $width/3,($height-1) ]; 
  }

  if ($shape eq "pent_d") {  #Pentagon (down)
    $no_points=5; 
    @points_list=[0,$height*2/3, $width/6,0, $width*5/6,0, $width,$height*2/3,  $width/2,$height ]; 
  }

  if ($shape eq "pent_l") {  #Pentagon (left)
    $no_points=5; 
    @points_list=[0,$height/2, $width/3,0, $width,$height/6, $width,$height*5/6,  $width/3,$height ];  
  }

  if ($shape eq "pent_r") {  #Pentagon (right)
    $no_points=5; 
    @points_list=[0,$height*5/6, 0,$height/6, $width*2/3,0, $width,$height/2,  $width*2/3,$height ]; 
  }

  if ($shape eq "pent_u") {  #Pentagon (up)
    $no_points=5; 
    @points_list=[0,$height/3, $width/2,0, $width,$height/3, $width*5/6,($height-1),  $width/6,($height-1) ]; 
  }

  if ($shape eq "rect") {  #Rectangle
    $no_points=4; 
    @points_list=[0,0, $width,0, $width,($height-1), 0,($height-1)]; 
  }

  if ($shape eq "rect_r") {  #Rectangle (rounded)
    $no_points=4; 
    @points_list=[0,0, $width,0, $width,($height-2), 0,($height-2)]; 
    $corner_radius=1; 
  }

  if ($shape eq "tri_r") {  #Triangle (right)
    $no_points=3;
    @points_list=[0,0, $width,$height/2, 0,$height]; 
    $tx1= int($width/8); $ty1=int($height/3); 
    $tx2= int($width*2/3); $ty2= int($height*2/3);
    $ty2=$ty1+12 if (($ty2 - $ty1) > 12);
    $offset_x= -((($width  - ($tx2-$tx1))/2)-$width/8); #for forcing left alignment of text
 }

  if ($shape eq "tri_l") {  #Triangle (left)
    $no_points=3;
    @points_list=[0,$height/2, $width,0, $width,$height]; 
    $tx1= int($width/3); $ty1=int($height/3); 
    $tx2= int($width*7/8); $ty2= int($height*2/3);
    $ty2=$ty1+12 if (($ty2 - $ty1) > 12);
    $offset_x= ((($width  - ($tx2-$tx1))/2));#-$width/8); #for forcing left alignment of text
  }

 if ($shape eq "elli") {$gradient_type=2;$ellipse=1;} #Ellipse

  #parse the text for NEWLINE

  if ($text =~ /(.*)NEWLINE(.)/) { #because perl-gimp cannot pass newlines in strings
    $text =~ s/NEWLINE/\n/;
    $ty2= $ty2 + ($ty2 - $ty1) unless ($width <= 30); #double the height
  }

  
  # Create a new image 
  my $img = gimp_image_new($width, $height, RGB);
  
  # Create a new layer for the background and add it to the image
  my $background = gimp_layer_new($img, $width, $height,
				  RGB, "Background", 100,
				  NORMAL_MODE);
  gimp_layer_add_alpha($background);
  gimp_image_add_layer($background, 1);
  gimp_palette_set_background("#ffffff");
  gimp_edit_fill($background, WHITE_IMAGE_FILL);

  #make a selection depending on the shape of the bubble
  if ($ellipse eq "0") {
    $img->free_select($no_points, @points_list, 0, 1, 0, 0);
    round_corners ($img, $width, $height) if ($corner_radius == 1);
  } else {
    $img->ellipse_select(0,0, $width, $height-1, 0, 1, 0, 0);
  }

  #fill selection with color gradient, shrink the selection and fill again in other direction
  gimp_palette_set_background($bgcolor);
  gimp_palette_set_foreground($bgcolor2);
  gimp_blend($background,0,0,$gradient_type, 100,0,0,0,0,0, 0,$height/2, $width,$height/2);
  my $big_sel=gimp_selection_save($img);
  gimp_selection_shrink($img, $height/15);
  gimp_blend($background,0,0,$gradient_type, 100,0,0,0,0,0,  $width,$height/2, 0,$height/2);

  #add the border line
  gimp_selection_load($big_sel);
  gimp_palette_set_foreground($bgcolor);
  gimp_brushes_set_brush("Circle (01)");
  gimp_edit_stroke($background);

  #blur
  gimp_selection_none($img);
  plug_in_gauss_rle(0, $background, 1,1,1);

  #add text
  #first find font size that fits best into text box then add the text
  gimp_palette_set_foreground($text_color);
  @list=(1, $tx1,$ty1,$tx2,$ty2);
  $font= fit_font($img,$background,$font,$text,$rotate,$antialias, $text_color, @list);
  @fontdesc = split /-/, $fontname;
  $size=int($fontdesc[7]);
  $fontdesc[7]=$size;
  $font = join "-", @fontdesc;
  $tmplay=$img->plug_in_dynamic_text(0 ,$text, $antialias, 1, $rotate, -1, $text_color, 0, $font);

  #center the text in the bubble
  $tmplay->set_offsets((($img->width  - $tmplay->width )/2+$offset_x), 
			 (($img->height - $tmplay->height)/2)+$offset_y);

  #create color gradient on the text
  gimp_palette_set_background($text_color2);
  $tmplay->blend(0,0,0, 100,0,0,0,0,0, 0,$tmplay->height/2, $tmplay->width, $tmplay->height/2);  

  #add text glow
  $img->perl_fu_add_glow( $tmplay, "#ffffff", 2);

  $img->image_flatten();
  
  return $img;
}


#used for rounding corners of a selection (for rect_r shape)
sub round_corners {
  my ($img, $width, $height)=@_;
  my ($y_round, $x_round);
  $y_round=$height/2;
  $x_round=$width/2;
  eval { $img->undo_push_group_start };
  my @bounds = $img->selection_bounds;

  # recreate the selection
  $img->rect_select($bounds[1], $bounds[2], $bounds[3]-$bounds[1], $bounds[4]-$bounds[2], 0, 0, 0.5);

  # cut out the corners
  $img->rect_select($bounds[1], $bounds[2], $x_round/2, $y_round/2, 1, 0, 0.5); 
  $img->rect_select($bounds[3]-int($x_round/2), $bounds[2], int($x_round/2), $y_round/2, 1, 0, 0.5); 
  $img->rect_select($bounds[3]-int($x_round/2), $bounds[4]-int($y_round/2), int($x_round/2), int($y_round/2), 1, 0, 0.5); 
  $img->rect_select($bounds[1], $bounds[4]-int($y_round/2), $x_round/2, int($y_round/2), 1, 0, 0.5); 

  # add them back as elipses
  
  $img->ellipse_select($bounds[1], $bounds[2], $x_round, $y_round, 0, 1, 0, 0.5); 
  $img->ellipse_select($bounds[3]-$x_round, $bounds[2], $x_round, $y_round, 0, 1, 0, 0.5); 
  $img->ellipse_select($bounds[3]-$x_round, $bounds[4]-$y_round, $x_round, $y_round, 0, 1, 0, 0.5); 
  $img->ellipse_select($bounds[1], $bounds[4]-$y_round, $x_round, $y_round, 0, 1, 0, 0.5); 
  eval { $img->undo_push_group_end };
}

sub growfont {
        ($fontname, $plussize) = @_;
        @fontdesc = split /-/, $fontname;
        $fontdesc[8] eq "*" ?  ($fontdesc[7] += $plussize/72) : ($fontdesc[8]+= $plussize);
        $outname = join "-", @fontdesc;
        return $outname;
        }  

#find font size that fits the selection
sub fit_font {
  my($img,$layer,$xlfd,$string,$rotate, $antialias, $text_color, @list) =@_;
  ($_,$x1,$y1,$x2,$y2) = @list;
  $width = $x2-$x1;
  $height = $y2-$y1;   
#  my @extents=Gimp->text_get_extents_fontname($string,xlfd_size($xlfd),$xlfd);
#  $growsize = ($extents[0]<$width && $extents[1]<$height) ? 288 : -288;
#  if ($growsize > 0 ) {
#    while ($extents[0]<$width && $extents[1]<$height) {
#      $xlfd = growfont($xlfd,$growsize);
#      @extents=Gimp->text_get_extents_fontname($string,xlfd_size($xlfd),$xlfd);
#    }
#    $xlfd = growfont($xlfd, -$growsize);
#  } else {
#    while ($extents[0]>$width || $extents[1]>$height) {
#      $xlfd = growfont($xlfd,$growsize);
#      @extents=Gimp->text_get_extents_fontname($string,xlfd_size($xlfd),$xlfd);
#    }
#  } while ($extents[0]<$width && $extents[1]<$height) {
#    $xlfd = growfont($xlfd,144); # precision for the last bit
#    @extents=Gimp->text_get_extents_fontname($string,xlfd_size($xlfd),$xlfd);
#  }
  $xlfd = growfont($xlfd, -144); $tmplay = $layer->text_fontname($x1,$y1,$string,0,1,xlfd_size($xlfd), $xlfd);
  $width2=$tmplay->width;
  $height2=$tmplay->height;	# X returns crap, so fine tune it here.
  #print "$width2, $height2:$width, $height\n";
  while ($width2<$width && $height2<$height) {
    $tmplay->remove;
    $xlfd = growfont($xlfd,288);
    $tmplay=$layer->text_fontname($x1,$y1,$string,0,1,xlfd_size($xlfd), $xlfd);
    $width2=$tmplay->width;
    $height2=$tmplay->height;
  }
  
  while ($width2>$width || $height2>$height) {
    $tmplay->remove;
    $xlfd = growfont($xlfd,-72);
    $tmplay=$layer->text_fontname($x1,$y1,$string,0,$antialias,xlfd_size($xlfd), $xlfd);
    $width2=$tmplay->width;
    $height2=$tmplay->height;
  }
  while ($width2<$width && $height2<$height) {
    $tmplay->remove;
    $xlfd = growfont($xlfd,12);
    $tmplay=$layer->text_fontname($x1,$y1,$string,0,$antialias,xlfd_size($xlfd), $xlfd);
    $width2=$tmplay->width;
    $height2=$tmplay->height;
  }
  $tmplay->remove;
  $xlfd = growfont($xlfd,-12);
  return $xlfd;
  
}

# register the script
register "bubbles", "Draw SMART bubbles", "Draw SMART bubbles",
    "Ivica Letunic", "Ivica Letunic",
    "2001-03-21",
    "<Toolbox>/Xtns/Perl-Fu/SMART Bubbles", 
    "*",
    [
    
    
     [PF_STRING, "shape",     "Bubble Shape", "hex_h"],
     [PF_COLOR,  "bgcolor", "Background color", [40,180,160]],
     [PF_COLOR,  "text_color", "Text color", [255,00,00]],
     [PF_COLOR,  "text_color2", "Text color2", [00,00,255]],
     [PF_STRING, "text",     "Bubble Text", "0"],
     [PF_INT, "width",     "Bubble Width", 150],
     [PF_INT, "height",     "Bubble Height", 60],
     [PF_FONT, "font",     "Font", "-*-*-*-*-*-*-16-*-*-*-p-*-iso8859-9"],
    ],
    \&start;

# Handle over control to gimp
exit main();
