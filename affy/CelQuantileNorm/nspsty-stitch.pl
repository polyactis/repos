#!/usr/bin/perl -w

# 2005-2006 (c) Hin-Tak Leung
#
# This perl script merges NSP/STY runs with annotation info.

use Getopt::Std;
use IO::File;

my %opts;

getopts( 'i:o:m:h', \%opts );

if ( $opts{'h'} ) {
	print << 'USAGE';
Usage: 
    nspsty-stitch.pl -o stem -i annotAB.padded+sorted NSPdir1 STYdir1 \
                 NSPdir2 STYdir2 NSPdir3 STYdir3 ...

Options:
    -o <stem>      output files are written as stem_01.txt.gz, stem_02.txt.gz, etc
    -i <annot>     annotation file extract, listing the merged order, 
                       plus extra info of the snps.
    -m <manifest>  use manifest to map the plate/well to sample id. 
                       (not recommended for duplicate sample id's)
    -h             displaying this help message
USAGE
	exit 0;
}

unless ( $opts{'i'} ) { die "You have to specify a concatenated list file!\n" }
unless ( $opts{'o'} ) { die "You have to specify an output file stem name!\n" }

my %pw2sample = ();

if ( $opts{'m'} ) {

	# print "opening ", $opts{'m'}, "\n";
	open( MAPFILE, $opts{'m'} );
	while (<MAPFILE>) {
		my ( $sample, $platewell );
		( $sample, undef, undef, undef, $platewell, undef ) = split;
		$platewell                      = uc $platewell;
		$pw2sample{ $platewell . "_A" } = $sample . "_A";
		$pw2sample{ $platewell . "_B" } = $sample . "_B";
	}
}

#foreach my $key (keys %pw2sample) {
#    print $key, " ", $pw2sample{$key}, "\n";
#}

my $current_chrom = "";

#03/16/08 yh: open the annotation file
my $full_list = IO::File->new( $opts{'i'}, "r" );

my @NSPdirs = ();
my @STYdirs = ();

while (@ARGV) {
	push @NSPdirs, ( shift @ARGV );
	push @STYdirs, ( shift @ARGV );
}

if ( ( scalar @NSPdirs ) != ( scalar @STYdirs ) ) {
	die "You don't have the files in pairs!\n";
}

my $num_pairs = scalar @NSPdirs;

my @nspfiles = ();
my @styfiles = ();

my $outfile = undef;

my $snp          = undef;
my $nsp_snp      = undef;
my $sty_snp      = undef;
my $nsp_the_rest = "";
my $sty_the_rest = "";

while ( my $current_line = <$full_list> ) {
	chomp $current_line;
	( $snp, undef, $chrom_pos, undef, undef ) = split /\s+/, $current_line, 4;

	#print $snp, " ", $chrom_pos, "\n";
	( my $chrom = $chrom_pos ) =~ s/_.+//;

	my $read_nsp        = 0;
	my $read_sty        = 0;
	my $current_outname = "";

	if ( $chrom ne $current_chrom ) {
		if ($outfile) {
			$outfile->close();
			foreach my $file ( @nspfiles, @nspfiles ) {
				$file->close();
			}
		}
		$current_outname = $opts{'o'} . "_${chrom}.txt.gz";
		$outfile         = IO::File->new("|gzip -c > $current_outname");

		print $outfile
		  join( "\t", ( "AFFYID", "RSID", "pos", "AlleleA", "AlleleB" ) );

		@nspfiles = ();
		@styfiles = ();

		foreach my $nspdir (@NSPdirs) {
			push @nspfiles,
			  ( IO::File->new( $nspdir . "/NSP_chromosome_${chrom}.txt", "r" )
			  );
		}
		foreach my $stydir (@STYdirs) {
			push @styfiles,
			  ( IO::File->new( $stydir . "/STY_chromosome_${chrom}.txt", "r" )
			  );
		}

		for ( my $i = 0 ; $i < $num_pairs ; $i++ ) {
			my $nspheader = ${ nspfiles [$i] }->getline();
			my $styheader = ${ styfiles [$i] }->getline();
			if ( $nspheader ne $styheader ) {
				print $nspheader;
				print $styheader;
				die "header don't match\n";
			}
			else {
				chomp $nspheader;
				if ( $opts{'m'} ) {
					$nspheader =~ s/^\s+//;    #remove tab at beginning
					my @a = split /\s+/, $nspheader;
					my @b = map { $_ =~ s/^WTCCC01-//; $pw2sample{$_} } @a;
					print $outfile "\t", join( "\t", @b );
				}
				else {
					print $outfile $nspheader;
				}
			}
		}
		print $outfile "\n";
		&read_nsp();
		&read_sty();
		$current_chrom = $chrom;
	}

	$current_line =~ s/[012]._0*//;    # remove the chromosome part
	print $outfile $current_line;

	if ( ($nsp_snp) && ( $snp eq $nsp_snp ) ) {
		print $outfile $nsp_the_rest;
		&read_nsp();
	}
	elsif ( ($sty_snp) && ( $snp eq $sty_snp ) ) {
		print $outfile $sty_the_rest;
		&read_sty();
	}
	else {
		if ( $nsp_snp && $sty_snp ) {
			print "next snps $nsp_snp $sty_snp\n";
			die "sty/nsp not in sequence\n";
		}
	}
}

#close the last filehandle
if ($outfile) {
	$outfile->close();
}

sub read_nsp {
	$nsp_the_rest = "";
	$nsp_snp      = undef;
	foreach my $nspfile (@nspfiles) {
		my $line = $nspfile->getline();
		if ($line) {
			my ( $snp_here, $the_rest ) = split /\s+/, $line, 2;
			chomp $the_rest;
			if ( ($nsp_snp) && ( $snp_here ne $nsp_snp ) ) {
				die "not all nsp files are the same\n";
			}
			else {
				$nsp_the_rest .= ( "\t" . $the_rest );
				unless ($nsp_snp) {
					$nsp_snp = $snp_here;
				}
			}
		}
	}
	$nsp_the_rest .= "\n";
}

sub read_sty {
	$sty_the_rest = "";
	$sty_snp      = undef;
	foreach my $styfile (@styfiles) {
		my $line = $styfile->getline();
		if ($line) {
			my ( $snp_here, $the_rest ) = split /\s+/, $line, 2;
			chomp $the_rest;
			if ( ($sty_snp) && ( $snp_here ne $sty_snp ) ) {
				die "not all sty files are the same\n";
			}
			else {
				$sty_the_rest .= ( "\t" . $the_rest );
				unless ($sty_snp) {
					$sty_snp = $snp_here;
				}
			}
		}
	}
	$sty_the_rest .= "\n";
}


