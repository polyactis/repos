#!/usr/bin/env perl
# Usage example: load.pl hostname dbname db_user db_passwd gff_file

use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;


# Open the sequence database
my $db  = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                        -dsn     => 'dbi:mysql:database='.$ARGV[1].';host='.$ARGV[0],
                        -user => $ARGV[2],
					    -pass => $ARGV[3],
                        -write   => 1 );
#$db->init_database();

my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(-store    => $db,
							   -verbose  => 1);

$loader->load($ARGV[4]);
