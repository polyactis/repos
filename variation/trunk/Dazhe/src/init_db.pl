#!/usr/bin/env perl
# 2010-2-18 Usage example: init_db.pl hostname dbname db_user db_passwd

use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;

# Open the sequence database
my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                                 -dsn     => 'dbi:mysql:database='.$ARGV[1].';host='.$ARGV[0],
                            -user => $ARGV[2],
            	            -pass => $ARGV[3],
		                    -create => 1,
                                                 -write   => 1 );
