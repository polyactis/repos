use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;

# Open the sequence database
my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                                 -dsn     => 'dbi:mysql:gbrowse',
					       -user => 'gbrowse_loader',
					       -pass => 'loader',
					       -create => 1,
                                                 -write   => 1 );
