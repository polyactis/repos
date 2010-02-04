use Bio::DB::SeqFeature::Store;
use Bio::DB::SeqFeature::Store::GFF3Loader;

# Open the sequence database
my $db      = Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::mysql',
                                                 -dsn     => 'dbi:mysql:gbrowse',
					       -user => 'gbrowse_loader',
					       -pass => 'loader',
                                                 -write   => 1 );
#$db->init_database();

my $loader = Bio::DB::SeqFeature::Store::GFF3Loader->new(-store    => $db,
							   -verbose  => 1);

$loader->load($ARGV[0]);
