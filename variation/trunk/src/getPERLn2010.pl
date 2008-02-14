#!/usr/bin/perl -w
#**************************
# formatCalls.pl
# by Tina Hu 
# created 12.03.07
#**************************
### - reformat data from Yan into a 96 x 250 matrix
### - also at how well things match up by annotation
### - comparisons to make:
###   1) yan vs. yan (replicability)
###   2) yan vs. PERL
###   3) yan vs. 2010
###   4) yan vs. CHLA
###   5) CHLA vs. PERL
###   5) CHLA vs. 2010
#**************************
use strict;
use tinaDBI;

### SETUP
#****************************************************
my $dbh_chip = tinaDBI::connectDB("chip");
my $dbh_at = tinaDBI::connectDB("at");

my %accessions;
open (FILE,'aroma_calls.121408.csv');
#open (FILE,'accname2file.csv');
while (<FILE>) { chomp $_; next if ($_ =~ /file_name,accession,file/);
   my ($file_acc,$db_acc,$file) = split(/,/,$_);
   $accessions{$file_acc} = $db_acc;
} close FILE;


### get the 2010 and PERL calls for each position
my @calls;
open (FILE,"250K_snp_pos.txt") || die "ERROR: Can't open 250K_snp_pos.txt: $!\n";
foreach my $snp_pos (<FILE>) { next if ($snp_pos =~ /^chr/); chomp $snp_pos;
   my ($chr,$pos) = split(/ /,$snp_pos);

   my $snp_id = $chr . ('0' x (8-length($pos))) . $pos;
   my ($in_final) = tinaDBI::count($dbh_chip,'orthologous.rc_for_allele_freq_analysis_2','id',$snp_id);
   my ($annotation) = tinaDBI::select($dbh_chip,'type','chip.snp_combined_may_9_06_no_van',"id=$snp_id AND ecotype='col-0'");

   ##### 2010 --- using database at
   ### select allele.base from locus,allele,genotype,accession
   ### where locus.chromosome=$chrom AND locus.position=$pos AND accession.name='Lov-5'
   ###       AND locus.id=allele.locus AND genotype.accession=accession.id AND genotype.allele=allele.id
   my %calls2010 = tinaDBI::select($dbh_at,'accession.name,allele.base','locus,allele,genotype,accession',
      "genotype.accession=accession.id AND genotype.allele=allele.id AND locus.id=allele.locus AND locus.chromosome=$chr AND locus.position=$pos");

   ##### PERL --- using database chip
   ### SELECT mbml98 FROM chip.snp_combined_may_9_06_no_van s where ecotype='lov-5' and chromosome=1 and position=29291
   ### SELECT specific_type FROM orthologous.rc_for_allele_freq_analysis_2 where id=$id --- also get annotation
   my %callsPERL = tinaDBI::select($dbh_chip,'ecotype,mbml98','chip.snp_combined_may_9_06_no_van',"chromosome=$chr AND position=$pos");
   $callsPERL{'Bay-0'} = $callsPERL{'bay-0'}; delete $callsPERL{'bay-0'};
   $callsPERL{'Bor-4'} = $callsPERL{'bor-4'}; delete $callsPERL{'bor-4'};
   $callsPERL{'Br-0'} = $callsPERL{'br-0'}; delete $callsPERL{'br-0'};
   $callsPERL{'Bur-0'} = $callsPERL{'bur-0'}; delete $callsPERL{'bur-0'};
   $callsPERL{'C24'} = $callsPERL{'c24'}; delete $callsPERL{'c24'};
   $callsPERL{'Col-0'} = $callsPERL{'col-0'}; delete $callsPERL{'col-0'};
   $callsPERL{'Cvi-0'} = $callsPERL{'cvi-0'}; delete $callsPERL{'cvi-0'};
   $callsPERL{'Est-1'} = $callsPERL{'est-1'}; delete $callsPERL{'est-1'};
   $callsPERL{'Fei-0'} = $callsPERL{'fei-0'}; delete $callsPERL{'fei-0'};
   $callsPERL{'Got-7'} = $callsPERL{'got-7'}; delete $callsPERL{'got-7'};
   $callsPERL{'Ler-1'} = $callsPERL{'ler-1'}; delete $callsPERL{'ler-1'};
   $callsPERL{'Lov-5'} = $callsPERL{'lov-5'}; delete $callsPERL{'lov-5'};
   $callsPERL{'NFA-8'} = $callsPERL{'nfa-8'}; delete $callsPERL{'nfa-8'};
   $callsPERL{'RRS-10'} = $callsPERL{'rrs-10'}; delete $callsPERL{'rrs-10'};
   $callsPERL{'RRS-7'} = $callsPERL{'rrs-7'}; delete $callsPERL{'rrs-7'};
   $callsPERL{'Shahdara'} = $callsPERL{'sha'}; delete $callsPERL{'sha'};
   $callsPERL{'Tamm-2'} = $callsPERL{'tamm-2'}; delete $callsPERL{'tamm-2'};
   $callsPERL{'Ts-1'} = $callsPERL{'ts-1'}; delete $callsPERL{'ts-1'};
   $callsPERL{'Tsu-1'} = $callsPERL{'tsu-1'}; delete $callsPERL{'tsu-1'};
   $callsPERL{'Van-0'} = $callsPERL{'van-0'}; delete $callsPERL{'van-0'};

   my $calls_from_2010 = '';
   my $calls_from_PERL = '';
   foreach my $file_acc (sort {$a cmp $b} keys %accessions) {
      my $db_acc = $accessions{$file_acc};
      my ($call_2010,$call_PERL) = ('?','?');

      if (exists $calls2010{$db_acc}) { $call_2010 = $calls2010{$db_acc} if ($calls2010{$db_acc} =~ /[ATGC]/); }
      if (exists $callsPERL{$db_acc}) { $call_PERL = $callsPERL{$db_acc} if ($callsPERL{$db_acc} =~ /[ATGC]/); }
      $calls_from_2010 .= ',' . $call_2010;
      $calls_from_PERL .= ',' . $call_PERL;
   }

   push @calls,"$snp_id,$in_final,$annotation,$chr,$pos" . $calls_from_2010 . $calls_from_PERL . "\n";
} close FILE;


my @acc_names = sort {$a cmp $b} keys %accessions;
open (OUT,">250K_PERL_2010.aroma_calls.121408.csv");
print OUT "snp_id,published,annotation,chr,pos,";
print OUT join('2010,',@acc_names) . '_2010,';
print OUT join('PERL,',@acc_names) . '_PERL' . "\n";
print OUT @calls;
close OUT;


tinaDBI::disconnectDB($dbh_chip);
tinaDBI::disconnectDB($dbh_at);
exit;
