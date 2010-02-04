"""
Parses the link into snp files
"""

fi = open("/etc/apache2/gbrowse.conf/arabidopsis1.conf")
fo = open("/etc/apache2/gbrowse.conf/arabidopsis.conf","w")

insertedtext="""link            = sub {
                    my $f = shift;
                    my $score = $f->score();
                    my @met = split("_",$f->method());
                    my @snp = split("_",$f->load_id());
                    return sprintf("http://banyan.usc.edu:5000/SNP/?chromosome=%d&position=%d&call_method_id=%d&phenotype_method_id=%d&analysis_method_id=%d&score=%.6f",$snp[0],$snp[1],$met[0],$met[1],$met[2],$score);
                }
"""

for line in fi:
	fo.write(line)
	if line[0]=="[": # start of a new track
		if line[len(line)-6:len(line)-1]=="_SNP]":
			fo.write(insertedtext)

fo.close()
fi.close()
