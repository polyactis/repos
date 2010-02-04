#!/usr/bin/env python
import os, math

def _run_():   
    user_250k=None
    passwd_250k=None
    host_250k='papaya'
    call_method_ids=["32"]
    analysis_method_ids=["1"]
    
    user_gb='gbrowse_loader'
    passwd_gb='loader'
    host_gb='mahogany'
    source_name='Nordborg'
    
    import MySQLdb
    print "Connecting to stock_250k, host="+host_250k
    if not user_250k:
        import sys
        sys.stdout.write("Username: ")
        user_250k = sys.stdin.readline().rstrip()
    if not passwd_250k:
        import getpass
        passwd_250k = getpass.getpass()
    try:
        conn_250k = MySQLdb.connect (host = host_250k, user = user_250k, passwd = passwd_250k, db = "stock_250k")
    except MySQLdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit (1)
    cursor_250k = conn_250k.cursor ()
    
    # Create the Phenotype and Association Method naming dictionaries
    analysis_method_dict={"1":"KW","7":"Emma"}
    phenotype_method_dict={}
    sqlcommand="select id, short_name from phenotype_method"
    cursor_250k.execute(sqlcommand)
    while(1):
        row = cursor_250k.fetchone()
        if not row:
            break
        phenotype_method_dict[str(row[0])]=str(row[1])
    
    # Retrieve the filenames
    print "Fetching data"
    sqlcommand="select call_method_id, phenotype_method_id, analysis_method_id, filename from results_method where ("
    sqlcommand+=" OR ".join(["call_method_id="+i for i in call_method_ids])+") AND ("
    sqlcommand+=" OR ".join(["analysis_method_id="+i for i in analysis_method_ids])+")"
    numRows = int(cursor_250k.execute(sqlcommand))
    
    # Store them in a file first? or use them straight away?
    # Try straight away
    
    # Connect to the gbrowse database
    try:
        conn_gb = MySQLdb.connect (host = host_gb, user = user_gb, passwd = passwd_gb, db = "gbrowse")
    except MySQLdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])
        sys.exit (1)
    cursor_gb = conn_gb.cursor()    
        
    while(1):
        rawrow = cursor_250k.fetchone()
        if not rawrow:
            break
        row=[str(i) for i in rawrow]
        sql_gb = "select tag from typelist where tag='"+"_".join(row[:3])+":"+source_name+"'"
        exist_or_not = int(cursor_gb.execute(sql_gb))
        if exist_or_not == 0: # Data not in database yet -)
            min_score, max_score = parse_and_load_gff3(row[3],"_".join(row[:3]),source_name)
            mk_conf("_".join(row[:3]), phenotype_method_dict[row[1]]+" Using "+analysis_method_dict[row[2]]+" (Call method "+row[0]+")", max_score, min_score)
    cursor_250k.close ()
    conn_250k.close ()

def parse_and_load_gff3(filename, feature_name, source_name):
    """
    Obtain a result tsv file, parse it into gff3 format then load it into database
    """
    percentage_threshold=10
    pfile=open(filename)
    
    lines = pfile.readlines()
    pfile.close()

    pvalues = []
    for line in lines[1:]:
    	pval = float(line.split("\t")[2])
	if pval > 0:
	    pvalues.append(-math.log10(pval))
	else:
	    pvalues.append(0)	
    pvalues.sort()
    threshold = pvalues[(100-percentage_threshold)*len(pvalues)/100]
    max = pvalues[len(pvalues)-1]

    pvoutfile = open(feature_name+".gff", "w")

    for line in lines[1:]:
        n=line.strip().split("\t")
	pval = float(n[2])
	if pval == 0:
	    continue
        elif -math.log10(pval)>threshold:
            pvoutfile.write("Chr"+n[0]+"\t"+source_name+"\t"+feature_name+"\t"+n[1]+"\t"+n[1]+"\t"+str(-math.log10(pval))+"\t.\t.\tName="+feature_name+";ID="+n[0]+"_"+n[1]+"\n")
    
    pvoutfile.close()
    
    # load the data... with a perl script
    os.system("perl load.pl "+feature_name+".gff")

    # and then destroy the file
    os.system("rm -f "+feature_name+".gff")
    return (threshold, max) # returns the minimum and maximum numbers

def mk_conf(feature_name, track_name, max, min):
    """
    Makes an entry in the arabidopsis gbrowse configuration file
    """
    f=open("/etc/apache2/gbrowse.conf/arabidopsis.conf","a")
    f.write("\n")
    f.write("["+feature_name+"]\n")
    f.write("feature         = "+feature_name)
    f.write(
"""
glyph           = xyplot
graph_type      = points
scale           = both
fgcolor         = black
bgcolor         = green
height          = 150
point_symbol    = filled_disc
point_radius    = 5
group_on        = display_name
category        = Association Mapping Results
""")
    f.write("min_score       = "+str(int(min/2))+"\n")
    f.write("max_score       = "+str(int(max+1))+"\n")
    f.write("key             = "+track_name+"\n")
    f.write("\n")
    f.write("["+feature_name+"_SNP]\n")
    f.write("feature         = "+feature_name)
    f.write(
"""
glyph           = diamond
fgcolor         = black
bgcolor         = green
category        = Association Mapping Results
label           = sub {
                   my $feature = shift;
                   my ($pos) = $feature->start();
                   return $pos;
                }
balloon hover   = sub {
                   my $f = shift;
                   my $score = $f->score();
                   my $position = $f->start();
                   return sprintf("Position: %d Score: %.6f", $position, $score);
                }
link            = sub {
                   my $f = shift;
                   my $score = $f->score();
                   my @met = split("_",$f->method());
                   my @snp = split("_",$f->load_id());
                   return sprintf("http://banyan.usc.edu:5000/SNP/?chromosome=%d&position=%d&call_method_id=%d&phenotype_method_id=%d&analysis_method_id=%d&score=%.6f",$snp[0],$snp[1],$met[0],$met[1],$met[2],$score);
		}
""")
    f.write("key             = "+track_name+" Clickable\n")
    f.close()
    
if __name__=="__main__":
    _run_()
