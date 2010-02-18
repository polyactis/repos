#!/usr/bin/env python
"""
Examples:
	# a different conf file, call method 32
	~/script//variation/Dazhe/src/get_data_from_db.py -u yh -p STOCK_250k_PASSWD -q GBROWSE_DB_PASSWD -w ./arabidopsis1.conf -l 32
	
	# a different conf file, call method 32, a different pylons server, a different track db ID
	~/script//variation/Dazhe/src/get_data_from_db.py -u yh -p STOCK_250k_PASSWD -q GBROWSE_DB_PASSWD -y arabidopsis.usc.edu -f cypress_pgdb -w ./arabidopsis1.conf -l 32
	
Description:
	This program reads association results from files whose filenames are stored in stock_250k, turns them into gff3 format,
		and call 'load.pl' to load them into the gbrowse db.
		It also modifies the gbrowse_conf_file to which more track confs are appended.
	
	The gbrowse db setting (gb_host, gb_dbname, gb_db_user, gb_db_passwd) in this program is passed to load.pl.
	
	So this program could be run from a machine where the gbrowse db doesn't reside.
	But afterwards the gbrowse_conf_file has to appended to the appropriate gbrowse configuration file on the server where the GBrowse runs. 
	
	
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, numpy, getopt
import traceback, gc
from pymodule import process_function_arguments, getListOutOfStr, PassingData
import numpy
from variation.src import Stock_250kDB

class LoadStock250kAssociationIntoGBrowseDB(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database stock_250k is? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the stock_250k db', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, '', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username for stock_250k', ],\
							('db_passwd', 1, ):[None, 'p', 1, 'database password for stock_250k', ],\
							
							('gb_host', 1, ): ['cypress.usc.edu', 's', 1, 'hostname of the GBrowse db', ],\
							('gb_dbname', 1, ): ['gbrowse', 'n', 1, 'database name for the GBrowse db', ],\
							('gb_db_user', 1, ): ['gbrowse_loader', 'g', 1, 'database username for the GBrowse db', ],\
							('gb_db_passwd', 1, ):[None, 'q', 1, 'database password for the GBrowse db', ],\
							
							('pylons_server', 1, ):['banyan.usc.edu:5000', 'y', 1, 'pylons server to which clicking the GBrowse SNP links', ],\
							('gbrowse_conf_fname', 1, ):['/etc/gbrowse2/arabidopsis1.conf', 'w', 1, ''],\
							('default_track_db', 1, ):['cypress_db', 'f', 1, 'database ID in gbrowse_conf_fname, in which the track data is stored.'], \
							
							('tmp_dir', 0, ): ['/tmp/', 'o', 1, 'The temporary directory hold contain track gff files before them being loaded into db'],\
							('call_method_id', 0, int):[32, 'l', 1, 'Restrict arrays included in this call_method. Default is no such restriction.'],\
							('phenotype_method_id_ls', 0, ): [None, 'e', 1, 'comma or dash-separated list of phenotype method ids, like 61-70,81. Not specifying this means no restriction.'],\
							('analysis_method_id_ls', 0, ):[None, 'a', 1, 'comma or dash-separated analysis method id lists'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-04-08
		"""
		from pymodule import ProcessOptions
		self.ad=ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		if self.analysis_method_id_ls:
			self.analysis_method_id_ls = getListOutOfStr(self.analysis_method_id_ls, data_type=int)
		if self.phenotype_method_id_ls:
			self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
	
	def parse_and_load_gff3(self, filename, feature_name, output_dir, source_name, gb_host, gb_dbname, gb_db_user, gb_db_passwd):
		"""
		2010-2-4
			- a new file whose name is '%s.gff'%feature_name would be created but then removed after it's loaded into db
			- will call load.pl
			
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
				pvalues.append(30)
		pvalues.sort()
		threshold = pvalues[(100-percentage_threshold)*len(pvalues)/100]	#only the top portion 
		max = pvalues[-1]
		
		feature_fname = os.path.join(output_dir, '%s.gff'%feature_name)
		pvoutfile = open(feature_fname, "w")
	
		for line in lines[1:]:
			n=line.strip().split("\t")
			pval = float(n[2])
			if pval == 0:
				continue
			elif -math.log10(pval)>threshold:
				pvoutfile.write("Chr"+n[0]+"\t"+source_name+"\t"+feature_name+"\t"+n[1]+"\t"+n[1]+"\t"+str(-math.log10(pval))+"\t.\t.\tName="+feature_name+";ID="+n[0]+"_"+n[1]+"\n")
			
		pvoutfile.close()
		
		# load the data... with a perl script
		program_path = os.path.dirname(sys.argv[0])
		load_pl_path = os.path.join(program_path, 'load.pl')
		os.system("perl %s %s %s %s %s %s"%(load_pl_path, gb_host, gb_dbname, gb_db_user, gb_db_passwd, feature_fname))
	
		# and then destroy the file
		os.system("rm -f "+feature_fname)
		return (threshold, max) # returns the minimum and maximum numbers
	
	def mk_conf(self, feature_name, track_name, max, min, pylons_server = 'arabidopsis.usc.edu', \
			gbrowse_conf_fname = '/etc/apache2/gbrowse.conf/arabidopsis.conf',\
			default_track_db = ''):
		"""
		2010-2-4
			add argument pylons_server, gbrowse_conf_fname, default_track_db to make it 
		Makes an entry in the arabidopsis gbrowse configuration file
		"""
		f=open(gbrowse_conf_fname,"a")
		f.write("\n")
		f.write("["+feature_name+"]\n")
		f.write("feature		 = %s\n"%feature_name)
		if default_track_db:
			f.write("database	= %s\n"%(default_track_db))
		f.write(
"""glyph		= xyplot
graph_type		= points
scale			= both
fgcolor			= black
bgcolor			= green
height			= 150
point_symbol	= filled_disc
point_radius	= 5
group_on		= display_name
category		= Association Mapping Results
""")
		f.write("min_score	   = "+str(int(min/2))+"\n")
		f.write("max_score	   = "+str(int(max+1))+"\n")
		f.write("key			 = "+track_name+"\n")
		f.write("\n")
		f.write("["+feature_name+"_SNP]\n")
		f.write("feature		 = %s\n"%feature_name)
		if default_track_db:
			f.write("database	= %s\n"%(default_track_db))
		f.write(
"""glyph		= diamond
fgcolor			= black
bgcolor			= green
category		= Association Mapping Results
label			= sub {
					my $feature = shift;
					my ($pos) = $feature->start();
					return $pos;
				}
balloon hover	= sub {
					my $f = shift;
					my $score = $f->score();
					my $position = $f->start();
					return sprintf("Position: %d Score: %.6f", $position, $score);
				}
link			= sub {
					my $f = shift;
					my $score = $f->score();
					my @met = split("_",$f->method());
					my @snp = split("_",$f->load_id());
					return sprintf("http://"""+ pylons_server + """/SNP/?chromosome=%d&position=%d&call_method_id=%d&phenotype_method_id=%d&analysis_method_id=%d&score=%.6f",$snp[0],$snp[1],$met[0],$met[1],$met[2],$score);
				}
""")
		f.write("key			 = "+track_name+" Clickable\n")
		f.close()
	
	def run(self):
		"""
		2010-2-4
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		source_name='Nordborg'
		
		phenotypes_in_mahogany = range(1,49)+range(57,83)+[158,159]+range(161,180)+range(182,187)+[272,273,274]+range(277,284)   
		# 2010-2-4 temporary solution, not to reload the features already in mahogany db (gbrowse2 allows multiple slaves) 
		feature_set_in_mahogany = set()
		for call_method_id in ["32"]:
			for phenotype_method_id in phenotypes_in_mahogany:
				for analysis_method_id in ["1", "7"]:
					feature_name = "%s_%s_%s"%(call_method_id, phenotype_method_id, analysis_method_id)
					feature_set_in_mahogany.add(feature_name)
		
		
		import MySQLdb
		print "Connecting to stock_250k, host="+self.hostname
		try:
			conn_250k = MySQLdb.connect (host=self.hostname, user=self.db_user, passwd=self.db_passwd, db=self.dbname)
		except MySQLdb.Error, e:
			print "Error %d: %s" % (e.args[0], e.args[1])
			sys.stderr.write("Can't connect to db %s on host %s\n"%(self.dbname, self.hostname))
			sys.exit (1)
		cursor_250k = conn_250k.cursor()
		
		# Create the Phenotype and Association Method naming dictionaries		
		phenotype_method_dict={}
		sqlcommand="select id, short_name from phenotype_method"
		cursor_250k.execute(sqlcommand)
		while(1):
			row = cursor_250k.fetchone()
			if not row:
				break
			phenotype_method_dict[row[0]]=row[1]
		
		# Retrieve the filenames
		print "Fetching data ..."
		sqlcommand="select r.call_method_id, r.phenotype_method_id, r.analysis_method_id, r.filename, a.short_name from results_method r,\
			analysis_method a where a.id = r.analysis_method_id "
		if self.call_method_id>0:
			sqlcommand += " AND r.call_method_id=%s "%self.call_method_id
		numRows = int(cursor_250k.execute(sqlcommand))
		
		# Store them in a file first? or use them straight away?
		# Try straight away
		
		# Connect to the gbrowse database
		try:
			conn_gb = MySQLdb.connect(host=self.gb_host, user=self.gb_db_user, passwd=self.gb_db_passwd, db=self.gb_dbname)
		except MySQLdb.Error, e:
			print "Error %d: %s" % (e.args[0], e.args[1])
			sys.stderr.write("Can't connect to db %s on host %s\n"%(self.gb_dbname, self.gb_host))
			sys.exit (1)
		cursor_gb = conn_gb.cursor()	
			
		while(1):
			rawrow = cursor_250k.fetchone()
			if not rawrow:
				break
			row=[str(i) for i in rawrow]
			
			phenotype_method_id = int(row[1])
			analysis_method_id = int(row[2])
			if self.phenotype_method_id_ls and phenotype_method_id not in self.phenotype_method_id_ls: 
				continue
			if self.analysis_method_id_ls and analysis_method_id not in self.analysis_method_id_ls:
				continue
			feature_name = "_".join(row[:3])
			
			# 2010-2-4 temporary solution to avoid overlap with mahogany gbrowse db
			if feature_name in feature_set_in_mahogany:
				sys.stderr.write("%s already in mahogany. ignored.\n"%feature_name)
				continue
			
			sql_gb = "select tag from typelist where tag='%s:%s'"%(feature_name, source_name)
			exist_or_not = int(cursor_gb.execute(sql_gb))
			if exist_or_not == 0: # Data not in database yet -)
				min_score, max_score = self.parse_and_load_gff3(row[3], feature_name, self.tmp_dir, source_name, self.gb_host, \
															self.gb_dbname, self.gb_db_user, self.gb_db_passwd)
				analysis_method_name = row[4]
				self.mk_conf(feature_name, phenotype_method_dict[phenotype_method_id]+" Using "+analysis_method_name, max_score, \
					min_score, pylons_server=self.pylons_server, gbrowse_conf_fname=self.gbrowse_conf_fname, \
					default_track_db=self.default_track_db)
		cursor_250k.close()
		conn_250k.close()
	
if __name__=="__main__":
	from pymodule import ProcessOptions
	main_class = LoadStock250kAssociationIntoGBrowseDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()