#!/usr/bin/env python
"""
Usage: FillStrainInfo.py [OPTIONS] -i INPUT_FILE

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	input file
	-s ...,	strain_info table, 'strain_info'(default)
	-t ...,	type of input, 1(Diane.csv GPS default), 2 (850 natural accessions.csv)
	-c, --commit	commit this database transaction
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	FillStrainInfo.py -i Diane.csv

Description:
	Program to expand the strain_info table with GPS and other strain info
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import getopt, csv
import Numeric
from sets import Set
from codense.common import db_connect


class FillStrainInfo:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema='dbsnp', input_fname=None, \
		strain_info_table='strain_info', input_type=1, commit=0, debug=0, report=0):
		"""
		2007-03-21
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.strain_info_table = strain_info_table
		self.input_type = int(input_type)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_strain_id_via_db_matching(self, curs, partial_acc, strain_info_table):
		"""
		2007-03-21
		2007-03-22
			1. RIV-xx causes ambiguous matching with MNF-Riv-xx
				so if len(rows)>1, try matching from beginning
			2. in strain_info_table, YNG47 should be YNG-47
				so modify the partial_acc a bit
			
		"""
		if partial_acc=='Yng-47':
			partial_acc = 'Yng47'
		curs.execute("select id from %s where acc~*'%s\\\D'"%(strain_info_table, partial_acc))
		rows = curs.fetchall()
		if len(rows)>1:
			curs.execute("select id from %s where acc~*'^%s\\\D'"%(strain_info_table, partial_acc))
			rows = curs.fetchall()
		return rows
	
	def ProcessDianeGPSInfo(self, curs, input_fname, strain_info_table, report=0):
		sys.stderr.write("Processing Diane GPS info ...\n")
		reader = csv.reader(open(input_fname), delimiter='\t')
		header = reader.next()
		old_pop_loc = None
		old_latitude = None
		old_longitude = None
		old_comments = None
		failed_barcode_ls = []
		counter = 0
		for row in reader:
			pop_location, pop, line, barcode, distance, latitude, longitude, comments = row
			if old_pop_loc == None:
				old_pop_loc = pop_location
				old_latitude = latitude
				old_longitude = longitude
				old_comments = comments
			elif old_pop_loc != pop_location:
				old_pop_loc = pop_location
				old_latitude = latitude
				old_longitude = longitude
				old_comments = comments
			elif comments:
				old_comments += comments
			if self.debug:
				import pdb
				pdb.set_trace()
			strain_id_ls = self.get_strain_id_via_db_matching(curs, barcode, strain_info_table)
			if len(strain_id_ls)==1:
				pop_location = pop_location.replace("'", ' ')
				if old_latitude and old_longitude:
					curs.execute("update %s set pop_location='%s', pop='%s', line=%s, barcode='%s',\
						distance='%s', latitude=%s, longitude=%s, comments='%s' where id=%s"%\
						(strain_info_table, pop_location, pop, line, barcode, distance, old_latitude, \
						old_longitude, old_comments, strain_id_ls[0][0]))
				else:	#the last block doesn't have old_latitude and old_longitude
					curs.execute("update %s set pop_location='%s', pop='%s', line=%s, barcode='%s',\
						distance='%s', comments='%s' where id=%s"%\
						(strain_info_table, pop_location, pop, line, barcode, distance, \
						old_comments, strain_id_ls[0][0]))
			else:
				failed_barcode_ls.append(barcode)
			counter += 1
			if report and counter%200==0:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
		if report:
			sys.stderr.write("%s%s\n"%('\x08'*20, counter))
		sys.stderr.write("Done.\n")
		return failed_barcode_ls
	
	def Process850NaturalAccessions(self, curs, input_fname, strain_info_table, report=0):
		sys.stderr.write("Processing 850 Natural Accessions ... \n")
		reader = csv.reader(open(input_fname), delimiter='\t')
		header = reader.next()
		failed_parental_stock_acc_ls = []
		counter = 0
		for row in reader:
			abrc_stock_acc, parental_stock_acc, abbr_name, full_name, country, days_to_flower = row
			if not days_to_flower:
				days_to_flower = -75
			else:
				days_to_flower = int(days_to_flower)
			if self.debug:
				import pdb
				pdb.set_trace()
			strain_id_ls = self.get_strain_id_via_db_matching(curs, parental_stock_acc, strain_info_table)
			if len(strain_id_ls)==1:
				full_name = full_name.replace("'", ' ')
				curs.execute("update %s set abrc_stock_acc='%s', parental_stock_acc='%s', abbr_name='%s', full_name='%s',\
					country='%s', days_to_flower=%s where id=%s"%\
					(strain_info_table, abrc_stock_acc, parental_stock_acc, abbr_name, full_name, country, \
					days_to_flower, strain_id_ls[0][0]))
			else:
				failed_parental_stock_acc_ls.append(parental_stock_acc)
			counter += 1
			if report and counter%200==0:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
		if report:
			sys.stderr.write("%s%s\n"%('\x08'*20, counter))
		sys.stderr.write("Done.\n")
		return failed_parental_stock_acc_ls
	
	def run(self):
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		if self.input_type==1:
			failed_ls = self.ProcessDianeGPSInfo(curs, self.input_fname, self.strain_info_table, self.report)
		elif self.input_type==2:
			failed_ls = self.Process850NaturalAccessions(curs, self.input_fname, self.strain_info_table, self.report)
		print "%s failures"%len(failed_ls)
		print failed_ls
		if self.commit:
			curs.execute("end")
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:s:t:brch", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'dbsnp'
	input_fname = None
	strain_info_table = 'strain_info'
	input_type = 1
	debug = 0
	report = 0
	commit = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-i",):
			input_fname = arg
		elif opt in ("-s",):
			strain_info_table = arg
		elif opt in ("-t",):
			input_type = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1

	if input_fname and hostname and dbname and schema:
		instance = FillStrainInfo(hostname, dbname, schema, input_fname, \
			strain_info_table, input_type, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
