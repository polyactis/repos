#!/usr/bin/env python
"""

Examples:
	data2dbSNP.py -a "384-SNP Illumina multiplex" -i /tmp/2ndPlate-data.csv -l /tmp/080501_accession2ecotype_id_384SNP.csv -s "384-SNP-Illumina-multiplex" -m "93 accessions from the so-called 2nd batch (or 2nd 96)"  -f /tmp/chosenR2.csv  -c

Description:
	put 384-SNP Illumina-multiplex data into db.
	
	all the data files are attached in ticket 17 in Trac.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import sys, getopt, csv, re
import psycopg2 as psycopg
from annot.bin.codense.common import db_connect, org_short2long, org2tax_id

class data2dbSNP_2007_03_06:
	"""
	2008-05-02 renamed to data2dbSNP_2007_03_06. a new data2dbSNP class is generated.
	
Usage: data2dbSNP.py [OPTIONS] -i INPUT_FILE -o OUTPUT_TABLE

Option:
	-z ..., --hostname=...	the hostname, dl324b-1(default)
	-d ..., --dbname=...	the database name, yhdb(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	INPUT_FILE
	-o ...,	SNP output table
	-s ...,	strain_info table, 'strain_info'(default)
	-n ...,	snp_locus_table, 'snp_locus'(default)
	-g ...	the two-letter abbreviation of organism, at(default)
		only for type =1
	-y ...	parser type(1, default)
	-b, --debug	just test running the program, no daytime restriction
	-r, --report	report the progress (time before each query)
	-c, --commit	commit this database transaction
	-h, --help	show this help

Examples:
	data2dbSNP.py -i UCHIC_HO_011107.csv -o justin_data -r -c

Description:
	put SNP data into database schema
	
	2007-02-17
	
	"""
	def __init__(self, hostname='dl324b-1', dbname='yhdb', schema='dbsnp', input_fname=None, \
		output_table=None, strain_info_table='strain_info', snp_locus_table='snp_locus', \
		organism='hs', type=1, debug=0, report=0, commit=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.output_table = output_table
		self.strain_info_table = strain_info_table
		self.snp_locus_table = snp_locus_table
		self.tax_id = org2tax_id(org_short2long(organism))
		self.type = int(type)
		self.debug = int(debug)
		self.report = int(report)
		self.commit = int(commit)
		
		self.snp_acc_category_pattern = re.compile("([a-zA-Z]*[\-]*[a-zA-Z]+)[\-_ ]*[\w]+")
		#- (followed by number, 2nd - followed by alphanumeric), _ or space is separator
	
	def get_strain_acc2id(self, curs, strain_info_table, tax_id):
		sys.stderr.write("Getting mapping data from %s ..."%strain_info_table)
		strain_acc2id = {}
		curs.execute("DECLARE crs CURSOR FOR select acc, id from %s where tax_id=%s"%(strain_info_table, tax_id))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				acc, id = row
				strain_acc2id[acc] = id
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write("Done.\n")
		return strain_acc2id
	
	def get_snp_acc2id(self, curs, snp_locus_table, tax_id):
		snp_acc2id = self.get_strain_acc2id(curs, snp_locus_table, tax_id)
		return snp_acc2id
	
	def submit_to_strain_info_table(self, curs, strain_acc, category, new_id, tax_id, strain_info_table):
		curs.execute("insert into %s(acc, category, id, tax_id) values('%s', '%s', %s, %s)"%\
			(strain_info_table, strain_acc, category, new_id, tax_id))
	
	def expand_snp_locus_table(self, curs, snp_locus_table, snp_acc_list, snp_acc2id, tax_id):
		sys.stderr.write("Expanding snp_locus_table ...")
		for snp_acc in snp_acc_list:
			if snp_acc not in snp_acc2id:
				new_id = len(snp_acc2id)+1
				curs.execute("insert into %s(acc, id, tax_id) values('%s', %s, %s)"%(snp_locus_table, snp_acc, new_id, tax_id))
				snp_acc2id[snp_acc] = new_id
		sys.stderr.write("Done.\n")
		return snp_acc2id
	
	def submit_to_snp_table(self, curs, strain_id, snp_index2acc, snp_acc2id, call_list, output_table):
		snp_counter = 0
		for i in range(len(call_list)):
			call = call_list[i]
			if call!='NA' and call!='0':	#0 or NA means Not Available
				snp_id = snp_acc2id[snp_index2acc[i]]
				curs.execute("insert into %s(strain_id, snp_id, call) values(%s, %s, '%s')"%(output_table, strain_id, snp_id, call))
				snp_counter += 1
		return snp_counter
	
	def parse_file(self, input_fname, curs, output_table, tax_id, snp_acc_category_pattern, strain_acc2id, snp_acc2id, snp_locus_table):
		"""
		2007-03-06
			if category extraction using 're' fails, just use the 1st 3 character of strain_acc
		"""
		sys.stderr.write("Starting to parse %s ...\n"%input_fname)
		reader = csv.reader(open(input_fname), delimiter='\t')
		snp_acc_list = reader.next()[1:]
		snp_index2acc = dict(zip(range(len(snp_acc_list)), snp_acc_list))
		snp_acc2id = self.expand_snp_locus_table(curs, snp_locus_table, snp_acc_list, snp_acc2id, tax_id)
		counter = 0
		snp_counter = 0
		for row in reader:
			strain_acc = row[0]
			call_list = row[1:]
			if strain_acc not in strain_acc2id:
				new_id = len(strain_acc2id)+1
				try:
					category = snp_acc_category_pattern.match(strain_acc).groups()[0]
				except:
					print 'category extraction error'
					print row
					print counter
					print strain_acc
					category = strain_acc[:3]
				self.submit_to_strain_info_table(curs, strain_acc, category, new_id, tax_id, strain_info_table)
				strain_acc2id[strain_acc] = new_id
			snp_counter += self.submit_to_snp_table(curs, new_id, snp_index2acc, snp_acc2id, call_list, output_table)
			counter += 1
			if self.report and counter%100==0:
				sys.stderr.write('%s%s\t%s'%('\x08'*20, counter, snp_counter))
		if self.report:
			sys.stderr.write('%s%s\t%s'%('\x08'*20, counter, snp_counter))
		del reader
		sys.stderr.write("Done.\n")
		
	def run(self):
		"""
		2007-02-17
		
		-db_connect()
		-get_strain_acc2id()
		-get_snp_acc2id()
		-parse_file()
			-expand_snp_locus_table()
			-submit_to_strain_info_table()
			-submit_to_snp_table()
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		strain_acc2id = self.get_strain_acc2id(curs, self.strain_info_table, self.tax_id)
		snp_acc2id = self.get_snp_acc2id(curs, self.snp_locus_table, self.tax_id)
		self.parse_file(self.input_fname, curs, self.output_table, self.tax_id, self.snp_acc_category_pattern, strain_acc2id, snp_acc2id, self.snp_locus_table)
		if self.commit:
			curs.execute("end")

from dbsnp import DBSNP, SNPs, SNPs2SNPset, SNPset, Calls, CallMethod, Accession

class data2dbSNP(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ):['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ):['dbsnp', 'd', 1, '',],\
							('user', 1, ):[None, 'u', 1, 'database username',],\
							('passwd',1, ):[None, 'p', 1, 'database password', ],\
							('input_fname',1, ): [None, 'i', 1, 'File containing 384-SNP illumina data. csv format.'],\
							('name2ecotype_id_fname',1, ): [None, 'l', 1, 'File containing the linking between accession name in input_fname and ecotype_id. csv format.'],\
							('snp_probe_fname',1, ): [None, 'f', 1, 'File containing snp name, probe sequence, chromosome, position. chosenR.csv sent from chris.'],\
							('snpset_acc',1, ): [None, 's'],\
							('snpset_description', 0, ): [None, 'n' ],\
							('call_method_short_name', 1, ): [None, 'a'],\
							('call_method_description', 0, ): [None,  'e'],\
							('call_method_data_description', 0, ): [None, 'm'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def get_name_duplicate2accession(self, name2ecotype_id_fname):
		"""
		2008-05-04
			consider name duplicates
		2008-05-02
		"""
		sys.stderr.write("Getting name_duplicate2accession ... ")
		reader = csv.reader(open(name2ecotype_id_fname))
		name_duplicate2accession = {}
		name2duplicate = {}	#to keep track of duplicates
		for row in reader:
			name, ecotype_id = row
			if name not in name2duplicate:
				name2duplicate[name] = 0
			name2duplicate[name] += 1
			name_duplicate = (name, name2duplicate[name])
			name_duplicate2accession[name_duplicate] = Accession(name=name, ecotype_id=ecotype_id, duplicate=name2duplicate[name])
		sys.stderr.write("Done.\n")
		return name_duplicate2accession
	
	def get_chromosome_pos2snp_obj(self, snp_probe_fname, session, snpset):
		"""
		2008-05-05
			snp_probe_fname is a csv file, a square matrix extracted from chosenR2.csv from chris.
			4-column: name, probe sequence, chromosome, position
		"""
		sys.stderr.write("Getting chromosome_pos2snp_obj ...")
		reader = csv.reader(open(snp_probe_fname))
		chromosome_pos2snp_obj = {}
		import re
		offset_pattern = re.compile(r'o(\d+)')
		allele_pattern = re.compile(r'\[([a-zA-Z])/([a-zA-Z])\]')
		for row in reader:
			snp_acc, probe_sequence, random_stuff, chromosome, position = row[:5]
			offset = None
			search_result = offset_pattern.search(snp_acc)
			if search_result:
				offset = int(search_result.groups()[0])
			search_result1 = allele_pattern.search(probe_sequence)
			allele1, allele2 = None, None
			if search_result1:
				allele1, allele2 = search_result1.groups()
			chromosome=int(chromosome)
			position=int(position)
			s = SNPs(name=snp_acc, chromosome=chromosome, position=position, offset=offset, \
					probe_sequence=probe_sequence, allele1=allele1, allele2=allele2)
			snpset.snps.append(s)
			s.snpset.append(snpset)
			chromosome_pos2snp_obj[(chromosome, position)] = s
			session.save(s)
		sys.stderr.write("Done.\n")
		return chromosome_pos2snp_obj
	
	def readin_calls(self, input_fname, name_duplicate2accession, session, callmethod, chromosome_pos2snp_obj):
		"""
		2008-05-04
			consider name duplicates
		2008-05-02
		"""
		sys.stderr.write("Reading in calls ...")
		reader = csv.reader(open(input_fname))
		#handle the SNPs first
		chr_ls = reader.next()[1:]
		position_ls = reader.next()[1:]
		chr_pos_ls = []
		no_of_snps = len(chr_ls)
		for i in range(no_of_snps):
			chromosome = int(chr_ls[i])
			position = int(position_ls[i])
			chr_pos_ls.append((chromosome, position))
		
		reader.next()	#skip this weird line
		name2duplicate = {}	#to keep track of duplicates
		for row in reader:
			name = row[0]
			if name not in name2duplicate:
				name2duplicate[name] = 0
			name2duplicate[name] += 1
			name_duplicate = (name, name2duplicate[name])
			accession = name_duplicate2accession[name_duplicate]
			for i in range(no_of_snps):
				chr_pos = chr_pos_ls[i]
				if chr_pos in chromosome_pos2snp_obj:
					s = chromosome_pos2snp_obj[chr_pos]
				else:
					sys.stderr.write("Error: %s not in chromosome_pos2snp_obj. Ignore. Won't commit db transaction.\n"%(repr(chr_pos)))
					self.commit = 0
					continue
				calls_obj = Calls(genotype=row[i+1])
				calls_obj.snps = s
				calls_obj.call_method = callmethod
				calls_obj.accession = accession
				session.save(calls_obj)
		sys.stderr.write("Done.\n")
	
	def main(self):
		"""
		2008-08-11
			the database interface changed in variation.src.dbsnp
		"""
		db = DBSNP(username=self.user,
				   password=self.passwd, host=self.hostname, database=self.dbname)
		session = db.session
		session.begin()
		#transaction = session.create_transaction()
		if self.debug:
			import pdb
			pdb.set_trace()
		snpset = SNPset(name=self.snpset_acc, description=self.snpset_description)
		callmethod = CallMethod(short_name=self.call_method_short_name, method_description=self.call_method_description, data_description=self.call_method_data_description)
		
		chromosome_pos2snp_obj = self.get_chromosome_pos2snp_obj(self.snp_probe_fname, session, snpset)
		
		name_duplicate2accession = self.get_name_duplicate2accession(self.name2ecotype_id_fname)
		
		self.readin_calls(self.input_fname, name_duplicate2accession, session, callmethod, chromosome_pos2snp_obj)
		session.flush()
		if self.commit:
			session.commit()
		else:	#default is also rollback(). to demonstrate good programming
			session.rollback()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = data2dbSNP
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.main()
	"""
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:o:s:n:g:y:brch", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'dl324b-1'
	dbname = 'yhdb'
	schema = 'dbsnp'
	input_fname = None
	output_table = None
	strain_info_table = 'strain_info'
	snp_locus_table = 'snp_locus'
	organism = 'at'
	type = 1
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
		elif opt in ("-o",):
			output_table = arg
		elif opt in ("-s",):
			strain_info_table = arg
		elif opt in ("-n",):
			snp_locus_table = arg
		elif opt in ("-g",):
			organism = arg
		elif opt in ("-y",):
			type = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1

	if input_fname and output_table and hostname and dbname and schema:
		instance = data2dbSNP_2007_03_06(hostname, dbname, schema, input_fname, output_table, \
			strain_info_table, snp_locus_table, organism, type, debug, report, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
	"""
