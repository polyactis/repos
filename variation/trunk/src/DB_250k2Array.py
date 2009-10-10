#!/usr/bin/env python
"""

Examples:
	#put SNP intensity matrices into designated file-system storage
	DB_250k2Array.py -o /Network/Data/250k/db/intensity/
	
	#output SNP intensity matrix in some temporary directory
	DB_250k2Array.py -o /tmp/arrays -a 616-710,811
	
	#output the intensity of CNV probes instead. a file called 'call_method_%s_CNV_intensity.tsv'%(call_method_id) would in the output dir.
	DB_250k2Array.py -o /tmp/CNV/ -l 17 -t 2 -u yh
	
	#calculate median_intensity for array 498,499 and store the value into db. '-o' is ju
	~/script/variation/src/DB_250k2Array.py -a 498,499 -t 3 -c
	
Description:
	2008-12-09 This program outputs the intensity of two types of probes on the array in run_type:
		1. SNP probes. each array (according to array_info_table) corresponds to one intensity matrix in output_dir.
		2. CNV probes. one giant ProbeXArray matrix. Columns are: probe_id, array1_id, array2_id, ..., chromosome, position
	2009-3-11 in run-type=3:
		3: calculate intensity medium of all probes in the array and store the value in db

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
from pymodule import process_function_arguments, getListOutOfStr
import numpy
import Stock_250kDB

class probe:
	def __init__(self, probes_id, snps_id, xpos, ypos, allele, strand):
		self.xpos = int(xpos)
		self.ypos = int(ypos)
		self.probes_id = probes_id
		self.snps_id = snps_id
		self.allele = allele
		self.strand = strand

class probes_class:
	def __init__(self):
		self.probes_id2probes_info = {}
	
	def add_one_probe(self, probes_id, snps_id, xpos, ypos, allele, strand):
		if probes_id in self.probes_id2probes_info:
			sys.stderr.write("Error: probes_id %s already in probes_id2probes_info.\n"%(probes_id))
			sys.exit(3)
		probe_ins = probe(probes_id, snps_id, xpos, ypos, allele, strand)
		self.probes_id2probes_info[probes_id] = probe_ins
	
	def get_one_probe(self, probes_id):
		return self.probes_id2probes_info.get(probes_id)

class snp:
	def __init__(self, snps_id, snpid):
		self.snps_id = snps_id
		self.snpid = snpid
		self.probes_id_ls = [-1]*4	#probes_id of [sense1, sense2, antisense1, antisense2]
		self.allele2index = {}

class snps_class:
	def __init__(self):
		self.snps_id2snps_info = {}
		self.snps_id_ls = []
	
	def add_one_snp(self, snps_id, snpid):
		if snps_id in self.snps_id2snps_info:
			sys.stderr.write("Error: snps_id %s already in snps_id2snps_info.\n"%(snps_id))
			sys.exit(3)
		snp_ins = snp(snps_id, snpid)
		self.snps_id2snps_info[snps_id] = snp_ins
		self.snps_id_ls.append(snps_id)
	
	def add_one_allele2snp(self, snps_id, allele):
		if snps_id not in self.snps_id2snps_info:
			sys.stderr.write("snps_id: %s not in snps_id2snps_info yet. no allele added."%snps_id)
			return
		allele2index = self.snps_id2snps_info[snps_id].allele2index
		self.snps_id2snps_info[snps_id].allele2index[allele] = len(allele2index)

	def add_one_probes_id2snp(self, snps_id, probes_id, allele, strand):
		allele_index = self.snps_id2snps_info[snps_id].allele2index[allele]
		if strand=='sense':
			self.snps_id2snps_info[snps_id].probes_id_ls[allele_index] = probes_id
		elif strand=='antisense':
			self.snps_id2snps_info[snps_id].probes_id_ls[allele_index+2] = probes_id
	
	def get_one_snp(self, snps_id):
		return self.snps_id2snps_info[snps_id]

class DB_250k2Array(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, '', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ):[None, 'p', 1, 'database password', ],\
							('output_dir', 0, ): [None, 'o', 1, 'directory to contain output files for run_type=1/2'],\
							('snps_table', 1, ): ['snps', 's', 1],\
							('probes_table', 1, ): ['probes', 'e'],\
							('array_info_table', 1, ):['array_info', 'y'],\
							('array_id_ls', 0, ): [None, 'a', 1, 'comma or dash-separated array id list, like 61-70,81. Not specifying this means all arrays.'],\
							('call_method_id', 0, int):[0, 'l', 1, 'Restrict arrays included in this call_method. Default is no such restriction.'],\
							('array_file_directory', 0, ):[None, 'f', 1, 'The results directory. Default is None. use the one given by db.'],\
							('run_type', 1, int):[1, 't', 1, '1: output SNP probe intensity, 2: output the intensity of CNV probes, 3: calculate intensity medium of all probes in the array and put into db'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	debug = 0
	def __init__(self, **keywords):
		"""
		2008-04-08
		"""
		from pymodule import ProcessOptions
		self.ad=ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		if self.array_id_ls:
			self.array_id_ls = getListOutOfStr(self.array_id_ls, data_type=str)
	
	def get_snps(self, curs, snps_table):
		"""
		2008-07-13
			add condition that alleles must be not null, to select SNPs in snps_table
		"""
		sys.stderr.write("Getting snps ... ")
		snps = snps_class()
		curs.execute("select id, name, allele1, allele2 from %s where allele1 is not null and allele2 is not null order by chromosome, position"%snps_table)
		rows = curs.fetchall()
		for row in rows:
			snps_id, snpid, allele1, allele2 = row
			snps.add_one_snp(snps_id, snpid)
			snps.add_one_allele2snp(snps_id, allele1)
			snps.add_one_allele2snp(snps_id, allele2)
		del rows
		sys.stderr.write("Done.\n")
		return snps
	
	def get_probes(cls, curs, probes_table, snps=None, run_type=1):
		"""
		2009-2-12
			become a classmethod
		2008-12-09
			add option run_type
		"""
		sys.stderr.write("Getting probes ... ")
		probes = probes_class()
		if run_type==2:
			curs.execute("select id, chromosome, position, xpos, ypos, allele, strand from %s where direction is not null order by chromosome, position"%(probes_table))
		else:
			curs.execute("select id, snps_id, xpos, ypos, allele, strand from %s where snps_id is not null"%(probes_table))
		rows = curs.fetchall()
		xy_ls = []
		chr_pos_ls = []
		probes_id_ls = []
		counter = 0
		for row in rows:
			if run_type==2:
				probes_id, chromosome, position, xpos, ypos, allele, strand = row
				xy_ls.append((xpos, ypos))
				chr_pos_ls.append((chromosome, position))
				probes_id_ls.append(probes_id)
			else:
				probes_id, snps_id, xpos, ypos, allele, strand = row
				probes.add_one_probe(probes_id, snps_id, xpos, ypos, allele, strand)
				snps.add_one_probes_id2snp(snps_id, probes_id, allele, strand)
			counter += 1
			if cls.debug and counter%1000==0:
				break
			
		del rows
		sys.stderr.write("Done.\n")
		return probes, xy_ls, chr_pos_ls, probes_id_ls
	
	get_probes = classmethod(get_probes)
	
	def outputArray(self, session, curs, output_dir, array_info_table, snps, probes, array_id_ls, xy_ls, chr_pos_ls, probes_id_ls,\
				call_method_id=0, run_type=1, array_file_directory=None):
		"""
		2009-10-9
			add argument array_file_directory.
		2009-3-11
			add run_type=3
				calculate intensity medium of all probes in the array and store the value in db
			array_id_ls is a list of array_ids in str type
		2009-3-5
			skip if no probes (if one_snp.probes_id_ls == [-1]*4:) for that SNP (fake SNP in the SNP table)
		2008-12-09
			add option run_type
		2008-07-12
			add option array_id
		2008-04-08
		"""
		sys.stderr.write("Outputting arrays ... \n")
		import rpy
		rpy.r.library('affy')
		array_size = None
		if run_type!=3 and not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		if array_id_ls:
			curs.execute("select id, filename, maternal_ecotype_id from %s where id in (%s)"%(array_info_table, ','.join(array_id_ls)))
		elif call_method_id:
			curs.execute("select v.array_id, a.filename, a.maternal_ecotype_id from view_call v, %s a where v.array_id=a.id and v.call_method_id=%s order by array_id"%\
						(array_info_table, call_method_id))
		else:
			curs.execute("select id, filename, maternal_ecotype_id from %s"%(array_info_table))
		rows = curs.fetchall()
		no_of_objects = len(rows)
		
		if run_type==2:	#2008-12-09 don't initialize the data_matrix if run_type is not 2 (CNV probe).
			data_matrix = numpy.zeros([len(probes_id_ls), no_of_objects], numpy.float)
		array_id_avail_ls = []
		ecotype_id_ls = []
		for i in range(no_of_objects):
			array_id, filename, ecotype_id = rows[i][:3]
			array_id_avail_ls.append(array_id)
			ecotype_id_ls.append(ecotype_id)
			
			if array_file_directory and os.path.isdir(array_file_directory):
				filename = os.path.join(array_file_directory, os.path.split(filename)[1])
			
			sys.stderr.write("\t%d/%d: Extracting intensity from %s ... \n"%(i+1, no_of_objects, filename))
			
			if run_type==1:	#output SNP probe intensity within the loop
				output_fname = os.path.join(output_dir, '%s_array_intensity.tsv'%(array_id))
				if os.path.isfile(output_fname):
					sys.stderr.write("\tFile %s already exists. Ignore.\n"%(output_fname))
					continue
			
			#read array by calling R
			array = rpy.r.read_affybatch(filenames=filename)
			intensity_array = rpy.r.intensity(array)	#return a lengthX1 2-Dimensional array.
			intensity_array_size = len(intensity_array)
			if array_size == None:
				array_size = int(math.sqrt(intensity_array_size))	#assume it's square array
			
			if run_type==2:	#CNV probe
				for j in range(len(xy_ls)):
					xpos, ypos = xy_ls[j]
					#chromosome, position = chr_pos_ls[j]
					intensity_array_index = array_size*(array_size - xpos - 1) + ypos
					#output_row = [chromosome, position]
					intensity = math.log10(intensity_array[intensity_array_index][0])
					#output_row.append(intensity)
					#writer.writerow(output_row)
					data_matrix[j][i] = intensity
			elif run_type==1:	#SNP probe intensity
				writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
				header = [ 'sense1', 'sense2', 'antisense1', 'antisense2']
				
				func = lambda x: '%s_%s'%(array_id, x)
				header = map(func, header)
				header = ['SNP_ID'] + header
				writer.writerow(header)
				for snps_id in snps.snps_id_ls:
					one_snp = snps.get_one_snp(snps_id)
					output_row = [one_snp.snpid]
					if one_snp.probes_id_ls == [-1]*4:	#2009-3-5 skip if no probes for that SNP (fake SNP in the SNP table)
						continue
					for probes_id in one_snp.probes_id_ls:
						one_probe = probes.get_one_probe(probes_id)
						intensity_array_index = array_size*(array_size - one_probe.xpos - 1) + one_probe.ypos
						output_row.append(intensity_array[intensity_array_index][0])
					writer.writerow(output_row)
				del writer
			elif run_type==3:	#calculate the intensity medium of all probes and store into db
				median_intensity = numpy.median(intensity_array)[0]
				array_info_entry = Stock_250kDB.ArrayInfo.get(array_id)
				array_info_entry.median_intensity = median_intensity
				session.save_or_update(array_info_entry)
			else:
				sys.stderr.write("Error: run_type %s is not supported.\n"%run_type)
				sys.exit(3)

			del intensity_array, array
		
		
		if run_type==2:
			#2008-11-13 output in Roger's multi-sample format
			header =['probes_id'] + array_id_avail_ls + ['chromosome', 'position']
			output_fname = os.path.join(output_dir, 'call_method_%s_CNV_intensity.tsv'%(call_method_id))
			
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			writer.writerow(header)
			for i in range(data_matrix.shape[0]):
				data_row = [probes_id_ls[i]] + list(data_matrix[i]) + list(chr_pos_ls[i])
				writer.writerow(data_row)
			del writer
		
		sys.stderr.write("Done.\n")
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
									password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		session.begin()
		
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.db_user, passwd = self.db_passwd)
		curs = conn.cursor()
		
		if self.run_type==1 or self.run_type==2:
			if not self.output_dir:
				sys.stderr.write("Run_Type 1 or 2 requires output_dir (-o).\n")
				sys.exit(2)
		
		if self.run_type==1:
			snps = self.get_snps(curs, self.snps_table)
		else:
			snps = None
		if self.run_type==1 or self.run_type==2:
			probes, xy_ls, chr_pos_ls, probes_id_ls = self.get_probes(curs, self.probes_table, snps, \
																			run_type=self.run_type)
		else:
			probes, xy_ls, chr_pos_ls, probes_id_ls = None, None, None, None
		
		self.outputArray(session, curs, self.output_dir, self.array_info_table, snps, probes, self.array_id_ls, xy_ls, \
						chr_pos_ls, probes_id_ls, call_method_id=self.call_method_id, run_type=self.run_type, \
						array_file_directory=self.array_file_directory)
		
		if self.commit:
			session.flush()
			session.commit()
			session.clear()
		else:
			session.rollback()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DB_250k2Array
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	
	"""
	if len(sys.argv) == 1:
		print DB_250k2Array.__doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "user=", "passwd=", "help", "type=", "debug", "report"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:u:p:i:o:y:a:e:f:g:j:cbr", long_options_list)
	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	hostname = None
	dbname = None
	user = None
	passwd = None
	input_fname = None
	output_fname = None
	type = None
	argument1 = None
	argument2 = None
	argument3 = None
	argument4 = None
	argument5 = None
	help = 0
	commit = 0
	debug = 0
	report = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-u", "--user"):
			user = arg
		elif opt in ("-p", "--passwd"):
			passwd = arg
		elif opt in ("-i",):
			input_fname = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-y", "--type"):
			type = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-a",):
			argument1 = arg
		elif opt in ("-e",):
			argument2 = arg
		elif opt in ("-f",):
			argument3 = arg
		elif opt in ("-g",):
			argument4 = arg
		elif opt in ("-j",):
			argument5 = arg
	
	ins = DB_250k2Array(hostname=hostname, dbname=dbname, user=user, passwd=passwd, output_dir=output_fname, snps_table=argument1, \
				probes_table=argument2, array_info_table=argument3,\
				debug=debug, report=report)
	ins.run()
	"""
