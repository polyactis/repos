#!/usr/bin/env python
"""

Examples:
	#put intensity matrices into designated file-system storage
	DB_250k2Array.py -o /Network/Data/250k/db/intensity/
	
	#put them in some temporary directory
	DB_250k2Array.py -o /tmp/arrays
Description:
	output all .cel array files (according to array_info_table) as intensity matrices in output_dir

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
from pymodule import process_function_arguments

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
		return self.probes_id2probes_info[probes_id]

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
	option_default_dict = {('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, '', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ):[None, 'p', 1, 'database password', ],\
						('output_dir', 1, ): [None, 'o', 1, ],\
						('snps_table', 1, ): ['snps', 's', 1],\
						('probes_table', 1, ): ['probes', 'e'],\
						('array_info_table', 1, ):['array_info', 'y'],\
						('array_id_ls', 0, ): [None, 'a', 1, 'comma-separated array id list, used to choose some arrays. Not specifying this means all arrays.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		"""
		2008-04-08
		"""
		from pymodule import ProcessOptions
		self.ad=ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
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
	
	def get_probes(self, curs, probes_table, snps):
		"""
		"""
		sys.stderr.write("Getting probes ... ")
		probes = probes_class()
		curs.execute("select id, snps_id, xpos, ypos, allele, strand from %s where snps_id is not null"%(probes_table))
		rows = curs.fetchall()
		for row in rows:
			probes_id, snps_id, xpos, ypos, allele, strand = row
			probes.add_one_probe(probes_id, snps_id, xpos, ypos, allele, strand)
			snps.add_one_probes_id2snp(snps_id, probes_id, allele, strand)
		del rows
		sys.stderr.write("Done.\n")
		return probes
	
	def outputArray(self, curs, output_dir, array_info_table, snps, probes, array_id_ls):
		"""
		2008-07-12
			add option array_id
		2008-04-08
		"""
		sys.stderr.write("Outputting arrays ... \n")
		import rpy
		rpy.r.library('affy')
		array_size = None
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		if array_id_ls:
			curs.execute("select id, filename from %s where id in (%s)"%(array_info_table, array_id_ls))
		else:
			curs.execute("select id, filename from %s"%(array_info_table))
		rows = curs.fetchall()
		no_of_objects = len(rows)
		for i in range(no_of_objects):
			array_id, filename = rows[i]
			sys.stderr.write("\t%d/%d: Extracting intensity from %s ... \n"%(i+1, no_of_objects, filename))
			
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
			
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			header = [ 'sense1', 'sense2', 'antisense1', 'antisense2']
			func = lambda x: '%s_%s'%(array_id, x)
			header = map(func, header)
			header = ['SNP_ID'] + header
			writer.writerow(header)
			for snps_id in snps.snps_id_ls:
				one_snp = snps.get_one_snp(snps_id)
				output_row = [one_snp.snpid]
				for probes_id in one_snp.probes_id_ls:
					one_probe = probes.get_one_probe(probes_id)
					intensity_array_index = array_size*(array_size - one_probe.xpos - 1) + one_probe.ypos
					output_row.append(intensity_array[intensity_array_index][0])
				writer.writerow(output_row)
			del writer, intensity_array, array
		sys.stderr.write("Done.\n")
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.db_user, passwd = self.db_passwd)
		curs = conn.cursor()
		snps = self.get_snps(curs, self.snps_table)
		probes = self.get_probes(curs, self.probes_table, snps)
		self.outputArray(curs, self.output_dir, self.array_info_table, snps, probes, self.array_id_ls)

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
