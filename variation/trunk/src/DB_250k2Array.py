#!/usr/bin/env python
"""
Usage: DB_250k2Array.py [OPTIONS] -o OUTPUT_DIR

Argument list:
	-z ..., --hostname=...	the hostname, papaya.usc.edu(default)
	-d ..., --dbname=...	the database name, stock_250k(default)
	-u ..., --user=...	the db username, (otherwise it will ask for it).
	-p ..., --passwd=...	the db password, (otherwise it will ask for it).
	-o ...,	output_dir*
	-a ...,	snps_table, 'stock_250k.snps'(default)
	-e ...,	probes_table, 'stock_250k.probes'(default)
	-g ...,	array_info_table, 'stock_250k.array_info'(default)
	-b,	toggle debug
	-r, toggle report
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
	def __init__(self, **keywords):
		"""
		2008-04-08
		"""
		argument_default_dict = {('hostname',1, ):'papaya.usc.edu',\
								('dbname',1, ):'stock_250k',\
								('user',1, ):None,\
								('passwd',1, ):None,\
								('output_dir',1, ):None,\
								('snps_table',1, ):'stock_250k.snps',\
								('probes_table',1, ):'stock_250k.probes',\
								('array_info_table',1, ):'stock_250k.array_info',\
								('commit',0, int):0,\
								('debug',0, int):0,\
								('report',0, int):0}
		"""
		2008-02-28
			argument_default_dict is a dictionary of default arguments, the key is a tuple, ('argument_name', is_argument_required, argument_type)
			argument_type is optional
		"""
		#argument dictionary
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def get_snps(self, curs, snps_table):
		"""
		"""
		sys.stderr.write("Getting snps ... ")
		snps = snps_class()
		curs.execute("select id, snpid, allele1, allele2 from %s order by chromosome, position"%snps_table)
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
	
	def outputArray(self, curs, output_dir, array_info_table, snps, probes):
		"""
		2008-04-08
		"""
		sys.stderr.write("Outputting arrays ... \n")
		import rpy
		rpy.r.library('affy')
		array_size = None
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
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
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.user, passwd = self.passwd)
		curs = conn.cursor()
		
		snps = self.get_snps(curs, self.snps_table)
		probes = self.get_probes(curs, self.probes_table, snps)
		self.outputArray(curs, self.output_dir, self.array_info_table, snps, probes)

if __name__ == '__main__':
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
		
