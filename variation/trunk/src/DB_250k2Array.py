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
	
	# 2009-10-09 specify a different directory for all the cel files.
	~/script/variation/src/DB_250k2Array.py -l 48 -t 2 -u yh  -o panfs/250k/CNV/ -f panfs/db/raw_data/
	
	#2009-11-24 calculate DLRSpread and 7-number summary & outliers, etc. into db
	~/script/variation/src/DB_250k2Array.py -t 4 -u yh -c
	
Description:
	2008-12-09 This program outputs the intensity of two types of probes on the array in run_type:
		1. SNP probes. each array (according to array_info_table) corresponds to one intensity matrix in output_dir.
		2. CNV probes. one giant ProbeXArray matrix. Columns are: probe_id, array1_id, array2_id, ..., chromosome, position
	2009-3-11 add run-type=3:
		3: calculate intensity medium of all probes in the array and store the value in db
	2009-11-24 add run-type=4:
		4: calculate DLRSpread, 7-number summary for both raw intensity from each chromosome and corresponding dLR.
			store outlier probes (defined below) in db as well.
			
			DLRSpread = IQR(dLR) / 4*erfinv(0.5), where dLR is an array of
		differences between log ratios of adjacent probes along the chromosome, erfinv is the Inverse Error
		Function and IQR is Inter Quartile Range. store the value in db.
			
			7-number summary = (minimum, first_decile, lower_quartile, median, upper_quartile, last_decile, maximum).
			
			outliers (<lower_quartile-1.5*IQR or >upper_quartile+1.5*IQR) are recorded in ArrayQuartileOutlier.
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
import Stock_250kDB
from scipy import stats	# for scoreatpercentile/percentileatscore to get quartiles
import scipy.special
from common import calculate7NumberSummaryForOneList

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
							('run_type', 1, int):[1, 't', 1, '1: output SNP probe intensity, 2: output the intensity of CNV probes, 3: calculate median intensity of all probes, 4: calculate DLRSpread & 7-number summary'],\
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
	
	@classmethod
	def get_probes(cls, curs, probes_table, snps=None, run_type=1):
		"""
		2009-11-21
			curs could either be elixirdb.metadata.bind or MySQLdb.connect().cursor
		2009-2-12
			become a classmethod
		2008-12-09
			add option run_type
				1: SNP probes
				2: CNV probes 	(argument snps is not used)
		"""
		sys.stderr.write("Getting probes ... ")
		probes = probes_class()
		if run_type==2 or run_type==4:
			rows = curs.execute("select id, chromosome, position, xpos, ypos, allele, strand from %s where direction is not null order by chromosome, position"%(probes_table))
		else:
			rows = curs.execute("select id, snps_id, xpos, ypos, allele, strand from %s where snps_id is not null"%(probes_table))
		
		is_elixirdb = 1	# 2009-11-20 By default, assume curs is elixirdb.metadata.bind
		if hasattr(curs, 'fetchall'):	# 2009-11-20 curs is MySQLdb.connect
			rows = curs.fetchall()
			is_elixirdb = 0
		xy_ls = []
		chr_pos_ls = []
		probes_id_ls = []
		counter = 0
		for row in rows:
			if run_type==2 or run_type==4:
				if is_elixirdb:
					probes_id = row.id
					chromosome = row.chromosome
					position = row.position
					xpos = row.xpos
					ypos = row.ypos
					allele = row.allele
					strand = row.strand
				else:
					probes_id, chromosome, position, xpos, ypos, allele, strand = row
				xy_ls.append((xpos, ypos))
				chr_pos_ls.append((chromosome, position))
				probes_id_ls.append(probes_id)
			else:
				if is_elixirdb:
					probes_id = row.id
					snps_id = row.snps_id
					xpos = row.xpos
					ypos = row.ypos
					allele = row.allele
					strand = row.strand
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
	
	@classmethod
	def getArrayWidth(cls, array_filename):
		"""
		2009-11-20
			get the width (number of probes on X-axis) of an array assuming it's a square array.
		"""
		import rpy
		rpy.r.library('affy')
		array = rpy.r.read_affybatch(filenames=array_filename)
		intensity_array = rpy.r.intensity(array)	#return a lengthX1 2-Dimensional array.
		intensity_array_size = len(intensity_array)
		array_width = int(math.sqrt(intensity_array_size))	#assume it's square array
		returnData = PassingData(array=array, intensity_array=intensity_array, array_width=array_width)
		return returnData
	
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
		array_width = None
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
		array_label_ls = []
		for i in range(no_of_objects):
			array_id, filename, ecotype_id = rows[i][:3]
			array_id_avail_ls.append(array_id)
			array_label_ls.append('%s_%s'%(array_id, ecotype_id))
			
			if array_file_directory and os.path.isdir(array_file_directory):
				filename = os.path.join(array_file_directory, os.path.split(filename)[1])
			
			sys.stderr.write("\t%d/%d: Extracting intensity from %s ... \n"%(i+1, no_of_objects, filename))
			
			if run_type==1:	#output SNP probe intensity within the loop
				output_fname = os.path.join(output_dir, '%s_array_intensity.tsv'%(array_id))
				if os.path.isfile(output_fname):
					sys.stderr.write("\tFile %s already exists. Ignore.\n"%(output_fname))
					continue
			
			#read array by calling R
			if array_width == None:
				returnData = self.getArrayWidth(filename)
				intensity_array = returnData.intensity_array
				array = returnData.array
				array_width = returnData.array_width
			else:
				array = rpy.r.read_affybatch(filenames=filename)
				intensity_array = rpy.r.intensity(array)	#return a lengthX1 2-Dimensional array.
			
			if run_type==2:	#CNV probe
				for j in range(len(xy_ls)):
					xpos, ypos = xy_ls[j]
					#chromosome, position = chr_pos_ls[j]
					intensity_array_index = array_width*(array_width - xpos - 1) + ypos
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
						intensity_array_index = array_width*(array_width - one_probe.xpos - 1) + one_probe.ypos
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
	
	@classmethod
	def organizeProbesIntoChromosome(cls, xy_ls, chr_pos_ls, probes_id_ls):
		"""
		2009-11-24
			split out of calculateProbeQuartilePerChromosome()
		"""
		sys.stderr.write("Getting probes into each chromosome ...")
		chr2xy_ls = {}
		chr2probe_id_ls = {}
		for i in range(len(xy_ls)):
			chr,pos = chr_pos_ls[i]
			if chr not in chr2xy_ls:
				chr2xy_ls[chr] = []
				chr2probe_id_ls[chr] = []	#initialize with the start_probe_id
			chr2xy_ls[chr].append(xy_ls[i])
			chr2probe_id_ls[chr].append(probes_id_ls[i])
		sys.stderr.write("Done.\n")
		return chr2xy_ls, chr2probe_id_ls
	
	@classmethod
	def calculateProbeQuartilePerChromosome(cls, session, xy_ls, chr_pos_ls, probes_id_ls):
		"""
		2009-11-24
			calculate DLRSpread = IQR(dLR) / 4*erfinv(0.5), where dLR is an array of
			differences between log ratios of adjacent probes along the chromosome, erfinv is the Inverse Error
			Function and IQR is Inter Quartile Range. store the value in db.
		2009-11-20
			calculate CNV probe intensity quartile (1st, median, 3rd) for each chromosome and store them in database
		"""
		import numpy, math
		import rpy
		rpy.r.library('affy')
		
		chr2xy_ls, chr2probe_id_ls = cls.organizeProbesIntoChromosome(xy_ls, chr_pos_ls, probes_id_ls)
		
		rows = Stock_250kDB.ArrayInfo.query()
		no_of_arrays = rows.count()
		array_width = None
		counter = 0
		for row in rows:
			array_id = row.id
			filename = row.filename
			ecotype_id = row.maternal_ecotype_id
						
			sys.stderr.write("\t%d/%d: Extracting intensity from %s ... "%(counter+1, no_of_arrays, array_id))
			
			#read array by calling R
			array = rpy.r.read_affybatch(filenames=filename)
			intensity_array = rpy.r.intensity(array)	#return a lengthX1 2-Dimensional array.
			if array_width == None:
				intensity_array_size = len(intensity_array)
				array_width = int(math.sqrt(intensity_array_size))	#assume it's square array
			
			chr2intensity_ls = {}
			chr2dlR_ls = {}
			no_of_outliers = 0
			no_of_outliers_saved = 0
			for chr, chr_xy_ls in chr2xy_ls.iteritems():
				chr2intensity_ls[chr] = []
				intensity_ls = chr2intensity_ls[chr]
				chr2dlR_ls[chr] = []	# dLR is an array of differences between log ratios of adjacent probes along the chromosome
				dlR_ls = chr2dlR_ls[chr]
				for xpos, ypos in chr_xy_ls:
					intensity_array_index = array_width*(array_width - xpos - 1) + ypos
					intensity = math.log10(intensity_array[intensity_array_index][0])
					intensity_ls.append(intensity)
					current_probe_index = len(intensity_ls)-1
					if current_probe_index>=1:
						dlR_ls.append(intensity_ls[current_probe_index]-intensity_ls[current_probe_index-1])
				start_probe_id=chr2probe_id_ls[chr][0]
				stop_probe_id=chr2probe_id_ls[chr][-1]
				# check if this array_quartile is already in db.
				array_quartile = Stock_250kDB.ArrayQuartile.query.filter_by(array_id=array_id).filter_by(start_probe_id=start_probe_id).\
						filter_by(stop_probe_id=stop_probe_id).first()
				if not array_quartile:	# not in db. create one.
					array_quartile = Stock_250kDB.ArrayQuartile(array_id=array_id, start_probe_id=chr2probe_id_ls[chr][0],\
														stop_probe_id=chr2probe_id_ls[chr][-1])
				array_quartile.no_of_probes = len(intensity_ls)
				calculate7NumberSummaryForOneList(intensity_ls, array_quartile)
				
				IQR_dlR = stats.scoreatpercentile(dlR_ls, 75) - stats.scoreatpercentile(dlR_ls, 25)
				array_quartile.dlrspread = IQR_dlR/(4*scipy.special.erfinv(0.5))
				dlR_7NumberSummary = calculate7NumberSummaryForOneList(dlR_ls)
				for summary_name in  ["minimum", "first_decile", "lower_quartile", "median", "upper_quartile", "last_decile", "maximum"]:
					summary_value = getattr(dlR_7NumberSummary, summary_name)
					attribute_name = 'dlr_%s'%summary_name
					setattr(array_quartile, attribute_name, summary_value)
				
				session.save_or_update(array_quartile)
				
				
				# find and store the outliers
				IQR = array_quartile.upper_quartile - array_quartile.lower_quartile
				lower_whisker = array_quartile.lower_quartile - 1.5*IQR
				upper_whisker = array_quartile.upper_quartile + 1.5*IQR
				for i in range(len(intensity_ls)):
					intensity = intensity_ls[i]
					if intensity<lower_whisker or intensity>upper_whisker:
						probe_id = chr2probe_id_ls[chr][i]
						no_of_outliers += 1
						if array_quartile.id:
							array_quartile_outlier = Stock_250kDB.ArrayQuartileOutlier.query.filter_by(probe_id=probe_id).\
								filter_by(array_quartile_id=array_quartile.id).first()
						else:
							array_quartile_outlier = None
						if not array_quartile_outlier:	# not in db. add it to db
							array_quartile_outlier = Stock_250kDB.ArrayQuartileOutlier(probe_id=probe_id, intensity=intensity)
							array_quartile_outlier.array_quartile = array_quartile
							session.save_or_update(array_quartile_outlier)
							no_of_outliers_saved += 1
				
				session.flush()
			counter += 1
			sys.stderr.write("%s outliers (%s saved). Done.\n"%(no_of_outliers, no_of_outliers_saved))
	
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
		if self.run_type!=3:
			probes, xy_ls, chr_pos_ls, probes_id_ls = self.get_probes(curs, self.probes_table, snps, \
																			run_type=self.run_type)
		else:
			probes, xy_ls, chr_pos_ls, probes_id_ls = None, None, None, None
		
		if self.run_type in (1,2,3):
			self.outputArray(session, curs, self.output_dir, self.array_info_table, snps, probes, self.array_id_ls, xy_ls, \
						chr_pos_ls, probes_id_ls, call_method_id=self.call_method_id, run_type=self.run_type, \
						array_file_directory=self.array_file_directory)
		elif self.run_type==4:
			self.calculateProbeQuartilePerChromosome(session, xy_ls, chr_pos_ls, probes_id_ls)
		else:
			sys.stderr.write("Run type %s is not supported.\n"%self.run_type)
			
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
