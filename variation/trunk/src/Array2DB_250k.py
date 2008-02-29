"""
2008-02-26
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, numpy, stat, rpy
import warnings, traceback, gc
from pymodule import process_function_arguments

class Array2DB_250k:
	"""
	2008-02-26
		class to process phenotype data from at.phenotype and at.experiment to get flowering time.
		Flowering time is "time of first flower open" - "date counted as germination"
		
	Argument list:
		-z ..., --hostname=...	the hostname, localhost(default)
		-d ..., --dbname=...	the database name, stock20071008(default)
		-k ..., --schema=...	which schema in the database, (IGNORE)
		-i ...,	input_dir*
		-a ...,	array_data_table, array_data(default)
		-e ...,	probes_table, 'probes'(default)
		-f ...,	strain_info_table, 'strain_info'(default)
		-g ...,	array_info_table, 'array_info'(default)
		-c,	commit db transaction
		-b,	toggle debug
		-r, toggle report
	Examples:
		main.py -y 3 -i /Network/Data/250k/raw_data/ -d stock_250k -c
	Description:
		Dump .cel array files in input_dir (files in up to 2-level directories) into db.
	"""
	def __init__(self, **keywords):
		"""
		2008-02-28
		"""
		argument_default_dict = {('hostname',1, ):'localhost',\
								('dbname',1, ):'stock_250k',\
								('schema',1, ):'',\
								('input_dir',1, ):None,\
								('array_data_table',1, ):'array_data',\
								('probes_table',1, ):'probes',\
								('strain_info_table',1, ):'strain_info',\
								('array_info_table',1, ):'array_info',\
								('commit',0, int):0,\
								('debug',0, int):0,\
								('report',0, int):0}
		"""
		2008-02-28
			argument_default_dict is a dictionary of default arguments, the key is a tuple, ('argument_name', is_argument_required, argument_type)
			argument_type is optional
		"""
		#argument dictionary
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__)
		self.debug = self.ad['debug']
		self.report = self.ad['report']
		
	def get_xypos2probes_id(self, curs, probes_table):
		"""
		2008-02-28
		"""
		sys.stderr.write("Getting xypos2probes_id ... ")
		curs.execute("select id, xpos, ypos from %s"%(probes_table))
		rows = curs.fetchall()
		xypos2probes_id = {}
		for row in rows:
			probes_id, xpos, ypos = row
			xypos2probes_id[(xpos, ypos)] = probes_id
		sys.stderr.write("Done.\n")
		return xypos2probes_id
	
	def get_filename2array_id_in_db(self, curs, array_info_table):
		"""
		2008-02-28
		"""
		sys.stderr.write("Getting filename2array_id_in_db ... ")
		filename2array_id_in_db = {}
		curs.execute("select id, filename from %s"%(array_info_table))
		rows = curs.fetchall()
		for row in rows:
			array_id, filename = row
			filename2array_id_in_db[filename] = array_id
		sys.stderr.write("Done.\n")
		return filename2array_id_in_db
	
	def get_filename2array_id(self, input_dir, filename2array_id_in_db):
		"""
		2008-02-29
			new_array_id starts from 1 + maximum avaible array_id in db
		2008-02-28
		"""
		sys.stderr.write("Getting filename2array_id ... ")
		file_dir_ls = os.listdir(input_dir)
		if self.debug:
			sys.stderr.write("\n\tTotally, %d 1st-order file/dir's to be processed.\n"%len(file_dir_ls))
		filename_ls = []
		for i in range(len(file_dir_ls)):
			file_dir = file_dir_ls[i]
			if self.debug:
				sys.stderr.write("%d/%d:\t%s\n"%(i+1,len(file_dir_ls),file_dir))
			pathname = os.path.join(input_dir, file_dir)
			if os.path.isdir(pathname):
				sub_file_dir_ls = os.listdir(pathname)
				for sub_file_dir in sub_file_dir_ls:
					sub_pathname = os.path.join(pathname, sub_file_dir)
					if os.path.isfile(sub_pathname):
						filename_ls.append(sub_pathname)
					else:
						warnings.warn("%s is not file. (might be 2nd directory). Ignored.\n"%(sub_pathname))
			elif os.path.isfile(pathname):
				filename_ls.append(pathname)
			else:
				warnings.warn("%s is neither directory nor file. Ignored.\n"%(pathname))
		
		filename2array_id = {}
		if filename2array_id_in_db:	#there are arrays already existing in db.
			new_array_id = max(filename2array_id_in_db.values())+1
		else:
			new_array_id = 1	#auto increment primary key in sql starts from 1. although mistake is already committed.
		for pathname in filename_ls:
			file_ext = os.path.splitext(pathname)[1].lower()	#the extension of a filename and lower case
			if file_ext=='.cel' and pathname not in filename2array_id_in_db:	#make sure it's .cel file and not already in db
				filename2array_id[pathname] = new_array_id
				new_array_id += 1
			else:
				warnings.warn("%s is either not .cel or already in db. Ignored.\n"%(pathname))
		sys.stderr.write("Done.\n")
		return filename2array_id
	
	def submit_filename2array_id(self, curs, filename2array_id, array_info_table):
		"""
		2008-02-28
		"""
		sys.stderr.write("Submitting filename2array_id ... ")
		for filename, array_id in filename2array_id.iteritems():
			curs.execute("insert into %s(filename, id) values ('%s', %s)"%(array_info_table, filename, array_id))
		sys.stderr.write("Done.\n")
	
	def submit_one_array(self, curs, array_data_table, array_id, intensity_array, array_data_with_xypos):
		"""
		2008-02-28
			intensity_array is 2D array although the 2nd dimension is only of size 1.
		"""
		sys.stderr.write("Submitting one array data ... ")
		for i in range(len(intensity_array)):
			array_data_with_xypos_row = array_data_with_xypos[i]
			if len(array_data_with_xypos_row)==2:
				xpos, ypos = array_data_with_xypos_row
				curs.execute("insert into %s(array_id, intensity, xpos, ypos) values (%s, %s, %s, %s)"%\
							(array_data_table, array_id, intensity_array[i][0], xpos, ypos))
			elif len(array_data_with_xypos_row)==3:
				probes_id, xpos, ypos = array_data_with_xypos_row
				curs.execute("insert into %s(array_id, probes_id, intensity, xpos, ypos) values (%s, %s, %s, %s, %s)"%\
							(array_data_table, array_id, probes_id, intensity_array[i][0], xpos, ypos))
			else:
				warnings.warn("%s is neither 2 nor 3 columns. Ignored."%(array_data_with_xypos_row))
		sys.stderr.write("Done.\n")
	
	def submit_all_array_data(self, filename2array_id, xypos2probes_id, curs, array_data_table):
		"""
		2008-02-28
		"""
		sys.stderr.write("Submitting array data ... \n")
		rpy.r.library('affy')
		files = filename2array_id.keys()
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		array_size = None
		for i in range(len(files)):
			f = files[i]
			array_id = filename2array_id.get(f)
			sys.stderr.write("%d/%d:\t%s\n"%(i+1,len(files),f))
			array = rpy.r.read_affybatch(filenames=f)
			intensity_array = rpy.r.intensity(array)	#return a lengthX1 2-Dimensional array.
			intensity_array_size = len(intensity_array)
			if array_size == None:
				array_size = int(math.sqrt(intensity_array_size))	#assume it's square array
			array_data_with_xypos = []
			for i in range(intensity_array_size):
				"""
				2008-02-28 mapping between array index and (xpos, ypos) is inferred from Xu Zhang's readcel.R
				basically, the array is stored by row and the last row (xpos is reversed) comes up first.
				"""
				xpos = array_size - i/array_size -1
				ypos = i%array_size
				xypos = (xpos, ypos)
				probes_id = xypos2probes_id.get(xypos)
				if probes_id!=None:
					array_data_with_xypos.append([probes_id, xpos, ypos])
				else:
					array_data_with_xypos.append([xpos, ypos])
			self.submit_one_array(curs, array_data_table, array_id, intensity_array, array_data_with_xypos)
			if self.debug:
				sys.stderr.write("Current thresholds: %s.\n"%repr(gc.get_threshold()) )
				gc.collect()
				sys.stderr.write("%s uncollectable objects are deleted.\n"%(len(gc.garbage)) )
				del gc.garbage[:]
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		2008-02-28
		"""
		import MySQLdb
		conn = MySQLdb.connect(db=self.ad['dbname'], host=self.ad['hostname'])
		curs = conn.cursor()
		filename2array_id_in_db = self.get_filename2array_id_in_db(curs, self.ad['array_info_table'])
		filename2array_id = self.get_filename2array_id(self.ad['input_dir'], filename2array_id_in_db)
		if self.ad['commit']:
			self.submit_filename2array_id(curs, filename2array_id, self.ad['array_info_table'])
		xypos2probes_id = self.get_xypos2probes_id(curs, self.ad['probes_table'])
		if self.ad['commit']:
			self.submit_all_array_data(filename2array_id, xypos2probes_id, curs, self.ad['array_data_table'])
			curs.execute("commit")
