#!/usr/bin/env python
"""
Usage:	Array2DB_250k.py [OPTIONS] -i INPUT_DIR

Argument list:
	-z ..., --hostname=...	the hostname, papaya.usc.edu(default)
	-d ..., --dbname=...	the database name, stock_250k(default)
	-u ..., --user=...	the db username, (otherwise it will ask for it).
	-p ..., --passwd=...	the db password, (otherwise it will ask for it).
	-i ...,	input_dir*	where all the .cel files sit
	-o ...,	output_dir	to store the renamed .cel files. '/Network/Data/250k/db/raw_data/'(default)
	-g ...,	array_info_table, 'array_info'(default)
	-c,	commit db transaction
	-b,	toggle debug
	-r, toggle report
Examples:
	Array2DB_250k.py -i /Network/Data/250k/raw_data/ -d stock_250k -c

Description:
	Dump .cel array files from input_dir (files in up to 2-level directories) into db and associated file-system storage (output-dir).

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, stat, getopt
import traceback, gc, subprocess

class ArrayInfo(object):
	def __init__(self, **keywords):
		from pymodule import process_function_arguments
		argument_default_dict = {('curs',1, ): None,\
								('array_info_table',1, ):None,\
								('user',1, ):None,\
								('debug',0, int):0,\
								('report',0, int):0}
		"""
		2008-02-28
			argument_default_dict is a dictionary of default arguments, the key is a tuple, ('argument_name', is_argument_required, argument_type)
			argument_type is optional
		"""
		#argument dictionary
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self, howto_deal_with_required_none=2)
		self.md5sum2array_id, self.max_array_id = self.get_md5sum2array_id(self.curs, self.array_info_table)
		
	def get_md5sum2array_id(self, curs, array_info_table):
		"""
		2008-04-10
		"""
		sys.stderr.write("Getting filename2array_id_in_db ... ")
		md5sum2array_id = {}
		curs.execute("select id, md5sum, filename from %s"%(array_info_table))
		rows = curs.fetchall()
		max_array_id = 0
		for row in rows:
			array_id, md5sum, filename = row
			if md5sum in md5sum2array_id:
				sys.stderr.write("Error: two array_id (%s and %s) have same md5sum %s in db. Check db.\n"%(array_id, md5sum2array_id[md5sum], md5sum))
				sys.exit(3)
			md5sum2array_id[md5sum] = array_id
			if array_id>max_array_id:
				max_array_id = array_id
		sys.stderr.write("Done.\n")
		return md5sum2array_id, max_array_id
	
	get_md5sum2array_id = classmethod(get_md5sum2array_id)
	
	def get_md5sum(self, filename):
		"""
		"""
		md5sum_command = 'md5sum'
		md5sum_p = subprocess.Popen([md5sum_command, filename], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		md5sum_stdout_out = md5sum_p.stdout.read()
		md5sum_stderr_out = md5sum_p.stderr.read()
		if md5sum_stderr_out:
			sys.stderr.write("%s %s failed with stderr: %s.\n"%(md5sum_command, filename, md5sum_stderr_out))
			sys.exit(4)
		else:
			return md5sum_stdout_out.split()[0]
	
	get_md5sum = classmethod(get_md5sum)
	
	def assignNewIdToThisArray(self, array_filename, output_dir):
		md5sum = self.get_md5sum(array_filename)
		if md5sum in self.md5sum2array_id:
			sys.stderr.write("%s already exists in db with md5sum=%s and array_id=%s. ignored.\n"%\
							(array_filename, md5sum, self.md5sum2array_id[md5sum]))
			return -1
		else:
			new_array_id = self.max_array_id + 1
			output_fname = os.path.join(output_dir, '%s_raw_data.cel'%new_array_id)
			
			cp_p = subprocess.Popen(['cp', array_filename, output_fname], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
			cp_p_stdout_out = cp_p.stdout.read()
			cp_p_stderr_out = cp_p.stderr.read()
			if cp_p_stderr_out:
				sys.stderr.write("copy error: %s\n"%cp_p_stderr_out)
				return -1
			else:
				self.curs.execute("insert into %s(id, filename, original_filename, md5sum, created_by) values (%s, '%s', '%s', '%s', '%s')"%\
							(self.array_info_table, new_array_id, output_fname, array_filename, md5sum, self.user))
				self.md5sum2array_id[md5sum] = new_array_id
				self.max_array_id += 1
				return new_array_id

class Array2DB_250k(object):
	__doc__ = __doc__	#use documentation in the beginning of the file as this class's doc
	def __init__(self, **keywords):
		"""
		2008-04-09
			array_data_table, probes_table, strain_info_table are no longer required. but leave them in the argument_default_dict.
		2008-02-28
		"""
		from pymodule import process_function_arguments
		argument_default_dict = {('hostname',1, ):'papaya.usc.edu',\
								('dbname',1, ):'stock_250k',\
								('user',1, ):None,\
								('passwd',1, ):None,\
								('input_dir',1, ):None,\
								('output_dir',1, ):'/Network/Data/250k/db/raw_data/',\
								('array_data_table',0, ):'array_data',\
								('probes_table',0, ):'probes',\
								('strain_info_table',0, ):'strain_info',\
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
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
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
	
	def get_all_files_in_input_dir(self, input_dir):
		"""
		2008-04-10
		"""
		sys.stderr.write("Getting all files from %s ... \n"%input_dir)
		file_dir_ls = os.listdir(input_dir)
		sys.stderr.write("\n\tTotally, %d 1st-order file/dir's to be processed.\n"%len(file_dir_ls))
		filename_ls = []
		no_of_objects = len(file_dir_ls)
		for i in range(no_of_objects):
			file_dir = file_dir_ls[i]
			sys.stderr.write("%d/%d:\t%s\n"%(i+1, no_of_objects, file_dir))
			pathname = os.path.join(input_dir, file_dir)
			if os.path.isdir(pathname):
				sub_file_dir_ls = os.listdir(pathname)
				for sub_file_dir in sub_file_dir_ls:
					sub_pathname = os.path.join(pathname, sub_file_dir)
					if os.path.isfile(sub_pathname):
						filename_ls.append(sub_pathname)
					else:
						sys.stderr.write("%s is not file. (might be 2nd directory). Ignored.\n"%(sub_pathname))
			elif os.path.isfile(pathname):
				filename_ls.append(pathname)
			else:
				sys.stderr.write("%s is neither directory nor file. Ignored.\n"%(pathname))
		sys.stderr.write("Done.\n")
		return filename_ls
			
	def get_filename2array_id(self, input_dir, filename2array_id_in_db):
		"""
		2008-02-29
			new_array_id starts from 1 + maximum avaible array_id in db
		2008-02-28
		"""
		sys.stderr.write("Getting filename2array_id ... \n")
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
						sys.stderr.write("%s is not file. (might be 2nd directory). Ignored.\n"%(sub_pathname))
			elif os.path.isfile(pathname):
				filename_ls.append(pathname)
			else:
				sys.stderr.write("%s is neither directory nor file. Ignored.\n"%(pathname))
		
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
				sys.stderr.write("%s is either not .cel or already in db. Ignored.\n"%(pathname))
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
		2008-02-29
			xpos and ypos are no longer inserted into db if probes_id is avaible.
		2008-02-28
			intensity_array is 2D array although the 2nd dimension is only of size 1.
		"""
		sys.stderr.write("Submitting one array data ... ")
		no_of_points_without_probes_id = 0
		count = 0
		for i in range(len(intensity_array)):
			array_data_with_xypos_row = array_data_with_xypos[i]
			if len(array_data_with_xypos_row)==2:
				xpos, ypos = array_data_with_xypos_row
				curs.execute("insert into %s(array_id, intensity, xpos, ypos) values (%s, %s, %s, %s)"%\
							(array_data_table, array_id, intensity_array[i][0], xpos, ypos))
				no_of_points_without_probes_id += 1
			elif len(array_data_with_xypos_row)==1:
				probes_id = array_data_with_xypos_row[0]
				curs.execute("insert into %s(array_id, probes_id, intensity) values (%s, %s, %s)"%\
							(array_data_table, array_id, probes_id, intensity_array[i][0]))
			else:
				sys.stderr.write("%s is neither 2 nor 3 columns. Ignored.\n"%(array_data_with_xypos_row))
			count += 1
		sys.stderr.write("%s/%s have no probes_id. Done.\n"%(no_of_points_without_probes_id, count))
	
	def submit_all_array_data(self, filename2array_id, xypos2probes_id, curs, array_data_table):
		"""
		2008-02-29
			xpos and ypos are no longer inserted into db if probes_id is avaible.
		2008-02-28
		"""
		sys.stderr.write("Submitting array data ... \n")
		import rpy
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
					array_data_with_xypos.append([probes_id])
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
		2008-04-11
			new way of handling raw cel files:
			
			md5sum each array, check if each array is in db already.
			if yes:
				skip it
			else:
				assign a new id and insert an entry into array_info_table
				copy the original cel file to output_dir and put array_id in the beginning of the filename.
		2008-02-28
		"""
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.user, passwd = self.passwd)
		curs = conn.cursor()
		"""
		filename2array_id_in_db = self.get_filename2array_id_in_db(curs, self.array_info_table)
		filename2array_id = self.get_filename2array_id(self.input_dir, filename2array_id_in_db)
		if self.commit:
			self.submit_filename2array_id(curs, filename2array_id, self.array_info_table)
		xypos2probes_id = self.get_xypos2probes_id(curs, self.probes_table)
		self.submit_all_array_data(filename2array_id, xypos2probes_id, curs, self.array_data_table)
		"""
		arrayInfo = ArrayInfo(curs=curs, array_info_table=self.array_info_table, user=self.user)
		input_fname_ls = self.get_all_files_in_input_dir(self.input_dir)
		for filename in input_fname_ls:
			sys.stderr.write("Assigning new id to %s ... "%filename)
			return_value = arrayInfo.assignNewIdToThisArray(filename, self.output_dir)
			if return_value==-1:
				sys.stderr.write("Failed.\n")
			else:
				sys.stderr.write("\n")
		if self.commit:
			curs.execute("commit")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "user=", "passwd=", "help", "type=", "debug", "report"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:u:p:i:o:g:cbr", long_options_list)
	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	hostname = None
	dbname = None
	user = None
	passwd = None
	input_dir = None
	output_dir = None
	array_info_table = None
	commit = None
	debug = None
	report = None
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-u", "--user"):
			user = arg
		elif opt in ("-p", "--passwd"):
			passwd = arg
		elif opt in ("-i",):
			input_dir = arg
		elif opt in ("-o",):
			output_dir = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-g",):
			array_info_table = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	instance = Array2DB_250k(hostname=hostname, dbname=dbname, user=user, passwd=passwd, input_dir=input_dir, output_dir=output_dir,
					array_info_table=array_info_table, \
					commit=commit, debug=debug, report=report)
	instance.run()
