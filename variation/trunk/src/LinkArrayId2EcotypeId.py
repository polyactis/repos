"""
2008-03-11
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
import traceback, gc
from pymodule import process_function_arguments

class LinkArrayId2EcotypeId:
	"""
	Argument list:
		-z ..., --hostname=...	the hostname, localhost(default)
		-d ..., --dbname=...	the database name, stock20071008(default)
		-k ..., --schema=...	which schema in the database, (IGNORE)
		-i ...,	array_info_table, stock_250k.array_info(default)*
		-a ...,	calls_comment_table, stock20071008.calls_250k_duplicate_comment(default)(argument1)*
		-e ...,	ecotype_table, 'stock20071008.ecotype'(default)(argument2)*
		-c,	commit db transaction
		-b,	toggle debug
		-r, toggle report
	Examples:
		main.py -y 4 -i stock_250k.array_info -c
	Description:
		2008-03-11
		Link 250k's array id to ecotype id. Either by matching strains already linked to ecotype id or de novo db querying.
	"""
	def __init__(self, **keywords):
		"""
		2008-03-11
		"""
		argument_default_dict = {('hostname',1, ):'localhost',\
								('dbname',1, ):'stock20071008',\
								('schema',0, ):'',\
								('ecotype_table',1, ):'stock20071008.ecotype',\
								('array_info_table',1, ):'stock_250k.array_info',\
								('calls_comment_table',1, ):'stock20071008.calls_250k_duplicate_comment',\
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
	
	def get_array_id_filename_ls(self, curs, array_info_table):
		"""
		2008-03-11
			
		"""
		sys.stderr.write("Getting array_id_filename_ls from %s ... "%(array_info_table))
		curs.execute("select id, filename from %s where ecotype_id is NULL and maternal_ecotype_id is NULL and paternal_ecotype_id is NULL"%\
					(array_info_table))
		rows = curs.fetchall()
		array_id_filename_ls = []
		for row in rows:
			array_id_filename_ls.append(row)
		sys.stderr.write("Done.\n")
		return array_id_filename_ls
	
	def get_known_filename2ecotype_id(self, curs, calls_comment_table):
		"""
		2008-03-11
			get the old mapping, which was done by Calls2DB_250k.py 2 months ago.
			the key is (last-level-dir, filename)
		"""
		sys.stderr.write("Getting known_filename2ecotype_id from %s ... "%(calls_comment_table))
		curs.execute("select ecotypeid, comment from %s"%(calls_comment_table))
		rows = curs.fetchall()
		known_filename2ecotype_id = {}
		yanli8_8_07_map = {'LY-1':'Col-0A',	#2008-03-11 yanli8_8_07_map is a map of 6 old strain names, later modified to its real strain name
						'LY-2':'Col-0B',
						'LY-3':'Ler-1A',
						'LY-4':'Ler-1B',
						'LY-5':'Van-0A',
						'LY-6':'Van-0B'}
		import re
		true_filename_without_index_rp = re.compile(r'[\d]+__(.*)')	#2008-03-11 chichi prepends something like 2__ to filename. i.e. 2__Omo2_3_base-calls.txt
		for row in rows:
			ecotypeid, pathname = row
			dir, filename = os.path.split(pathname)
			rp_search_result = true_filename_without_index_rp.search(filename)
			if rp_search_result:
				filename = rp_search_result.groups()[0]	#take the first group out as true filename
				strain_name = filename[:-15]		#2008-03-11 chichi not only prepends 2__ to filename, also replaced '-' with '_' in the strain names
				strain_name = strain_name.replace('_', '-')
				filename = '%s%s'%(strain_name, filename[-15:])
			if filename[:-15] in yanli8_8_07_map:	#index from -15  to the last corresponds to '_base-calls.txt'
				filename = '%s%s'%(yanli8_8_07_map[filename[:-15]], filename[-15:])
			last_dir = os.path.split(dir)[1]
			pathname = os.path.join(last_dir, filename)
			if pathname in known_filename2ecotype_id:
				sys.stderr.write("Error: %s already exists in known_filename2ecotype_id with ecotypeid=%s.\n"%\
								(pathname, known_filename2ecotype_id[pathname]))
				sys.exit(2)
			known_filename2ecotype_id[pathname] = ecotypeid
		sys.stderr.write("Done.\n")
		return known_filename2ecotype_id
	
	def link_array_id_using_old_map(self, array_id_filename_ls, known_filename2ecotype_id):
		"""
		2008-03-11
		"""
		sys.stderr.write("Linking array id to ecotype id given known_filename2ecotype_id ... ")
		array_id2ecotype_id_pair = {}
		array_id_filename_left_ls = []
		for array_id, pathname in array_id_filename_ls:
			dir, filename = os.path.split(pathname)
			last_dir = os.path.split(dir)[1]
			filename_pre, filename_ext = os.path.splitext(filename)
			#conform to the key in known_filename2ecotype_id
			filename_pre = filename_pre.replace(' ', '')	#all spaces were deleted in Chichi's base-call output filename
			base_call_filename = '%s_base-calls.txt'%filename_pre
			new_pathname = os.path.join(last_dir, base_call_filename)
			if new_pathname in known_filename2ecotype_id:
				ecotype_id = known_filename2ecotype_id[new_pathname]
				array_id2ecotype_id_pair[array_id] = [ecotype_id, ecotype_id]
			else:
				array_id_filename_left_ls.append([array_id, pathname])
		sys.stderr.write("%s matched and %s left.\n"%(len(array_id2ecotype_id_pair), len(array_id_filename_left_ls)))
		return array_id2ecotype_id_pair, array_id_filename_left_ls
	
	def denovo_link_array_id(self, curs, array_id_filename_ls, ecotype_table, array_id2ecotype_id_pair):
		"""
		2008-03-11
			modelled after find_out_ecotypeid_given_strain_name() from Calls2DB_250k.py
		"""
		sys.stderr.write("Trying to find ecotype id by querying database ...\n")
		for array_id, pathname in array_id_filename_ls:
			strain_name_pair = raw_input("Enter one or two strain names, separated by ',' for %s (q to stop manual linking and submit already linked):"%(pathname))
			if strain_name_pair=='q':
				break
			strain_name_pair = strain_name_pair.split(',')
			ecotype_id_pair = []
			for strain_name in strain_name_pair:
				ecotype_id = self.find_out_ecotypeid_given_strain_name(curs, strain_name, ecotype_table)
				if ecotype_id!=None:
					ecotype_id_pair.append(ecotype_id)
			if ecotype_id_pair:
				if len(ecotype_id_pair)==1:	#maternal/paternal are same.
					ecotype_id_pair.append(ecotype_id_pair[0])
				array_id2ecotype_id_pair[array_id] = ecotype_id_pair
		sys.stderr.write("Done.\n")
	
	def find_out_ecotypeid_given_strain_name(self, curs, strain_name, ecotype_table='ecotype'):
		"""
		2008-03-11
			modelled after find_out_ecotypeid_given_strain_name() from Calls2DB_250k.py
		"""
		while 1:
			curs.execute("select e.id, e.name, e.stockparent, e.nativename, s.name, c.abbr from %s e, site s, address a, country c where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id and (e.nativename rlike '%s'or e.name rlike '%s')"%(ecotype_table, strain_name, strain_name))
			rows = curs.fetchall()
			choice_id2ecotypeid = {}
			choice_id = 0
			header_ls = ['', 'ecotypeid', 'name', 'stockparent', 'nativename', 'site_name', 'country']
			sys.stderr.write("%s\n"%('\t'.join(header_ls)))
			for row in rows:
				ecotypeid, name, stockparent, nativename, site_name, country = row
				choice_id += 1
				choice_id2ecotypeid[repr(choice_id)] = ecotypeid
				ls = ['%s:'%choice_id]
				ls.append(list(row))
				ls = map(repr, ls)
				sys.stderr.write("%s\n"%'\t'.join(ls))
			if len(choice_id2ecotypeid)==0:
				strain_name = raw_input("Found no match for %s. Enter a new strain name to retry (q to give up):"%strain_name)
				if strain_name == 'q':
					return None
			else:
				break
		choice_instructions = "Enter the row number corresponding to the ecotypeid(c to continue on next strain):"
		choice_id = raw_input(choice_instructions)
		while 1:
			if choice_id =='c':
				break
			elif choice_id in choice_id2ecotypeid:
				break
			else:
				sys.stderr.write("%s is not a choice.\n"%(choice_id))
				choice_id = raw_input(choice_instructions)
		if choice_id=='c':
			return None
		else:
			return choice_id2ecotypeid[choice_id]
	
	def update_array_info_table(self, curs, array_info_table, array_id2ecotype_id_pair):
		"""
		2008-03-11
		"""
		sys.stderr.write("Updating %s entries in %s ... "%(len(array_id2ecotype_id_pair), array_info_table))
		for array_id, ecotype_id_pair in array_id2ecotype_id_pair.iteritems():
			curs.execute("update %s set maternal_ecotype_id=%s, paternal_ecotype_id=%s where id=%s"%(array_info_table, ecotype_id_pair[0], ecotype_id_pair[1], array_id))
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		2008-03-11
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		import MySQLdb
		conn = MySQLdb.connect(db=self.ad['dbname'], host=self.ad['hostname'])
		curs = conn.cursor()
		array_id_filename_ls = self.get_array_id_filename_ls(curs, self.ad['array_info_table'])
		known_filename2ecotype_id = self.get_known_filename2ecotype_id(curs, self.ad['calls_comment_table'])
		array_id2ecotype_id_pair, array_id_filename_left_ls = self.link_array_id_using_old_map(array_id_filename_ls, known_filename2ecotype_id)
		self.denovo_link_array_id(curs, array_id_filename_left_ls, self.ad['ecotype_table'], array_id2ecotype_id_pair)
		if self.ad['commit']:
			self.update_array_info_table(curs, self.ad['array_info_table'], array_id2ecotype_id_pair)
			curs.execute("commit")