#!/usr/bin/env python
"""

Examples:
	LinkArrayId2EcotypeId.py -i stock_250k.array_info -c
	
Description:
	2008-03-11
	Link 250k's array id to ecotype id. Either by matching strains already linked to ecotype id or de novo db querying.
	
	calls_comment_table has two columns, ecotypeid and comment (array filename). If not supplied, it's skipped.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
import traceback, gc

class LinkArrayId2EcotypeId(object):
	__doc__ = __doc__
	option_default_dict = {('z', 'hostname', 1, 'hostname of the db server', 1, ): 'papaya.usc.edu',\
							('d', 'dbname', 1, '', 1, ): 'stock_250k',\
							('u', 'user', 1, 'database username', 1, ):None,\
							('p', 'passwd', 1, 'database password', 1, ):None,\
							('e', 'ecotype_table', 1, 'Query this table to get ecotypeid', 1, ): 'stock.ecotype',\
							('m', 'mapping_file',1, 'a file mapping filename in input_dir to ecotypeid', 0, ): None,\
							('a', 'calls_comment_table', 1, 'table from which to retrieve existing array-ecotypeid mapping.', 0, ): None,\
							('i', 'array_info_table', 1, 'Table where all arrays are', 1, ): 'array_info',\
							('s', 'stock_149SNP_db', 1, 'database name for 149SNP data', 1, ): 'stock',\
							('c', 'commit', 0, 'commit db transaction', 0, int):0,\
							('b', 'debug', 0, 'toggle debug mode', 0, int):0,\
							('r', 'report', 0, 'toggle report, more verbose stdout/stderr.', 0, int):0}
	"""
	2008-04-40
		option_default_dict is a dictionary for option handling, including argument_default_dict info
		the key is a tuple, ('short_option', 'long_option', has_argument, description_for_option, is_option_required, argument_type)
		argument_type is optional
	"""
	def __init__(self, **keywords):
		"""
		2008-03-11
		"""
		from pymodule import process_function_arguments, turn_option_default_dict2argument_default_dict
		argument_default_dict = turn_option_default_dict2argument_default_dict(self.option_default_dict)
		#argument dictionary
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
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
	
	def link_array_id_using_original_filename(self, array_id_filename_ls, curs, array_info_table, ecotype_table, stock_149SNP_db):
		"""
		2008-04-24
			query ecotype's nativename, name, alias, stockparent via the original filename of the array. only valid when the db returns only 1 entry.
		"""
		sys.stderr.write("Linking array id to ecotype id using original filename ... ")
		array_id2ecotype_id_pair = {}
		array_id_filename_left_ls = []
		for array_id, pathname in array_id_filename_ls:
			curs.execute("select id, original_filename from %s where id=%s"%(array_info_table, array_id))
			rows = curs.fetchall()
			original_filename = rows[0][1]
			dir, filename = os.path.split(original_filename)
			filename_pre, filename_ext = os.path.splitext(filename)
			strain_name = filename_pre
			curs.execute("select e.id, e.name, e.stockparent, e.nativename, e.alias, s.name, c.abbr from %s e, %s.site s, %s.address a, %s.country c where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id and (e.nativename rlike '%s'or e.name rlike '%s' or e.alias rlike '%s' or e.stockparent rlike '%s')"%\
						(ecotype_table, stock_149SNP_db, stock_149SNP_db, stock_149SNP_db, strain_name, strain_name, strain_name, strain_name))
			rows = curs.fetchall()
			if len(rows)==1:	#exactly one entry, not less, not more
				ecotype_id = rows[0][0]
				if self.debug:
					print original_filename
					print rows[0]
					yes_or_no = raw_input("confirm?")
					yes_or_no = yes_or_no.lower()
					if yes_or_no=='y' or yes_or_no =='yes':
						array_id2ecotype_id_pair[array_id] = [ecotype_id, ecotype_id]
				else:
					array_id2ecotype_id_pair[array_id] = [ecotype_id, ecotype_id]
			else:
				array_id_filename_left_ls.append([array_id, pathname])
		sys.stderr.write("%s matched and %s left.\n"%(len(array_id2ecotype_id_pair), len(array_id_filename_left_ls)))
		return array_id2ecotype_id_pair, array_id_filename_left_ls
	
	def denovo_link_array_id(self, curs, array_id_filename_ls, ecotype_table, array_id2ecotype_id_pair, stock_149SNP_db):
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
				ecotype_id = self.find_out_ecotypeid_given_strain_name(curs, strain_name, ecotype_table, stock_149SNP_db)
				if ecotype_id!=None:
					ecotype_id_pair.append(ecotype_id)
			if ecotype_id_pair:
				if len(ecotype_id_pair)==1:	#maternal/paternal are same.
					ecotype_id_pair.append(ecotype_id_pair[0])
				array_id2ecotype_id_pair[array_id] = ecotype_id_pair
		sys.stderr.write("Done.\n")
	
	def find_out_ecotypeid_given_strain_name(self, curs, strain_name, ecotype_table='ecotype', stock_149SNP_db='stock'):
		"""
		2008-03-11
			modelled after find_out_ecotypeid_given_strain_name() from Calls2DB_250k.py
		"""
		while 1:
			curs.execute("select e.id, e.name, e.stockparent, e.nativename, e.alias, s.name, c.abbr from %s e, %s.site s, %s.address a, %s.country c where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id and (e.nativename rlike '%s'or e.name rlike '%s' or e.alias rlike '%s' or e.stockparent rlike '%s')"%\
						(ecotype_table, stock_149SNP_db, stock_149SNP_db, stock_149SNP_db, strain_name, strain_name, strain_name, strain_name))
			rows = curs.fetchall()
			choice_id2ecotypeid = {}
			choice_id = 0
			header_ls = ['', 'ecotypeid', 'name', 'stockparent', 'nativename', 'alias', 'site_name', 'country']
			sys.stderr.write("%s\n"%('\t'.join(header_ls)))
			for row in rows:
				ecotypeid, name, stockparent, nativename, alias, site_name, country = row
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
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.user, passwd = self.passwd)
		curs = conn.cursor()
		array_id_filename_ls = self.get_array_id_filename_ls(curs, self.array_info_table)
		if self.calls_comment_table:
			known_filename2ecotype_id = self.get_known_filename2ecotype_id(curs, self.calls_comment_table)
		else:	#no calls_comment_table, empty dictionary
			known_filename2ecotype_id = {}
		#1st round automatic id linking via old mapping
		array_id2ecotype_id_pair, array_id_filename_left_ls = self.link_array_id_using_old_map(array_id_filename_ls, known_filename2ecotype_id)
		#2nd round automatic id linking
		array_id2ecotype_id_pair2, array_id_filename_left_ls = self.link_array_id_using_original_filename(array_id_filename_left_ls, curs, self.array_info_table, self.ecotype_table, self.stock_149SNP_db)
		array_id2ecotype_id_pair.update(array_id2ecotype_id_pair2)
		self.denovo_link_array_id(curs, array_id_filename_left_ls, self.ecotype_table, array_id2ecotype_id_pair, self.stock_149SNP_db)
		if self.commit:
			self.update_array_info_table(curs, self.array_info_table, array_id2ecotype_id_pair)
			curs.execute("commit")

if __name__ == '__main__':
	from pymodule import process_options, generate_program_doc
	main_class = LinkArrayId2EcotypeId
	opts_dict = process_options(sys.argv, main_class.option_default_dict, error_doc=generate_program_doc(sys.argv[0], main_class.option_default_dict)+main_class.__doc__)
	
	instance = main_class(**opts_dict)
	instance.run()
