#!/usr/bin/env python
"""

Examples:
	#test run without commiting database (no records in database in the end)
	./src/Results2DB_250k.py -i /home/nordborglab/pvalue.log -a 1 -e 1 -f kw_test_96_LD -m kw -n best_96_by_tina -u yh
	
	#commit transaction
	./src/Results2DB_250k.py -i /home/nordborglab/pvalue.log -a 1 -e 1 -f kw_test_96_LD -m kw -n best_96_by_tina -u yh -c
	
	#omit short_name
	Results2DB_250k.py -a 17 -e 186 -i /Network/KW_newDataset_186_Bact_titer.pvals -l 1 -u yh -c
	
Description:
	This program would submit simple association results into database.

	The input file format could be 3-column, tab-delimited or 4-column with MAF (minor allele frequency)
		chromosome	position	score/pvalue	MAF
		
	If analysis_method_id is 13 (boolean SNP pair), format is totally different. This program just copies it over.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, stat, getopt, re
import traceback, gc, subprocess
from Stock_250kDB import Stock_250kDB, Results, ResultsMethod, PhenotypeMethod, CallMethod, ResultsMethodType, AnalysisMethod
from Stock_250kDB import Snps as SNPs
#from db import Results, ResultsMethod, Stock_250kDatabase, PhenotypeMethod, CallMethod, SNPs, ResultsMethodType
from pymodule import figureOutDelimiter

import sqlalchemy as sql
"""
2008-04-16 temporarily put here
	-s ...,	short_name*	give a short name of what you did. try to incorporate method, genotype data and phenotype data
	-m ...,	method_description*	a longer description of your method
	-a ...,	data_description*	which data you used
"""

class Results2DB_250k(object):
	__doc__ = __doc__	#use documentation in the beginning of the file as this class's doc
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, '', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname',1, ): [None, 'i', 1, 'File containing association results'],\
							('output_dir',1, ): ['/Network/Data/250k/db/results/', 'o', 1, 'file system storage for the results files. results_method.filename'],\
							('short_name', 0, ): [None, 'f', 1, 'short name for this result. Must be unique from previous ones. combining phenotype, data, method is a good one. If not given, will be automatically generated.' ],\
							('phenotype_method_id',1,int): [None, 'e', 1, 'which phenotype you used, check table phenotype_method'],\
							('call_method_id', 1, int ): [None, 'a', 1, 'data from which call_method, field id in table call_method'],\
							('data_description', 0, ): [None, 'n', 1, 'Describe how your data is derived from that call method. like non-redundant set, 1st 96, etc.'],\
							('method_description', 0, ): [None, 'm', 1, 'Describe your method and what type of score, association (-log or not), recombination etc.'],\
							('results_method_type_id', 1, int): [1, 's', 1, 'which type of method. field id in table results_method_type. 1="association"',],\
							('analysis_method_id', 1, int): [None, 'l', 1, ''],\
							('comment',0, ): [None, 't', 1, 'Anything more worth for other people to know?'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	pa_has_characters = re.compile(r'[a-zA-Z_]')	#2008-05-30 check if a string has character in it, used to judge whether the 1st line is header or not.
	"""
	04/28/08 no longer needed
							('results_table',1, ): 'results',\
							('results_method_table',1, ):'results_method',\
							('phenotype_method_table',1, ):'phenotype_method',\	
	"""
	def __init__(self, **keywords):
		"""
		2008-04-28
			use ProcessOptions, newer option handling class
		2008-04-16
		"""
		
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def check_if_phenotype_method_id_in_db(self, curs, phenotype_method_table, phenotype_method_id):
		"""
		"""
		curs.execute("select id from %s where id=%s"%(phenotype_method_table, phenotype_method_id))
		rows = curs.fetchall()
		if len(rows)>0:
			return 1
		else:
			return 0
	
	def submit_results_method(self, curs, results_method_table, short_name, method_description, data_description):
		"""
		2008-04-16
			submit the method part into db and return the id
		"""
		sys.stderr.write("Submitting results method ...")
		curs.execute("insert into %s(short_name, method_description, data_description) values ('%s', '%s', '%s')"%\
					(results_method_table, short_name, method_description, data_description))
		curs.execute("select id from %s where short_name='%s'"%\
					(results_method_table, short_name))
		rows = curs.fetchall()
		sys.stderr.write("Done.\n")
		return rows[0][0]
		
	marker_pos2snp_id = None
	is_new_marker_added = False	#2008-05-26 flag whether new markers were generated. to check after commit/rollback
	def reset_marker_pos2snp_id(cls):
		"""
		2008-05-27
			set "cls.is_new_marker_added = False"
		2008-05-26
			after commit or rollback in plone, session is closed and those new marker objects are gone. need to reset everything.
		"""
		if cls.is_new_marker_added:
			del cls.marker_pos2snp_id
			cls.marker_pos2snp_id = cls.get_marker_pos2snp_id(db)
			cls.is_new_marker_added = False
	reset_marker_pos2snp_id = classmethod(reset_marker_pos2snp_id)
	
	def get_marker_pos2snp_id(cls, db):
		"""
		2008-05-24
		"""
		sys.stderr.write("Getting marker_pos2snp_id ...")
		marker_pos2snp_id = {}
		snps_table = db.tables['snps'].alias()
		conn = db.connection
		results = conn.execute(sql.select([snps_table.c.id, snps_table.c.chromosome, snps_table.c.position, snps_table.c.end_position]))
		for row in results:
			key = (row.chromosome, row.position, row.end_position)
			marker_pos2snp_id[key] = row.id
		sys.stderr.write("Done.\n")
		return marker_pos2snp_id
	get_marker_pos2snp_id = classmethod(get_marker_pos2snp_id)
	
	def submit_results(cls, db, input_fname, rm, user, output_fname=None):
		"""
		2009-1-7
			insert float into the middle below
				column_5th=int(float(row[4]))	#int('89.0') would raise an exception
		2008-11-12
			parse lines with column_6(genotype_var_perc) and more (comment)
		2008-09-30
			deal with 5-column file. The 5-th column is minor allele count.
			also return True in the end. return False if error in the middle.
		2008-08-19
			add original_filename to ResultsMethod
		2008-07-16
			if input_fname is neither file name nor file object, exit the program
			better handling of the column_4th and its header
		2008-07-16
			if it's 4-column, the last one is MAF.
			can't deal with segment score anymore.
		2008-05-30
			merged with store_file()
				dump the file onto file system storage if output_fname is given
				db submission is too slow
		2008-05-26
			input_fname from plone is not file object although it has file object interface.
		2008-05-26
			csv.Sniffer() can't figure out delimiter if '\n' is in the string, use own dumb function figureOutDelimiter()
		2008-05-25
			save marker(snps) in database if it's not there.
			use marker id in results table
		2008-05-24
			figure out delimiter automatically
			input_fname could be a file object (from plone)
			phenotype method doesn't go with results anymore. it goes with results_method
		2008-04-28
			changed to use Stock_250kDatabase (SQLAlchemy) to do db submission
		"""
		if isinstance(input_fname, str) and os.path.isfile(input_fname):
			sys.stderr.write("Submitting results from %s ..."%(os.path.basename(input_fname)))
			delimiter = figureOutDelimiter(input_fname)
			reader = csv.reader(open(input_fname), delimiter=delimiter)
			rm.original_filename = input_fname
		elif hasattr(input_fname, 'readline') or hasattr(input_fname, 'read'):	#input_fname is not a file name, but direct file object. it could also be <ZPublisher.HTTPRequest.FileUpload instance at 0xa1774f4c>
			sys.stderr.write("Submitting results from %s on plone ..."%input_fname.filename)
			cs = csv.Sniffer()
			input_fname.seek(0)	#it's already read by plone to put int data['input_fname'], check results2db_250k.py
			if getattr(input_fname, 'readline', None) is not None:
				test_line = input_fname.readline()
				delimiter = cs.sniff(test_line).delimiter
			else:
				test_line = input_fname.read(200)
				delimiter = figureOutDelimiter(test_line)	#counting is a safer solution. if test_line include '\n', cs.sniff() won't figure it out.
			input_fname.seek(0)
			reader = csv.reader(input_fname, delimiter=delimiter)
			if getattr(input_fname, 'filename', None):
				rm.original_filename = getattr(input_fname, 'filename', None)
			else:
				rm.original_filename = getattr(input_fname, 'name', None)
		else:
			sys.stderr.write("Error: %s is neither a file name nor a file object.\n"%input_fname)
			sys.exit(4)
		
		if output_fname:
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		elif cls.marker_pos2snp_id is None:
			cls.marker_pos2snp_id = cls.get_marker_pos2snp_id(db)
		
		header_outputted = 0
		no_of_lines = 0
		
		session = db.session
		for row in reader:
			#check if 1st line is header or not
			if no_of_lines ==0 and cls.pa_has_characters.search(row[1]):	#check the 2nd one, which is strict digits. while the 1st column, chromosome could be 'X' or something
				continue
			chr = int(row[0])
			start_pos = int(row[1])
			score = row[2]
			stop_pos = None
			column_4th = None
			column_5th = None
			column_6 = None
			rest_of_row = []
			rest_of_header = []
			
			marker_name = '%s_%s'%(chr, start_pos)
			if len(row)>=4:
				column_4th=row[3]
				#stop_pos = int(row[2])
				#score = row[3]
			if len(row)>=5:
				#column_4th=row[3]
				column_5th=int(float(row[4]))	#2009-1-7 int('89.0') would raise an exception
			if len(row)>=6:
				column_6 = row[5]
			if len(row)>=7:
				rest_of_row = row[6:]
				rest_of_header = ['beta%s'%i for i in range(len(rest_of_row))]
				#sys.stderr.write("ERROR: Found %s columns.\n"%(len(row)))
				#return False
			
			if output_fname:	#go to file system
				if not header_outputted:	#3-column or 4-column header
					if stop_pos is not None:
						position_header = ['start_position', 'stop_position']
					else:
						position_header = ['position']
					header = ['chromosome'] + position_header + ['score']
					if column_4th is not None:
						header.append('MAF')
					if column_5th is not None:
						header.append('MAC')	#Minor Allele Count
					if column_6 is not None:
						header.append('genotype_var_perc')	#genotype variance percentage
					if rest_of_row:
						header += rest_of_header
					writer.writerow(header)
					header_outputted = 1
				data_row = [chr, start_pos]
				if stop_pos is not None:
					data_row.append(stop_pos)
				data_row.append(score)
				if column_4th is not None:
					data_row.append(column_4th)
				if column_5th is not None:
					data_row.append(column_5th)
				if column_6 is not None:
					data_row.append(column_6)
				if rest_of_row:
					data_row += rest_of_row
				writer.writerow(data_row)
			else:
				key = (chr, start_pos, stop_pos)
				if key in cls.marker_pos2snp_id:
					snps_id = cls.marker_pos2snp_id[key]
					if isinstance(snps_id, SNPs):	#it's a new marker object
						r = Results(score=score)
						r.snps = snps_id
					else:	#others are all integer ids
						r = Results(snps_id=snps_id, score=score)
				else:
					#construct a new marker
					marker = SNPs(name=marker_name, chromosome=chr, position=start_pos, end_position=stop_pos, created_by=user)
					#save it in database to get id
					session.save(marker)
					cls.marker_pos2snp_id[key] = marker	#for the next time to encounter same marker
					cls.is_new_marker_added = True	#set this flag as new marker was inputted into the dict
					r = Results(score=score)
					r.snps = marker
					del marker
				r.results_method = rm
				session.save(r)
				del r
			no_of_lines += 1
		
		del reader
		if output_fname:
			del writer
		sys.stderr.write("Done.\n")
		return True
	
	submit_results = classmethod(submit_results)
	
	def come_up_new_results_filename(cls, output_dir, results_method_id, results_method_type_id):
		"""
		2008-05-30
			to save the filename into db
		"""
		output_dir = os.path.join(output_dir, 'type_%s'%results_method_type_id)
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		output_fname = os.path.join(output_dir, '%s_results.tsv'%results_method_id)
		return output_fname
	
	come_up_new_results_filename = classmethod(come_up_new_results_filename)
	
	def copyResultsFile(cls, db, input_fname, rm, user, output_fname=None):
		"""
		2008-09-30
			return True
		2008-09-09
			similar task to submit_results, but not look into the file, just copy the file
		"""
		sys.stderr.write("Copying results from %s ..."%(os.path.basename(input_fname)))
		rm.original_filename = input_fname
		pipe_f = os.popen('cp %s %s'%(input_fname, output_fname))
		pipe_f_out = pipe_f.read()
		if pipe_f_out:
			sys.stderr.write("\tcp output: %s\n"%pipe_f_out)
		sys.stderr.write("Done.\n")
		return True
		
	copyResultsFile = classmethod(copyResultsFile)
	
	def plone_run(cls, db, short_name, phenotype_method_id, call_method_id, data_description, \
				method_description, comment, input_fname, user, results_method_type_id=None, \
				analysis_method_id=None, results_method_type_short_name=None, output_dir=None, commit=0):
		"""
		2008-09-30
			don't save results_method into database if bad thing happend when getting data out of the file.
		2008-09-09
			directly copy the result file if analysis_method_id==13
		2008-08-19
			automatically generate short_name if it's NULL
		2008-07-16
			adjust to new Elixir-based db api.
			new analysis_method_id is added to results_method.
		2008-05-30
			go to output_dir
			drop submit_results()
			use store_file()
		2008-05-26
			add results_method_type_id and results_method_type_short_name
		2008-05-24
			to conveniently wrap up all codes so that both this program and plone can call
		"""
		session = db.session
		session.begin()
		
		if not output_dir:
			output_dir = cls.option_default_dict[('output_dir',1)][0]	#to get default from option_default_dict
		
		rmt = session.query(ResultsMethodType).get(results_method_type_id)
		if not rmt and results_method_type_short_name is not None:	#create a new results method type
			rmt = ResultsMethodType(short_name=results_method_type_short_name)
			session.save(rmt)
		
		if not rmt:
			sys.stderr.write("No results method type available for results_method_type_id=%s.\n"%results_method_type_id)
			sys.exit(3)
		
		pm = PhenotypeMethod.query.get(phenotype_method_id)
		if not pm:
			sys.stderr.write("No phenotype method available for phenotype_method_id=%s.\n"%phenotype_method_id)
			sys.exit(3)
		
		cm = CallMethod.query.get(call_method_id)
		if not cm:
			sys.stderr.write("No call method available for call_method_id=%s.\n"%call_method_id)
			sys.exit(3)
		
		am = AnalysisMethod.query.get(analysis_method_id)
		if not am:
			sys.stderr.write("No analysis method available for analysis_method_id=%s.\n"%analysis_method_id)
			sys.exit(3)
		
		rm = ResultsMethod.query.filter_by(call_method_id=cm.id).filter_by(phenotype_method_id=pm.id).\
			filter_by(analysis_method_id=am.id).filter_by(results_method_type_id=rmt.id)
		if rm.count()>0:
			rm = rm.first()
			sys.stderr.write("There is already an entry in results_method (id=%s) with same (call_method_id, phenotype_method_id, analysis_method_id, results_method_type_id)=(%s, %s, %s, %s).\n"\
							%(rm.id, call_method_id, phenotype_method_id, analysis_method_id, results_method_type_id))
			sys.exit(3)
		
		if not short_name:
			short_name = '%s_%s_%s'%(am.short_name, pm.short_name, cm.id)
		
		rm = ResultsMethod(short_name=short_name, method_description=method_description, \
						data_description=data_description, comment=comment, created_by=user)
		rm.phenotype_method = pm
		rm.call_method = cm
		rm.analysis_method = am
		
		session.save(rm)
		if rmt:
			rm.results_method_type = rmt
		
		#2008-05-30 no submit_results() to database
		#cls.submit_results(db, input_fname, rm, user)
		
		session.flush()	#not necessary as no immediate query on the new results after this and commit() would execute this.
		if commit:
			rm.filename = cls.come_up_new_results_filename(output_dir, rm.id, rm.results_method_type.id)
			if rm.analysis_method_id==13:
				return_value = cls.copyResultsFile(db, input_fname, rm, user, rm.filename)
			else:
				return_value = cls.submit_results(db, input_fname, rm, user, rm.filename)
			if return_value:
				session.save_or_update(rm)
			else:	#bad thing happend when getting data out of the file. don't save this results_method.
				session.delete(rm)
			session.flush()
			session.commit()
			session.clear()
			cls.reset_marker_pos2snp_id()
		else:	#default is also rollback(). to demonstrate good programming
			session.rollback()
			cls.reset_marker_pos2snp_id()
	plone_run = classmethod(plone_run)
	
	def run(self):
		"""
		2008-07-15
			adjust to new Elixir-based db api.
		2008-04-28
			use Stock_250kDatabase to do database stuff
		2008-04-16
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup()
		#db = Stock_250kDatabase(username=self.user,
		#		   password=self.passwd, hostname=self.hostname, database=self.dbname)
		if not os.path.isfile(self.input_fname):
			sys.stderr.write("Error: file, %s,  is not a file.\n"%(self.input_fname))
			sys.exit(3)
		self.plone_run(db, self.short_name, self.phenotype_method_id, self.call_method_id, self.data_description, \
				self.method_description, self.comment, self.input_fname, self.db_user, results_method_type_id=self.results_method_type_id,\
				analysis_method_id=self.analysis_method_id, output_dir=self.output_dir, commit=self.commit)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Results2DB_250k
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()