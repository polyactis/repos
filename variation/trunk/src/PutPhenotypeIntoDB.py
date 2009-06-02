#!/usr/bin/env python
"""

Examples:
	PutPhenotypeIntoDB.py -i /Network/Data/250k/anthocyanin.csv -A 4 -u yh -c
	
	PutPhenotypeIntoDB.py -i /Network/Data/250k/hyaloperonaspora2.csv -s Hyaloperonaspora2 -u yh -c
	
	#2009-5-29 put raw lesioning phenotype into db.
	PutPhenotypeIntoDB.py -i ./acc2lesioning.csv -s LES_raw -u yh  -t -c
	
Description:
	2009-5-28
		Put phenotype data into table phenotype_avg, phenotype_method.
		Input format is two column (either delimiter), accession name and phenotype value. If accession name cannot be uniquely or not at all
			linked to an ecotype id. Program would stop and ask user.

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback

from Stock_250kDB import Stock_250kDB, PhenotypeAvg, PhenotypeMethod
from pymodule import figureOutDelimiter
from pymodule.utils import getGeneIDSetGivenAccVer
from common import getNativename2TgEcotypeIDSet

class PutPhenotypeIntoDB(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("input_fname", 1, ): [None, 'i', 1, ''],\
							("phenotype_id", 0, int): [-1, 'A', 1, 'ID of phenotype method.'],\
							("phenotype_name", 0, ): [None, 's', 1, 'phenotype_method short name. if given, phenotype_id would be ignored.'],\
							('skip_1st_line', 0, int):[0, 't', 0, 'skip the first line in the input file'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2009-5-28
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
		if self.phenotype_id==-1 and self.phenotype_name is None:
			sys.stderr.write("Error: None of phenotype_id and phenotype_name is specified.\n")
			sys.exit(3)
		
	def putPhenotypeIntoDB(self, input_fname, phenotype_id, phenotype_name, nativename2tg_ecotypeid_set, db, skip_1st_line=False):
		"""
		2009-5-28
		"""
		import csv, sys, os
		session = db.session
		if phenotype_name:	#if the short name is given, forget about phenotype_id
			pm = PhenotypeMethod.query.filter_by(short_name=phenotype_name).first()	#try search the db first.
			if not pm:
				pm = PhenotypeMethod(short_name=phenotype_name)
				session.save(pm)
				session.flush()
		else:	#use the list_type_id to get it
			pm = PhenotypeMethod.get(phenotype_id)
		session.save_or_update(pm)
		
		
		delimiter=figureOutDelimiter(input_fname)
		
		reader = csv.reader(open(input_fname), delimiter=delimiter)
		if skip_1st_line:
			reader.next()	#skips the 1st line
		counter = 0
		success_counter = 0
		for row in reader:
			counter += 1
			if not row:	#skip empty lines
				continue
			original_name = row[0].strip()	#2008-12-11 remove spaces/tabs in the beginning/end
			if original_name in nativename2tg_ecotypeid_set:
				tg_ecotypeid_set = nativename2tg_ecotypeid_set.get(original_name)
				if len(tg_ecotypeid_set)>1:
					while 1:
						ecotype_id = raw_input("Pick one ecotype id from "+repr(list(tg_ecotypeid_set)) + " for %s: "%original_name)
						yes_or_no = raw_input("Sure about ecotype id %s?(Y/n)"%ecotype_id)
						if yes_or_no =='n' or yes_or_no=='N' or yes_or_no == 'No' or yes_or_no =='no':
							pass
						else:
							break
						
				else:
					ecotype_id = list(tg_ecotypeid_set)[0]
			else:
				
				choice_instructions = "Enter ecotypeid for accession %s: "%original_name
				while 1:
					ecotype_id = raw_input(choice_instructions)
					yes_or_no = raw_input("Sure about ecotype id %s?(y/N)"%ecotype_id)
					if yes_or_no =='n' or yes_or_no=='N' or yes_or_no == 'No' or yes_or_no =='no':
						pass
					else:
						break
			ecotype_id = int(ecotype_id)
			phenotype_value = float(row[1])
			pa = PhenotypeAvg(ecotype_id=ecotype_id, value=phenotype_value)
			pa.phenotype_method = pm
			session.save(pa)
			success_counter += 1
		del reader
		sys.stderr.write("%s/%s linked successfully.\n"%(success_counter, counter))
		

	def run(self):
		"""
		2009-5-28
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		
		nativename2tg_ecotypeid_set = getNativename2TgEcotypeIDSet(db.metadata.bind)
		
		session = db.session
		session.begin()
		self.putPhenotypeIntoDB(self.input_fname, self.phenotype_id, self.phenotype_name, nativename2tg_ecotypeid_set, db, self.skip_1st_line)
		if self.commit:
			session.commit()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PutPhenotypeIntoDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()