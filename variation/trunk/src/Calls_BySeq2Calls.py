#!/usr/bin/env python
"""

Examples:
	Calls_BySeq2Calls.py -u yh -c

Description:
	program to split table calls_byseq in db stock into table strain and calls.
	
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
from StockDB import StockDB, SNPs, Strain, Calls, Calls_BySeq
from pymodule import figureOutDelimiter, nt2number, number2nt

class Calls_BySeq2Calls(object):
	__doc__ = __doc__	#use documentation in the beginning of the file as this class's doc
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock', 'd', 1, '', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	pa_has_characters = re.compile(r'[a-zA-Z_]')
	def __init__(self, **keywords):
		"""
		2008-08-11
		"""
		
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def splitCalls_BySeq(self, session):
		"""
		2008-08-11
		"""
		sys.stderr.write("Splitting calls_byseq ....\n")
		offset_index = 0
		block_size = 5000
		rows = Calls_BySeq.query.offset(offset_index).limit(block_size)
		counter = 0
		while rows.count()!=0:
			for row in rows:
				#search it first
				strain = Strain.query.filter_by(ecotypeid=row.ecotypeid).filter_by(plateid=row.plateid).first()
				if not strain:
					strain = Strain(ecotypeid=row.ecotypeid, plateid=row.plateid, extractionid=row.extractionid, \
								seqinfoid=row.seqinfoid, wellid=row.wellid, replicate=row.replicate)
					session.save(strain)
					session.flush()
				call = row.call1
				if call:
					call = call.upper()
				if row.call2:
					call2 = row.call2.upper()
					call = call+call2
				call = number2nt[nt2number[call]]
				calls = Calls(snpid=row.snpid, allele=call)
				calls.strain = strain
				session.save(calls)
				counter += 1
				offset_index += 1
			sys.stderr.write("%s%s\t%s"%('\x08'*40, offset_index, counter))
			if self.debug and offset_index > 1000:
				break
			rows = Calls_BySeq.query.offset(offset_index).limit(block_size)
			session.flush()
		sys.stderr.write("Done.\n")
		
	def run(self):
		"""
		"""
		db = StockDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		session = db.session
		if not self.debug:
			session.begin()
		self.splitCalls_BySeq(session)
		
		if not self.debug:
			if self.commit:
				session.flush()
				session.commit()
				session.clear()
			else:	#default is also rollback(). to demonstrate good programming
				session.rollback()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Calls_BySeq2Calls
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()