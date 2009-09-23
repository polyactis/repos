#!/usr/bin/env python
"""

Examples:
	Calls_BySeq2Calls.py -u yh -c

Description:
	program to split table calls_byseq in db stock into table strain and calls.
	2009-9-22 Tweaked to be able to deal with incremental changes in calls_byseq.
	
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
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode, test-run 5000 entries from calls_byseq.'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	pa_has_characters = re.compile(r'[a-zA-Z_]')
	def __init__(self, **keywords):
		"""
		2008-08-11
		"""
		
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def return_seqinfoid_name_given_snpid(self, snpid):
		"""
		2008-09-12
			the 149 SNPs are processed in 4 plates by sequenom. The order is preserved in table snps.
			snps.id from 1-38 is plate 1. 39-75 is plate 2. 76-112 is plate 3. 113-149 is plate 4.
		"""
		if snpid>=1 and snpid<39:
			return 'seqinfoid1'
		elif snpid>=39 and snpid<76:
			return 'seqinfoid2'
		elif snpid>=76 and snpid<113:
			return 'seqinfoid3'
		elif snpid>=113 and snpid<=149:
			return 'seqinfoid4'
	
	def get_strain_unique_key2strain(self, session):
		"""
		2009-9-22
			get strain_unique_key2strain from db for usage in splitCalls_BySeq()
		"""
		sys.stderr.write("Getting strain_unique_key2strain ... ")
		strain_unique_key2strain = {}
		rows = Strain.query()
		counter = 0
		for row in rows:
			strain_unique_key = (row.ecotypeid, row.plateid, row.wellid, row.replicate)
			strain_unique_key2strain[strain_unique_key] = row
		sys.stderr.write("%s strains. Done.\n"%(len(strain_unique_key2strain)))
		return strain_unique_key2strain
	
	def splitCalls_BySeq(self, session, strain_unique_key2strain={}):
		"""
		2009-9-22
			add argument strain_unique_key2strain
			identify a strain by (ecotypeid, plateid, wellid, replicate)
			check if (strainid, snpid) is already in table Calls
		2008-09-12
			assign 4 seqinfoids if they all exist
		2008-08-11
			when checking out Calls_BySeq, order by ecotypeid, plateid
		2008-08-11
		"""
		sys.stderr.write("Splitting calls_byseq ....\n")
		offset_index = 0
		block_size = 5000
		rows = Calls_BySeq.query.order_by(Calls_BySeq.ecotypeid).offset(offset_index).limit(block_size)
		counter = 0
		while rows.count()!=0:
			for row in rows:
				#search it first
				strain_unique_key = (row.ecotypeid, row.plateid, row.wellid, row.replicate)
				if strain_unique_key in strain_unique_key2strain:
					strain = strain_unique_key2strain[strain_unique_key]
				else:
					strain = Strain(ecotypeid=row.ecotypeid, plateid=row.plateid, extractionid=row.extractionid, \
								wellid=row.wellid, replicate=row.replicate)
					strain_unique_key2strain[strain_unique_key] = strain	#add this to strain_unique_key2strain
				seqinfoid_name = self.return_seqinfoid_name_given_snpid(row.snpid)
				seqinfoid = getattr(strain, seqinfoid_name, None)
				if seqinfoid is None:
					setattr(strain, seqinfoid_name, row.seqinfoid)
				else:
					if seqinfoid!=row.seqinfoid:
						sys.stderr.write("Strain with ecotypeid=%s, plateid=%s, wellid=%s has two different seqinfoids in this SNP block: %s, %s(snpid=%s).\n"%\
										(strain.ecotypeid, strain.plateid, strain.wellid, seqinfoid, row.seqinfoid, row.snpid))
				session.save_or_update(strain)
				"""
				strain_in_db = Strain.query.filter_by(ecotypeid=row.ecotypeid).filter_by(plateid=row.plateid).filter_by(wellid=row.wellid).filter_by(extractionid=row.extractionid)
				if strain_in_db.count()>1:
					sys.stderr.write("Warning: has more than 1 strain with ecotypeid=%s, plateid=%s, wellid=%s, extractionid=%s.\n"%(row.ecotypeid, row.plateid, row.wellid, row.extractionid))
				
				strain = strain_in_db.first()
				if not strain:
				"""
				offset_index += 1
				if strain.id is not None:
					calls_in_db = Calls.query.filter_by(strainid=strain.id).filter_by(snpid=row.snpid)
					if calls_in_db.count()>0:	# in db, skip
						continue
				
				call = row.call1
				if call:
					call = call.upper()
				if row.call2:
					call2 = row.call2.upper()
					call = call+call2
				call = number2nt[nt2number[call]]
				calls = Calls(snpid=row.snpid, allele=call)
				calls.strain = strain
				session.save_or_update(calls)
				counter += 1
			sys.stderr.write("%s%s out of %s moved"%('\x08'*60, counter, offset_index))
			#if self.debug and offset_index > 1000:
			#	break
			session.flush()	#flush now and then to avoid possible memory error due to too many objects residing in memory.
			rows = Calls_BySeq.query.order_by(Calls_BySeq.ecotypeid).offset(offset_index).limit(block_size)
		"""
		#no need to explicityly save/flush them, it'll be all handled during 'commit' stage.
		#save strains in ecotypeid order
		strain_unique_key_ls = strain_unique_key2strain.keys()
		strain_unique_key_ls.sort()
		for strain_unique_key in strain_unique_key_ls:
			session.save(strain_unique_key2strain[strain_unique_key])
			session.flush()
		#save calls
		for calls in calls_ls:
			session.save(calls)
		"""
		sys.stderr.write("Done.\n")
		
	def run(self):
		"""
		2009-9-22
			add the debug section
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db = StockDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup()
		session = db.session
		session.begin()
		strain_unique_key2strain = self.get_strain_unique_key2strain(session)
		self.splitCalls_BySeq(session, strain_unique_key2strain)
		
		if self.commit:
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