#!/usr/bin/env python
"""

Examples:
	PutLDIntoDB.py -i  ~/panfs/250k/call_method_17_LD_m0.2.tsv -j 17 -c
	
Description:
	2008-10-15 Program to put LD data into database. (very slow, not sure if there is a better way to store >7 millions entries.)
	
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math
from pymodule import PassingData, getColName2IndexFromHeader
import Stock_250kDB

class PutLDIntoDB(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("input_fname", 1, ): [None, 'i', 1, 'the LD file outputted by MpiLD.py'],\
							('call_method_id', 0, int):[0, 'j', 1, 'which call_method_id is the LD calculated based on.'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-10-15
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def save_LD(self, session, LD_fname, call_method_id, commit=0):
		"""
		2008-10-15
			adapted from DrawSNPRegion.get_LD() 
		"""
		sys.stderr.write("Reading in LD info from %s ...\n"%(LD_fname))
		reader = csv.reader(open(LD_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		for row in reader:
			snp1 = row[col_name2index['snp1']].split('_')
			snp1 = map(int, snp1)
			snp2 = row[col_name2index['snp2']].split('_')
			snp2 = map(int, snp2)
			allele1_freq = float(row[col_name2index['allele1_freq']])
			allele2_freq = float(row[col_name2index['allele2_freq']])
			r2 = float(row[col_name2index['r2']])
			D_prime = float(row[col_name2index['D_prime']])
			D = float(row[col_name2index['D']])
			if snp1<snp2:
				snp_pair = (snp1[0], snp1[1], snp2[0], snp2[1])
			else:
				snp_pair = (snp2[0], snp2[1], snp1[0], snp1[1])
			no_of_pairs = int(float(row[col_name2index['no_of_pairs']])/2)	#MpiLD.py outputs this double (due to haploid regarded as diploid)
			ld = Stock_250kDB.LD(snp1_maf=allele1_freq, snp2_maf=allele2_freq, d=D, d_prime=D_prime, r2=r2, no_of_pairs=no_of_pairs)
			ld.chr1 = snp_pair[0]
			ld.pos1 = snp_pair[1]
			ld.chr2 = snp_pair[2]
			ld.pos2 = snp_pair[3]
			ld.call_method_id = call_method_id
			if commit:
				session.save(ld)
				session.flush()
			counter += 1
			if counter%100000==0:
				sys.stderr.write('%s\t%s'%('\x08'*100, counter))
			if counter%1000==0 and self.debug>0:
				break
				pass
		sys.stderr.write("%s entries. Done.\n"%counter)
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		#session.begin()	#transaction hoards huge memory in the end!
		
		self.save_LD(session, self.input_fname, self.call_method_id, self.commit)
		
		"""
		if self.commit:
			session.flush()
			session.commit()
			session.clear()
		else:	#default is also rollback(). to demonstrate good programming
			session.rollback()
		"""

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PutLDIntoDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()