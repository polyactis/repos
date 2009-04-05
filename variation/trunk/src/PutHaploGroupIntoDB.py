#!/usr/bin/env python
"""

Examples:
	PutHaploGroupIntoDB.py -i /tmp/cnd.fN.data.hg.trust2.proposal.filtered2 -u yh -r -c
	
Description:
	program to put haplotype groups constructed by alex into database.
	
	Input file is tab or coma-delimited with the 1st line as the header. Alternatively saved from Alex's file.
	
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
import StockDB
from pymodule import figureOutDelimiter, nt2number, number2nt, getColName2IndexFromHeader
from common import get_ecotypeid2tg_ecotypeid

class PutHaploGroupIntoDB(object):
	__doc__ = __doc__	#use documentation in the beginning of the file as this class's doc
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock', 'd', 1, '', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname', 1, ):[None, 'i', 1, 'the spreadsheet from alex'],\
							('max_snp_typing_error_rate', 1, float): [0.005, '', 1, 'maximum SNP-typing error rate allowed in constructing haplotype groups'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode, test-run 5000 entries from calls_byseq.'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	pa_has_characters = re.compile(r'[a-zA-Z_]')
	def __init__(self, **keywords):
		"""
		2009-3-31
		"""
		
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def getSNPIDLs(self):
		"""
		2009-3-31
		"""
		sys.stderr.write("Getting snp_id_ls ...")
		rows = StockDB.SNPs.query.order_by(StockDB.SNPs.chromosome).order_by(StockDB.SNPs.position).all()
		snp_id_ls = []
		for row in rows:
			snp_id_ls.append(row.id)
		sys.stderr.write("Done.\n")
		return snp_id_ls
	
	def putHaplotypeGroupIntoDB(self, session, input_fname, tg_ecotypeid2row, max_snp_typing_error_rate, snp_id_ls):
		"""
		2009-3-31
		2009-4-4
			add argument tg_ecotypeid2row
		"""
		sys.stderr.write("Constructing haplotype groups ...\n")
		pattern_ecotypeid = re.compile(r'(?<=\))\d+')
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		col_name2col_index = getColName2IndexFromHeader(reader.next())
		ecotypeid_idx = col_name2col_index['ecotypeid']
		haplo_name_idx = col_name2col_index['haplogroup']
		geographic_integrity_idx = col_name2col_index['geographic_integrity']
		filtered_SNPs_idx = col_name2col_index['filtered_SNPs']
		counter = 0
		for tg_ecotypeid, row in tg_ecotypeid2row.iteritems():
			ecotypeid = int(row[ecotypeid_idx])
			ecotypeid = tg_ecotypeid	#2009-4-4 use tg_ecotypeid instead
			haplo_name = row[haplo_name_idx]
			geographic_integrity_name = row[geographic_integrity_idx]
			filtered_SNPs = row[filtered_SNPs_idx]
			ref_ecotypeid = int(pattern_ecotypeid.search(haplo_name).group(0))
			haplo_group = StockDB.HaploGroup.query.filter_by(short_name=haplo_name).first()
			if not haplo_group:
				haplo_group = StockDB.HaploGroup(short_name=haplo_name, ref_ecotypeid=ref_ecotypeid, max_snp_typing_error_rate=max_snp_typing_error_rate)
				session.save(haplo_group)
				session.flush()
			
			ecotype = StockDB.Ecotype.get(ecotypeid)
			haplo_group.ecotypes.append(ecotype)
			geographic_integrity = StockDB.GeographicIntegrity.query.filter_by(short_name=geographic_integrity_name).first()
			if not geographic_integrity:
				geographic_integrity = StockDB.GeographicIntegrity(short_name=geographic_integrity_name)
				session.save(geographic_integrity)
				session.flush()
			ecotype.geographic_integrity = geographic_integrity
			session.save_or_update(ecotype)
			#one bit of ecotype: link the ecotypeid to tg_ecotype_id
			
			
			#deal with filtered SNPs
			for i in range(len(filtered_SNPs)):
				allele = filtered_SNPs[i]
				if allele=='_':
					continue
				fc = StockDB.FilteredCalls(ecotypeid=ecotypeid, snpid=snp_id_ls[i], allele=allele)
				session.save(fc)
				session.flush()
			counter += 1
			if counter%500==0 and self.report:
				sys.stderr.write('%s%s'%('\x08'*80, counter))
		session.flush()
		sys.stderr.write("Done.\n")
	
	def dropRedundantEcotypes(self, input_fname, ecotypeid2tg_ecotypeid):
		"""
		2009-4-4
			retain only one row out of duplicated ecotype rows based on ecotypeid2tg_ecotypeid.
				it's not random. usually the one with same ecotype id as tg_ecotypeid unless tg_ecotypeid doesn't appear.
			if duplicated ecotypes belong to different haplotype group, choose the one with tg_ecotypeid otherwise random.
		"""
		sys.stderr.write("Dropping redundant ecotypes ...\n")
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		col_name2col_index = getColName2IndexFromHeader(reader.next())
		ecotypeid_idx = col_name2col_index['ecotypeid']
		haplo_name_idx = col_name2col_index['haplogroup']
		nativename_idx = col_name2col_index['nativename']
		tg_ecotypeid2row = {}
		no_of_duplicates = 0
		no_of_duplicates_with_different_haplogroups = 0
		counter = 0
		for row in reader:
			ecotypeid = int(row[ecotypeid_idx])
			haplo_name = row[haplo_name_idx]
			nativename = row[nativename_idx]
			if ecotypeid in ecotypeid2tg_ecotypeid:
				tg_ecotypeid = ecotypeid2tg_ecotypeid[ecotypeid]
				if tg_ecotypeid not in tg_ecotypeid2row:
					tg_ecotypeid2row[tg_ecotypeid] = row
				else:
					no_of_duplicates += 1
					old_row = tg_ecotypeid2row[tg_ecotypeid]
					old_ecotypeid = int(old_row[ecotypeid_idx])
					old_haplo_name = old_row[haplo_name_idx]
					old_nativename = row[nativename_idx]
					if old_haplo_name!=haplo_name:
						sys.stderr.write("ecotype %s(%s) in haplotype group %s, while duplicate %s(%s) in haplotype group %s.\n"%\
										 (ecotypeid, nativename, haplo_name, old_ecotypeid, old_nativename, old_haplo_name))
						no_of_duplicates_with_different_haplogroups += 1
					if ecotypeid==tg_ecotypeid:	#replace if the new ecotypeid matching the tg_ecotypeid whether the haplotype group is same or not.
						tg_ecotypeid2row[tg_ecotypeid] = row
			else:
				sys.stderr.write("Warning: ecotype %s not in ecotypeid2tg_ecotypeid.\n"%(ecotypeid))
			counter += 1
		sys.stderr.write("no_of_duplicates: %s, out of which %s encompass different haplotype groups. %s accessions in total. Done.\n"%\
						 (no_of_duplicates, no_of_duplicates_with_different_haplogroups, counter))
		return tg_ecotypeid2row
	
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
			
		db = StockDB.StockDB(drivername=self.drivername, username=self.db_user,
				 		password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup()
		session = db.session
		session.begin()
		ecotypeid2tg_ecotypeid = get_ecotypeid2tg_ecotypeid(db.metadata.bind, debug=self.debug)
		tg_ecotypeid2row = self.dropRedundantEcotypes(self.input_fname, ecotypeid2tg_ecotypeid)
		snp_id_ls = self.getSNPIDLs()
		self.putHaplotypeGroupIntoDB(session, self.input_fname, tg_ecotypeid2row, self.max_snp_typing_error_rate, snp_id_ls)
		
		if self.commit:
			session.commit()
			session.clear()
		else:	#default is rollback(). to demonstrate good programming
			session.rollback()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PutHaploGroupIntoDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()