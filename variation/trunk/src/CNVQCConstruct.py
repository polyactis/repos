#!/usr/bin/env python
"""

Examples:
	# Clark et al. Science 2007's data
	CNVQCConstruct.py -i ~/script/variation/data/CNV/PERL_deletions_in_PRPs.txt -s Clark2007a -t deletion -m PCR -u yh -n 1 -c
	
	# Korbinian Schneeberger & Stephan Ossowski Paired-end
	CNVQCConstruct.py -i ~/script/variation/data/CNV/PERL_deletions_in_PRPs.txt -s SchneebergerOssowski2009 -t deletion -m PairedEndSolexa -u yh -y 2 -c
	
	# Bob schmitz data
	CNVQCConstruct.py -i ~/script/variation/data/CNV/cvi_chr1_indels -s BobSchmitz2009 -t deletion -m PairedEndSolexa -o Cvi-0 -u yh -y 3 -c
	
Description:
	2009-10-28 Fill up relevant db tables with CNV QC data from different sources.
	
	Argument -t "cnv_type" is useless for run_type 3. it depends on whether size_affected is positive or negative.
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv, numpy, traceback, subprocess
from pymodule import figureOutDelimiter, getListOutOfStr
import Stock_250kDB

class CNVQCConstruct(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname', 1, ): ['', 'i', 1, 'input file.', ],\
							('data_source', 1, ): ['', 's', 1, 'a short name for data source'],\
							('cnv_type', 1, ): ['', 't', 1, 'which type of CNV? insertion, deletion, etc. in table cnv_type or will be inserted.'],\
							('cnv_method', 0, ): ['', 'm', 1, 'method used to derive CNV. If not provided, it is ignored.', ],\
							('no_of_lines_to_skip', 1, int): [0, 'n', 1, 'Number of lines to skip in the beginning of the input file.', ],\
							('original_id', 0, ): ['', 'o', 1, 'only used for run_type 3. cuz data is from one accession and not specifying which one.', ],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.'],\
							('run_type', 0, int):[1, 'y', 1, 'Run type 1: Clark et al. 2007, 2: Korbinian Schneeberger & Stephan Ossowski Paired-end, 3: Bob Schmitz']}
	
	def __init__(self, **keywords):
		"""
		2009-10-28
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def getDBObj(self, session, table_class, short_name):
		"""
		2009-10-28
		"""
		sys.stderr.write("Getting/Creating a %s object ..."%(table_class.table.name))
		db_obj = table_class.query.filter_by(short_name=short_name).first()
		if not db_obj:
			db_obj = table_class(short_name=short_name)
			session.save(db_obj)
			session.flush()
		sys.stderr.write("Done.\n")
		return db_obj
	
	def getCNVQCAccessionObj(self, session, original_id, data_source_obj):
		"""
		2009-10-28
		"""
		db_obj = Stock_250kDB.CNVQCAccession.query.filter_by(original_id=original_id).filter_by(data_source_id=data_source_obj.id).first()
		if not db_obj:
			db_obj = Stock_250kDB.CNVQCAccession(original_id=original_id)
			db_obj.data_source = data_source_obj
			session.save(db_obj)
			session.flush()
		return db_obj
	
	def generateCNVQCCallObjType1(self, session, row, data_source_obj, cnv_type_obj, cnv_method_obj=None, original_id=None):
		"""
		2009-10-29
			run_type 1 (clark et al 2007) format:
				Accession
				Chromosome
				Start deletion position
				End deletion position
				Deletion length (bp)
				Forward primer for validation
				Reverse primer for validation
		"""
		original_id, chromosome, start, stop, size_affected = row[:5]
		score = None
		chromosome = int(chromosome)
		start = int(start)
		stop = int(stop)
		size_affected = int(size_affected)
		acc_key = (original_id, data_source_obj.short_name)
		if acc_key in self.original_id_data_source2acc_obj:
			acc_obj = self.original_id_data_source2acc_obj.get(acc_key)
		else:
			acc_obj = self.getCNVQCAccessionObj(session, original_id, data_source_obj)
			self.original_id_data_source2acc_obj[acc_key] = acc_obj
		cnv_qc_call = Stock_250kDB.CNVQCCalls(chromosome=chromosome, start=start, stop=stop, \
											size_affected=size_affected, score=score)
		cnv_qc_call.accession = acc_obj
		cnv_qc_call.cnv_type = cnv_type_obj
		cnv_qc_call.cnv_method = cnv_method_obj
		return cnv_qc_call
	
	def generateCNVQCCallObjType2(self, session, row, data_source_obj, cnv_type_obj, cnv_method_obj=None, original_id=None):
		"""
		2009-10-29
			run_type 2 (Schneeberger & Ossowski) format:
				SampleID
				Number of read pairs supporting deletion
				chromosome
				begin (it is not the beginning of the deletion but the beginning of a region in which the deletion is located)
				end (dito)
				length (as supported by the variation of the distances of the alignments)
				length (as supported by the reference)
				number of positions without core-alignments between begin and end (only looking at non-repetitive)
				number of positions without core-alignments between begin and end
				p-value (likeliness that the given insert size distribution generated the stretched mate pairs)
		"""
		original_id, no_of_read_pairs, chromosome, start, stop, length_supported_by_var_of_dist_of_align, size_affected = row[:7]
		score = float(row[-1])	# Data from Korbinian and Stephan has pvalues.
		chromosome = int(chromosome)
		start = int(start)
		stop = int(stop)
		size_affected = int(size_affected)
		acc_key = (original_id, data_source_obj.short_name)
		if acc_key in self.original_id_data_source2acc_obj:
			acc_obj = self.original_id_data_source2acc_obj.get(acc_key)
		else:
			acc_obj = self.getCNVQCAccessionObj(session, original_id, data_source_obj)
			self.original_id_data_source2acc_obj[acc_key] = acc_obj
		
		comment = '#read pairs: %s, length by dist var: %s, #positions with no non-repeat-reads: %s, #positions with no reads: %s'%\
				(no_of_read_pairs, length_supported_by_var_of_dist_of_align, row[7], row[8])
		cnv_qc_call = Stock_250kDB.CNVQCCalls(chromosome=chromosome, start=start, stop=stop, \
											size_affected=size_affected, score=score, comment=comment)
		cnv_qc_call.accession = acc_obj
		cnv_qc_call.cnv_type = cnv_type_obj
		cnv_qc_call.cnv_method = cnv_method_obj
		return cnv_qc_call
	
	cnv_type2cnv_type_obj = {}
	def generateCNVQCCallObjFromBobSchmitzData(self, session, row, data_source_obj, cnv_type_obj, cnv_method_obj=None, original_id=None):
		"""
		2009-10-29
			cnv_type_obj is useless. it depends on whether size_affected is positive or negative.
			
		"""
		chromosome, start_stop, contig_id, contig_start_stop, size_affected = row[:5]
		chromosome = int(chromosome)
		start, stop = start_stop.split(' - ')
		start = int(start)
		stop = int(stop)
		size_affected = int(size_affected)
		if size_affected<0:
			cnv_type = "deletion"
		elif size_affected>0:
			cnv_type = "insertion"
		else:
			cnv_type = "normal"
		if cnv_type in self.cnv_type2cnv_type_obj:
			cnv_type_obj = self.cnv_type2cnv_type_obj.get(cnv_type)
		else:
			cnv_type_obj = self.getDBObj(session, Stock_250kDB.CNVType, cnv_type)
			self.cnv_type2cnv_type_obj[cnv_type] = cnv_type_obj
		
		comment = '%s %s'%(contig_id, contig_start_stop)
		acc_key = (original_id, data_source_obj.short_name)
		if acc_key in self.original_id_data_source2acc_obj:
			acc_obj = self.original_id_data_source2acc_obj.get(acc_key)
		else:
			acc_obj = self.getCNVQCAccessionObj(session, original_id, data_source_obj)
			self.original_id_data_source2acc_obj[acc_key] = acc_obj
		cnv_qc_call = Stock_250kDB.CNVQCCalls(chromosome=chromosome, start=start, stop=stop, \
											size_affected=abs(size_affected), comment=comment)
		cnv_qc_call.accession = acc_obj
		cnv_qc_call.cnv_type = cnv_type_obj
		cnv_qc_call.cnv_method = cnv_method_obj
		return cnv_qc_call
	
	original_id_data_source2acc_obj = {}
	def putQCIntoDB(self, session, input_fname, no_of_lines_to_skip, data_source_obj, cnv_type_obj, cnv_method_obj=None, \
				run_type=1, original_id=None):
		"""
		2009-10-28
		"""
		sys.stderr.write("Putting QC data into database ... \n")
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		for i in range(no_of_lines_to_skip):
			reader.next()
		
		counter = 0
		for row in reader:
			if run_type ==1:
				cnv_qc_call = self.generateCNVQCCallObjType1(session, row, data_source_obj, cnv_type_obj, cnv_method_obj)
			elif run_type ==2:
				cnv_qc_call = self.generateCNVQCCallObjType2(session, row, data_source_obj, cnv_type_obj, cnv_method_obj)
			elif run_type==3:
				cnv_qc_call = self.generateCNVQCCallObjFromBobSchmitzData(session, row, data_source_obj, cnv_type_obj,\
																		cnv_method_obj, original_id=original_id)
			else:
				sys.stderr.write("Run type %s not supported.\n"%run_type)
			session.save(cnv_qc_call)
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
		session.flush()
		sys.stderr.write("%s records. Done.\n"%counter)
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		if self.run_type==3 and  not self.original_id:
			sys.stderr.write("Need original_id for run_type %s.\n"%self.run_type)
			sys.exit(2)
		
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				  			password=self.db_passwd, hostname=self.hostname, database=self.dbname, 
				   			schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		session.begin()
		
		data_source_obj = self.getDBObj(session, Stock_250kDB.DataSource, self.data_source)
		cnv_type_obj = self.getDBObj(session, Stock_250kDB.CNVType, self.cnv_type)
		if self.cnv_method:
			cnv_method_obj = self.getDBObj(session, Stock_250kDB.CNVMethod, self.cnv_method)
		else:
			cnv_method_obj = None
		self.putQCIntoDB(session, self.input_fname, self.no_of_lines_to_skip, data_source_obj, cnv_type_obj, \
						cnv_method_obj=cnv_method_obj, run_type=self.run_type, original_id=self.original_id)
		
		if self.commit:
			session.commit()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CNVQCConstruct
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()