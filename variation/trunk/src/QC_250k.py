#!/usr/bin/env python
"""

Examples:
	#test run (without -c) QC between 250k and 2010
	QC_250k.py -i /home/crocea/script/variation/data/2010/data_2010_x_250k_y0001.tsv -m 1
	
	#test run (without -c) QC between 250k and perlegen
	QC_250k.py -i /home/crocea/script/variation/data/perlegen/data_perlegen_ecotype_id_x_250k_y0101.tsv -m 2
	
	#test run (without -c) QC between 250k and 149SNP
	QC_250k.py -i /home/crocea/script/variation/stock20080403/data_y10001101.tsv -m 3
	
	#QC for the call files (input for Calls2DB_250k.py) against 149SNP
	QC_250k.py -i /home/crocea/script/variation/stock20080403/data_y10001101.tsv -m 3 -n /Network/Data/250k/tmp-jr/calls-oligo/ -o /tmp/tmp-jr_calls-oligo_149SNP_QC.tsv
	
	QC_250k.py -i /home/crocea/script/variation/data/2010/data_2010_x_250k_y0001.tsv -m 5 -l 3 -y 0.85
	QC_250k.py -i /home/crocea/script/variation/data/perlegen/data_perlegen_ecotype_id_x_250k_y0101.tsv -m 6 -l 3 -y 0.85 -c
	QC_250k.py -i /home/crocea/script/variation/stock20080403/data_y00001101.tsv -m 7 -l 3 -y 0.85 -c
	
	#calculate missing rate for call_info entries (call_method_id=3), minimum calling probability=0.85
	QC_250k.py -m 0 -l 3 -c -y 0.85
	
	#accession-wise QC between 250k (call_method_id=3) and 149SNP
	QC_250k.py -m 3 -l 3 -y 0.85 -c
	#accession-wise QC between 250k (call_method_id=3) and 384-illumina, one by one
	QC_250k.py -m 8 -l 3 -y 0.85 -c -a
	
	#do snp-wise QC between 250k (call_method_id=3, excluding arrays with > 20% mismatch rate, according to the QC with maximum no_of_non_NA_pairs) and Perlegen
	QC_250k.py -m 2 -l 3 -y 0.85 -e 2 -x 0.20 -c
	
	#QC on an NPUTE output file.
	QC_250k.py -n /mnt/nfs/NPUTE_data/250k_l3_y0..6_w0.2_x0.2_h170_w10.npute -m 3
	QC_250k.py -n /mnt/nfs/NPUTE_data/NPUTE_output/250k_l3_y0.70_v0.6_w0.2_x0.2.tsv_w10.npute -m 8 -o /tmp/250k_l3_y0.70_v0.6_w0.2_x0.2.tsv_w10.npute_m8_QC.tsv -l 3
	
Description:
	QC for 250k call data from call_info_table against 2010, perlegen, 149SNP data.
	It will select the call_info entries that haven't been QCed for a particular QC method.
	The output is in the call_QC_table with NA_rate, mismatch_rate and etc.
	
	QC_method_id:
		All but 0 correspond to field/column id in table QC_method.
		0 is not a real QC method id. It calculates NA rate for call_info entries whose NA_rates haven't been calculated.
	

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
from pymodule import ProcessOptions, PassingData, SNPData, TwoSNPData, read_data
from variation.src.common import number2nt, nt2number
from variation.src import Stock_250kDB
from variation.src.Stock_250kDB import Results, ResultsMethod, PhenotypeMethod, QCMethod, CallQC, Snps, SnpsQC, CallInfo, README
from pymodule.db import formReadmeObj
import sqlalchemy, numpy

class QC_250k(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('user', 1, ): [None, 'u', 1, 'database username', ],\
							('passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_dir', 0, ): [None, 'n', 1, 'If given, call data would be read from this directory/file rather than from call info table.\
								Does not work for QC_method_id=0. The data is sorted according to array_id. No call_info_id available. No db submission.', ],\
							('cmp_data_filename', 0,): [None, 'i', 1, 'the data file to be compared with. if not given, it gets figured out by QC_method_id.'],\
							('min_probability', 0, float, ):[-1, 'y', 1, 'minimum probability for a call to be non-NA if there is a 3rd column for probability.'],\
							('output_fname', 0, ): [None, 'o', 1, 'if given, QC results will be outputed into it.'],\
							('call_QC_table', 1, ): ['call_QC', 'q', 1, ''],\
							('QC_method_id', 1, int): [None, 'm', 1, 'id in table QC_method'],\
							('call_method_id', 1, int): [None, 'l', 1, 'id in table call_method'],\
							('run_type', 1, int): [1, 'e', 1, 'QC on 1=accession-wise or 2=snp-wise'],\
							('one_by_one', 0, int): [0, 'a', 0, 'whether to do QC of call_info entries one by one or as a whole. only for run_type=1.'],\
							('max_call_info_mismatch_rate', 0, float,):[-1, 'x', 1, 'maximum mismatch rate of an array call_info entry. used to exclude bad arrays to calculate snp-wise QC.'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-07-01
			use the new-version option_default_dict, and ProcessOptions
		2008-05-04
			if cmp_data_filename not specified, try to find in the data_description column in table QC_method.
		2008-04-20
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		#QC_method_id2cmp_data_filename = {1:"/home/crocea/script/variation/data/2010/data_2010_x_250k_y0001.tsv",\
		#		2:'/home/crocea/script/variation/data/perlegen/data_perlegen_ecotype_id_x_250k_y0101.tsv',\
		#		3:'/home/crocea/script/variation/stock20080403/data_y00001101.tsv'}
		#if self.cmp_data_filename == None:
		#	self.cmp_data_filename = QC_method_id2cmp_data_filename[self.QC_method_id]

	def get_call_info_id2fname(cls, db, QC_method_id, call_method_id, filter_calls_QCed=1, max_call_info_mismatch_rate=1, debug=0):
		"""
		2008-07-01
			add max_call_info_mismatch_rate <1 in the condition.
		2008-05-18
			add array_id to the call_info_id2fname structure
		2008-05-06
			use sqlalchemy connection
		2008-05-05
			add option filter_calls_QCed
		2008-04-20
			call_info entries that QC haven't been done.
			and its corresponding maternal_ecotype_id=paternal_ecotype_id (no crosses).
			
		"""
		sys.stderr.write("Getting call_info_id2fname ... ")
		#curs.execute("select distinct c.id, c.filename, a.maternal_ecotype_id, a.paternal_ecotype_id, q.id as qc_id, q.QC_method_id \
		#		from %s a, %s c left join %s q on c.id=q.call_info_id where c.array_id=a.id and a.maternal_ecotype_id=a.paternal_ecotype_id and c.method_id=%s"%\
		#		(array_info_table, call_info_table, call_QC_table, call_method_id))
		#rows = curs.fetchall()
		call_info_ls = CallInfo.query.filter_by(method_id=call_method_id).all()
		
		call_info_id2fname = {}
		call_info_ls_to_return = []
		for call_info in call_info_ls:
			if call_info.array.maternal_ecotype_id!=call_info.array.paternal_ecotype_id:	#ignore crosses
				continue
			if not call_info.array.maternal_ecotype_id:	#not linked to ecotypeid yet
				continue
			ignore_this = 0
			if filter_calls_QCed:
				for call_QC in call_info.call_qc_ls:
					if call_QC.qc_method_id==QC_method_id:	#same QC method has been done on this
						ignore_this = 1
						break
			
			#choose the call_QC with maximum no of non-NA pairs to get mismatch_rate
			if call_info.call_qc_ls and max_call_info_mismatch_rate<1:	#2008-07-01
				call_QC_with_max_no_of_non_NA_pairs = call_info.call_qc_ls[0]
				for call_QC in call_info.call_qc_ls:
					if call_QC.no_of_non_NA_pairs>call_QC_with_max_no_of_non_NA_pairs.no_of_non_NA_pairs:
						call_QC_with_max_no_of_non_NA_pairs = call_QC
				if call_QC_with_max_no_of_non_NA_pairs.mismatch_rate>max_call_info_mismatch_rate:
					ignore_this = 1
				call_info.call_QC_with_max_no_of_non_NA_pairs = call_QC_with_max_no_of_non_NA_pairs				
			else:
				call_info.call_QC_with_max_no_of_non_NA_pairs = None
			
			if ignore_this:
				continue
			call_info_id2fname[call_info.id] = [call_info.array.maternal_ecotype_id, call_info.filename, call_info.array_id]
			call_info_ls_to_return.append(call_info)
			#if debug and len(call_info_id2fname)>40:
			#	break
		
		sys.stderr.write("%s call files. Done.\n"%(len(call_info_id2fname)))
		return call_info_id2fname, call_info_ls_to_return
	
	get_call_info_id2fname = classmethod(get_call_info_id2fname)
	
	def get_array_id2fname(cls, curs, input_dir, array_info_table='array_info'):
		"""
		2008-05-18
			add array_id to the array_id2fname  structure
		2008-04-21
			the call files in input_dir are named like 'array_id'_call.tsv
		"""
		sys.stderr.write("Getting array_id2fname ... \n")
		file_dir_ls = os.listdir(input_dir)

		array_id2fname = {}
		for  i in range(len(file_dir_ls)):
			file_dir = file_dir_ls[i]
			#sys.stderr.write("%d/%d:\t%s\n"%(i+1,len(file_dir_ls),file_dir))
			pathname = os.path.join(input_dir, file_dir)
			array_id = int(file_dir.split('_')[0])
			curs.execute("select maternal_ecotype_id, paternal_ecotype_id from %s where id=%s"%(array_info_table, array_id))
			rows = curs.fetchall()
			if len(rows)>0:	#assuming array_id is unique and return only one entry if there's one
				maternal_ecotype_id, paternal_ecotype_id = rows[0]
				if maternal_ecotype_id==None and paternal_ecotype_id==None:
					sys.stderr.write("%s doesn't have corresponding ecotypeids. Ignored.\n"%(pathname))
					continue
				if maternal_ecotype_id!=paternal_ecotype_id:
					sys.stderr.write("%s is a cross with maternal_ecotype_id=%s and paternal_ecotype_id=%s. No QC now. Ignored.\n"%(pathname, maternal_ecotype_id, paternal_ecotype_id))
					continue
				array_id2fname[array_id] = [maternal_ecotype_id, pathname, array_id]
			else:
				sys.stderr.write("%s doesn't have entries in %s. Ignored.\n"%(pathname, array_info_table))
				continue
		sys.stderr.write("%s call files. Done.\n"%(len(array_id2fname)))
		return array_id2fname
	
	get_array_id2fname = classmethod(get_array_id2fname)
	
	def read_call_matrix(cls, call_info_id2fname, min_probability=-1, snps_name_set=None):
		"""
		2008-05-20
			return PassingData to wrap everything
		2008-05-19
			the header is now 'chromosome_position'. no candidate alleles behind them.
			header.append('%s_%s'%(SNP_id_ls[0], SNP_id_ls[1]))
		2008-05-18
			rename strain_acc_list, category_list to call_info_id_ls and ecotype_id_ls
			add array_id_ls
			return array_id_ls
		2008-05-06 add snps_name_set
		2008-05-03 add min_probability
		2008-04-20
			fake header, strain_acc_list, category_list, data_matrix
		"""
		sys.stderr.write("Creating call matrix ... \n")
		header = ['', '']	#1st and 2nd is header for 1st two columns.
		call_matrix = []
		call_info_id_ls = []
		ecotype_id_ls = []
		array_id_ls = []
		counter = 0
		call_info_id_ls = call_info_id2fname.keys()
		no_of_entries = len(call_info_id_ls)
		for i in range(no_of_entries):
			call_info_id = call_info_id_ls[i]
			ecotype_id, fname, array_id = call_info_id2fname[call_info_id]
			sys.stderr.write("%s%d/%d:\t\t%s"%('\x08'*100, i+1, no_of_entries, fname))
			call_info_id_ls.append(call_info_id)
			ecotype_id_ls.append(ecotype_id)
			array_id_ls.append(array_id)
			reader = csv.reader(open(fname), delimiter='\t')
			reader.next()	#throw away the first line
			data_row = []
			
			for row in reader:
				SNP_id, call = row[:2]
				if snps_name_set and SNP_id not in snps_name_set:
					continue
				if counter==0:	#first file
					SNP_id_ls = SNP_id.split('_')
					header.append('%s_%s'%(SNP_id_ls[0], SNP_id_ls[1]))
				if len(row)==3:
					probability = float(row[2])
					if probability < min_probability:
						call = 'NA'
				data_row.append(nt2number[call])
			del reader
			call_matrix.append(data_row)
			counter += 1
		pdata = PassingData(header=header, call_info_id_ls=call_info_id_ls, array_id_ls=array_id_ls,\
			ecotype_id_ls=ecotype_id_ls, data_matrix=call_matrix)
		sys.stderr.write("Done.\n")
		return pdata
	
	read_call_matrix = classmethod(read_call_matrix)
	
	def submit_to_call_QC(cls, session, row_id2NA_mismatch_rate, QC_method_id, user, min_probability, row_id12row_id2, call_method_id, readme):
		"""
		2008-05-21
			ecotype_id, call_info_id = row_id	#bug here, order changed.
		2008-05-19
			NA_mismatch_ls was expanded
		2008-05-06
			add readme
		2008-05-05
			add ecotype_id, min_probability, tg_ecotype_id
		"""
		sys.stderr.write("Submitting row_id2NA_mismatch_rate to database ...")
		row_id_ls = row_id2NA_mismatch_rate.keys()
		row_id_ls.sort()	#try to keep them in call_info_id order
		for row_id in row_id_ls:
			NA_mismatch_ls = row_id2NA_mismatch_rate[row_id]
			ecotype_id, call_info_id = row_id	#bug here, order changed.
			tg_ecotype_id = row_id12row_id2[row_id]
			NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs, relative_NA_rate, relative_no_of_NAs, relative_no_of_totals = NA_mismatch_ls
			#call_QC stores the relative NA rate. call_info already stores the independent NA rate
			NA_rate, no_of_NAs, no_of_totals = relative_NA_rate, relative_no_of_NAs, relative_no_of_totals
			callqc = CallQC(call_info_id=call_info_id, min_probability=min_probability, ecotype_id=ecotype_id, tg_ecotype_id=tg_ecotype_id,\
						qc_method_id=QC_method_id, call_method_id=call_method_id, NA_rate=NA_rate, mismatch_rate=mismatch_rate,\
						no_of_NAs=no_of_NAs, no_of_totals=no_of_totals, no_of_mismatches=no_of_mismatches, no_of_non_NA_pairs=no_of_non_NA_pairs,\
						created_by=user)
			callqc.readme = readme
			session.save(callqc)
			"""
			data_insert_ls = [row_id[0]] + NA_mismatch_ls + [QC_method_id, user]	#row_id is (call_info_id, ecotypeid)
			curs.execute("insert into " + call_QC_table + " (call_info_id, NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs, QC_method_id, created_by)\
				values(%s, %s, %s, %s, %s, %s, %s, %s, %s)", data_insert_ls)
			"""
		sys.stderr.write("Done.\n")
	
	submit_to_call_QC = classmethod(submit_to_call_QC)
	
	def cal_independent_NA_rate(cls, db, min_probability, readme):
		"""
		2008-07-13
			save_or_update(call_info)
		2008-05-06
			use sqlalchemy connection
		2008-04-20
			calculate indepent (no data to be compared) NA rates.
			update it in the db.
		"""
		sys.stderr.write("Calculating indepent NA rate ... \n")
		#curs.execute("select c.id, c.filename from %s c where c.NA_rate is null order by id"%\
		#			(call_info_table))
		#rows = curs.fetchall()
		#call_info_ls = db.session.query(CallInfo).filter_by(NA_rate=None).all()
		call_info_ls = CallInfo.query.filter_by(NA_rate=None).all()
		no_of_rows = len(call_info_ls)
		sys.stderr.write("\tTotally, %d call_info entries to be processed.\n"%no_of_rows)
		for i in range(no_of_rows):
			call_info = call_info_ls[i]
			sys.stderr.write("%d/%d:\t%s\n"%(i+1, no_of_rows, call_info.filename))
			reader = csv.reader(open(call_info.filename), delimiter='\t')
			reader.next()	#throw away the first line
			no_of_totals = 0
			no_of_NAs = 0
			for row in reader:
				SNP_id, call = row[:2]
				no_of_totals += 1
				if len(row)==3:
					probability = float(row[2])
					if probability < min_probability:
						call = 'NA'
				if call=='NA':
					no_of_NAs += 1
			if no_of_totals!=0:
				call_info.NA_rate = float(no_of_NAs)/no_of_totals
				sys.stderr.write("%s\n"%call_info.NA_rate)
				call_info.readme = readme
				db.session.save_or_update(call_info)
			del reader
			#curs.execute("update " + call_info_table + " set NA_rate=%s where id=%s",\
			#		(NA_rate, call_info_id))
		sys.stderr.write("Done.\n")
	
	def output_row_id2NA_mismatch_rate(cls, row_id2NA_mismatch_rate, output_fname, file_1st_open=1):
		"""
		2008-05-12
			smart way to figure out how to deal with row_id and its labels
		2008-05-12
			tsv => csv fromat
		2008-05-11
			add file_1st_open argument
		2008-04-22
		"""
		sys.stderr.write("Outputting row_id2NA_mismatch_rate to %s ..."%(output_fname))
		if file_1st_open:
			open_flag = 'w'
		else:
			open_flag = 'a'
		writer = csv.writer(open(output_fname, open_flag))
		NA_mismatch_ls_header = ['NA_rate', 'mismatch_rate', 'no_of_NAs', 'no_of_totals', \
				'no_of_mismatches', 'no_of_non_NA_pairs', 'relative_NA_rate', 'relative_no_of_NAs', 'relative_no_of_totals']
		row_id_ls = row_id2NA_mismatch_rate.keys()
		row_id_ls.sort()	#try to keep them in call_info_id order
		if len(row_id_ls)>0:
			row_id0 = row_id_ls[0]
			if not isinstance(row_id0, str) and hasattr(row_id0, '__len__'):
				header = ['']*len(row_id0)
			else:
				header = ['']
			header += NA_mismatch_ls_header
			writer.writerow(header)
			for row_id in row_id_ls:
				NA_mismatch_ls = row_id2NA_mismatch_rate[row_id]
				if isinstance(row_id, tuple):
					row_id_ls = list(row_id)
				elif isinstance(row_id, list):
					row_id_ls = row_id
				else:
					row_id_ls = [row_id]
				writer.writerow(row_id_ls + NA_mismatch_ls)
		del writer
		sys.stderr.write("Done.\n")
	
	output_row_id2NA_mismatch_rate = classmethod(output_row_id2NA_mismatch_rate)
	
	def get_snps_name2snps_id(cls, db):
		"""
		2008-08-17
			new db interface from Stock_250kDB.py
		2008-05-05
		"""
		sys.stderr.write("Getting snps_name2snps_id ...")
		#from sqlalchemy.sql import select
		#s = select([db.tables['snps'].c.id, db.tables['snps'].c.name])
		#result = db.connection.execute(s)
		result = db.session.execute("select id, name from %s"%(Snps.table.name))	#2008-08-17
		snps_name2snps_id = {}
		for row in result:
			snps_id, snps_name = row
			snps_name2snps_id[snps_name] = snps_id
		sys.stderr.write("Done.\n")
		return snps_name2snps_id
	
	get_snps_name2snps_id = classmethod(get_snps_name2snps_id)
	
	QC_method_id2snps_table = {1:'at.locus',\
									2:'',\
									3:'stock.snps',\
									4:'',\
									5:'at.locus',\
									6:'',\
									7:'stock.snps',\
									8:'dbsnp.snps'}
	
	
	def qcDataMatrixVSsnpData(self, pdata, snps_name2snps_id, snpData2, curs, session, readme):
		"""
		2008-08-16
			split from run() to enable one_by_one option
		"""
		#swap the ecotype_id_ls and call_info_id_ls when passing them to SNPData. now strain_acc_list=ecotype_id_ls
		snpData1 = SNPData(header=pdata.header, strain_acc_list=pdata.ecotype_id_ls, category_list=pdata.call_info_id_ls, data_matrix=pdata.data_matrix, \
						min_probability=self.min_probability, call_method_id=self.call_method_id, col_id2id=snps_name2snps_id,\
						max_call_info_mismatch_rate=self.max_call_info_mismatch_rate, snps_table='stock_250k.snps')	#snps_table is set to the stock_250k snps_table
		
		twoSNPData = TwoSNPData(SNPData1=snpData1, SNPData2=snpData2, curs=curs, \
							QC_method_id=self.QC_method_id, user=self.user, row_matching_by_which_value=0, debug=self.debug)
		
		if self.run_type==1:
			row_id2NA_mismatch_rate = twoSNPData.cmp_row_wise()
		elif self.run_type==2:
			twoSNPData.save_col_wise(session, readme)
			row_id2NA_mismatch_rate = {}
		else:
			sys.stderr.write("run_type=%s is not supported.\n"%self.run_type)
			sys.exit(5)
		passingdata = PassingData()
		passingdata.row_id2NA_mismatch_rate = row_id2NA_mismatch_rate
		passingdata.row_id12row_id2 = twoSNPData.row_id12row_id2
		return passingdata
	
	def findOutCmpDataFilename(cls, cmp_data_filename, QC_method_id, QCMethod_class):
		"""
		2008-08-26
			add QCMethod_class
		2008-08-16
			split from run() to let QC_149.py to call it
		"""
		# if cmp_data_filename not specified, try to find in the data_description column in table QC_method.
		if not cmp_data_filename and QC_method_id!=0:
			qm = QCMethod_class.query.get(QC_method_id)
			if qm and qm.data_description:
				data_description_ls = qm.data_description.split('=')
				if len(data_description_ls)>1:
					cmp_data_filename = qm.data_description.split('=')[1].strip()
					return cmp_data_filename
		#after db query, cmp_data_filename is still nothing, exit program.
		if not cmp_data_filename and QC_method_id!=0:
			sys.stderr.write("cmp_data_filename is still nothing even after db query. please specify it on the commandline.\n")
			sys.exit(3)
	
	findOutCmpDataFilename = classmethod(findOutCmpDataFilename)
	
	def run(self):
		"""
		2008-04-25
			return None if QC_method_id==0
		2008-04-20
			for plone to call it just to get row_id2NA_mismatch_rate
		"""
		#database connection and etc
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.user,
				   password=self.passwd, hostname=self.hostname, database=self.dbname)
		db.setup()
		session = db.session
		session.begin()
		#transaction = session.create_transaction()
		
		self.cmp_data_filename = self.findOutCmpDataFilename(self.cmp_data_filename, self.QC_method_id)
		
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.user, passwd = self.passwd)
		curs = conn.cursor()
		self.curs = curs
		if self.debug:
			import pdb
			pdb.set_trace()
		
		readme = formReadmeObj(sys.argv, self.ad, README)
		session.save(readme)
		
		QC_method_id2snps_table = self.QC_method_id2snps_table
		
		if self.QC_method_id==0:
			self.cal_independent_NA_rate(db, self.min_probability, readme)
			row_id2NA_mismatch_rate = None
		else:
			#from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
			header, strain_acc_list, category_list, data_matrix = read_data(self.cmp_data_filename)
			strain_acc_list = map(int, strain_acc_list)	#it's ecotypeid, cast it to integer to be compatible to the later ecotype_id_ls from db
			snpData2 = SNPData(header=header, strain_acc_list=strain_acc_list, \
							data_matrix=data_matrix, snps_table=QC_method_id2snps_table.get(self.QC_method_id))	#category_list is not used.
			
			if self.input_dir and os.path.isdir(self.input_dir):
				#04/22/08 Watch: call_info_id2fname here is fake, it's actually keyed by (array_id, ecotypeid)
				#no submission to db
				call_info_id2fname = self.get_array_id2fname(curs, self.input_dir)
			elif self.input_dir and os.path.isfile(self.input_dir):	#it's file
				call_info_id2fname = None
			else:
				if self.run_type==2:	#no filtering on call_info entries that have been QCed.
					filter_calls_QCed=0
				elif self.run_type==1:
					filter_calls_QCed = 1
					self.max_call_info_mismatch_rate = 1	#don't use this when doing accession-wise QC
				else:
					sys.stderr.write("run_type=%s is not supported.\n"%self.run_type)
					sys.exit(5)
				call_info_id2fname, call_info_ls_to_return = self.get_call_info_id2fname(db, self.QC_method_id, self.call_method_id, \
																						filter_calls_QCed, self.max_call_info_mismatch_rate, self.debug)
			if self.run_type==2:
				snps_name2snps_id = self.get_snps_name2snps_id(db)
			else:
				snps_name2snps_id = None
			
			if call_info_id2fname:
				if self.one_by_one and self.run_type==1:	#one_by_one only for QC by accession
					row_id2NA_mismatch_rate = {}
					row_id12row_id2 = {}
					counter = 0
					for call_info_id, value in call_info_id2fname.iteritems():
						counter += 1
						print "No", counter
						tmp_dict = {}
						tmp_dict[call_info_id] = value
						pdata = self.read_call_matrix(tmp_dict, self.min_probability)
						passingdata = self.qcDataMatrixVSsnpData(pdata, snps_name2snps_id, snpData2, curs, session, readme)
						row_id2NA_mismatch_rate.update(passingdata.row_id2NA_mismatch_rate)
						row_id12row_id2.update(passingdata.row_id12row_id2)
						del pdata
						
						if self.debug and counter==10:
							break
				else:
					pdata = self.read_call_matrix(call_info_id2fname, self.min_probability)
					passingdata = self.qcDataMatrixVSsnpData(pdata, snps_name2snps_id, snpData2, curs, session, readme)
					row_id2NA_mismatch_rate = passingdata.row_id2NA_mismatch_rate
					row_id12row_id2 = passingdata.row_id12row_id2
					del pdata
			else:
				#input file is SNP by strain format. double header (1st two lines)
				header, snps_name_ls, category_list, data_matrix = FilterStrainSNPMatrix.read_data(self.input_dir, double_header=1)
				pdata = PassingData()
				pdata.ecotype_id_ls = header[0][2:]
				pdata.call_info_id_ls = header[1][2:]
				data_matrix = numpy.array(data_matrix)
				pdata.data_matrix = data_matrix.transpose()
				pdata.header = ['', ''] + snps_name_ls	#fake a header for SNPData
				passingdata = self.qcDataMatrixVSsnpData(pdata, snps_name2snps_id, snpData2, curs, session, readme)
				row_id2NA_mismatch_rate = passingdata.row_id2NA_mismatch_rate
				row_id12row_id2 = passingdata.row_id12row_id2
				del pdata
		
		if self.output_fname and self.run_type==1 and row_id2NA_mismatch_rate:
			self.output_row_id2NA_mismatch_rate(row_id2NA_mismatch_rate, self.output_fname)
		
		if self.run_type==1 and self.commit and not self.input_dir and row_id2NA_mismatch_rate:
			#if self.input_dir is given, no db submission. call_info_id2fname here is fake, it's actually keyed by (array_id, ecotypeid)
			#row_id2NA_mismatch_rate might be None if it's method 0.
			self.submit_to_call_QC(session, row_id2NA_mismatch_rate, self.QC_method_id, self.user, self.min_probability, \
								row_id12row_id2, self.call_method_id, readme)
		if self.commit:
			curs.execute("commit")
			session.commit()
		else:
			session.rollback()
		
		self.row_id2NA_mismatch_rate = row_id2NA_mismatch_rate	#for plone to get the data structure

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = QC_250k
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	
	"""
	from pymodule import process_options, generate_program_doc
	opts_dict = process_options(sys.argv, QC_250k.option_default_dict, error_doc=generate_program_doc(sys.argv[0], QC_250k.option_default_dict)+QC_250k.__doc__)
	
	instance = QC_250k(**opts_dict)
	instance.run()
	"""