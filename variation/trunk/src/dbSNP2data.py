#!/usr/bin/env python
"""

Examples:
	dbSNP2data.py -i justin_data -o justin_data.csv -r
	
	dbSNP2data.py -i justin_data -o justin_data.csv -r -t
	
	dbSNP2data.py -i calls -o /tmp/chicago.data.y -s ecotype -n snps
	
	dbSNP2data.py -z localhost -d stock20071008 -i calls -o stock20071008/data.tsv -s ecotype -n snps
	
	dbSNP2data.py -z localhost -d stock20071008 -i calls -o /tmp/data.tsv -s ecotype -n snps
	
	#to see how many all-NA strains were discarded
	dbSNP2data.py -z localhost -d stock20071008 -i calls -o /tmp/data10101100.tsv -s ecotype -n snps -y 10101100
	
	#all strains, but resolve duplicated (nativename, stockparent)s
	dbSNP2data.py -z localhost -d stock20071008 -i calls -o /tmp/data00101100.tsv -s ecotype -n snps -y 00101100
	
	#output data for all ecotypeid. however duplicated calls with same ecotypeid are imputed randomly
	dbSNP2data.py -z localhost -d stock20071008 -i calls -o /tmp/data00001100.tsv -s ecotype -n snps -y 00001100

	#output 250k SNP data (not to resolve duplicated calls, all SNPs. uncomment another sql line to get only SNPs shared with 149SNP)
	dbSNP2data.py -z localhost -d stock20071008 -i calls_250k -o /tmp/data_250k.tsv -s ecotype -n snps_250k -y 10001101 -r
	
	#output 149SNP (table calls_byseq)
	dbSNP2data.py -z localhost -d stock20071227 -i calls_byseq -o stock20071227/data_y10001101.tsv -s ecotype -n snps -y 10001101 -r
	
	#output all 149SNP data
	dbSNP2data.py -o /tmp/stock_149SNP.tsv
	
	#output only 149SNP data with GPS info
	dbSNP2data.py -o /tmp/stock_149SNP_y10001111.tsv -y 1 -r
	
	#output 384-illumina data
	dbSNP2data.py -o /tmp/384-illumina_y00001101.tsv -n dbsnp.snps -s dbsnp.accession -i dbsnp.calls -y 00001101
	
	#output 384-illumina data in csv format. output matrix transposed, in SNP by strain.
	dbSNP2data.py  -o /tmp/384-illumina_y0000111111.csv -n dbsnp.snps -s dbsnp.accession -i dbsnp.calls -y 0000111111
	
Description:
	output SNP data from database schema
	
	Turning on each bit in processing_bits (0=off, 1=on):
	1. 0: everything in strain_info_table, 1: only include strains with GPS info, 2: north american strains only, 3: 2010's 192 strains
	2. include columns of other strain info (latitude, longitude, nativename, stockparent, site, country)
	3. resolve duplicated calls (unique constraint on (nativename, stockparent))
	4. toss out rows to make distance matrix NA free
	5. need heterozygous call
	6. with header line
	7. use alphabet to represent nucleotide, not number
	8. discard strains with all-NA data
	9. output matrix type (in terms of row X column). 0=(strain X SNP, default), 1=(SNP X strain)
	10. delimiter. 0=tab default, 1=comma.
	
	you can specify the bits up to the one you want to change and omit the rest. i.e.
	-y 11
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/')))
import psycopg2 as psycopg
import sys, getopt, csv, re
from annot.bin.codense.common import db_connect, org_short2long, org2tax_id
from variation.src.common import nt2number, number2nt, ab2number, number2ab
import Numeric as num
from sets import Set
from pymodule import write_data_matrix

class dbSNP2data(object):
	__doc__ = __doc__
	"""
	2007-02-19
	"""
	option_default_dict = {('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock', 'd', 1, '', ],\
							('user', 1, ): [None, 'u', 1, 'database username', ],\
							('passwd', 1, ):[None, 'p', 1, 'database password', ],\
							('input_table', 1, ): ['calls', 'i', 1, 'table containing calls for each accession & snp.'],\
							('strain_info_table', 1, ): ['ecotype', 's', 1, 'Table with info about each accession/ecotype. could be "strain_info" or "ecotype"'],\
							('output_fname', 1, ): [None, 'o', 1, 'Output Filename'],\
							('snp_locus_table', 1, ): ['snps', 'n', 1, 'Table with info about snps. could be "snp_locus" or "snps"'],\
							('processing_bits', 1, ): ['0000111100', 'y', 1, 'processing bits to control which processing step should be turned on.\
								default is 10101101. for what each bit stands, see Description.' ],\
							('db_connection_type', 1, int): [1, 'm', 1, 'which type of database. 1=MySQL. 2=PostgreSQL.',],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	"""
	2008-04-28
	option_default_dict is a dictionary for option handling, including argument_default_dict info
		the key is a tuple, ('argument_name', is_argument_required, argument_type) and argument_type is optional.
		
		the value could just be the default_value or a list = [default_value, 'short_option', has_argument, description_for_option]
		'short_option', has_argument, description_for_option are all orderly optional. You can just supply 'short_option' or 'short_option', has_argument.
	"""
	def __init__(self, **keywords):
		"""
		2008-05-08
			add self.output_matrix_type, self.delimiter
		2008-04-28
			use ProcessOptions, newer option handling class
		2008-04-25 use new option handling
		2007-02-25
			add argument toss_out_rows
		2007-07-11
			add mysql_connection
		2007-09-23
			use processing_bits to control processing steps
		"""
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		
		#below are all default values
		processing_bits_ls = [0,0,0,0,1,1,1,1, 0, 0]
		
		for i in range(len(self.processing_bits)):
			processing_bits_ls[i] = int(self.processing_bits[i])
		#now pass all values
		self.only_include_strains_with_GPS,\
		self.include_other_strain_info,\
		self.resolve_duplicated_calls,\
		self.toss_out_rows,\
		self.need_heterozygous_call,\
		self.with_header_line,\
		self.nt_alphabet,\
		self.discard_all_NA_strain,\
		self.output_matrix_type, \
		self.delimiter_type = processing_bits_ls
		
		delimiter_dict = {0: '\t', \
						1: ','}
		
		self.delimiter = delimiter_dict[self.delimiter_type]
	
	def get_snp_id2index(self, curs, input_table, snp_locus_table):
		"""
		2007-03-05
			add snp_locus_table to sort snp based on chromosome and position
		"""
		sys.stderr.write("Getting snp_id2index ...")
		snp_id2index = {}
		snp_id_list = []
		curs.execute("select distinct i.snp_id, s.chromosome, s.position from %s i, %s s where i.snp_id=s.id order by chromosome, position"%(input_table, snp_locus_table))
		rows = curs.fetchall()
		for row in rows:
			snp_id = row[0]
			snp_id_list.append(snp_id)
			snp_id2index[snp_id] = len(snp_id2index)
		sys.stderr.write("Done.\n")
		return snp_id2index, snp_id_list
	
	def get_snp_id2index_m(self, curs, input_table, snp_locus_table):
		"""
		2008-05-05
			deal with tables in schema dbsnp
		2007-07-11
			mysql version of get_snp_id2index
		2007-12-13
			for 'snps': take only SNPs appearing in input_table
			for 'snps_250k': take only the SNPs overlapping with 'snps'/149SNP
		"""
		sys.stderr.write("Getting snp_id2index ..m.")
		snp_id2index = {}
		snp_id_list = []
		if snp_locus_table == 'snps':
			curs.execute("select distinct i.snpid, s.chromosome, s.position from %s i, %s s where i.snpid=s.id order by chromosome, position"%(input_table, snp_locus_table))
		elif snp_locus_table == 'snps_250k':
			#2007-12-13 only the snps overlapping with 149SNP
			#curs.execute("select distinct s1.id, s1.chromosome, s1.position from %s s2, %s s1 where s1.chromosome=s2.chromosome and s1.position=s2.position order by chromosome, position"%('snps', snp_locus_table))
			#2007-12-13 all 250k SNPs
			curs.execute("select distinct s1.id, s1.chromosome, s1.position from %s s1"%(snp_locus_table))
		elif snp_locus_table == 'dbsnp.snps':
			curs.execute("select distinct s.id, s.chromosome, s.position from %s i, %s s, %s m where s.id=m.snps_id and i.snps_id=s.id order by chromosome, position"%(input_table, snp_locus_table, 'dbsnp.snps_ab_allele_mapping'))
		rows = curs.fetchall()
		for row in rows:
			snp_id = row[0]
			snp_id_list.append(snp_id)
			snp_id2index[snp_id] = len(snp_id2index)
		sys.stderr.write("Done.\n")
		return snp_id2index, snp_id_list
	
	def get_strain_id2index(self, curs, input_table):
		sys.stderr.write("Getting strain_id2index ...")
		strain_id2index = {}
		strain_id_list = []
		curs.execute("select distinct strain_id from %s order by strain_id"%(input_table))
		rows = curs.fetchall()
		for row in rows:
			strain_id = row[0]
			strain_id_list.append(strain_id)
			strain_id2index[strain_id] = len(strain_id2index)
		sys.stderr.write("Done.\n")
		return strain_id2index, strain_id_list
	
	def get_strain_id2index_m(self, curs, input_table, strain_info_table, only_include_strains_with_GPS=0, resolve_duplicated_calls=0):
		"""
		2008-05-05
			deal with tables in schema dbsnp
		2007-07-11
			mysql version of get_strain_id2index
		2007-09-12
			strain_info_table
			only strains whose latitude and longitude not null
		2007-09-22
			solve the call-duplication problem.
		2007-09-23
			add only_include_strains_with_GPS and resolve_duplicated_calls
		2007-10-09
			only_include_strains_with_GPS has more meanings
		2007-12-13
			add duplicate
			strain_id = (strain_id, duplicate)
		2007-12-16
			add strain_id2acc, strain_id2category
			abandon get_strain_id_info_m()
		2008-01-02
			replicate replaces duplicate for table calls_byseq
		"""
		sys.stderr.write("Getting strain_id2index ..m.")
		strain_id2index = {}
		strain_id_list = []
		nativename2strain_id = {}
		strain_id2acc = {}
		strain_id2category = {}
		if input_table=='calls_byseq':
			common_sql_string = "select distinct d.ecotypeid, d.replicate, s.nativename, s.stockparent from %s d, %s s"%(input_table, strain_info_table)
		elif input_table=='dbsnp.calls':
			common_sql_string = "select distinct s.ecotype_id, s.duplicate, s.id, s.id from %s d, %s s"%\
				(input_table, strain_info_table)
		else:
			common_sql_string = "select distinct d.ecotypeid, d.replicate, s.nativename, s.stockparent from %s d, %s s"%(input_table, strain_info_table)
			
		if input_table=='dbsnp.calls':
			curs.execute("%s where s.id=d.accession_id order by ecotype_id, duplicate"%(common_sql_string))
		else:
			if only_include_strains_with_GPS==1:
				curs.execute("%s where d.ecotypeid=s.id and s.latitude is not null and s.longitude is not null  order by ecotypeid, nativename, stockparent"%(common_sql_string))
			elif only_include_strains_with_GPS==2:	#2007-10-01 north american samples
				curs.execute("%s where d.ecotypeid=s.id and s.latitude is not null and s.longitude is not null and s.longitude<-60 and s.longitude>-130 order by ecotypeid, nativename, stockparent"%(common_sql_string))
			elif only_include_strains_with_GPS==3:
				curs.execute("%s, batch_ecotype be, batch b where b.batchname='192' and b.id=be.batchid and s.id=be.ecotypeid and d.ecotypeid=s.id and s.latitude is not null and s.longitude is not null  order by ecotypeid, nativename, stockparent"%(common_sql_string))
			else:
				curs.execute("%s where d.ecotypeid=s.id order by ecotypeid, nativename, stockparent"%(common_sql_string))
		rows = curs.fetchall()
		for row in rows:
			ecotypeid, duplicate, nativename, stockparent = row
			if input_table=='dbsnp.calls':
				strain_id = row[2]	#accession_id
			else:
				strain_id = (ecotypeid, duplicate)
				nativename = nativename.upper()
			if resolve_duplicated_calls:
				key_pair = (nativename, stockparent)
				if key_pair not in nativename2strain_id:
					nativename2strain_id[key_pair] = strain_id
					strain_id_list.append(strain_id)
					strain_id2index[strain_id] = len(strain_id2index)
					strain_id2acc[strain_id] = ecotypeid
					strain_id2category[strain_id] = duplicate
			else:
				strain_id_list.append(strain_id)
				strain_id2index[strain_id] = len(strain_id2index)
				strain_id2acc[strain_id] = ecotypeid
				strain_id2category[strain_id] = duplicate
		sys.stderr.write("Done.\n")
		return strain_id2index, strain_id_list, nativename2strain_id, strain_id2acc, strain_id2category
	
	def get_strain_id_info(self, curs, strain_id_list, strain_info_table):
		sys.stderr.write("Getting strain_id_info ...")
		strain_id2acc = {}
		strain_id2category = {}
		for strain_id in strain_id_list:
			curs.execute("select acc, category from %s where id=%s"%(strain_info_table, strain_id))
			rows = curs.fetchall()
			acc, category = rows[0]
			strain_id2acc[strain_id] = acc
			strain_id2category[strain_id] = category
		sys.stderr.write("Done.\n")
		return strain_id2acc, strain_id2category
	
	def get_strain_id_info_m(self, curs, strain_id_list, strain_info_table):
		"""
		2008-05-08
			deprecated
		2007-07-11
			mysql version of get_strain_id_info
		2007-09-13
			replace "name, nativename" with "id, name"
		2007-09-20
			replace 'name' with 'nativename'
		2007-12-13
			changes following strain_id = (strain_id, duplicate)
		2007-12-16
			abandoned, get_strain_id2index_m() supercedes this function. 
		"""
		sys.stderr.write("Getting strain_id_info ..m.")
		strain_id2acc = {}
		strain_id2category = {}
		for strain_id in strain_id_list:
			curs.execute("select id, nativename from %s where id=%s"%(strain_info_table, strain_id[0]))
			rows = curs.fetchall()
			acc, category = rows[0]
			strain_id2acc[strain_id] = acc
			strain_id2category[strain_id] = category
		sys.stderr.write("Done.\n")
		return strain_id2acc, strain_id2category
	
	def get_snp_id_info(self, curs, snp_id_list, snp_locus_table):
		sys.stderr.write("Getting snp_id_info ...")
		snp_id2acc = {}
		for snp_id in snp_id_list:
			curs.execute("select acc from %s where id=%s"%(snp_locus_table, snp_id))
			rows = curs.fetchall()
			acc = rows[0][0]
			snp_id2acc[snp_id] = acc
		sys.stderr.write("Done.\n")
		return snp_id2acc
	
	def get_snp_id_info_m(self, curs, snp_id_list, snp_locus_table):
		"""
		2007-07-11
			mysql version of get_snp_id_info
		2007-12-13
			snp_locus_table could be 'snps' or 'snps_250k'
		2007-12-18
			table 'snps' and 'snps_250k' have same field names, no need to differentiate two
		"""
		sys.stderr.write("Getting snp_id_info ..m.")
		snp_id2info = {}
		if snp_locus_table=='dbsnp.snps':
			snp_acc_column = 'name'
		else:
			snp_acc_column = 'snpid'
		#sql_sentence = "select %s, chromosome, position " + " from %s where id=%s"%(snp_locus_table, snp_id)
		for snp_id in snp_id_list:
			curs.execute("select %s, chromosome, position from %s where id=%s"%(snp_acc_column, snp_locus_table, snp_id))
			rows = curs.fetchall()
			acc, chromosome, position = rows[0]
			snp_id2info[snp_id] = (acc, chromosome, position)
		sys.stderr.write("Done.\n")
		return snp_id2info
	
	def get_data_matrix(self, curs, strain_id2index, snp_id2index, nt2number, input_table, need_heterozygous_call):
		"""
		2007-03-20
			check whether strain_id and snp_id are in strain_id2index, snp_id2index
		"""
		sys.stderr.write("Getting data_matrix ...\n")
		data_matrix = num.zeros([len(strain_id2index), len(snp_id2index)])
		curs.execute("DECLARE crs CURSOR FOR select strain_id, snp_id, call from %s"%(input_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				strain_id, snp_id, call = row
				if strain_id in  strain_id2index and snp_id in snp_id2index:	#2007-03-20
					call_number = nt2number[call]
					if need_heterozygous_call:
						data_matrix[strain_id2index[strain_id], snp_id2index[snp_id]] = call_number
					elif call_number<=4:	#single letter or NA
						data_matrix[strain_id2index[strain_id], snp_id2index[snp_id]] = call_number
				counter += 1
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
		curs.execute("close crs")
		sys.stderr.write("Done.\n")
		return data_matrix
	
	def get_data_matrix_m(self, curs, strain_id2index, snp_id2index, nt2number, input_table, need_heterozygous_call, snps_id2mapping=None):
		"""
		2008-05-05
			deal with tables in schema dbsnp
		2007-07-11
			mysql version of get_data_matrix
			
			callhet tells whether it's heterozygous or not.
		2007-09-20
			upper case for callhet
		2007-12-13
			input_table: 'calls', 'calls_250k'
			changes following strain_id = (strain_id, duplicate)
		2007-12-14
			after an index on calls_250k(snpid) was created in db, output data matrix one snpid after another.
			otherwise, the program blew up the memory by sifting through all data.
		2008-01-01
			add code to deal with table calls_byseq
		"""
		sys.stderr.write("Getting data_matrix ..m.\n")
		data_matrix = num.zeros([len(strain_id2index), len(snp_id2index)])
		if input_table == 'calls':
			common_sql_string = "select ecotypeid, replicate, snpid, call1, call2 from %s"%(input_table)
		elif input_table == 'calls_byseq':
			common_sql_string = "select ecotypeid, replicate, snpid, call1, call2 from %s"%(input_table)
		elif input_table == 'calls_250k':
			common_sql_string = "select ecotypeid, duplicate, snpid, snpcall from %s"%(input_table)
		elif input_table == 'dbsnp.calls':
			common_sql_string = "select accession_id, snps_id, genotype from %s"%(input_table)
		else:
			sys.stderr.write("Error: SNP call table '%s' not supported.\n"%(input_table))
			sys.exit(3)
		snp_counter = 0
		for snp_id in snp_id2index:
			snp_counter += 1
			counter = 0
			if self.report:
				sys.stderr.write("%s\tSNP %s=%s"%('\x08'*80, snp_counter, snp_id))
			if input_table == 'dbsnp.calls':
				curs.execute("%s where snps_id=%s"%(common_sql_string, snp_id))
			else:
				curs.execute("%s where snpid=%s"%(common_sql_string, snp_id))
			rows = curs.fetchmany(5000)
			while rows:
				for row in rows:
					if input_table == 'calls' or input_table=='calls_byseq':
						strain_id, duplicate, snp_id, call, callhet = row
					elif input_table == 'calls_250k':
						strain_id, duplicate, snp_id, call = row
						callhet = None
					elif input_table == 'dbsnp.calls':
						strain_id, snp_id, call = row
						callhet = None
					if input_table != 'dbsnp.calls':
						strain_id = (strain_id, duplicate)
					call = call.upper()
					if callhet:
						callhet.upper()	#2007-09-20	just in case
						call = call+callhet
					if strain_id in  strain_id2index and snp_id in snp_id2index:	#2007-03-20
						if input_table == 'dbsnp.calls' and not snps_id2mapping:
							call_number = ab2number[call]
						else:
							if snps_id2mapping:
								call = snps_id2mapping[snp_id][call]
							call_number = nt2number[call]
						if need_heterozygous_call:
							data_matrix[strain_id2index[strain_id], snp_id2index[snp_id]] = call_number
						elif call_number<=4:	#single letter or NA
							data_matrix[strain_id2index[strain_id], snp_id2index[snp_id]] = call_number
					counter += 1
				rows = curs.fetchmany(5000)
				if self.report:
					sys.stderr.write("%s\tSNP %s=%s\t\t%s"%('\x08'*80, snp_counter, snp_id, counter))
		sys.stderr.write("Done.\n")
		return data_matrix
	
	def get_majority_call_number(cls, call_counter_ls):
		"""
		2008-04-04
			make it class method
		2007-09-22
			return the call with maximum vote or NA
		"""
		index_with_max_value = num.argmax(call_counter_ls)
		if index_with_max_value>0:	#if 0(NA) is maximum or shared-maximum (call_counter_ls=0), no need to pursue
			for i in range(1, len(call_counter_ls)):
				if i!=index_with_max_value and call_counter_ls[i]==call_counter_ls[index_with_max_value]:	#a shared-maximum call, can't be resolved. =>NA
					index_with_max_value = 0
					break
		return index_with_max_value
	get_majority_call_number = classmethod(get_majority_call_number)
	
	def get_nativename_snpid2call_m(self, curs, ecotype_table, calls_table):
		"""
		2007-09-22
			only for ((nativename, stockparent),snpid)s with >1 calls (duplicated calls).
		"""
		sys.stderr.write("Checking inconsistent duplicate calls ...")
		old_strain_snp_pair = None
		call_counter_ls = [0]*11
		strain_snp_pair2call_number = {}
		duplicated_times = 0
		curs.execute("select e.nativename, e.stockparent, c.snpid, c.call1, c.call2 from %s e, %s c where e.id=c.ecotypeid order by nativename, stockparent, snpid"%(ecotype_table, calls_table))
		rows = curs.fetchall()
		no_of_distinct_pairs = 0
		no_of_duplicated_pairs = 0
		counter = 0
		for row in rows:
			nativename, stockparent, snpid, call, callhet = row
			nativename = nativename.upper()	#bug here. same name appears in >1 forms differing in cases
			call = call.upper()
			if callhet:
				callhet.upper()
				call = call+callhet
			call_number = nt2number[call]
			strain_snp_pair = ((nativename, stockparent), snpid)
			if old_strain_snp_pair == None:	#1st time
				old_strain_snp_pair = strain_snp_pair
				duplicated_times += 1
			elif strain_snp_pair != old_strain_snp_pair:
				if duplicated_times > 1:	#only for duplicated calls
					no_of_duplicated_pairs += 1
					majority_call_number = self.get_majority_call_number(call_counter_ls)
					strain_snp_pair2call_number[old_strain_snp_pair] = majority_call_number
				no_of_distinct_pairs += 1
				old_strain_snp_pair = strain_snp_pair
				call_counter_ls = [0]*11
				duplicated_times = 1
			else:
				duplicated_times += 1
			counter += 1
			if call_number !=0:	#dont' need NA
				call_counter_ls[call_number] += 1	#have to be put last cuz in the 2nd condition, call_counter_ls might need to be cleared up.
		#don't miss out the last one
		if duplicated_times>1:	#only for duplicated calls
			no_of_duplicated_pairs += 1
			majority_call_number = self.get_majority_call_number(call_counter_ls)
			strain_snp_pair2call_number[old_strain_snp_pair] = majority_call_number
		sys.stderr.write("%s total calls. %s/%s(duplicated distinct pairs/distinct pairs). Done.\n"%(counter, no_of_duplicated_pairs, no_of_distinct_pairs))
		return strain_snp_pair2call_number
	
	def fill_in_resolved_duplicated_calls(self, data_matrix, strain_id2index, snp_id2index, nativename2strain_id, nativename_snpid2call):
		"""
		2007-09-22
			replace the calls in data_matrix with resolved duplicated calls
		2007-10-01
			to make sure key in nativename_snpid2call appear in nativename2strain_id
		"""
		sys.stderr.write("Filling in resolved duplicated calls ...")
		for strain_snp_pair, call_number in nativename_snpid2call.iteritems():
			if strain_snp_pair[0] in nativename2strain_id:
				strain_id = nativename2strain_id[strain_snp_pair[0]]
				snpid = strain_snp_pair[1]
				data_matrix[strain_id2index[strain_id]][snp_id2index[snpid]] = call_number
		sys.stderr.write("Done.\n")
		return data_matrix
	
	def get_strain_id2other_info(self, curs, strain_id_list, strain_info_table):
		"""
		2007-09-22
		2007-12-13
			changes following strain_id = (strain_id, duplicate)
		2007-12-16
			add nativename
		"""
		sys.stderr.write("Getting strain_id2other_info ..m.")
		strain_id2other_info = {}
		for strain_id in strain_id_list:
			curs.execute("select e.id, e.latitude, e.longitude, e.nativename, e.stockparent, s.name, c.abbr from %s e, address a, site s, country c where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id and e.id=%s"%(strain_info_table, strain_id[0]))
			rows = curs.fetchall()
			id, latitude, longitude, nativename, stockparent, site_name, abbr = rows[0]
			strain_id2other_info[strain_id] = [latitude, longitude, nativename, stockparent, site_name, abbr]
		sys.stderr.write("Done.\n")
		return strain_id2other_info
	
	def write_data_matrix(self, data_matrix, output_fname, strain_id_list, snp_id_list, snp_id2acc, with_header_line, nt_alphabet, strain_id2acc=None, strain_id2category=None, rows_to_be_tossed_out=Set(), strain_id2other_info=None, discard_all_NA_strain=0, predefined_header_row=['strain', 'duplicate', 'latitude', 'longitude', 'nativename', 'stockparent', 'site', 'country']):
		"""
		2008-05-08
			defunct use write_data_matrix from pymodule
		2007-02-19
			if strain_id2acc is available, translate strain_id into strain_acc,
			if strain_id2category is available, add 'category'
		2007-02-25
			if one strain's SNP row is all NA, it'll be skipped
		2007-02-25
			add argument rows_to_be_tossed_out
		2007-09-23
			add discard_all_NA_strain
		2007-10-22
			add no_of_all_NA_rows
		2007-12-13
			add predefined_header_row
		2007-12-16
			add 'duplicate' into predefined_header_row
		"""
		sys.stderr.write("Writing data_matrix ...")
		no_of_all_NA_rows = 0
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		if with_header_line:
			header_row = [predefined_header_row[0]]
			if strain_id2category:
				header_row.append(predefined_header_row[1])
			if strain_id2other_info:
				no_of_fields = len(strain_id2other_info.values()[0])	#2007-12-13
				for i in range(no_of_fields):
					header_row.append(predefined_header_row[2+i])
			for snp_id in snp_id_list:
				header_row.append(snp_id2acc[snp_id])
			writer.writerow(header_row)
		for i in range(len(data_matrix)):
			if strain_id2acc:
				new_row = [strain_id2acc[strain_id_list[i]]]
			else:
				new_row = [strain_id_list[i]]
			if strain_id2category:
				new_row.append(strain_id2category[strain_id_list[i]])
			if strain_id2other_info:
				new_row += strain_id2other_info[strain_id_list[i]]
			if discard_all_NA_strain and sum(data_matrix[i]==0)==data_matrix.shape[1]:
				no_of_all_NA_rows += 1
				continue
			elif i not in rows_to_be_tossed_out:	#2007-02-25
				for j in data_matrix[i]:
					if nt_alphabet:
						j = number2nt[j]
					new_row.append(j)
				writer.writerow(new_row)
		del writer
		sys.stderr.write("%s all NA rows ."%no_of_all_NA_rows)
		sys.stderr.write("Done.\n")
	
	def sort_file(self, ofname):
		"""
		2007-02-19
		1. sort file1 into file2
		2. move file2 to file1
		"""
		sys.stderr.write("\tSorting %s"%ofname)
		commandline = 'sort %s > %s.post'%(ofname, ofname)
		exit_code = system_call(commandline)
		commandline = 'mv %s.post %s'%(ofname, ofname)
		exit_code = system_call(commandline)
		sys.stderr.write(".\n")
	
	def toss_rows_to_make_distance_matrix_NA_free(self, data_matrix):
		"""
		2007-02-25
		"""
		sys.stderr.write("Removing rows to make distance matrix NA free ...")
		NA_pair_list = []
		no_of_rows, no_of_columns = data_matrix.shape
		for i in range(no_of_rows):
			for j in range(i+1, no_of_rows):
				no_of_NA_pairs = 0
				for k in range(no_of_columns):
					if data_matrix[i,k]!=0 and data_matrix[j,k]!=0:
						break
					else:
						no_of_NA_pairs += 1
				if no_of_NA_pairs==no_of_columns:
					NA_pair_list.append([i,j])
		
		import networkx as nx
		g = nx.Graph()
		g.add_edges_from(NA_pair_list)
		vertex_list_to_be_deleted = self.find_smallest_vertex_set_to_remove_all_edges(g)
		if self.debug:
			print
			print "NA_pair_list:"
			print NA_pair_list
			print "vertex_list_to_be_deleted:"
			print vertex_list_to_be_deleted
		sys.stderr.write(" %s removed, done.\n"%(len(vertex_list_to_be_deleted)))
		return vertex_list_to_be_deleted
	
	def find_smallest_vertex_set_to_remove_all_edges(self, g):
		"""
		2007-02-25
		2007-04-17
			no recursive function, cuz it could reach the maximum recursion depth (whatever that is)
			
		"""
		vertex_with_max_degree = -1
		max_degree = 0
		for v in g:
			degree_of_v = g.degree(v)
			if degree_of_v > max_degree:
				max_degree = degree_of_v
				vertex_with_max_degree = v
		vertex_list_to_be_deleted = []
		while max_degree>0:	#to avoid empty-edge graph
			g.delete_node(vertex_with_max_degree)
			vertex_list_to_be_deleted.append(vertex_with_max_degree)
			max_degree = 0
			for v in g:
				degree_of_v = g.degree(v)
				if degree_of_v > max_degree:
					max_degree = degree_of_v
					vertex_with_max_degree = v
		return vertex_list_to_be_deleted
	
	def run(self):
		"""
		2008-05-08
			transpose everything if output_matrix_type=1 (bjarni's SNP matrix format)
		2007-02-19
			--db_connect
			--get_snp_id2index()
			--get_strain_id2index()
			--get_strain_id_info()
			--get_snp_id_info()
			--get_data_matrix()
			if self.toss_out_rows:
				--toss_rows_to_make_distance_matrix_NA_free()
					--find_smallest_vertex_set_to_remove_all_edges()
			--write_data_matrix()
			#--sort_file()
		2007-09-22
			for mysql_connection
				add get_nativename_snpid2call_m()
				add fill_in_resolved_duplicated_calls()
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		if self.db_connection_type==1:
			import MySQLdb
			#conn = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
			conn = MySQLdb.connect(db=self.dbname,host=self.hostname, user=self.user, passwd = self.passwd)
			curs = conn.cursor()
			snp_id2index, snp_id_list = self.get_snp_id2index_m(curs, self.input_table, self.snp_locus_table)
			strain_id2index, strain_id_list, nativename2strain_id, strain_id2acc, strain_id2category = self.get_strain_id2index_m(curs, self.input_table, self.strain_info_table, self.only_include_strains_with_GPS, self.resolve_duplicated_calls)
			
			#strain_id2acc, strain_id2category = self.get_strain_id_info_m(curs, strain_id_list, self.strain_info_table)
			snp_id2info = self.get_snp_id_info_m(curs, snp_id_list, self.snp_locus_table)
			if self.input_table == 'dbsnp.calls':
				from variation.src.FigureOut384IlluminaABMapping import get_snps_id2mapping
				snps_id2mapping = get_snps_id2mapping(self.hostname, dbname='dbsnp', user=self.user, passwd=self.passwd)
			else:
				snps_id2mapping = None
			data_matrix = self.get_data_matrix_m(curs, strain_id2index, snp_id2index, nt2number, self.input_table, self.need_heterozygous_call, snps_id2mapping)
			if self.resolve_duplicated_calls:
				nativename_snpid2call = self.get_nativename_snpid2call_m(curs, self.strain_info_table, self.input_table)
				data_matrix = self.fill_in_resolved_duplicated_calls(data_matrix, strain_id2index, snp_id2index, nativename2strain_id, nativename_snpid2call)
			if self.include_other_strain_info:
				strain_id2other_info = self.get_strain_id2other_info(curs, strain_id_list, self.strain_info_table)
			else:
				strain_id2other_info = None
		elif self.db_connection_type==2:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			snp_id2index, snp_id_list = self.get_snp_id2index(curs, self.input_table, self.snp_locus_table)
			strain_id2index, strain_id_list = self.get_strain_id2index(curs, self.input_table)
			
			strain_id2acc, strain_id2category = self.get_strain_id_info(curs, strain_id_list, self.strain_info_table)
			snp_id2info = self.get_snp_id_info(curs, snp_id_list, self.snp_locus_table)
			data_matrix = self.get_data_matrix(curs, strain_id2index, snp_id2index, nt2number, self.input_table, self.need_heterozygous_call)
			strain_id2other_info = None
		
		if self.toss_out_rows:
			rows_to_be_tossed_out = self.toss_rows_to_make_distance_matrix_NA_free(data_matrix)
			rows_to_be_tossed_out = Set(rows_to_be_tossed_out)
		else:
			rows_to_be_tossed_out = Set()
		
		#05/08/08
		if self.discard_all_NA_strain:
			from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
			rows_with_too_many_NAs_set, row_index2no_of_NAs = FilterStrainSNPMatrix.remove_rows_with_too_many_NAs(data_matrix, row_cutoff=1)
			rows_to_be_tossed_out.update(rows_with_too_many_NAs_set)
		
		strain_acc_list = [strain_id2acc[strain_id] for strain_id in strain_id_list]
		category_list = [strain_id2category[strain_id] for strain_id in strain_id_list]
		if self.output_matrix_type==1:
			#transpose everything
			data_matrix = num.array(data_matrix)
			data_matrix = num.transpose(data_matrix)
			
			header = ['Chromosomes', 'Positions'] + strain_acc_list
			chromosome_ls = []
			position_ls = []
			for snp_id, info in snp_id2info.iteritems():
				snp_name, chromosome, position = info
				chromosome_ls.append(chromosome)
				position_ls.append(position) 
			
			strain_acc_list = chromosome_ls
			category_list = position_ls
			cols_to_be_tossed_out = rows_to_be_tossed_out
			rows_to_be_tossed_out = None
			strain_id2other_info = None	#make up one
		else:
			header = ['strain', 'category']
			for snp_id, info in snp_id2info.iteritems():
				snp_name, chromosome, position = info
				header.append(snp_name)
			cols_to_be_tossed_out = None
		
		write_data_matrix(data_matrix, self.output_fname, header, strain_acc_list, category_list, rows_to_be_tossed_out=rows_to_be_tossed_out, \
					cols_to_be_tossed_out=cols_to_be_tossed_out, nt_alphabet=self.nt_alphabet,\
					strain_acc2other_info=strain_id2other_info, delimiter=self.delimiter)
		
		#self.write_data_matrix(data_matrix, self.output_fname, strain_id_list, snp_id_list, snp_id2acc, self.with_header_line,\
		#	self.nt_alphabet, strain_id2acc, strain_id2category, rows_to_be_tossed_out, strain_id2other_info, self.discard_all_NA_strain)

		#self.sort_file(self.output_fname)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = dbSNP2data
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	"""
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:o:s:n:g:y:mbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'dl324b-1'
	dbname = 'yhdb'
	schema = 'dbsnp'
	input_table = None
	output_fname = None
	strain_info_table = 'strain_info'
	snp_locus_table = 'snp_locus'
	organism = 'at'
	processing_bits = '10101101'
	mysql_connection = 0
	debug = 0
	report = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-i",):
			input_table = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-s",):
			strain_info_table = arg
		elif opt in ("-n",):
			snp_locus_table = arg
		elif opt in ("-g",):
			organism = arg
		elif opt in ("-y",):
			processing_bits = arg
		elif opt in ("-m",):
			mysql_connection = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if input_table and output_fname and hostname and dbname and schema:
		instance = dbSNP2data(hostname, dbname, schema, input_table, output_fname, \
			strain_info_table, snp_locus_table, organism, processing_bits, \
			mysql_connection, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
	"""