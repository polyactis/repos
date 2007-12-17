#!/usr/bin/env python
"""
Usage: dbSNP2data.py [OPTIONS] -i INPUT_TABLE -o OUTPUT_FILE

Option:
	-z ..., --hostname=...	the hostname, dl324b-1(default)
	-d ..., --dbname=...	the database name, yhdb(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	input table
	-o ...,	output file
	-s ...,	strain_info table, 'strain_info'(default), 'ecotype'
	-n ...,	snp_locus_table, 'snp_locus'(default), 'snps'
	-y ...,	processing bits to control which processing step should be turned on.
		default is 10101101. for what each bit stands, see Description.
	-m,	mysql connection, change dbname and hostname
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	dbSNP2data.py -i justin_data -o justin_data.csv -r
	
	dbSNP2data.py -i justin_data -o justin_data.csv -r -t
	
	dbSNP2data.py -i calls -o /tmp/chicago.data.y -s ecotype -n snps -m
	
	dbSNP2data.py -z localhost -d stock20071008 -i calls -o stock20071008/data.tsv -s ecotype -n snps -m
	
	dbSNP2data.py -z localhost -d stock20071008 -i calls -o /tmp/data.tsv -s ecotype -n snps -m
	
	#to see how many all-NA strains were discarded
	dbSNP2data.py -z localhost -d stock20071008 -i calls -o /tmp/data10101100.tsv -s ecotype -n snps -m -y 10101100
	
	#all strains, but resolve duplicated (nativename, stockparent)s
	dbSNP2data.py -z localhost -d stock20071008 -i calls -o /tmp/data00101100.tsv -s ecotype -n snps -m -y 00101100
	
	#output data for all ecotypeid. however duplicated calls with same ecotypeid are imputed randomly
	dbSNP2data.py -z localhost -d stock20071008 -i calls -o /tmp/data00001100.tsv -s ecotype -n snps -m -y 00001100

	#output 250k SNP data (not to resolve duplicated calls, only SNPs shared with 149SNP)
	dbSNP2data.py -z localhost -d stock20071008 -i calls_250k -o /tmp/data_250k.tsv -s ecotype -n snps_250k -m -y 10001101 -r

Description:
	output SNP data from database schema
	
	definition of each bit in processing_bits (0=off, 1=on), default is 10101101.
	1. 1: only include strains with GPS info, 2: north american strains only, 3: 2010's 192 strains
	2. include columns of other strain info (latitude, longitude, stockparent, site, country)
	3. resolve duplicated calls (unique constraint on (nativename, stockparent))
	4. toss out rows to make distance matrix NA free
	5. need heterozygous call
	6. with header line
	7. use alphabet to represent nucleotide, not number
	8. discard strains with all-NA data
	
	you can specify the bits up to the one you want to change and omit the rest. i.e.
	-y 11
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import psycopg, sys, getopt, csv, re
from codense.common import db_connect, org_short2long, org2tax_id
from common import nt2number, number2nt
import Numeric as num
from sets import Set


class dbSNP2data:
	"""
	2007-02-19
	"""
	def __init__(self, hostname='dl324b-1', dbname='yhdb', schema='dbsnp', input_table=None, \
		output_fname=None, strain_info_table='strain_info', snp_locus_table='snp_locus', \
		organism='hs', processing_bits='10101101', mysql_connection=0, debug=0, report=0):
		"""
		2007-02-25
			add argument toss_out_rows
		2007-07-11
			add mysql_connection
		2007-09-23
			use processing_bits to control processing steps
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_table = input_table
		self.output_fname = output_fname
		self.strain_info_table = strain_info_table
		self.snp_locus_table = snp_locus_table
		#self.tax_id = org2tax_id(org_short2long(organism))
		self.processing_bits = processing_bits
		
		#below are all default values
		processing_bits_ls = [1,0,1,0,1,1,0,1]
		
		for i in range(len(processing_bits)):
			processing_bits_ls[i] = int(processing_bits[i])
		#now pass all values
		self.only_include_strains_with_GPS,\
		self.include_other_strain_info,\
		self.resolve_duplicated_calls,\
		self.toss_out_rows,\
		self.need_heterozygous_call,\
		self.with_header_line,\
		self.nt_alphabet,\
		self.discard_all_NA_strain = processing_bits_ls
		
		self.mysql_connection = int(mysql_connection)
		self.debug = int(debug)
		self.report = int(report)
	
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
			#curs.execute("select distinct s1.id, s1.chromosome, s1.position from %s s2, %s s1 where s1.chromosome=s2.chromosome and s1.position=s2.position order by chromosome, position"%('snps', snp_locus_table))
			curs.execute("select distinct s1.id, s1.chromosome, s1.position from %s s1"%(snp_locus_table))
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
		"""
		sys.stderr.write("Getting strain_id2index ..m.")
		strain_id2index = {}
		strain_id_list = []
		nativename2strain_id = {}
		common_sql_string = "select distinct d.ecotypeid, d.duplicate, s.nativename, s.stockparent from %s d, %s s"%(input_table, strain_info_table)
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
			strain_id, duplicate, nativename, stockparent = row
			strain_id = (strain_id, duplicate)
			nativename = nativename.upper()
			if resolve_duplicated_calls:
				key_pair = (nativename, stockparent)
				if key_pair not in nativename2strain_id:
					nativename2strain_id[key_pair] = strain_id
					strain_id_list.append(strain_id)
					strain_id2index[strain_id] = len(strain_id2index)
			else:
				strain_id_list.append(strain_id)
				strain_id2index[strain_id] = len(strain_id2index)
		sys.stderr.write("Done.\n")
		return strain_id2index, strain_id_list, nativename2strain_id
	
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
		2007-07-11
			mysql version of get_strain_id_info
		2007-09-13
			replace "name, nativename" with "id, name"
		2007-09-20
			replace 'name' with 'nativename'
		2007-12-13
			changes following strain_id = (strain_id, duplicate)
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
		"""
		sys.stderr.write("Getting snp_id_info ..m.")
		snp_id2acc = {}
		for snp_id in snp_id_list:
			if snp_locus_table == 'snps':
				curs.execute("select snpid from %s where id=%s"%(snp_locus_table, snp_id))
			elif snp_locus_table == 'snps_250k':
				curs.execute("select snpacc from %s where id=%s"%(snp_locus_table, snp_id))
			else:
				sys.stderr.write("Error: SNP table '%s' not supported.\n"%(snp_locus_table))
				sys.exit(3)
			rows = curs.fetchall()
			acc = rows[0][0]
			snp_id2acc[snp_id] = acc
		sys.stderr.write("Done.\n")
		return snp_id2acc
	
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
	
	def get_data_matrix_m(self, curs, strain_id2index, snp_id2index, nt2number, input_table, need_heterozygous_call):
		"""
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

		"""
		sys.stderr.write("Getting data_matrix ..m.\n")
		data_matrix = num.zeros([len(strain_id2index), len(snp_id2index)])
		if input_table == 'calls':
			common_sql_string = "select ecotypeid, duplicate, snpid, call1, callhet from %s"%(input_table)
		elif input_table == 'calls_250k':
			common_sql_string = "select ecotypeid, duplicate, snpid, snpcall from %s"%(input_table)
		else:
			sys.stderr.write("Error: SNP call table '%s' not supported.\n"%(input_table))
			sys.exit(3)
		snp_counter = 0
		for snp_id in snp_id2index:
			snp_counter += 1
			counter = 0
			if self.report:
				sys.stderr.write("%s\tSNP %s=%s"%('\x08'*80, snp_counter, snp_id))
			curs.execute("%s where snpid=%s"%(common_sql_string, snp_id))
			rows = curs.fetchmany(5000)
			while rows:
				for row in rows:
					if input_table == 'calls':
						strain_id, duplicate, snp_id, call, callhet = row
					elif input_table == 'calls_250k':
						strain_id, duplicate, snp_id, call = row
						callhet = None
					strain_id = (strain_id, duplicate)
					call = call.upper()
					if callhet:
						callhet.upper()	#2007-09-20	just in case
						call = call+callhet
					if strain_id in  strain_id2index and snp_id in snp_id2index:	#2007-03-20
						call_number = nt2number[call]
						if need_heterozygous_call:
							data_matrix[strain_id2index[strain_id], snp_id2index[snp_id]] = call_number
						elif call_number<=4:	#single letter or NA
							data_matrix[strain_id2index[strain_id], snp_id2index[snp_id]] = call_number
					counter += 1
				rows = curs.fetchmany(5000)
				if self.report:
					sys.stderr.write("%s\tSNP %s=%s\t%s"%('\x08'*80, snp_counter, snp_id, counter))
		sys.stderr.write("Done.\n")
		return data_matrix
	
	def get_majority_call_number(self, call_counter_ls):
		"""
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
		curs.execute("select e.nativename, e.stockparent, c.snpid, c.call1, c.callhet from %s e, %s c where e.id=c.ecotypeid order by nativename, stockparent, snpid"%(ecotype_table, calls_table))
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
		"""
		sys.stderr.write("Getting strain_id2other_info ..m.")
		strain_id2other_info = {}
		for strain_id in strain_id_list:
			curs.execute("select e.id, e.latitude, e.longitude, e.stockparent, s.name, c.abbr from %s e, address a, site s, country c where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id and e.id=%s"%(strain_info_table, strain_id[0]))
			rows = curs.fetchall()
			id, latitude, longitude, stockparent, site_name, abbr = rows[0]
			strain_id2other_info[strain_id] = [latitude, longitude, stockparent, site_name, abbr]
		sys.stderr.write("Done.\n")
		return strain_id2other_info
	
	def write_data_matrix(self, data_matrix, output_fname, strain_id_list, snp_id_list, snp_id2acc, with_header_line, nt_alphabet, strain_id2acc=None, strain_id2category=None, rows_to_be_tossed_out=Set(), strain_id2other_info=None, discard_all_NA_strain=0, predefined_header_row=['strain', 'nativename', 'latitude', 'longitude', 'stockparent', 'site', 'country']):
		"""
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
		if self.mysql_connection:
			import MySQLdb
			#conn = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
			conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
			curs = conn.cursor()
			snp_id2index, snp_id_list = self.get_snp_id2index_m(curs, self.input_table, self.snp_locus_table)
			strain_id2index, strain_id_list, nativename2strain_id = self.get_strain_id2index_m(curs, self.input_table, self.strain_info_table, self.only_include_strains_with_GPS, self.resolve_duplicated_calls)
			
			strain_id2acc, strain_id2category = self.get_strain_id_info_m(curs, strain_id_list, self.strain_info_table)
			snp_id2acc = self.get_snp_id_info_m(curs, snp_id_list, self.snp_locus_table)
			data_matrix = self.get_data_matrix_m(curs, strain_id2index, snp_id2index, nt2number, self.input_table, self.need_heterozygous_call)
			if self.resolve_duplicated_calls:
				nativename_snpid2call = self.get_nativename_snpid2call_m(curs, self.strain_info_table, self.input_table)
				data_matrix = self.fill_in_resolved_duplicated_calls(data_matrix, strain_id2index, snp_id2index, nativename2strain_id, nativename_snpid2call)
			if self.include_other_strain_info:
				strain_id2other_info = self.get_strain_id2other_info(curs, strain_id_list, self.strain_info_table)
			else:
				strain_id2other_info = None
		else:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			snp_id2index, snp_id_list = self.get_snp_id2index(curs, self.input_table, self.snp_locus_table)
			strain_id2index, strain_id_list = self.get_strain_id2index(curs, self.input_table)
			
			strain_id2acc, strain_id2category = self.get_strain_id_info(curs, strain_id_list, self.strain_info_table)
			snp_id2acc = self.get_snp_id_info(curs, snp_id_list, self.snp_locus_table)
			data_matrix = self.get_data_matrix(curs, strain_id2index, snp_id2index, nt2number, self.input_table, self.need_heterozygous_call)
			strain_id2other_info = None
		
		if self.toss_out_rows:
			rows_to_be_tossed_out = self.toss_rows_to_make_distance_matrix_NA_free(data_matrix)
			rows_to_be_tossed_out = Set(rows_to_be_tossed_out)
		else:
			rows_to_be_tossed_out = Set()
		if not self.resolve_duplicated_calls:	#if not to resolve duplicated calls, use strain_id_list to print 1st column
			strain_id2acc = None
		self.write_data_matrix(data_matrix, self.output_fname, strain_id_list, snp_id_list, snp_id2acc, self.with_header_line,\
			self.nt_alphabet, strain_id2acc, strain_id2category, rows_to_be_tossed_out, strain_id2other_info, self.discard_all_NA_strain)

		#self.sort_file(self.output_fname)

if __name__ == '__main__':
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
