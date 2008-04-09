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
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/')))

import psycopg2, sys, getopt, csv, re
from annot.bin.codense.common import db_connect, org_short2long, org2tax_id
from variation.src.common import nt2number, number2nt
import Numeric as num
from sets import Set

def createEcotypeid2duplicate_view(curs, stock_db):
	"""
	2007-10-22
	"""
	curs.execute("create or replace view %s.ecotypeid2duplicate_view as select distinct ecotypeid, replicate from %s.calls order by ecotypeid, replicate"%(stock_db, stock_db))

def createNoGenotypingEcotypeView(curs, stock_db):
	"""
	2007-10-22
	"""
	curs.execute("create or replace view %s.no_genotyping_ecotype_view as select e.* from %s.ecotype e where not exists (select e1.id,c.id from %s.ecotype e1, %s.calls c where c.ecotypeid=e1.id and e1.id=e.id)"%(stock_db, stock_db, stock_db, stock_db))

def createNoGPSEcotypeView(curs, stock_db):
	"""
	2007-10-22
	"""
	curs.execute("create or replace view %s.no_gps_ecotype_view as select e.* from %s.ecotype e where latitude is null or longitude is null"%(stock_db, stock_db))

def createGenotypingAllNAEcotypeTable(curs, stock_db, table_name='genotyping_all_na_ecotype', commit=0):
	"""
	2007-10-22
		create a table to store all ecotypeid which have been genotyped (in table calls) but all results are NA.
	"""
	curs.execute("select distinct ecotypeid, replicate, call1 from %s.calls"%(stock_db))
	rows = curs.fetchall()
	genotype_run2call_ls = {}
	for row in rows:
		ecotypeid, duplicate, call1 = row
		key_pair = (ecotypeid, duplicate)
		if key_pair not in genotype_run2call_ls:
			genotype_run2call_ls[key_pair] = []
		genotype_run2call_ls[key_pair].append(call1)
	
	genotyping_all_na_ecotypeid_duplicate_ls = []
	for key_pair, call_ls in genotype_run2call_ls.iteritems():
		if len(call_ls)==1 and (call_ls[0]=='N' or call_ls[0]=='n'):
			genotyping_all_na_ecotypeid_duplicate_ls.append(key_pair)
	
	if commit:
		curs.execute("create table %s.%s(id	integer primary key auto_increment,\
			ecotypeid	integer,\
			duplicate	integer)"%(stock_db, table_name))
		for key_pair in genotyping_all_na_ecotypeid_duplicate_ls:
			ecotypeid, duplicate = key_pair
			curs.execute("insert into %s.%s(ecotypeid, duplicate) values (%s, %s)"%(stock_db, table_name, ecotypeid, duplicate))
	
	return genotyping_all_na_ecotypeid_duplicate_ls, genotype_run2call_ls


def get_no_genotyping_ecotypeid_set(curs, stock_db, no_genotyping_ecotype_view_name='no_genotyping_ecotype_view'):
	from sets import Set
	no_genotyping_ecotypeid_set = Set()
	curs.execute("select id from %s.%s"%(stock_db, no_genotyping_ecotype_view_name))
	rows = curs.fetchall()
	for row in rows:
		no_genotyping_ecotypeid_set.add(row[0])
	return no_genotyping_ecotypeid_set

def get_no_gps_ecotypeid_set(curs, stock_db, no_gps_ecotype_view_name='no_gps_ecotype_view'):
	from sets import Set
	no_gps_ecotypeid_set = Set()
	curs.execute("select id from %s.%s"%(stock_db, no_gps_ecotype_view_name))
	rows = curs.fetchall()
	for row in rows:
		no_gps_ecotypeid_set.add(row[0])
	return no_gps_ecotypeid_set
	

def createTableStructureToGroupEcotypeid(curs, stock_db, no_genotyping_ecotypeid_set, no_gps_ecotypeid_set,\
			genotyping_all_na_ecotypeid_duplicate_ls, ecotypeid2duplicate_view_name='ecotypeid2duplicate_view', \
			nativename_stkparent2tg_ecotypeid_table='nativename_stkparent2tg_ecotypeid', \
			ecotype_duplicate2tg_ecotypeid_table='ecotype_duplicate2tg_ecotypeid', commit=0):
	"""
	2007-10-22
		map (nativename, stockparent) to an ecotypeid (AKA tg_ecotypeid)
		map all (ecotype,duplicate) to that ecotypeid
	2007-12-16
		mysql is case-insensitive, python is case-sensitive.
		key_pair (nativename, stockparent) should be case-insensitive. like 'Kz-9' and 'KZ-9'
		uppercase nativename and stockparent
	"""
	sys.stderr.write("Getting nativename_stkparent2ecotypeid_duplicate_ls...")
	nativename_stkparent2ecotypeid_duplicate_ls = {}
	ecotypeid2nativename_stockparent = {}
	curs.execute("select e.nativename, e.stockparent, ed.ecotypeid, ed.replicate from %s.ecotype e, %s.%s ed where ed.ecotypeid=e.id order by nativename, stockparent, ecotypeid, replicate"%(stock_db, stock_db, ecotypeid2duplicate_view_name))
	rows = curs.fetchall()
	for row in rows:
		nativename, stockparent, ecotypeid, duplicate = row
		ecotypeid2nativename_stockparent[ecotypeid] = (nativename, stockparent)	#2007-12-16 before upper()
		nativename = nativename.upper()
		if stockparent:
			stockparent = stockparent.upper()
		key_pair = (nativename, stockparent)
		if key_pair not in nativename_stkparent2ecotypeid_duplicate_ls:
			nativename_stkparent2ecotypeid_duplicate_ls[key_pair] = []
		nativename_stkparent2ecotypeid_duplicate_ls[key_pair].append((ecotypeid, duplicate))
	sys.stderr.write("Done.\n")
	sys.stderr.write("Constructing nativename_stkparent2tg_ecotypeid ecotype_duplicate2tg_ecotypeid...\n")
	nativename_stkparent2tg_ecotypeid = {}
	ecotype_duplicate2tg_ecotypeid = {}
	from sets import Set
	genotyping_all_na_ecotypeid_duplicate_set = Set(genotyping_all_na_ecotypeid_duplicate_ls)
	no_of_solid_mappings = 0
	no_of_mappings_with_data_but_no_gps = 0
	no_of_mappings_with_gps_but_no_data = 0
	no_of_worst_random_mappings = 0
	for key_pair, ecotypeid_duplicate_ls in nativename_stkparent2ecotypeid_duplicate_ls.iteritems():
		tg_ecotypeid_quality_pair = None
		ecotypeid_ls_with_data_but_no_gps = []
		ecotypeid_ls_with_gps_but_no_data = []
		for pair in ecotypeid_duplicate_ls:
			if pair[0] not in no_genotyping_ecotypeid_set and pair[0] not in no_gps_ecotypeid_set and pair not in genotyping_all_na_ecotypeid_duplicate_set:
				tg_ecotypeid_quality_pair = (pair[0], 'solid')
				no_of_solid_mappings += 1
				break
			elif pair[0] not in no_genotyping_ecotypeid_set and pair not in genotyping_all_na_ecotypeid_duplicate_set:
				ecotypeid_ls_with_data_but_no_gps.append(pair[0])
			elif pair[0] not in no_gps_ecotypeid_set:
				ecotypeid_ls_with_gps_but_no_data.append(pair[0])
		if tg_ecotypeid_quality_pair==None:	#use ecotypeid with data
			if ecotypeid_ls_with_data_but_no_gps:
				tg_ecotypeid_quality_pair = (ecotypeid_ls_with_data_but_no_gps[0], 'with_data_but_no_gps')
				no_of_mappings_with_data_but_no_gps += 1
			elif ecotypeid_ls_with_gps_but_no_data:
				tg_ecotypeid_quality_pair = (ecotypeid_ls_with_gps_but_no_data[0], 'with_gps_but_no_data')
				no_of_mappings_with_gps_but_no_data += 1
			else:
				tg_ecotypeid_quality_pair = (ecotypeid_duplicate_ls[0][0], 'worst_random')
				no_of_worst_random_mappings += 1
		nativename_stkparent2tg_ecotypeid[key_pair] = tg_ecotypeid_quality_pair
		for pair in ecotypeid_duplicate_ls:
			ecotype_duplicate2tg_ecotypeid[pair] = tg_ecotypeid_quality_pair[0]
	no_of_total_mappings = float(len(nativename_stkparent2tg_ecotypeid))
	sys.stderr.write("\t%s(%s) solid mappings\n"%(no_of_solid_mappings, no_of_solid_mappings/no_of_total_mappings))
	sys.stderr.write("\t%s(%s) mappings_with_data_but_no_gps\n"%(no_of_mappings_with_data_but_no_gps, no_of_mappings_with_data_but_no_gps/no_of_total_mappings))
	sys.stderr.write("\t%s(%s) mappings_with_gps_but_no_data\n"%(no_of_mappings_with_gps_but_no_data, no_of_mappings_with_gps_but_no_data/no_of_total_mappings))
	sys.stderr.write("\t%s(%s) worst_random_mappings\n"%(no_of_worst_random_mappings, no_of_worst_random_mappings/no_of_total_mappings))
	sys.stderr.write("Done.\n")
	if commit:
		sys.stderr.write("Submitting to db...")
		curs.execute("create table %s.%s(id	integer primary key auto_increment,\
			nativename	varchar(50),\
			stockparent	varchar(10),\
			tg_ecotypeid	integer,\
			quality	varchar(50))"%(stock_db, nativename_stkparent2tg_ecotypeid_table))
		for key_pair, tg_ecotypeid_quality_pair in nativename_stkparent2tg_ecotypeid.iteritems():
			tg_ecotypeid, quality = tg_ecotypeid_quality_pair
			nativename, stockparent = ecotypeid2nativename_stockparent[tg_ecotypeid]
			curs.execute("insert into %s.%s(nativename, stockparent, tg_ecotypeid, quality) values ('%s', '%s', %s, '%s')"%(stock_db, nativename_stkparent2tg_ecotypeid_table, nativename, stockparent, tg_ecotypeid, quality))
		
		curs.execute("create table %s.%s(id	integer primary key auto_increment,\
			ecotypeid	integer,\
			duplicate	integer,\
			tg_ecotypeid	integer)"%(stock_db, ecotype_duplicate2tg_ecotypeid_table))
		
		for pair, tg_ecotypeid in ecotype_duplicate2tg_ecotypeid.iteritems():
			curs.execute("insert into %s.%s(ecotypeid, duplicate, tg_ecotypeid) values (%s, %s, %s)"%(stock_db, ecotype_duplicate2tg_ecotypeid_table, pair[0], pair[1], tg_ecotypeid))
		sys.stderr.write("Done.\n")
	return nativename_stkparent2ecotypeid_duplicate_ls, nativename_stkparent2tg_ecotypeid, ecotype_duplicate2tg_ecotypeid

def get_ecotypeid2duplicate_times(curs, stock_db, ecotypeid2duplicate_view_name='ecotypeid2duplicate_view'):
	"""
	2007-10-23
		not finished yet
	"""
	curs.execute("select ecotypeid, count(duplicate) from %s.%s group by ecotypeid"%(stock_db, ecotypeid2duplicate_view_name))


def get_nn_sp_duplicated_time2ecotype_duplicate_ls_ls(nativename_stkparent2ecotypeid_duplicate_ls):
	"""
	2007-10-23
	"""
	nn_sp_duplicated_time2ecotype_duplicate_ls_ls = {}
	for nativename_stkparent, ecotypeid_duplicate_ls in nativename_stkparent2ecotypeid_duplicate_ls.iteritems():
		nn_sp_duplicated_time = len(ecotypeid_duplicate_ls)
		if nn_sp_duplicated_time not in nn_sp_duplicated_time2ecotype_duplicate_ls_ls:
			nn_sp_duplicated_time2ecotype_duplicate_ls_ls[nn_sp_duplicated_time] = []
		nn_sp_duplicated_time2ecotype_duplicate_ls_ls[nn_sp_duplicated_time].append(ecotypeid_duplicate_ls)
	return nn_sp_duplicated_time2ecotype_duplicate_ls_ls

def convert_nn_sp_duplicated_time2ecotype_duplicate_ls_ls_2matrix(nn_sp_duplicated_time2ecotype_duplicate_ls_ls):
	"""
	2007-10-23
	"""
	nn_sp_duplicated_time_ls = nn_sp_duplicated_time2ecotype_duplicate_ls_ls.keys()
	nn_sp_duplicated_time_ls.sort()
	nn_sp_duplicated_time2index = {}
	import numpy
	m = numpy.zeros([len(nn_sp_duplicated_time_ls), 2], numpy.integer)
	for i in range(len(nn_sp_duplicated_time_ls)):
		nn_sp_duplicated_time = nn_sp_duplicated_time_ls[i]
		m[i,0] = nn_sp_duplicated_time
		m[i,1] = len(nn_sp_duplicated_time2ecotype_duplicate_ls_ls[nn_sp_duplicated_time])
		nn_sp_duplicated_time2index[nn_sp_duplicated_time] = i
	return m, nn_sp_duplicated_time2index

def get_nn_sp_duplicated_time_cross_another_set(nn_sp_duplicated_time2ecotype_duplicate_ls_ls, nn_sp_duplicated_time2index, set_to_be_intersected=None, is_ecotypeid_duplicate_pair_in_set=0):
	"""
	2007-10-23
	"""
	import numpy
	x_dim = len(nn_sp_duplicated_time2index)
	y_dim = max(nn_sp_duplicated_time2index.keys())	#the possible maximum number of occcurrences
	m = numpy.zeros([x_dim, y_dim], numpy.integer)
	new_ecotype_duplicate_ls_ls = []	#list of nn_sp_duplicated_time2ecotype_duplicate_ls_ls.values() - those in the set_to_be_intersected
	for nn_sp_duplicated_time, ecotype_duplicate_ls_ls in nn_sp_duplicated_time2ecotype_duplicate_ls_ls.iteritems():
		key2times = {}
		for ecotypeid_duplicate_ls in ecotype_duplicate_ls_ls:
			new_ecotype_duplicate_ls = []
			for ecotypeid_duplicate in ecotypeid_duplicate_ls:
				if is_ecotypeid_duplicate_pair_in_set:	#set_key might be (ecotypeid, duplicate) or just ecotypeid
					set_key = ecotypeid_duplicate
				else:
					set_key = ecotypeid_duplicate[0]
				if set_to_be_intersected:	#if intersected (either no gps, or all-NA), shall be removed
					if set_key not in set_to_be_intersected:
						new_ecotype_duplicate_ls.append(ecotypeid_duplicate)
						continue
				key = ecotypeid_duplicate[0]	#this is always ecotypeid
				if key not in key2times:
					key2times[key] = 0
				key2times[key] += 1
			new_ecotype_duplicate_ls_ls.append(new_ecotype_duplicate_ls)
		times2key_ls = {}
		for key, times in key2times.iteritems():
			if times not in times2key_ls:
				times2key_ls[times] = []
			times2key_ls[times].append(key)
		x_index = nn_sp_duplicated_time2index[nn_sp_duplicated_time]
		for times, key_ls in times2key_ls.iteritems():
			y_index = times-1
			m[x_index, y_index] = len(key_ls)
	reduced_nn_sp_duplicated_time2ecotype_duplicate_ls_ls = {}
	for ecotype_duplicate_ls in new_ecotype_duplicate_ls_ls:
		nn_sp_duplicated_time = len(ecotype_duplicate_ls)
		if nn_sp_duplicated_time not in reduced_nn_sp_duplicated_time2ecotype_duplicate_ls_ls:
			reduced_nn_sp_duplicated_time2ecotype_duplicate_ls_ls[nn_sp_duplicated_time] = []
		reduced_nn_sp_duplicated_time2ecotype_duplicate_ls_ls[nn_sp_duplicated_time].append(ecotype_duplicate_ls)
	
	return m, reduced_nn_sp_duplicated_time2ecotype_duplicate_ls_ls


"""
2007-10-22
functions to understand the duplication structure inside the database

hostname='localhost'
dbname='stock20071008'
dbname = 'stock'
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname)
curs = conn.cursor()
stock_db = dbname
createEcotypeid2duplicate_view(curs, stock_db)

createNoGenotypingEcotypeView(curs, stock_db)

createNoGPSEcotypeView(curs, stock_db)

genotyping_all_na_ecotypeid_duplicate_ls, genotype_run2call_ls = createGenotypingAllNAEcotypeTable(curs, stock_db, \
	table_name='genotyping_all_na_ecotype', commit=1)

no_genotyping_ecotypeid_set = get_no_genotyping_ecotypeid_set(curs, stock_db, no_genotyping_ecotype_view_name='no_genotyping_ecotype_view')

no_gps_ecotypeid_set = get_no_gps_ecotypeid_set(curs, stock_db, no_gps_ecotype_view_name='no_gps_ecotype_view')

nativename_stkparent2ecotypeid_duplicate_ls, nativename_stkparent2tg_ecotypeid, ecotype_duplicate2tg_ecotypeid = \
createTableStructureToGroupEcotypeid(curs, stock_db, no_genotyping_ecotypeid_set, no_gps_ecotypeid_set, \
genotyping_all_na_ecotypeid_duplicate_ls, ecotypeid2duplicate_view_name='ecotypeid2duplicate_view', \
nativename_stkparent2tg_ecotypeid_table='nativename_stkparent2tg_ecotypeid', \
ecotype_duplicate2tg_ecotypeid_table='ecotype_duplicate2tg_ecotypeid', commit=1)


##2007-10-23 to output tables
output_fname = 'script/variation/doc/paper/tables_figures.tex'
outf = open(output_fname, 'w')
nn_sp_duplicated_time2ecotype_duplicate_ls_ls = get_nn_sp_duplicated_time2ecotype_duplicate_ls_ls(nativename_stkparent2ecotypeid_duplicate_ls)

nn_sp_m, nn_sp_duplicated_time2index = convert_nn_sp_duplicated_time2ecotype_duplicate_ls_ls_2matrix(nn_sp_duplicated_time2ecotype_duplicate_ls_ls)

##cross nn_sp_duplicated_time2ecotype_duplicate_ls_ls to no set. just count how many duplicated ecotypeid's
m_X_ecotypeid_duplicated_times, useless_dict = get_nn_sp_duplicated_time_cross_another_set(nn_sp_duplicated_time2ecotype_duplicate_ls_ls, nn_sp_duplicated_time2index, set_to_be_intersected=None, is_ecotypeid_duplicate_pair_in_set=0)


##cross nn_sp_duplicated_time2ecotype_duplicate_ls_ls to sets of ecotypeid_duplicate with all-NA data
##reorganize nn_sp_duplicated_time2ecotype_duplicate_ls_ls
from sets import Set
m_X_genotyping_all_na_ecotypeid_duplicate_ls, reduced_1_nn_sp_duplicated_time2ecotype_duplicate_ls_ls= get_nn_sp_duplicated_time_cross_another_set(nn_sp_duplicated_time2ecotype_duplicate_ls_ls, nn_sp_duplicated_time2index, set_to_be_intersected=Set(genotyping_all_na_ecotypeid_duplicate_ls), is_ecotypeid_duplicate_pair_in_set=1)

#convert the reorganized nn_sp_duplicated_time2ecotype_duplicate_ls_ls into matrix
reduced_1_m, reduced_1_nn_sp_duplicated_time2index = convert_nn_sp_duplicated_time2ecotype_duplicate_ls_ls_2matrix(reduced_1_nn_sp_duplicated_time2ecotype_duplicate_ls_ls)


##cross nn_sp_duplicated_time2ecotype_duplicate_ls_ls to no-gps-ecotypeid set
m_X_no_gps_duplicated_times, reduced_2_nn_sp_duplicated_time2ecotype_duplicate_ls_ls = get_nn_sp_duplicated_time_cross_another_set(nn_sp_duplicated_time2ecotype_duplicate_ls_ls, nn_sp_duplicated_time2index, set_to_be_intersected=no_gps_ecotypeid_set, is_ecotypeid_duplicate_pair_in_set=0)

#convert the reorganized nn_sp_duplicated_time2ecotype_duplicate_ls_ls into matrix
reduced_2_m, reduced_2_nn_sp_duplicated_time2index = convert_nn_sp_duplicated_time2ecotype_duplicate_ls_ls_2matrix(reduced_2_nn_sp_duplicated_time2ecotype_duplicate_ls_ls)


from pymodule.latex import outputMatrixInLatexTable
import numpy
t1 = numpy.concatenate([nn_sp_m, m_X_ecotypeid_duplicated_times, m_X_genotyping_all_na_ecotypeid_duplicate_ls, m_X_no_gps_duplicated_times], 1)
caption = ' Cross (nativename,stkparent) to ecotypeid duplicated times, ecotypeid-with-all-NA,  ecotypeid-with-no-gps'
table_label = 't_nn_sp_e_1'
header_ls = [(2,'nativename stkparent'),(m_X_ecotypeid_duplicated_times.shape[1], 'ecotypeid duplicate'), (m_X_genotyping_all_na_ecotypeid_duplicate_ls.shape[1], 'ecotypeid-with-all-NA'), (m_X_no_gps_duplicated_times.shape[1], 'ecotypeid-with-no-gps')]
outf.write(outputMatrixInLatexTable(t1, caption, table_label, header_ls))
outf.flush()



caption = ' (nativename,stkparent) duplicates after ecotypeid-with-all-NA removal'
table_label = 't_nn_sp_e_2'
header_ls = [(2,'nativename stkparent')]
outf.write(outputMatrixInLatexTable(reduced_1_m, caption, table_label, header_ls))
outf.flush()

caption = ' (nativename,stkparent) duplicates after ecotypeid-with-no-gps removal'
table_label = 't_nn_sp_e_3'
header_ls = [(2,'nativename stkparent')]
outf.write(outputMatrixInLatexTable(reduced_2_m, caption, table_label, header_ls))
outf.flush()

####optional here. if 0 is not removed, codes in next block still run
##delete the values with 0 ecotypeid_duplicate associated
del reduced_1_nn_sp_duplicated_time2ecotype_duplicate_ls_ls[0]
##convert the reorganized nn_sp_duplicated_time2ecotype_duplicate_ls_ls into matrix
reduced_1_m, reduced_1_nn_sp_duplicated_time2index = convert_nn_sp_duplicated_time2ecotype_duplicate_ls_ls_2matrix(reduced_1_nn_sp_duplicated_time2ecotype_duplicate_ls_ls)


#cross reduced_1_nn_sp_duplicated_time2ecotype_duplicate_ls_ls to no set. just count how many duplicated ecotypeid's
reduced_1_m_X_ecotypeid_duplicated_times, useless_dict = get_nn_sp_duplicated_time_cross_another_set(reduced_1_nn_sp_duplicated_time2ecotype_duplicate_ls_ls, reduced_1_nn_sp_duplicated_time2index, set_to_be_intersected=None, is_ecotypeid_duplicate_pair_in_set=0)

##cross reduced_1_nn_sp_duplicated_time2ecotype_duplicate_ls_ls to no-gps-ecotypeid set
reduced_1_m_X_no_gps_duplicated_times, reduced_11_nn_sp_duplicated_time2ecotype_duplicate_ls_ls = get_nn_sp_duplicated_time_cross_another_set(reduced_1_nn_sp_duplicated_time2ecotype_duplicate_ls_ls, reduced_1_nn_sp_duplicated_time2index, set_to_be_intersected=no_gps_ecotypeid_set, is_ecotypeid_duplicate_pair_in_set=0)

t4 = numpy.concatenate([reduced_1_m, reduced_1_m_X_ecotypeid_duplicated_times, reduced_1_m_X_no_gps_duplicated_times], 1)
caption = ' (nativename,stkparent) duplicates after ecotypeid-with-all-NA removal, cross it to ecotypeid duplicated times, ecotypeid-with-no-gps'
table_label = 't_nn_sp_e_4'
header_ls = [(2,'nativename stkparent all-NA removal'), (reduced_1_m_X_ecotypeid_duplicated_times.shape[1], 'ecotypeid duplicate'), (reduced_1_m_X_no_gps_duplicated_times.shape[1], 'ecotypeid-with-no-gps')]
outf.write(outputMatrixInLatexTable(t4, caption, table_label, header_ls))
outf.flush()

#convert the reorganized nn_sp_duplicated_time2ecotype_duplicate_ls_ls into matrix
reduced_11_m, reduced_11_nn_sp_duplicated_time2index = convert_nn_sp_duplicated_time2ecotype_duplicate_ls_ls_2matrix(reduced_11_nn_sp_duplicated_time2ecotype_duplicate_ls_ls)


caption = ' (nativename,stkparent) duplicates after all-NA and no-gps removal'
table_label = 't_nn_sp_e_5'
header_ls = [(2,'nativename stkparent')]
outf.write(outputMatrixInLatexTable(reduced_11_m, caption, table_label, header_ls))
outf.flush()



####final table
#no_genotyping_ecotypeid_set
#no_gps_ecotypeid_set

genotyping_all_na_ecotypeid_duplicate_set = Set(genotyping_all_na_ecotypeid_duplicate_ls)
genotyping_all_na_ecotypeid_ls = [row[0] for row in genotyping_all_na_ecotypeid_duplicate_ls]
genotyping_all_na_ecotypeid_set = Set(genotyping_all_na_ecotypeid_ls)
m = numpy.zeros([1,8], numpy.integer)

header_ls = ['all-NA runs', 'all-NA ecotypeid', 'no-genotyping', 'no-gps', 'all-NA and no-genotyping', 'all-NA and no-gps', 'no-genotyping and no-gps', 'intersect all 3']
m[0,0] = len(genotyping_all_na_ecotypeid_duplicate_set)
m[0,1] = len(genotyping_all_na_ecotypeid_set)
m[0,2] = len(no_genotyping_ecotypeid_set)
m[0,3] = len(no_gps_ecotypeid_set)
m[0,4] = len(genotyping_all_na_ecotypeid_set&no_genotyping_ecotypeid_set)
m[0,5] = len(genotyping_all_na_ecotypeid_set&no_gps_ecotypeid_set)
m[0,6] = len(no_genotyping_ecotypeid_set&no_gps_ecotypeid_set)
m[0,7] = len(genotyping_all_na_ecotypeid_set&no_genotyping_ecotypeid_set&no_gps_ecotypeid_set)

caption = 'overlapping between all weirdos'
table_label = 't_nn_sp_e_6'
outf.write(outputMatrixInLatexTable(m, caption, table_label, header_ls))
outf.flush()

"""
