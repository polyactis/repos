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
	
	dbSNP2data.py -z localhost -d stock20070829 -i calls -o stock20070829/data.tsv -s ecotype -n snps -m
	
	dbSNP2data.py -z localhost -d stock20070919 -i calls -o /tmp/data.tsv -s ecotype -n snps -m
	
	#to see how many all-NA strains were discarded
	dbSNP2data.py -z localhost -d stock20070919 -i calls -o /tmp/data10101100.tsv -s ecotype -n snps -m -y 10101100
	
	#all strains, but resolve duplicated (nativename, stockparent)s
	dbSNP2data.py -z localhost -d stock20070919 -i calls -o /tmp/data00101100.tsv -s ecotype -n snps -m -y 00101100
	
	#output data for all ecotypeid. however duplicated calls with same ecotypeid are imputed randomly
	dbSNP2data.py -z localhost -d stock20070919 -i calls -o /tmp/data00001100.tsv -s ecotype -n snps -m -y 00001100

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

def createEcotypeid2duplicate_view(curs, stock_db):
	"""
	2007-10-22
	"""
	curs.execute("create or replace view %s.ecotypeid2duplicate_view as select distinct ecotypeid, duplicate from %s.calls order by ecotypeid, duplicate"%(stock_db, stock_db))

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
	curs.execute("select distinct ecotypeid, duplicate, call1 from %s.calls"%(stock_db))
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
	

def createTableStructureToGroupEcotypeid(curs, stock_db, no_genotyping_ecotypeid_set, no_gps_ecotypeid_set, genotyping_all_na_ecotypeid_duplicate_ls, ecotypeid2duplicate_view_name='ecotypeid2duplicate_view', nativename_stkparent2tg_ecotypeid_table='nativename_stkparent2tg_ecotypeid', ecotype_duplicate2tg_ecotypeid_table='ecotype_duplicate2tg_ecotypeid', commit=0):
	"""
	2007-10-22
		map (nativename, stockparent) to an ecotypeid (AKA tg_ecotypeid)
		map all (ecotype,duplicate) to that ecotypeid
	"""
	sys.stderr.write("Getting nativename_stkparent2ecotypeid_duplicate_ls...")
	nativename_stkparent2ecotypeid_duplicate_ls = {}
	curs.execute("select e.nativename, e.stockparent, ed.ecotypeid, ed.duplicate from %s.ecotype e, %s.%s ed where ed.ecotypeid=e.id order by nativename, stockparent, ecotypeid, duplicate"%(stock_db, stock_db, ecotypeid2duplicate_view_name))
	rows = curs.fetchall()
	for row in rows:
		nativename, stockparent, ecotypeid, duplicate = row
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
			nativename, stockparent = key_pair
			tg_ecotypeid, quality = tg_ecotypeid_quality_pair
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
hostname='localhost'
dbname='stock20071008'
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname)
curs = conn.cursor()
stock_db = dbname
createEcotypeid2duplicate_view(curs, stock_db)

createNoGenotypingEcotypeView(curs, stock_db)

createNoGPSEcotypeView(curs, stock_db)

genotyping_all_na_ecotypeid_duplicate_ls, genotype_run2call_ls = createGenotypingAllNAEcotypeTable(curs, stock_db, table_name='genotyping_all_na_ecotype', commit=1)

no_genotyping_ecotypeid_set = get_no_genotyping_ecotypeid_set(curs, stock_db, no_genotyping_ecotype_view_name='no_genotyping_ecotype_view')

no_gps_ecotypeid_set = get_no_gps_ecotypeid_set(curs, stock_db, no_gps_ecotype_view_name='no_gps_ecotype_view')

nativename_stkparent2ecotypeid_duplicate_ls, nativename_stkparent2tg_ecotypeid, ecotype_duplicate2tg_ecotypeid = createTableStructureToGroupEcotypeid(curs, stock_db, no_genotyping_ecotypeid_set, no_gps_ecotypeid_set, genotyping_all_na_ecotypeid_duplicate_ls, ecotypeid2duplicate_view_name='ecotypeid2duplicate_view', nativename_stkparent2tg_ecotypeid_table='nativename_stkparent2tg_ecotypeid', ecotype_duplicate2tg_ecotypeid_table='ecotype_duplicate2tg_ecotypeid', commit=0)


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
		"""
		sys.stderr.write("Getting snp_id2index ..m.")
		snp_id2index = {}
		snp_id_list = []
		curs.execute("select distinct i.snpid, s.chromosome, s.position from %s i, %s s where i.snpid=s.id order by chromosome, position"%(input_table, snp_locus_table))
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
		"""
		sys.stderr.write("Getting strain_id2index ..m.")
		strain_id2index = {}
		strain_id_list = []
		nativename2strain_id = {}
		if only_include_strains_with_GPS==1:
			curs.execute("select distinct d.ecotypeid, s.nativename, s.stockparent from %s d, %s s where d.ecotypeid=s.id and s.latitude is not null and s.longitude is not null  order by ecotypeid, nativename, stockparent"%(input_table, strain_info_table))
		elif only_include_strains_with_GPS==2:	#2007-10-01 north american samples
			curs.execute("select distinct d.ecotypeid, s.nativename, s.stockparent from %s d, %s s where d.ecotypeid=s.id and s.latitude is not null and s.longitude is not null and s.longitude<-60 and s.longitude>-130 order by ecotypeid, nativename, stockparent"%(input_table, strain_info_table))
		elif only_include_strains_with_GPS==3:
			curs.execute("select distinct d.ecotypeid, s.nativename, s.stockparent from %s d, %s s, batch_ecotype be, batch b where b.batchname='192' and b.id=be.batchid and s.id=be.ecotypeid and d.ecotypeid=s.id and s.latitude is not null and s.longitude is not null  order by ecotypeid, nativename, stockparent"%(input_table, strain_info_table))
		else:
			curs.execute("select distinct d.ecotypeid, s.nativename, s.stockparent from %s d, %s s where d.ecotypeid=s.id order by ecotypeid, nativename, stockparent"%(input_table, strain_info_table))
		rows = curs.fetchall()
		for row in rows:
			strain_id, nativename, stockparent = row
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
		"""
		sys.stderr.write("Getting strain_id_info ..m.")
		strain_id2acc = {}
		strain_id2category = {}
		for strain_id in strain_id_list:
			curs.execute("select id, nativename from %s where id=%s"%(strain_info_table, strain_id))
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
		"""
		sys.stderr.write("Getting snp_id_info ..m.")
		snp_id2acc = {}
		for snp_id in snp_id_list:
			curs.execute("select snpid from %s where id=%s"%(snp_locus_table, snp_id))
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
		"""
		sys.stderr.write("Getting data_matrix ..m.\n")
		data_matrix = num.zeros([len(strain_id2index), len(snp_id2index)])
		curs.execute("select ecotypeid, snpid, call1, callhet from %s"%(input_table))
		rows = curs.fetchmany(5000)
		counter = 0
		while rows:
			for row in rows:
				strain_id, snp_id, call, callhet = row
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
				sys.stderr.write('%s%s'%('\x08'*20, counter))
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
			
		"""
		sys.stderr.write("Getting strain_id2other_info ..m.")
		strain_id2other_info = {}
		for strain_id in strain_id_list:
			curs.execute("select e.id, e.latitude, e.longitude, e.stockparent, s.name, c.abbr from %s e, address a, site s, country c where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id and e.id=%s"%(strain_info_table, strain_id))
			rows = curs.fetchall()
			id, latitude, longitude, stockparent, site_name, abbr = rows[0]
			strain_id2other_info[id] = [latitude, longitude, stockparent, site_name, abbr]
		sys.stderr.write("Done.\n")
		return strain_id2other_info
	
	def write_data_matrix(self, data_matrix, output_fname, strain_id_list, snp_id_list, snp_id2acc, with_header_line, nt_alphabet, strain_id2acc=None, strain_id2category=None, rows_to_be_tossed_out=Set(), strain_id2other_info=None, discard_all_NA_strain=0):
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
		"""
		sys.stderr.write("Writing data_matrix ...")
		no_of_all_NA_rows = 0
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		if with_header_line:
			header_row = ['strain']
			if strain_id2category:
				header_row.append('nativename')
			if strain_id2other_info:
				header_row.append('latitude')
				header_row.append('longitude')
				header_row.append('stockparent')
				header_row.append('site')
				header_row.append('country')
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
