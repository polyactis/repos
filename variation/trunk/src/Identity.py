#!/usr/bin/env python
"""
Usage: Identity.py [OPTIONS] -i XXX -j XXX -p XXX

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock20071008(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	cpickle_fname, the cPickle file containing strain_pair_ls, geno_dist_ls, geo_dist_ls
	-j ...,	data_matrix_fname, the data matrix file where the cPickle data is based on
	-m ..., max_identity_geno_dist, 0.0(default)
	-g ..., geo_distance_window_size, 100 (km) (default)
		specify this latex_output_fname if you want output
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help
	

Examples:
	Identity.py -i ~/script/variation/stock20071008/data_d110_c0_5_strain_pair_geno_dist_geo_dist_ls -j ~/script/variation/stock20071008/data_d110_c0_5.tsv

	Identity.py -i ~/script/variation/stock20071008/data_d110_c0_5_strain_pair_geno_dist_0_0_0_4_geo_dist_ls -j ~/script/variation/stock20071008/data_d110_c0_5.tsv
	
	Identity.py -i ~/script/variation/genotyping/250ksnp/data/data_250k_strain_pair_geno_dist_0_0_0_4_geo_dist_ls -j ~/script/variation/genotyping/250ksnp/data/data_250k.tsv

Description:
	program to how identity pairs density are distributed along geographic distance in different regions.
"""
import getopt

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/test/python')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/test/python')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:j:m:g:o:brh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock20071008'
	schema = 'dbsnp'
	cpickle_fname = None
	data_matrix_fname = None
	max_identity_geno_dist = 0.0
	geo_distance_window_size = 100
	diff_type_to_be_outputted = 5
	latex_output_fname = None
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
			cpickle_fname = arg
		elif opt in ("-m",):
			max_identity_geno_dist = float(arg)
		elif opt in ("-g",):
			geo_distance_window_size = float(arg)
		elif opt in ("-j",):
			data_matrix_fname = arg
		elif opt in ("-o",):
			latex_output_fname = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	from annot.bin.codense.common import db_connect, form_schema_tables
	hostname='dl324b-1'
	dbname='yhdb'
	schema = 'dbsnp'
	hostname='zhoudb'
	dbname='graphdb'
	schema = 'dbsnp'
	#conn, curs = db_connect(hostname, dbname, schema)
	
	hostname='localhost'
	dbname='stock20070829'
	import MySQLdb
	conn0 = MySQLdb.connect(db=dbname,host=hostname)
	curs0 = conn0.cursor()
	
	hostname='localhost'
	dbname='stock20071008'
	import MySQLdb
	conn = MySQLdb.connect(db=dbname,host=hostname)
	curs = conn.cursor()
	
	from variation.src.misc import *
	
	if cpickle_fname and data_matrix_fname:
		import cPickle
		#inf = open('%s_strain_pair_geno_dist_0_0_0_3_geo_dist_ls'%os.path.splitext(cpickle_fname)[0], 'r')
		inf = open(cpickle_fname, 'r')
		strain_pair_ls, geno_dist_ls, geo_dist_ls = cPickle.load(inf)
		del inf
		
		europe_lon_span = [-12,50]
		norame_lon_span = [-130,-60]
		centralasia_lon_span = [60,90]
		japan_lon_span = [130, 150]
		
		span_ls = [europe_lon_span, norame_lon_span]
		span_label_ls = ['europe', 'nor america']
		
		ecotypeid2pos = get_ecotypeid2pos(curs, 'ecotype')
		
		span_pair_count_ls = get_span_pair_count_ls(data_matrix_fname, ecotypeid2pos, span_ls)
		print span_pair_count_ls
		matrix_of_counts_by_region_x_geo_dist, row_label_ls, col_label_ls, col_index_ls = group_identity_pair_according_to_region_distance(strain_pair_ls, geno_dist_ls, geo_dist_ls, ecotypeid2pos, span_ls, span_label_ls, geo_distance_window_size=geo_distance_window_size, max_identity_geno_dist=max_identity_geno_dist, span_pair_count_ls=span_pair_count_ls)
		output_fname_prefix = '%s_iden_%s_geo_window_%s'%(os.path.splitext(cpickle_fname)[0], max_identity_geno_dist, geo_distance_window_size)
		output_fname_prefix = output_fname_prefix.replace('.', '_')
		draw_table_as_bar_chart(matrix_of_counts_by_region_x_geo_dist, row_label_ls, col_label_ls, col_index_ls, output_fname_prefix)
	else:
		print __doc__
		sys.exit(2)
