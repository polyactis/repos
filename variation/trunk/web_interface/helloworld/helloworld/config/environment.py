"""Pylons environment configuration"""
import os

from pylons import config

import helloworld.lib.app_globals as app_globals
import helloworld.lib.helpers
from helloworld.config.routing import make_map
from mako.lookup import TemplateLookup
import helloworld.model as model

#import sys
#sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

def load_environment(global_conf, app_conf):
	"""Configure the Pylons environment via the ``pylons.config``
	object
	"""
	# Pylons paths
	root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	paths = dict(root=root,
				 controllers=os.path.join(root, 'controllers'),
				 static_files=os.path.join(root, 'public'),
				 templates=[os.path.join(root, 'templates')])

	# Initialize config with the basic options
	config.init_app(global_conf, app_conf, package='helloworld',
					template_engine='mako', paths=paths)

	config['routes.map'] = make_map()
	config['pylons.g'] = app_globals.Globals()
	config['pylons.h'] = helloworld.lib.helpers
	
	# Customize templating options via this variable
	#tmpl_options = config['buffet.template_options']

	# CONFIGURATION OPTIONS HERE (note: all config options will override
	# any Pylons config options)
	
	#2008-10-05 Use the strict behaviour of the template context object
	#config['pylons.strict_c'] = True
	
	#2008-10-05 Create the Mako TemplateLookup, with the default auto-escaping. it doesn't work though.
	config['pylons.g'].mako_lookup = TemplateLookup(
	directories=paths['templates'],
	module_directory=os.path.join(app_conf['cache_dir'], 'templates'),
	input_encoding='utf-8', output_encoding='utf-8',
	imports=['from webhelpers.html import escape'],
	default_filters=['escape'])
	
	#2008-10-05 setup the database connection
	drivername = config['app_conf']['drivername']
	hostname = config['app_conf']['hostname']
	dbname = config['app_conf']['dbname']
	schema = config['app_conf']['schema']
	db_user = config['app_conf']['db_user']
	db_passwd = config['app_conf']['db_passwd']
	pool_recycle = int(config['app_conf']['pool_recycle'])
	#model.setup()
	
	model.db = model.Stock_250kDB.Stock_250kDB(drivername=drivername, username=db_user, password=db_passwd, \
							hostname=hostname, database=dbname, schema=schema, pool_recycle=pool_recycle)
	model.db.setup(create_tables=False)

	model.genome_db = model.GenomeDB.GenomeDatabase(drivername=drivername, username=db_user, password=db_passwd, \
					hostname=hostname, database='genome', schema=schema, pool_recycle=pool_recycle)
	
	#from variation.src import dbsnp
	#snp_db = dbsnp.DBSNP(drivername=drivername, username=db_user, password=db_passwd, \
	#				hostname=hostname, database='dbsnp', schema=schema, pool_recycle=pool_recycle)
	
	#from variation.src import StockDB
	#stock_db = StockDB.StockDB(drivername=drivername, username=db_user, password=db_passwd, \
	#				hostname=hostname, database='stock', schema=schema, pool_recycle=pool_recycle)
	
	model.at_db = model.AtDB.AtDB(drivername=drivername, username=db_user, password=db_passwd, \
					hostname=hostname, database='at', schema=schema, pool_recycle=pool_recycle)
	"""
	for entity in entities:
		if entity.__module__==db.__module__:	#entity in the same module
			entity.metadata = metadata
			#using_table_options_handler(entity, schema=self.schema)
	"""
	model.genome_db.setup(create_tables=False)
	#snp_db.setup(create_tables=False)
	#stock_db.setup(create_tables=False)
	model.at_db.setup(create_tables=False)
	
	from variation.src.DrawSNPRegion import DrawSNPRegion
	def dealWithGeneAnnotation():
		gene_annotation_picklef = '/Network/Data/250k/tmp-yh/at_gene_model_pickelf'
		DrawSNPRegion_ins = DrawSNPRegion(db_user=db_user, db_passwd=db_passwd, hostname=hostname, database=dbname,\
										input_fname='/tmp/dumb', output_dir='/tmp', debug=0)
		gene_annotation = DrawSNPRegion_ins.dealWithGeneAnnotation(gene_annotation_picklef)
		return gene_annotation

	model.gene_annotation = dealWithGeneAnnotation()
	
	#2008-11-05 a dictionary to link two tables of types in order for cross-linking between pages of DisplayTopSNPTestRM and ScoreRankHistogram/DisplayResultsGene
	model.CandidateGeneTopSNPTestRMType_id_min_distance2ScoreRankHistogramType_id = {}
	ScoreRankHistogramType = model.Stock_250kDB.ScoreRankHistogramType
	CandidateGeneTopSNPTestRMType = model.Stock_250kDB.CandidateGeneTopSNPTestRMType
	
	rows = model.db.metadata.bind.execute("select s.id as sid, s.min_distance, s.call_method_id, c.id as cid from %s s, %s c where s.null_distribution_type_id=c.null_distribution_type_id and\
					s.results_type=c.results_type and s.get_closest=c.get_closest and s.min_MAF=c.min_MAF and \
					s.allow_two_sample_overlapping=c.allow_two_sample_overlapping"%\
					(ScoreRankHistogramType.table.name, CandidateGeneTopSNPTestRMType.table.name))	#2008-1-8 temporarily set call_method_id=17 cuz CandidateGeneTopSNPTestRMType doesn't include call_method_id
	
	for row in rows:
		key_tuple = (row.cid, row.min_distance, row.call_method_id)
		model.CandidateGeneTopSNPTestRMType_id_min_distance2ScoreRankHistogramType_id[key_tuple] = row.sid
	
	# 2009-4-10 takes too long in individual request, put here. used in Accession.py
	from variation.src.common import map_perlegen_ecotype_name2accession_id
	model.ecotype_name2accession_id = map_perlegen_ecotype_name2accession_id(model.db.metadata.bind)