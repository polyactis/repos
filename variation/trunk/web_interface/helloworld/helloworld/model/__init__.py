import os,sys
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
#import sqlalchemy
from elixir import setup_all, session, metadata, entities
from variation.src import Stock_250kDB
#from variation.src.Stock_250kDB import GeneListType
#from variation.src.Stock_250kDB import *
drivername = 'mysql'
hostname = 'banyan.usc.edu'
dbname = 'stock_250k'
schema = ""
db_user = 'yh'
db_passwd = 'yh324'
pool_recycle = 3600

db = Stock_250kDB.Stock_250kDB(drivername=drivername, username=db_user, password=db_passwd, \
				hostname=hostname, database=dbname, schema=schema, pool_recycle=pool_recycle)

from transfac.src import GenomeDB

genome_db = GenomeDB.GenomeDatabase(drivername=drivername, username=db_user, password=db_passwd, \
				hostname=hostname, database='genome', schema=schema, pool_recycle=pool_recycle)
"""
for entity in entities:
	if entity.__module__==db.__module__:	#entity in the same module
		entity.metadata = metadata
		#using_table_options_handler(entity, schema=self.schema)
"""
db.setup(create_tables=False)
genome_db.setup(create_tables=False)

import os,sys
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from variation.src.DrawSNPRegion import DrawSNPRegion
def dealWithGeneAnnotation():
	gene_annotation_picklef = '/Network/Data/250k/tmp-yh/at_gene_model_pickelf'
	DrawSNPRegion_ins = DrawSNPRegion(db_user=db_user, db_passwd=db_passwd, hostname=hostname, database=dbname,\
									input_fname='/tmp/dumb', output_dir='/tmp', debug=0)
	gene_annotation = DrawSNPRegion_ins.dealWithGeneAnnotation(gene_annotation_picklef)
	return gene_annotation

gene_annotation = dealWithGeneAnnotation()