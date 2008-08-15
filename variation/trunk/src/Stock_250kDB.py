#!/usr/bin/env python
"""
Examples:
	#setup database in postgresql
	Stock_250kDB.py -v postgres -u crocea -z localhost -d graphdb -k public
	
	#setup database in mysql
	Stock_250kDB.py -u yh
	
Description:
	2008-07-09
	This is a wrapper for the stock_250k database, build on top of elixir. supposed to supercede the table definitions in mysql.sql.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from sqlalchemy.engine.url import URL
from elixir import Unicode, DateTime, String, Integer, UnicodeText, Text, Boolean, Float
from elixir import Entity, Field, using_options, using_table_options
from elixir import OneToMany, ManyToOne, ManyToMany
from elixir import setup_all, session, metadata, entities
from elixir.options import using_table_options_handler	#using_table_options() can only work inside Entity-inherited class.
from sqlalchemy import UniqueConstraint
from sqlalchemy.schema import ThreadLocalMetaData

from datetime import datetime

from pymodule.db import ElixirDB

class README(Entity):
	#2008-08-07
	title = Field(String(2000))
	description = Field(String(60000))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='readme')
	using_table_options(mysql_engine='InnoDB')

class Phenotype(Entity):
	ecotype_id = Field(Integer, nullable=False)
	value = Field(Float)
	replicate = Field(Integer)
	phenotype_method = ManyToOne('Phenotype', colname='method_id', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='phenotype')
	using_table_options(mysql_engine='InnoDB')

class PhenotypeAvg(Entity):
	ecotype_id = Field(Integer, nullable=False)
	value = Field(Float)
	stdev = Field(Float)	
	sample_size = Field(Integer)
	ready_for_publication = Field(Integer, default=0)
	phenotype_method = ManyToOne('Phenotype', colname='method_id', ondelete='CASCADE', onupdate='CASCADE')
	readme = ManyToOne("README", colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	transformed_value = Field(Float)
	using_options(tablename='phenotype_avg')
	using_table_options(mysql_engine='InnoDB')

class GeneListType(Entity):
	short_name = Field(String(256), unique=True)
	original_filename = Field(String(760), unique=True)	#for unique constraint in mysql, max key length is 767 bytes
	description = Field(String(8192))
	gene_list = OneToMany('GeneList')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_list_type')
	using_table_options(mysql_engine='InnoDB')

class GeneList(Entity):
	gene_id = Field(Integer)
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE', inverse='gene_list')
	original_name = Field(String(128))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_list')
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('gene_id', 'list_type_id'))


class Snps(Entity):
	name = Field(String(200), unique=True, nullable = False)
	chromosome = Field(Integer)
	position = Field(Integer)
	end_position = Field(Integer)
	allele1 = Field(String(1))
	allele2 = Field(String(2))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snps')
	using_table_options(mysql_engine='InnoDB')

class SnpsContext(Entity):
	snp = ManyToOne('Snps', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	disp_pos = Field(Integer)
	gene_id = Field(Integer)
	gene_strand = Field(String(1))
	left_or_right = Field(String(200))
	disp_pos_comment = Field(String(2000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snps_context')
	using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('snps_id', 'gene_id'))

class CallMethod(Entity):
	short_name = Field(String(20))
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	comment = Field(String(8000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_method')
	using_table_options(mysql_engine='InnoDB')
	
class PhenotypeMethod(Entity):
	short_name = Field(String(20))
	only_first_96 = Field(Boolean, default=0)
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	comment = Field(String(8000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	data_type = Field(String(200))
	transformation_description = Field(String(8000))
	using_options(tablename='phenotype_method')
	using_table_options(mysql_engine='InnoDB')

class AnalysisMethod(Entity):
	short_name = Field(String(60))
	method_description = Field(String(8000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='analysis_method')
	using_table_options(mysql_engine='InnoDB')
	
class ResultsMethodType(Entity):
	short_name = Field(String(30), unique=True)
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='results_method_type')
	using_table_options(mysql_engine='InnoDB')

class ResultsMethod(Entity):
	short_name = Field(String(30), unique=True)
	filename = Field(String(1000), unique=True)
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	comment = Field(String(8000))
	phenotype_method = ManyToOne('PhenotypeMethod', colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	results_method_type = ManyToOne('ResultsMethodType', colname='results_method_type_id', ondelete='CASCADE', onupdate='CASCADE')
	analysis_method = ManyToOne('AnalysisMethod', colname='analysis_method_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='results_method')
	using_table_options(mysql_engine='InnoDB')

class Results(Entity):
	snp = ManyToOne('Snps', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	results_method = ManyToOne('ResultsMethod', colname='results_method_id', ondelete='CASCADE', onupdate='CASCADE')
	score = Field(Float)
	using_options(tablename='results')
	using_table_options(mysql_engine='InnoDB')

class CandidateGeneRankSumTestResult(Entity):
	"""
	2008-07-17
	"""
	results_method = ManyToOne('ResultsMethod', colname='results_method_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	statistic = Field(Float)
	pvalue = Field(Float)
	min_distance = Field(Integer)
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	max_pvalue_per_gene = Field(Integer)
	comment = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_rank_sum_test_result')
	using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('results_method_id', 'list_type_id'))

class ResultsByGene(Entity):
	"""
	2008-07-19
	"""
	results_method = ManyToOne('ResultsMethod', colname='results_method_id', ondelete='CASCADE', onupdate='CASCADE')
	snp = ManyToOne('Snps', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	gene_id = Field(Integer)
	disp_pos = Field(Integer)
	score = Field(Float)
	rank = Field(Float)
	using_options(tablename='results_by_gene')
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_method_id', 'snps_id', 'gene_id'))

class SnpsQC(Entity):
	"""
	2008-08-07
	"""
	snp = ManyToOne('Snps', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	min_probability = Field(Float)
	max_call_info_mismatch_rate = Field(Float)
	snps_name = Field(String(200))
	tg_snps_name = Field(String(200))
	NA_rate = Field(Float)
	no_of_NAs = Field(Integer)
	no_of_totals = Field(Integer)
	relative_NA_rate = Field(Float)
	relative_no_of_NAs = Field(Integer)
	relative_no_of_totals = Field(Integer)
	mismatch_rate = Field(Float)
	no_of_mismatches = Field(Integer)
	no_of_non_NA_pairs = Field(Integer)
	qc_method = ManyToOne("QCMethod", colname='qc_method_id', ondelete='CASCADE', onupdate='CASCADE')
	call_method = ManyToOne("CallMethod", colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	readme = ManyToOne("README", colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snps_qc')
	using_table_options(mysql_engine='InnoDB')

class QCMethod(Entity):
	short_name = Field(String(30), unique=True)
	data1_type = Field(String(30), nullable=False)
	data2_type = Field(String(30), nullable=False)
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	comment = Field(String(8000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='qc_method')
	using_table_options(mysql_engine='InnoDB')

class ContaminantType(Entity):
	"""
	2008-08-07
	"""
	short_name = Field(String(100), unique=True)
	description = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='contaminant_type')
	using_table_options(mysql_engine='InnoDB')

class ArrayInfo(Entity):
	name = Field(String(40))
	filename = Field(String(1000))
	original_filename = Field(String(1000))
	description = Field(String(2000))
	ecotype_id = Field(Integer)
	maternal_ecotype_id = Field(Integer)
	paternal_ecotype_id = Field(Integer)
	contaminant_type = ManyToOne("ContaminantType", colname='contaminant_type_id', ondelete='CASCADE', onupdate='CASCADE')
	strain_id = Field(Integer)
	md5sum = Field(String(100))
	experimenter = Field(String(200))
	samples = Field(String(20))
	dna_amount = Field(String(20))
	S260_280 = Field(Float)
	total_vol = Field(String(20))
	hyb_vol = Field(String(20))
	seed_source = Field(String(100))
	method_name = Field(String(250))
	readme = ManyToOne("README", colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='array_info')
	using_table_options(mysql_engine='InnoDB')
	
class CallInfo(Entity):
	filename = Field(String(1000))
	description = Field(String(2000))
	array = ManyToOne("ArrayInfo", colname='array_id', ondelete='CASCADE', onupdate='CASCADE')
	NA_rate = Field(Float)
	readme = ManyToOne("README", colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	call_method = ManyToOne('CallMethod', colname='method_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_info')
	using_table_options(mysql_engine='InnoDB')

class CallQC(Entity):
	call_info = ManyToOne("CallInfo", colname='call_info_id', ondelete='CASCADE', onupdate='CASCADE')
	min_probability = Field(Float)
	ecotype_id = Field(Integer)
	duplicate = Field(Integer)
	tg_ecotype_id = Field(Integer)
	tg_duplicate = Field(Integer)
	NA_rate = Field(Float)
	no_of_NAs = Field(Integer)
	no_of_totals = Field(Integer)
	mismatch_rate = Field(Float)
	no_of_mismatches = Field(Integer)
	no_of_non_NA_pairs = Field(Integer)
	call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	readme = ManyToOne("README", colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	qc_method = ManyToOne("QCMethod", colname='qc_method_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_qc')
	using_table_options(mysql_engine='InnoDB')

class Probes(Entity):
	snp = ManyToOne('Snps', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	seq = Field(String(25))
	chromosome = Field(Integer)
	position = Field(Integer)
	allele = Field(String(1))
	strand = Field(String(20))
	xpos = Field(Integer)
	ypos = Field(Integer)
	direction = Field(String(20))
	gene = Field(String(50))
	RNA = Field(String(50))
	tu = Field(String(50))
	flank = Field(String(50))
	expressedClones = Field(Float)
	totalClones = Field(Integer)
	multiTranscript = Field(String(50))
	LerDel = Field(String(50))
	LerCopy = Field(Integer)
	LerSNPdelL = Field(Integer)
	LerSNPdelR = Field(Integer)
	LerSNPpos = Field(Integer)
	promoter = Field(Boolean)
	utr5 = Field(Boolean)
	utr3 = Field(Boolean)
	intron = Field(Boolean)
	intergenic = Field(Boolean)
	downstream = Field(Boolean)
	cda = Field(Boolean)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='probes')
	using_table_options(mysql_engine='InnoDB')

class Stock_250kDB(ElixirDB):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ):['localhost', 'z', 1, 'hostname of the db server', ],\
							('database', 1, ):['stock_250k', 'd', 1, '',],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('username', 1, ):[None, 'u', 1, 'database username',],\
							('password', 1, ):[None, 'p', 1, 'database password', ],\
							('port', 0, ):[None, 'o', 1, 'database port number'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-07-09
		"""
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		if getattr(self, 'schema', None):	#for postgres
			for entity in entities:
				using_table_options_handler(entity, schema=self.schema)
		
		#metadata = ThreadLocalMetaData()
		metadata.bind = self._url
		setup_all(create_tables=True)	#create_tables=True causes setup_all to call elixir.create_all(), which in turn calls metadata.create_all()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Stock_250kDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	
	
	rows = GeneListType.query.all()
	for row in rows:
		print row.gene_list[0].list_type_id
	