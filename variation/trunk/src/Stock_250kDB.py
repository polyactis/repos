#!/usr/bin/env python
"""
Examples:
	#setup database in postgresql
	Stock_250kDB.py -v postgres -u crocea -z localhost -d graphdb -k public
	
	#setup database in mysql
	Stock_250kDB.py -u yh -z papaya.usc.edu
	
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
from elixir import Unicode, DateTime, String, Integer, UnicodeText, Text, Boolean, Float, Binary
from elixir import Entity, Field, using_options, using_table_options
from elixir import OneToMany, ManyToOne, ManyToMany
from elixir import setup_all, session, metadata, entities
from elixir.options import using_table_options_handler	#using_table_options() can only work inside Entity-inherited class.
from sqlalchemy import UniqueConstraint, create_engine
from sqlalchemy.schema import ThreadLocalMetaData, MetaData
from sqlalchemy.orm import scoped_session, sessionmaker

from datetime import datetime

from pymodule import SNPData
from pymodule.db import ElixirDB

__session__ = scoped_session(sessionmaker(autoflush=False, transactional=False))
#__metadata__ = ThreadLocalMetaData()	#2008-10 not good for pylon

__metadata__ = MetaData()

class README(Entity):
	#2008-08-07
	title = Field(String(2000))
	description = Field(String(60000))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='readme', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class Phenotype(Entity):
	ecotype_id = Field(Integer, nullable=False)
	value = Field(Float)
	replicate = Field(Integer)
	phenotype_method = ManyToOne('Phenotype', colname='method_id', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='phenotype', metadata=__metadata__, session=__session__)
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
	using_options(tablename='phenotype_avg', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class BiologyCategory(Entity):
	#2008-08-21
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	gene_list_type_ls = OneToMany("GeneListType")
	phenotype_method_ls = OneToMany("PhenotypeMethod")
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='biology_category', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class GeneListType(Entity):
	short_name = Field(String(256))
	biology_category = ManyToOne("BiologyCategory", colname='biology_category_id', ondelete='CASCADE', onupdate='CASCADE')
	type = ManyToOne("GeneListSuperType", colname='super_type_id', ondelete='CASCADE', onupdate='CASCADE')
	original_filename = Field(String(760), unique=True)	#for unique constraint in mysql, max key length is 767 bytes
	description = Field(String(8192))
	gene_list = OneToMany('GeneList')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_list_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('short_name', 'super_type_id'))

class GeneListSuperType(Entity):
	#2008-08-22
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	list_type_ls = OneToMany(GeneListType)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_list_super_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class GeneList(Entity):
	gene_id = Field(Integer)
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE', inverse='gene_list')
	original_name = Field(String(128))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_list', metadata=__metadata__, session=__session__)
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
	using_options(tablename='snps', metadata=__metadata__, session=__session__)
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
	using_options(tablename='snps_context', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('snps_id', 'gene_id'))

class CallMethod(Entity):
	"""
	2009-3-12
		add column filename
	2009-2-22
		add field unique_ecotype
	2009-2-6
		add eigen_value_ls
	2009-1-7
		add various columns to denote QC parameters
	"""
	short_name = Field(String(20))
	min_oligo_call_prob = Field(Float)
	max_array_mismatch_rate = Field(Float)
	max_array_NA_rate = Field(Float)
	max_snp_mismatch_rate = Field(Float)
	max_snp_NA_rate = Field(Float)
	npute_window_size = Field(Integer)
	imputed = Field(Boolean)
	unique_ecotype = Field(Boolean)
	avg_array_mismatch_rate = Field(Float)
	filename = Field(String(2024))
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	comment = Field(String(8000))
	call_info_ls = OneToMany("CallInfo")
	eigen_value_ls = OneToMany("CallMethodEigenValues", order_by='which_eigen')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class PhenotypeMethod(Entity):
	short_name = Field(String(20), unique=True)
	only_first_96 = Field(Boolean, default=0)
	biology_category = ManyToOne("BiologyCategory", colname='biology_category_id', ondelete='CASCADE', onupdate='CASCADE')
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	comment = Field(String(8000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	data_type = Field(String(200))
	transformation_description = Field(String(8000))
	using_options(tablename='phenotype_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class AnalysisMethod(Entity):
	"""
	2008-09-16
		add smaller_score_more_significant
	"""
	short_name = Field(String(120))
	method_description = Field(String(8000))
	smaller_score_more_significant = Field(Integer)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='analysis_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class ResultsMethodType(Entity):
	short_name = Field(String(30), unique=True)
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='results_method_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class ResultsMethod(Entity):
	short_name = Field(String(60), unique=True)
	filename = Field(String(1000), unique=True)
	original_filename = Field(Text)
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
	using_options(tablename='results_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class ResultsMethodJson(Entity):
	"""
	2009-5-4
		table to store json structures sent over to clients
	"""
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	no_of_top_snps = Field(Integer)
	min_MAF = Field(Float)
	json_data = Field(Binary(134217728), deferred=True)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='results_method_json', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_id', 'no_of_top_snps', 'min_MAF'))

class Results(Entity):
	"""
	2009-2-22
		change column name of "results_method" from results_method_id to results_id
	2009-2-6
		a table to store data derived from filenames stored in results_method
		add more fields
	"""
	snp = ManyToOne('Snps', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	score = Field(Float)
	rank = Field(Integer)
	beta = Field(Float)
	maf = Field(Float)
	mac = Field(Integer)
	genotype_var_perc = Field(Float)
	object = Field(Binary(134217728), deferred=True)	#a python dictionary to store other attributes
	using_options(tablename='results', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_id', 'snps_id'))

class CandidateGeneRankSumTestResult(Entity):
	"""
	2008-10-09
		add test_type = Field(Integer) to conform to CandidateGeneRankSumTestResultMethod
		rename results_by_gene to result
		rename results_by_gene_id to results_id
	2008-09-16
		table linked to results_by_gene, not results_method
	2008-07-17
	"""
	result = ManyToOne('ResultsByGene', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	statistic = Field(Float)
	pvalue = Field(Float)
	min_distance = Field(Integer)
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	max_pvalue_per_gene = Field(Integer)
	candidate_sample_size = Field(Integer)
	non_candidate_sample_size = Field(Integer)
	test_type = Field(Integer)
	comment = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_rank_sum_test_rbg', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	#using_table_options(UniqueConstraint('results_method_id', 'list_type_id'))

class CandidateGeneRankSumTestResultMethod(Entity):
	"""
	2008-10-09
		similar CandidateGeneRankSumTestResult linked to results_method

	"""
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	statistic = Field(Float)
	pvalue = Field(Float)
	min_distance = Field(Integer)
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	max_pvalue_per_gene = Field(Integer)
	candidate_sample_size = Field(Integer)
	non_candidate_sample_size = Field(Integer)
	test_type = Field(Integer)
	comment = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_rank_sum_test_rm', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_id', 'list_type_id', 'min_distance', 'get_closest', 'min_MAF', 'test_type'))


class ResultsByGene(Entity):
	"""
	2008-11-14
		more sensible unique constraints
	2008-09-15
		modify it to contain gene-based results(files) deriving from snp-based results in table ResultsMethod
	2008-08-27
		add readme_id
		modify unique constraint to include readme_id
	2008-07-19
	"""
	short_name = Field(String(60))
	filename = Field(String(767), unique=True)
	min_distance = Field(Integer)
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	comment = Field(Text)
	results_method = ManyToOne('ResultsMethod', colname='results_method_id', ondelete='CASCADE', onupdate='CASCADE')
	rank_test = OneToMany("CandidateGeneRankSumTestResult")
	top_snp_test = OneToMany("CandidateGeneTopSNPTest")
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='results_by_gene', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_method_id', 'min_distance', 'get_closest', 'min_MAF'))

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
	using_options(tablename='snps_qc', metadata=__metadata__, session=__session__)
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
	using_options(tablename='qc_method', metadata=__metadata__, session=__session__)
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
	using_options(tablename='contaminant_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class ArrayInfo(Entity):
	"""
	2009-3-25
		add median_intensity
	"""
	name = Field(String(40))
	filename = Field(String(1000))
	original_filename = Field(String(1000))
	description = Field(String(2000))
	ecotype_id = Field(Integer)
	maternal_ecotype_id = Field(Integer)
	paternal_ecotype_id = Field(Integer)
	contaminant_type = ManyToOne("ContaminantType", colname='contaminant_type_id', ondelete='CASCADE', onupdate='CASCADE')
	md5sum = Field(String(100))
	experimenter = Field(String(200))
	median_intensity = Field(Float)
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
	using_options(tablename='array_info', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class CallInfo(Entity):
	"""
	2009-2-2
		add pc_value_ls
	"""
	filename = Field(String(1000))
	description = Field(String(2000))
	array = ManyToOne("ArrayInfo", colname='array_id', ondelete='CASCADE', onupdate='CASCADE')
	NA_rate = Field(Float)
	readme = ManyToOne("README", colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	call_method = ManyToOne('CallMethod', colname='method_id', ondelete='CASCADE', onupdate='CASCADE')
	call_qc_ls = OneToMany("CallQC")
	pc_value_ls = OneToMany("CallInfoPCValues", order_by='which_pc')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_info', metadata=__metadata__, session=__session__)
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
	using_options(tablename='call_qc', metadata=__metadata__, session=__session__)
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
	using_options(tablename='probes', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class CandidateGeneTopSNPTest(Entity):
	"""
	2008-10-15
		add candidate_sample_size, non_candidate_sample_size, starting_rank, test_type to be compatible with CandidateGeneTopSNPTestRM
	2008-09-16
		table linked to results_by_gene, not results_method
	2008-08-20
	"""
	result = ManyToOne('ResultsByGene', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	pvalue = Field(Float)
	min_distance = Field(Integer)
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	candidate_sample_size = Field(Integer)
	non_candidate_sample_size = Field(Integer)
	no_of_top_candidate_genes = Field(Integer)
	no_of_top_genes = Field(Integer)
	no_of_top_snps = Field(Integer)
	starting_rank = Field(Integer)
	test_type = Field(Integer)
	comment = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_top_snp_test', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class CandidateGeneTopSNPTestRMType(Entity):
	"""
	2008-10-26
		move min_distance to CandidateGeneTopSNPTestRM
	2008-10-22
		hierarchical for CandidateGeneTopSNPTestRM
	"""
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	allow_two_sample_overlapping = Field(Integer)
	results_type = Field(Integer)
	test_type = ManyToOne('AnalysisMethod', colname='test_type_id', ondelete='CASCADE', onupdate='CASCADE')
	null_distribution_type = ManyToOne('NullDistributionType', colname='null_distribution_type_id', ondelete='CASCADE', onupdate='CASCADE')	#2008-10-30 useless
	comment = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_gene_top_snp_test_rm_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('get_closest', 'min_MAF', \
										'allow_two_sample_overlapping', 'results_type', 'test_type_id',\
										'null_distribution_type_id'))
	
class CandidateGeneTopSNPTestRM(Entity):
	"""
	2008-10-31
		add pvalue_gw_looping, pvalue_random_gene_list.
		null_distribution_type_id in CandidateGeneTopSNPTestRMType doesn't mean much.
		NULL data under different simulations stored in TopSNPTestRMNullData.
	2008-10-26
		add min_distance, no_of_tests_passed, no_of_tests
	2008-10-23
		restructure it for the new varieties of tests and link to CandidateGeneTopSNPTestRMType
	2008-10-15
		similar to CandidateGeneTopSNPTest but stores test results from ResultsMethod
	"""
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	pvalue = Field(Float)
	pvalue_gw_looping = Field(Float)	#2008-10-30	now CandidateGeneTopSNPTestRMNullData stores NULL data under different NULL distribution. all pvalues are merged into one entry.
	pvalue_random_gene_list = Field(Float)	#2008-10-30
	candidate_sample_size = Field(Integer)
	non_candidate_sample_size = Field(Integer)
	candidate_gw_size = Field(Integer)
	non_candidate_gw_size = Field(Integer)
	no_of_top_candidate_genes = Field(Integer)
	no_of_top_genes = Field(Integer)
	no_of_top_snps = Field(Integer)
	min_distance = Field(Integer)
	starting_rank = Field(Integer)
	no_of_tests_passed = Field(Integer)
	no_of_tests = Field(Integer)
	max_score = Field(Float)	#2008-10-22 keep record of max/min score in among these snps
	min_score = Field(Float)
	type = ManyToOne('CandidateGeneTopSNPTestRMType', colname='type_id', ondelete='CASCADE', onupdate='CASCADE')
	null_data_ls = OneToMany('TopSNPTestRMNullData')
	comment = Field(Text)
	using_options(tablename='candidate_gene_top_snp_test_rm', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_id', 'list_type_id', 'no_of_top_snps', \
										'min_distance', 'starting_rank', 'type_id'))

class CandidateGeneTopSNPTestRG(Entity):
	"""
	2008-11-11
		a table similar to CandidateGeneTopSNPTestRM. but candidate_sample_size, non_candidate_sample_size and etc refer to distinct genes instead of not SNPs.
		
	"""
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	pvalue = Field(Float)
	pvalue_gw_looping = Field(Float)
	pvalue_random_gene_list = Field(Float)
	candidate_sample_size = Field(Integer)
	non_candidate_sample_size = Field(Integer)
	candidate_gw_size = Field(Integer)
	non_candidate_gw_size = Field(Integer)
	no_of_top_candidate_genes = Field(Integer)
	no_of_top_genes = Field(Integer)
	no_of_top_snps = Field(Integer)
	min_distance = Field(Integer)
	starting_rank = Field(Integer)
	no_of_tests_passed = Field(Integer)
	no_of_tests = Field(Integer)
	max_score = Field(Float)
	min_score = Field(Float)
	type = ManyToOne('CandidateGeneTopSNPTestRMType', colname='type_id', ondelete='CASCADE', onupdate='CASCADE')
	comment = Field(Text)
	using_options(tablename='candidate_gene_top_snp_test_rg', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_id', 'list_type_id', 'no_of_top_snps', \
										'min_distance', 'starting_rank', 'type_id'))



class CandidateVsNonRatioPlot(Entity):
	"""
	2008-11-04
		table to store candidate-ratio/non-candidate-ratio curve plots in database
	"""
	type = ManyToOne('CandidateGeneTopSNPTestRMType', colname='type_id', ondelete='CASCADE', onupdate='CASCADE')
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	png_thumbnail = Field(Binary(64217728), deferred=True)
	png_data = Field(Binary(134217728), deferred=True)
	svg_data = Field(Binary(length=134217728), deferred=True)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='candidate_vs_non_ratio_plot', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('type_id', 'results_id', 'list_type_id'))

class TopSNPTestRMNullData(Entity):
	"""
	2008-10-31
		table to store NULL data of CandidateGeneTopSNPTestRM from different simulations under different NULL distributions
	"""
	observed = ManyToOne('CandidateGeneTopSNPTestRM', colname='observed_id', ondelete='CASCADE', onupdate='CASCADE')
	candidate_sample_size = Field(Integer)
	candidate_gw_size = Field(Integer)
	run_no = Field(Integer)
	null_distribution_type = ManyToOne('NullDistributionType', colname='null_distribution_type_id', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='top_snp_test_rm_null_data', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('observed_id', 'run_no', 'null_distribution_type_id'))
	

class SNPRegionPlotType(Entity):
	"""
	2008-10-06
	"""
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	snp_region_plot_ls = OneToMany('SNPRegionPlot')
	#ManyToMany('SNPRegionPlot', tablename='snp_region_plot2type')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snp_region_plot_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class SNPRegionPlot(Entity):
	"""
	2009-5-1
		add column call_method_id
	2008-10-21
		rename img_data to png_data and add svg_data
	2008-10-16
		fix an error in the definitions of img_data and original_filename. they were swapped
	2008-10-06
		table to store binary SNP region plots
	"""
	chromosome = Field(Integer)
	start = Field(Integer)
	stop = Field(Integer)
	png_data = Field(Binary(134217728), deferred=True)
	svg_data = Field(Binary(length=134217728), deferred=True)
	center_snp_position = Field(Integer)
	original_filename = Field(Text)	#no bigger than 128M
	call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	phenotype_method = ManyToOne('PhenotypeMethod', colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	plot_type = ManyToOne('SNPRegionPlotType', colname='plot_type_id', ondelete='CASCADE', onupdate='CASCADE')
	"""
	plot_types = ManyToMany('SNPRegionPlotType', colname='plot_type_id', ondelete='CASCADE', onupdate='CASCADE', \
						tablename='snp_region_plot2type')	#local_side, remote_side
	"""
	plot2gene_ls = OneToMany("SNPRegionPlotToGene")
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snp_region_plot', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('chromosome', 'start', 'stop', 'phenotype_method_id', 'plot_type_id'))

class SNPRegionPlotToGene(Entity):
	"""
	2008-10-06
	"""
	snp_region_plot = ManyToOne('SNPRegionPlot', colname='plot_id', ondelete='CASCADE', onupdate='CASCADE')
	gene_id = Field(Integer)
	using_options(tablename='snp_region_plot2gene', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('plot_id', 'gene_id'))

class LD(Entity):
	"""
	2008-10-15
		table to store Linkage Disequilibrium data (output of MpiLD.py)
	"""
	chr1 = Field(Integer)
	pos1 = Field(Integer)
	chr2 = Field(Integer)
	pos2 = Field(Integer)
	d = Field(Float)
	d_prime = Field(Float)
	r2 = Field(Float)
	snp1_maf = Field(Float)
	snp2_maf = Field(Float)
	no_of_pairs = Field(Integer)
	call_method_id = Field(Integer)
	#call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='ld', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class NullDistributionType(Entity):
	"""
	2008-10-21
		a table to store definition for different types null distributions
	"""
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='null_distribution_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class ScoreRankHistogramType(Entity):
	"""
	2008-10-30
		connect it to "ResultsGene" via ManyToMany()
	2008-10-21
		add null_distribution_type_id
	2008-10-15
		type of score/rank histograms in ScoreRankHistogram
	"""
	call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	min_distance = Field(Integer)
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	allow_two_sample_overlapping = Field(Integer)
	results_type = Field(Integer)
	null_distribution_type = ManyToOne('NullDistributionType', colname='null_distribution_type_id', ondelete='CASCADE', onupdate='CASCADE')
	score_rank_hist_ls = OneToMany('ScoreRankHistogram')
	results_gene_ls = ManyToMany("ResultsGene", tablename='results_gene2type', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='score_rank_histogram_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('call_method_id', 'min_distance', 'get_closest', 'min_MAF', \
										'allow_two_sample_overlapping', 'results_type', 'null_distribution_type_id'))

class ScoreRankHistogram(Entity):
	"""
	2008-10-15
		table to store score/rank histogram divided by candidate gene list. output of CheckCandidateGeneRank.py
	"""
	phenotype_method = ManyToOne('PhenotypeMethod', colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	hist_type = ManyToOne('ScoreRankHistogramType', colname='hist_type_id', ondelete='CASCADE', onupdate='CASCADE')
	score_hist = Field(Binary(length=134217728), deferred=True)
	score_hist_svg = Field(Binary(length=134217728), deferred=True)
	rank_hist = Field(Binary(length=134217728), deferred=True)
	rank_hist_svg = Field(Binary(length=134217728), deferred=True)
	original_filename = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='score_rank_histogram', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('phenotype_method_id', 'list_type_id', 'hist_type_id'))

class MAFVsScorePlot(Entity):
	"""
	2008-10-17
		table to store pvalue vs MAF plots
	"""
	results_method = ManyToOne('ResultsMethod', colname='results_method_id', ondelete='CASCADE', onupdate='CASCADE')
	png_data = Field(Binary(134217728), deferred=True)
	svg_data = Field(Binary(length=134217728), deferred=True)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='maf_vs_score_plot', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_method_id'))

class ResultsGene(Entity):
	"""
	2009-4-7
		add rank into the unique combo. rank could be different for the same SNP under different MAF scenarios.
	2008-10-30
		remove order_by='rank' in types = ManyToMany('ScoreRankHistogramType', ...)
	2008-10-28
		store genes associated with SNPs in results_method under certain association rule
	"""
	snp = ManyToOne('Snps', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	disp_pos = Field(Integer)
	gene_id = Field(Integer)
	score = Field(Float)
	rank = Field(Float)
	result = ManyToOne('ResultsMethod', colname='results_id', ondelete='CASCADE', onupdate='CASCADE')
	types = ManyToMany('ScoreRankHistogramType', \
						tablename='results_gene2type')	#local_side, remote_side
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='results_gene', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('snps_id', 'gene_id', 'results_id', 'rank'))


class FTPathway(Entity):
	"""
	2008-11-14
		table storing flowering time pathways
	"""
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	gene_ls = OneToMany('FTGene')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='ft_pathway', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')


class FTPathwayRelationship(Entity):
	"""
	2008-11-14
		table storing flowering time pathway relationships
	"""
	pathway1 = ManyToOne('FTPathway', colname='pathway1_id', ondelete='CASCADE', onupdate='CASCADE')
	pathway2 = ManyToOne('FTPathway', colname='pathway2_id', ondelete='CASCADE', onupdate='CASCADE')
	relationship_type_id = Field(Integer)
	arrow_start_point_loc = Field(String(512))
	arrow_end_point_loc = Field(String(512))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='ft_pathway_relationship', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('pathway1_id', 'pathway2_id', 'relationship_type_id'))


class FTGene(Entity):
	"""
	2008-11-14
		table storing flowering time pathway genes
	"""
	gene_id = Field(Integer)
	pathway = ManyToOne('FTPathway', colname='pathway_id', ondelete='CASCADE', onupdate='CASCADE')
	paper = Field(String(512))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='ft_gene', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('gene_id', 'pathway_id'))


class SNPAnnotation(Entity):
	"""
	2009-01-05
		table to store annotation of SNPs, like synonymous, non-synonymous, ...
		information finer than SnpsContext
	"""
	snp = ManyToOne('Snps', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	gene_id = Field(Integer)
	gene_commentary_id = Field(Integer)
	snp_annotation_type = ManyToOne('SNPAnnotationType', colname='snp_annotation_type_id', ondelete='CASCADE', onupdate='CASCADE')
	which_exon_or_intron = Field(Integer)
	pos_within_codon = Field(Integer)
	comment = Field(String(512))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snp_annotation', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('snps_id', 'gene_commentary_id', 'snp_annotation_type_id'))

class SNPAnnotationType(Entity):
	"""
	2009-01-05
		table to store types of SNP annotation, like synonymous, non-synonymous, ...
	"""
	short_name = Field(String(256), unique=True)
	description = Field(String(8192))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snp_annotation_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')


class CallInfoPCValues(Entity):
	"""
	2009-2-2
		table to store principle component values after PCA is applied on a certain call method dataset
	"""
	call_info = ManyToOne('CallInfo', colname='call_info_id', ondelete='CASCADE', onupdate='CASCADE')
	which_pc = Field(Integer)
	pc_value = Field(Float)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_info_pc_values', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('call_info_id', 'which_pc'))


class CallMethodEigenValues(Entity):
	"""
	2009-2-2
		table to store eigen values after PCA is applied on a certain call method dataset
	"""
	call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	which_eigen = Field(Integer)
	eigen_value = Field(Float)
	variance_perc = Field(Float)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_method_eigen_values', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('call_method_id', 'which_eigen'))

class PhenotypeHistPlots(Entity):
	"""
	2009-2-6
		phenotype histogram
	"""
	phenotype_method = ManyToOne('PhenotypeMethod', colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	hist_thumb = Field(Binary(64217728), deferred=True)
	hist_plot = Field(Binary(134217728), deferred=True)
	hist_plot_fname = Field(String(1000))
	hist_log_thumb = Field(Binary(64217728), deferred=True)
	hist_log_plot = Field(Binary(134217728), deferred=True)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='phenotype_hist_plots', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('phenotype_method_id'))

class CallPhenotypeQQPlots(Entity):
	"""
	2009-2-6
		qq plot
	"""
	call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	phenotype_method = ManyToOne('PhenotypeMethod', colname='phenotype_method_id', ondelete='CASCADE', onupdate='CASCADE')
	qq_thumb = Field(Binary(64217728), deferred=True)
	qq_plot = Field(Binary(134217728), deferred=True)
	qq_plot_fname = Field(String(1000))
	qq_log_thumb = Field(Binary(64217728), deferred=True)
	qq_log_plot = Field(Binary(134217728), deferred=True)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_phenotype_qq_plots', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('call_method_id', 'phenotype_method_id'))

class CmpEnrichmentOfTwoAnalysisMethods(Entity):
	"""
	2009-4-16 to be filled up by PlotCmpTwoAnalysisMethods.py
	"""
	result1 = ManyToOne('ResultsMethod', colname='results_id1', ondelete='CASCADE', onupdate='CASCADE')
	result2 = ManyToOne('ResultsMethod', colname='results_id2', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	type = ManyToOne('ScoreRankHistogramType', colname='type_id', ondelete='CASCADE', onupdate='CASCADE')
	r1_min_score = Field(Float)
	r1_max_score = Field(Float)
	r2_min_score = Field(Float)
	r2_max_score = Field(Float)
	threshold_type = ManyToOne('CmpEnrichmentOfTwoAnalysisMethodsType', colname='threshold_type_id', ondelete='CASCADE', onupdate='CASCADE')
	candidate_sample_size = Field(Integer)
	non_candidate_sample_size = Field(Integer)
	candidate_gw_size = Field(Integer)
	non_candidate_gw_size = Field(Integer)
	enrichment_ratio = Field(Float)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='cmp_enrichment_of_two_analysis_methods', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_id1', 'results_id2', 'list_type_id', 'type_id', 'r1_min_score', \
										'r1_max_score', 'r2_min_score', 'r2_max_score', 'threshold_type_id'))

class CmpEnrichmentOfTwoAnalysisMethodsType(Entity):
	"""
	2009-5-4 table to store threshold types of CmpEnrichmentOfTwoAnalysisMethods
	"""
	r1_threshold = Field(Float)
	r2_threshold = Field(Float)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='cmp_enrichment_of_two_analysis_methods_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('r1_threshold', 'r2_threshold'))


"""
#2008-10-29 'results_gene2type' automatically generated by ManyToMany is sub-optimal because it uses MyIAM engine and doesn't allow foreign key.
#two ways around this
#1. set "default-storage-engine=InnoDB" in mysqld section in my.cnf
#OR
#2. create the table manually with same columns + innodb engine + foreign keys 
class ResultsGene2Type(Entity):
	results_gene = ManyToOne('ResultsGene', colname='results_gene_id', ondelete='CASCADE', onupdate='CASCADE')
	score_rank_histogram_type = ManyToOne('ScoreRankHistogramType', colname='score_rank_histogram_type_id', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='results_gene2type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
"""

"""
class TopSNPTestRMPlotType(Entity):
	#2008-10-20
	#	type for TopSNPTestRMPlot
	no_of_snps = Field(Integer)
	test_type = Field(Integer)
	results_type = Field(Integer)
	min_distance = Field(Integer)
	get_closest = Field(Integer)
	min_MAF = Field(Float)
	
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='top_snp_test_rm_plot_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('no_of_snps', 'test_type', 'results_type', \
										'min_distance', 'get_closest', 'min_MAF'))

class TopSNPTestRMPlot(Entity):
	#2008-10-20
	#	store plots based on data from CandidateGeneTopSNPTestRM
	
	png_data = Field(Binary(134217728), deferred=True)
	svg_data = Field(Binary(length=134217728), deferred=True)
	results_method = ManyToOne('ResultsMethod', colname='results_method_id', ondelete='CASCADE', onupdate='CASCADE')
	list_type = ManyToOne('GeneListType', colname='list_type_id', ondelete='CASCADE', onupdate='CASCADE')
	plot_type = ManyToOne('TopSNPTestRMPlotType', colname='plot_type_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='top_snp_test_rm_plot', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('results_method_id', 'list_type_id', 'plot_type_id'))
"""

class WebHelp(Entity):
	"""
	2009-4-26 table to store help document of the web interface (accessible on the web thru a controller)
	"""
	short_name = Field(String(256), unique=True)
	title = Field(String(1024))
	content = Field(Text(8192*20))
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='web_help', metadata=__metadata__, session=__session__)
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
							('pool_recycle', 0, int):[3600, '', 1, 'the length of time to keep connections open before recycling them.'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-10-06
			add option 'pool_recycle' to recycle connection. MySQL typically close connections after 8 hours.
			__metadata__.bind = create_engine(self._url, pool_recycle=self.pool_recycle)
		2008-07-09
		"""
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		if getattr(self, 'schema', None):	#for postgres
			for entity in entities:
				if entity.__module__==self.__module__:	#entity in the same module
					using_table_options_handler(entity, schema=self.schema)
		
		#2008-10-05 MySQL typically close connections after 8 hours resulting in a "MySQL server has gone away" error.
		__metadata__.bind = create_engine(self._url, pool_recycle=self.pool_recycle)
		
		self.metadata = __metadata__
		self.session = __session__
	
	def setup(self, create_tables=True):
		"""
		2008-09-07
			expose option create_tables, default=True. assign it to False if no new table is to be created.
		"""
		setup_all(create_tables=create_tables)	#create_tables=True causes setup_all to call elixir.create_all(), which in turn calls metadata.create_all()
		#2008-08-26 setup_all() would setup other databases as well if they also appear in the program. Seperate this to be envoked after initialization
		# to ensure the metadata of other databases is setup properly.
	
	def getGWA(self, call_method_id, phenotype_method_id, analysis_method_id, results_directory=None, min_MAF=0.1, construct_chr_pos2index=False, pdata=None):
		"""
		2009-3-12
			convenient function to get genome-wide-association result from database
		"""
		rm = ResultsMethod.query.filter_by(call_method_id=call_method_id).filter_by(phenotype_method_id=phenotype_method_id).\
					filter_by(analysis_method_id=analysis_method_id).first()
		from GeneListRankTest import GeneListRankTest
		return GeneListRankTest.getResultMethodContent(rm, results_directory=results_directory, min_MAF=min_MAF, 
																construct_chr_pos2index=construct_chr_pos2index,\
																pdata=pdata)
	
	def getSNPMatrix(self, call_method_id, ignore_2nd_column=True):
		"""
		2009-3-12
			given a call_method_id, fetch a whole SNP dataset
		"""
		cm = CallMethod.get(call_method_id)
		snpData = SNPData(input_fname=cm.filename, turn_into_array=1, ignore_2nd_column=ignore_2nd_column)	#use 1st column (ecotype id) as main ID
		return snpData
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Stock_250kDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.setup()
	import pdb
	pdb.set_trace()
	rows = GeneListType.query.all()
	for row in rows[:10]:
		print row.gene_list[0].list_type_id
