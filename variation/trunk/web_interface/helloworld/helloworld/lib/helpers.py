"""Helper functions

Consists of functions to typically be used within templates, but also
available to Controllers. This module is available to both as 'h'.
"""
from webhelpers import *
import os,sys
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import helloworld.model as model
from variation.src.DrawSNPRegion import DrawSNPRegion
from pymodule import PassingData

def returnGeneDescLs(gene_annotation, gene_id_ls=[]):
	DrawSNPRegion_ins = DrawSNPRegion(db_user=model.db_user, db_passwd=model.db_passwd, hostname=model.hostname, database=model.dbname,\
							input_fname='/tmp/dumb', output_dir='/tmp', debug=0)
	
	matrix_of_gene_descriptions = []
	for gene_id in gene_id_ls:
		gene_model = gene_annotation.gene_id2model.get(gene_id)
		if gene_model:
			if len(gene_model.gene_commentaries)==0:
				gene_commentaries = [gene_model]	#fake one here
			else:
				gene_commentaries = gene_model.gene_commentaries
			for gene_commentary in gene_commentaries:
				gene_desc_ls = DrawSNPRegion_ins.returnGeneDescLs(DrawSNPRegion_ins.gene_desc_names, gene_model, gene_commentary)
				matrix_of_gene_descriptions.append(gene_desc_ls)
	return matrix_of_gene_descriptions

def getPhenotypeInfo(affiliated_table_name, extra_condition=None, extra_tables=None):
	"""
	2008-10-19
		add option extra_tables
	2008-10-16
		sort phenotype by biology_category_id and return other info as well
	"""
	table_str = '%s s, %s p'%(affiliated_table_name, model.Stock_250kDB.PhenotypeMethod.table.name)
	if extra_tables:
		table_str += ', %s'%extra_tables
	where_condition = 'p.id=s.phenotype_method_id'
	if extra_condition:
		where_condition += ' and %s'%extra_condition
	
	rows = model.db.metadata.bind.execute("select distinct p.id, p.biology_category_id, p.short_name from %s\
		where %s order by p.biology_category_id, p.id"\
		%(table_str, where_condition))
	phenotype_method_id_ls = []
	phenotype_method_id2index = {}
	phenotype_method_label_ls = []
	prev_biology_category_id = -1
	no_of_separators = 0
	for row in rows:
		if prev_biology_category_id == -1:
			prev_biology_category_id = row.biology_category_id
		elif row.biology_category_id!=prev_biology_category_id:
			prev_biology_category_id = row.biology_category_id
			#add a blank phenotype id as separator
			no_of_separators += 1
			phenotype_method_id2index[-no_of_separators] = len(phenotype_method_id_ls)
			phenotype_method_id_ls.append(-no_of_separators)
			phenotype_method_label_ls.append('=====')
		phenotype_method_id2index[row.id] = len(phenotype_method_id_ls)
		phenotype_method_id_ls.append(row.id)
		phenotype_method_label_ls.append('%s %s'%(row.id, row.short_name))
	phenotype_info = PassingData()
	phenotype_info.phenotype_method_id2index = phenotype_method_id2index
	phenotype_info.phenotype_method_id_ls = phenotype_method_id_ls
	phenotype_info.phenotype_method_label_ls = phenotype_method_label_ls
	return phenotype_info

def getListTypeInfo(affiliated_table_name, extra_condition=None, extra_tables=None):
	"""
	2008-10-19
		add option extra_tables
	2008-10-16
		sort gene list type by biology_category_id and return other info as well
		add -1 as a separator into list_type_id_ls
	"""
	table_str = '%s s, %s p'%(affiliated_table_name, model.Stock_250kDB.GeneListType.table.name)
	if extra_tables:
		table_str += ', %s'%extra_tables
	where_condition = 'p.id=s.list_type_id'
	if extra_condition:
		where_condition += ' and %s'%extra_condition
	rows = model.db.metadata.bind.execute("select distinct p.id, p.biology_category_id, p.short_name from %s \
		where %s order by p.biology_category_id, p.id"\
		%(table_str, where_condition))
	list_type_id_ls = []
	list_type_id2index = {}
	list_type_label_ls = []
	prev_biology_category_id = -1
	no_of_separators = 0
	for row in rows:
		if prev_biology_category_id == -1:
			prev_biology_category_id = row.biology_category_id
		elif row.biology_category_id!=prev_biology_category_id:
			prev_biology_category_id = row.biology_category_id
			no_of_separators += 1
			list_type_id2index[-no_of_separators] = len(list_type_id_ls)
			list_type_id_ls.append(-no_of_separators)
			list_type_label_ls.append('====\n====')
		list_type_id2index[row.id] = len(list_type_id_ls)
		list_type_id_ls.append(row.id)
		list_type_label_ls.append('%s %s'%(row.id, row.short_name))
	list_info = PassingData()
	list_info.list_type_id2index = list_type_id2index
	list_info.list_type_id_ls = list_type_id_ls
	list_info.list_type_label_ls = list_type_label_ls
	return list_info

def getAnalysisMethodInfo(affiliated_table_name, extra_condition=None, extra_tables=None):
	"""
	2008-10-19
	"""
	table_str = '%s s, %s a'%(affiliated_table_name, model.Stock_250kDB.AnalysisMethod.table.name)
	if extra_tables:
		table_str += ', %s'%extra_tables
	where_condition = 'a.id=s.analysis_method_id'
	if extra_condition:
		where_condition += ' and %s'%extra_condition
	rows = model.db.metadata.bind.execute("select distinct a.id, a.short_name from %s \
		where %s order by a.id"\
		%(table_str, where_condition))
	id_ls = []
	id2index = {}
	label_ls = []
	prev_biology_category_id = -1
	no_of_separators = 0
	for row in rows:
		id2index[row.id] = len(id_ls)
		id_ls.append(row.id)
		label_ls.append('%s %s'%(row.id, row.short_name))
	list_info = PassingData()
	list_info.id2index = id2index
	list_info.id_ls = id_ls
	list_info.label_ls = label_ls
	return list_info