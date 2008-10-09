"""Helper functions

Consists of functions to typically be used within templates, but also
available to Controllers. This module is available to both as 'h'.
"""
from webhelpers import *
import os,sys
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import helloworld.model as model
from variation.src.DrawSNPRegion import DrawSNPRegion


def returnGeneDescLs(gene_annotation, plot2gene_ls=[]):
	DrawSNPRegion_ins = DrawSNPRegion(db_user=model.db_user, db_passwd=model.db_passwd, hostname=model.hostname, database=model.dbname,\
							input_fname='/tmp/dumb', output_dir='/tmp', debug=0)
	
	matrix_of_gene_descriptions = []
	for plot2gene in plot2gene_ls:
		gene_model = gene_annotation.gene_id2model.get(plot2gene.gene_id)
		if gene_model:
			if len(gene_model.gene_commentaries)==0:
				gene_commentaries = [gene_model]	#fake one here
			else:
				gene_commentaries = gene_model.gene_commentaries
			for gene_commentary in gene_commentaries:
				gene_desc_ls = DrawSNPRegion_ins.returnGeneDescLs(DrawSNPRegion_ins.gene_desc_names, gene_model, gene_commentary)
				matrix_of_gene_descriptions.append(gene_desc_ls)
	return matrix_of_gene_descriptions