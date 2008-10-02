#!/usr/bin/env python
"""
2008-10-01
        example to call DrawSNPRegion as an API
"""
import sys, os, math
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from variation.src.DrawSNPRegion import DrawSNPRegion

instance = DrawSNPRegion(db_user='yh', db_passwd='', hostname='papaya.usc.edu', database='stock_250k',debug=2, input_fname='/tmp/dumb', output_dir='/tmp')
gene_annotation_picklef = '/Network/Data/250k/tmp-yh/at_gene_model_pickelf'
#LD_info_picklef = '/Network/Data/250k/tmp-yh/call_method_17_LD_m0.2_n0.1_m40000'
LD_info_picklef = None
LD_fname = '/Network/Data/250k/tmp-yh/call_method_17_LD_m0.2.tsv'
phenotype_method_id = 1
chromosome, start, stop  = 2, 9588685-40000, 9588685+40000
output_fname = '/tmp/SNP_2_9588685'
center_snp_position = 9588685
instance.drawSNPRegion(gene_annotation_picklef, LD_info_picklef, phenotype_method_id, chromosome, start, stop, \
		output_fname, center_snp_position=None, LD_fname=LD_fname,\
                which_LD_statistic=2, call_method_id=17, min_MAF=0.1)
