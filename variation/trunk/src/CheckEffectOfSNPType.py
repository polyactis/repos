#!/usr/bin/env python
"""

Examples:
	#check KW on LD
	CheckEffectOfSNPType.py -u yh -m 4 -o /Network/Data/250k/tmp-yh/EffectOfSNPType/kw_LD_m4_snp_type -a1 -q1 -j /Network/Data/250k/tmp-yh/snp_annotation_dstruc.pickle

Description:
	2008-1-8 program to check the effect of different types (intron, synonymous, non-syn ...) of SNPs.
		1. draws histograms of SNP pvalues from different types
		2. draws piecharts showing composition of types for SNPs above a value threshold

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback
from pymodule import ProcessOptions, PassingData, getGenomeWideResultFromFile, getListOutOfStr
import Stock_250kDB

import pylab, cPickle

"""
2008-01-06
"""
class SNPAnnotationDataStruc(object):
	def __init__(self):
		self.data_ls = []
		self.chr_pos2index_ls = {}
		self.table_field_ls = ['id', 'snps_id', 'gene_id', 'gene_commentary_id', 'snp_annotation_type_id', \
						'which_exon_or_intron', 'pos_within_codon', 'comment']
	def addOneEntry(self, row):
		from pymodule import PassingData
		data_obj = PassingData()
		
		for table_field in self.table_field_ls:
			setattr(data_obj, table_field, getattr(row, table_field, None))
		data_obj.chr = row.snp.chromosome
		data_obj.pos = row.snp.position
		chr_pos_key = (data_obj.chr, data_obj.pos)
		if chr_pos_key not in self.chr_pos2index_ls:
			self.chr_pos2index_ls[chr_pos_key] = []
		self.chr_pos2index_ls[chr_pos_key].append(len(self.data_ls))
		self.data_ls.append(data_obj)
	
	def returnAnnotationGivenChrPos(self, chr, pos):
		chr_pos_key = (chr, pos)
		annotation_ls = []
		obj_index_ls = self.chr_pos2index_ls.get(chr_pos_key)
		if obj_index_ls:
			for obj_index in obj_index_ls:
				annotation_ls.append(self.data_ls[obj_index])
		return annotation_ls

class CheckEffectOfSNPType(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('call_method_id', 0, int):[17, '', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
							('analysis_method_id_ls', 0, ):['7', 'a', 1, 'Restrict results based on this list of analysis_method. Default is no such restriction.'],\
							('phenotype_method_id_ls', 0, ):[None, 'q', 1, 'Restrict results based on this list of phenotype_method. Default is no such restriction.'],\
							('output_fname_prefix', 0, ): [None, 'o', 1, ''],\
							("snp_annotation_dstruc_picklef", 0, ): [None, 'j', 1, 'given the option, If the file does not exist yet, store a pickled data struc into it. If the file exists, load data out of it.'],\
							('min_value', 1, float): [4, 'm', 1, 'minimum value from genome_wide_result for a SNP to be counted in'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2009-1-5
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		self.analysis_method_id_ls = getListOutOfStr(self.analysis_method_id_ls, data_type=int)
		self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
		
	def getGenomeWideResult(self, call_method_id, phenotype_method_id, analysis_method_id):
		rows = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id).filter_by(analysis_method_id=analysis_method_id).\
					filter_by(phenotype_method_id=phenotype_method_id).filter_by(results_method_type_id=1)
		
		pdata = PassingData()
		if rows.count()==1:
			rm = rows.first()
		elif rows.count()==0:
			sys.stderr.write("No result fetched from db based on call_method_id=%s, analysis_method_id=%s, phenotype_method_id=%s.\n"%\
							(call_method_id, analysis_method_id, phenotype_method_id))
			rm = None
		else:
			sys.stderr.write("First result out of %s results fetched from db based on call_method_id=%s, analysis_method_id=%s, phenotype_method_id=%s.\n"%\
							(rows.count(), call_method_id, analysis_method_id, phenotype_method_id))
			rm = rows.first()
		if rm:
			input_fname = rm.filename
			pdata.gwr_name = '%s_%s_%s'%(rm.analysis_method.short_name, rm.phenotype_method_id, rm.phenotype_method.short_name)
		else:
			return
		
		genome_wide_result = getGenomeWideResultFromFile(input_fname, min_value_cutoff=None, do_log10_transformation=True, pdata=pdata)
		return genome_wide_result
	
	def getSNPAnnotationType2value_ls(self, genome_wide_result, snp_annotation_dstruc):
		
		snp_annotation_type_id2value_ls = {}
		no_of_snps_with_no_annotation = 0
		for data_obj in genome_wide_result.data_obj_ls:
			snp_annotation_ls = snp_annotation_dstruc.returnAnnotationGivenChrPos(data_obj.chromosome, data_obj.position)
			for snp_annotation in snp_annotation_ls:
				ty_id = snp_annotation.snp_annotation_type_id
				if ty_id not in snp_annotation_type_id2value_ls:
					snp_annotation_type_id2value_ls[ty_id] = []
				snp_annotation_type_id2value_ls[ty_id].append(data_obj.value)
			if len(snp_annotation_ls)==0:
				no_of_snps_with_no_annotation += 1
				print data_obj.chromosome, data_obj.position
		print "%s SNPs have no annotation."%no_of_snps_with_no_annotation
		return snp_annotation_type_id2value_ls
	
	def getSNPAnnotationDataStruc(self):
		sys.stderr.write("Getting snp_annotation_dstruc ...\n")
		snp_annotation_dstruc = SNPAnnotationDataStruc()
		i = 0
		block_size = 10000
		rows = Stock_250kDB.SNPAnnotation.query.offset(i).limit(block_size)
		while rows.count()!=0:
			for row in rows:
				i += 1
				snp_annotation_dstruc.addOneEntry(row)
			sys.stderr.write("%s%s"%('\x08'*40, i))
			rows = Stock_250kDB.SNPAnnotation.query.offset(i).limit(block_size)
		sys.stderr.write("Done.\n")
		return snp_annotation_dstruc	
	
	def get_snp_annotation_type_id2name(self):
		sys.stderr.write("Getting snp_annotation_type_id2name ...")
		snp_annotation_type_id2name = {}
		rows = Stock_250kDB.SNPAnnotationType.query.all()
		for row in rows:
			snp_annotation_type_id2name[row.id] = row.short_name
		sys.stderr.write("Done.\n")
		return snp_annotation_type_id2name
	
	def drawHistFromSNPAnnotationType2value_ls(self, snp_annotation_type_id2value_ls, snp_annotation_type_id2name, output_fname_prefix, \
											snp_annotation_type_id_ls=None, normed=False):
		sys.stderr.write("Drawing Histogram ...")
		pylab.clf()
		pylab.grid(True, alpha=0.3)
		if snp_annotation_type_id_ls is None:
			snp_annotation_type_id_ls = snp_annotation_type_id2value_ls.keys()
			snp_annotation_type_id_ls.sort()
		hist_patch_ls = []
		legend_ls = []
		facecolors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
		no_of_facecolors = len(facecolors)
		for i in range(len(snp_annotation_type_id_ls)):
			snp_annotation_type_id = snp_annotation_type_id_ls[i]
			snp_annotation_type_name = snp_annotation_type_id2name.get(snp_annotation_type_id)
			value_ls = snp_annotation_type_id2value_ls.get(snp_annotation_type_id)
			if snp_annotation_type_name is None or value_ls is None:
				continue
			if i>=no_of_facecolors:
				linewidth=0.3
			else:
				linewidth=0
			if len(value_ls)>30:
				no_of_bins = min(len(value_ls)/10, 40)
				h1 = pylab.hist(value_ls, no_of_bins, alpha=0.3, facecolor=facecolors[i%no_of_facecolors], linewidth=linewidth, normed=normed)
				hist_patch_ls.append(h1[2][0])
				legend_ls.append(snp_annotation_type_name)
		if hist_patch_ls and legend_ls:
			pylab.legend(hist_patch_ls, legend_ls, handlelen=0.02)
		if len(hist_patch_ls)>0:
			pylab.savefig('%s_hist.png'%output_fname_prefix, dpi=300)
			pylab.savefig('%s_hist.svg'%output_fname_prefix)
		sys.stderr.write("Done.\n")
	
	def drawPieChart(self, genome_wide_result, snp_annotation_dstruc, snp_annotation_type_id2name, output_fname_prefix, min_value=4):
		"""
		"""
		sys.stderr.write("Drawing Pie chart ...")
		snp_annotation_type_id2count = {}
		no_of_snps_with_no_annotation = 0
		for data_obj in genome_wide_result.data_obj_ls:
			snp_annotation_ls = snp_annotation_dstruc.returnAnnotationGivenChrPos(data_obj.chromosome, data_obj.position)
			for snp_annotation in snp_annotation_ls:
				if data_obj.value>=min_value:
					ty_id = snp_annotation.snp_annotation_type_id
					if ty_id not in snp_annotation_type_id2count:
						snp_annotation_type_id2count[ty_id] = 0
					snp_annotation_type_id2count[ty_id] += 1
			if len(snp_annotation_ls)==0:
				no_of_snps_with_no_annotation += 1
				print data_obj.chromosome, data_obj.position
		print "%s SNPs have no annotation."%no_of_snps_with_no_annotation
		
		#need to subtract redundant counting
		#intron - splice donor/acceptor
		#non-synonymous - premature-stop-codon/init-Met
		snp_annotation_type_id_pairs_to_subtract = [(4,9), (4,10), (6,7), (6,8)]
		for snp_annotation_type_id_pair in snp_annotation_type_id_pairs_to_subtract:
			ty_id1, ty_id2 = snp_annotation_type_id_pair
			if ty_id1 in snp_annotation_type_id2count and ty_id2 in snp_annotation_type_id2count:
				snp_annotation_type_id2count[ty_id1] = snp_annotation_type_id2count[ty_id1] - snp_annotation_type_id2count[ty_id2]
		#merge splice donor/acceptor
		if 10  in snp_annotation_type_id2count:
			if 9 not in snp_annotation_type_id2count:
				snp_annotation_type_id2count[9] = 0
			snp_annotation_type_id2count[9] += snp_annotation_type_id2count[10]
			del snp_annotation_type_id2count[10]
		
		#ty_id_value_ls = snp_annotation_type_id2count.items()
		#ty_id_value_ls.sort()
		count_ls = []
		label_ls = []
		snp_annotation_type_id_ls = [1,7,2,8,3,9,4,14,5,15,6,12,13,16]	#order of types in the pie char to minimize small pies next to each other
		for ty_id in snp_annotation_type_id_ls:
			if ty_id in snp_annotation_type_id2count:
				count = snp_annotation_type_id2count[ty_id]
				count_ls.append(count)
				label_ls.append(snp_annotation_type_id2name[ty_id])
		
		pylab.clf()
		if len(count_ls)>0:
			pylab.pie(count_ls, labels=label_ls, autopct='%.2f', pctdistance=0.75,  colors=('b', 'g', 'r', 'c', 'm', 'y', 'w'))
			pylab.savefig('%s_pie.png'%output_fname_prefix, dpi=300)
			pylab.savefig('%s_pie.svg'%output_fname_prefix)
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		"""		
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		
		snp_annotation_type_id_ls=[]
		if os.path.isfile(self.snp_annotation_dstruc_picklef):
			picklef = open(self.snp_annotation_dstruc_picklef)
			snp_annotation_dstruc = cPickle.load(picklef)
			del picklef
		else:
			snp_annotation_dstruc = self.getSNPAnnotationDataStruc()
			if self.snp_annotation_dstruc_picklef:
				picklef = open(self.snp_annotation_dstruc_picklef, 'w')
				cPickle.dump(snp_annotation_dstruc, picklef, -1)
				picklef.close()
		
		for analysis_method_id in self.analysis_method_id_ls:
			for phenotype_method_id in self.phenotype_method_id_ls:
				genome_wide_result = self.getGenomeWideResult(self.call_method_id, phenotype_method_id, analysis_method_id)
				if genome_wide_result:
					phenotype_method = Stock_250kDB.PhenotypeMethod.get(phenotype_method_id)
					analysis_method = Stock_250kDB.AnalysisMethod.get(analysis_method_id)
					output_fname_prefix = '%s_%s_%s_%s'%(self.output_fname_prefix, analysis_method.short_name, phenotype_method_id, phenotype_method.short_name)
					snp_annotation_type_id2value_ls = self.getSNPAnnotationType2value_ls(genome_wide_result, snp_annotation_dstruc)
					snp_annotation_type_id2name = self.get_snp_annotation_type_id2name()
					#self.drawHistFromSNPAnnotationType2value_ls(snp_annotation_type_id2value_ls, snp_annotation_type_id2name, output_fname_prefix)
					#for i in range(16):
					#	self.drawHistFromSNPAnnotationType2value_ls(snp_annotation_type_id2value_ls, snp_annotation_type_id2name, '%s_%s'%(output_fname_prefix,i+1), [i+1])
					self.drawPieChart(genome_wide_result, snp_annotation_dstruc, snp_annotation_type_id2name, output_fname_prefix, min_value=self.min_value)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CheckEffectOfSNPType
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()