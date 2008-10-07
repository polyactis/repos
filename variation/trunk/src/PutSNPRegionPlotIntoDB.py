#!/usr/bin/env python
"""

Examples:
	PutSNPRegionPlotIntoDB.py -i /Network/Data/250k/tmp-yh/snp_region_phenotype_32_c4_f200  -l Phenotype32Top200ShareBy4  -j /Network/Data/250k/tmp-yh/at_gene_model_pickelf  -c
	
Description:
	2008-10-06 Program to put plots outputted by DrawSNPRegion.py into database.
	
	After base64 encoding, the file size should not exceed 134217728(~128MB).
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math
import cPickle
from pymodule import PassingData, importNumericArray
import Stock_250kDB
from sets import Set
from DrawSNPRegion import DrawSNPRegion
import base64

class PutSNPRegionPlotIntoDB(DrawSNPRegion):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("input_dir", 1, ): [None, 'i', 1, 'directory containing the plot files'],\
							("plot_type_short_name", 1, ): [None, '', 1, 'A short name (<256 chars) to characterize the type of all the plots in input_dir'],\
							('distance_flanking_center_snp', 1, int): [20000, '', 0, 'The distance flanking the center SNP. used only if the filename has the center snp info.'],\
							("gene_annotation_picklef", 0, ): [None, 'j', 1, 'given the option, If the file does not exist yet, store a pickled gene_annotation into it. If the file exists, load gene_annotation out of it.'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-10-06
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def save_plot2gene(self, session, snp_region_plot, chromosome, start, stop, gene_annotation):
		"""
		2008-10-06
		"""
		no_of_genes = 0
		chromosome = str(chromosome)	#in gene_annotation, chromosome is in text format to include mitochondria etc.
		for gene_id in gene_annotation.chr_id2gene_id_ls[chromosome]:
			gene_model = gene_annotation.gene_id2model[gene_id]
			if gene_model.start!=None and gene_model.stop!=None and gene_model.stop>start:
				if gene_model.start>stop:	#totally out of range, skip it
						continue
				plot2gene = Stock_250kDB.SNPRegionPlotToGene(gene_id=gene_id)
				plot2gene.snp_region_plot = snp_region_plot
				session.save(plot2gene)
				no_of_genes += 1
		return no_of_genes
		
	def save_one_plot(self, session, filename, plot_type, gene_annotation, distance_flanking_center_snp):
		"""
		2008-10-06
		"""
		sys.stderr.write("Saving plot ...")
		file_basename = os.path.split(filename)[1]
		basename_ls = file_basename.split('_')
		chromosome = int(basename_ls[4])
		if basename_ls[6]=='phenotype':	#filename only contains the center snp position
			center_snp_position = int(basename_ls[5])
			start = center_snp_position-distance_flanking_center_snp
			stop = center_snp_position + distance_flanking_center_snp - 1	#follow DrawSNPRegion.getSNPsAroundThisSNP()
			phenotype_method_id = int(basename_ls[7])
		else:	#there's one more position, stop.
			center_snp_position = None
			start = int(basename_ls[5])
			stop = int(basename_ls[6])
			phenotype_method_id = int(basename_ls[8])
		
		snp_region_plot = Stock_250kDB.SNPRegionPlot(chromosome=chromosome, start=start, stop=stop,\
													center_snp_position=center_snp_position, phenotype_method_id=phenotype_method_id)
		snp_region_plot.original_filename = filename
		snp_region_plot.plot_type = plot_type
		f = open(filename, 'rb')
		
		snp_region_plot.img_data = base64.b64encode(f.read())
		session.save(snp_region_plot)
		
		no_of_genes = self.save_plot2gene(session, snp_region_plot, chromosome, start, stop, gene_annotation)
		sys.stderr.write("%s genes. Done.\n"%(no_of_genes))
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		session.begin()
		
		gene_annotation = self.dealWithGeneAnnotation(self.gene_annotation_picklef)
		
		plot_type = Stock_250kDB.SNPRegionPlotType.query.filter_by(short_name=self.plot_type_short_name).first()
		if not plot_type:
			plot_type = Stock_250kDB.SNPRegionPlotType(short_name=self.plot_type_short_name)
			session.save(plot_type)
		
		input_files = os.listdir(self.input_dir)
		no_of_files = len(input_files)
		sys.stderr.write("\tTotally, %d files.\n"%(no_of_files))
		i=0
		for fname in input_files:
			i += 1
			if os.path.splitext(fname)[1] =='.png':
				sys.stderr.write("%d/%d:\t%s\n"%(i, no_of_files, fname))
				fname_path = os.path.join(self.input_dir, fname)
				self.save_one_plot(session, fname_path, plot_type, gene_annotation, self.distance_flanking_center_snp)
		
		if self.commit:
			session.flush()
			session.commit()
			session.clear()
		else:	#default is also rollback(). to demonstrate good programming
			session.rollback()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PutSNPRegionPlotIntoDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()