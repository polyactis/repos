#!/usr/bin/env python

"""
Examples:
	PlotMAFVsScore.py -j 17 -o /tmp/maf_vs_score

Description:
	2008-10-19 program to draw maf(minor allele frequency) vs score plots
"""

import sys, os, math
"""
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
"""
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib as mpl; mpl.use("Agg")
import csv, stat, getopt
import traceback, gc, subprocess
import Stock_250kDB
from matplotlib import rcParams
rcParams['font.size'] = 6
rcParams['legend.fontsize'] = 6
#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
rcParams['axes.labelsize'] = 6
rcParams['axes.titlesize'] = 8
rcParams['xtick.labelsize'] = 6
rcParams['ytick.labelsize'] = 6
import pylab
import StringIO

from GeneListRankTest import GeneListRankTest

class PlotMAFVsScore(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('call_method_id', 0, int):[0, 'j', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							("output_dir", 0, ): [None, 'o', 1, 'directory to store output'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. no transaction. commit on every flush.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  **keywords):
		"""
		2008-10-19
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		if self.output_dir and not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
	
	def plot_maf_vs_score(self, rm, genome_wide_result, output_fname_prefix, commit=0):
		"""
		2008-10-19
		"""
		sys.stderr.write("MAF vs score plot for %s ..."%rm.id)
		pylab.clf()
		maf_ls = []
		score_ls = []
		for data_obj in genome_wide_result.data_obj_ls:
			maf_ls.append(data_obj.maf)
			score_ls.append(data_obj.value)
		
		pylab.plot(maf_ls, score_ls, '.', alpha=0.3, markersize=4)
		pylab.title('MAF vs score for %s on %s (rm.id=%s)'%(rm.analysis_method.short_name, rm.phenotype_method.short_name, rm.id))
		pylab.xlabel('MAF')
		pylab.ylabel('Score')
		
		png_data = None
		svg_data = None
		if commit:
			png_data = StringIO.StringIO()
			#svg_data = StringIO.StringIO()
			pylab.savefig(png_data, format='png', dpi=300)
			#pylab.savefig(svg_data, format='svg', dpi=300)
		elif output_fname_prefix:
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
			#pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		
		sys.stderr.write("Done.\n")
		return png_data, svg_data
		
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
		#session.begin()
		
		query = Stock_250kDB.ResultsMethod.query.filter_by(results_method_type_id=1)
		if self.call_method_id:
			query = query.filter_by(call_method_id=self.call_method_id)
		for rm in query.all():
			rows = Stock_250kDB.MAFVsScorePlot.query.filter_by(results_method_id=rm.id)
			if rows.count()>0:
				sys.stderr.write("Ignore. MAF vs score plot done for %s.\n"%rm.id)
				continue
			gwr = GeneListRankTest.getResultMethodContent(rm, results_directory=self.results_directory, min_MAF=0)
			if self.output_dir:
				output_fname_prefix = os.path.join(self.output_dir, '%s_maf_vs_score'%rm.id)
			else:
				output_fname_prefix = None
			png_data, svg_data = self.plot_maf_vs_score(rm, gwr, output_fname_prefix, commit=self.commit)
			if png_data or svg_data:
				maf_vs_score_plot = Stock_250kDB.MAFVsScorePlot(results_method_id=rm.id)
				if png_data:
					maf_vs_score_plot.png_data = png_data.getvalue()
				if svg_data:
					maf_vs_score_plot.svg_data = svg_data.getvalue()
				session.save(maf_vs_score_plot)
				del png_data, svg_data
				if self.commit:
					session.flush()
		"""
		if self.commit:
			session.flush()
			session.commit()
			session.clear()
		else:	#default is also rollback(). to demonstrate good programming
			session.rollback()
		"""

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PlotMAFVsScore
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()