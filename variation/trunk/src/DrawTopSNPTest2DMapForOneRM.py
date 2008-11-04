#!/usr/bin/env python
"""

Examples:
	#draw the 2D significance map of a result (results_id=2316 emma on LD) on candidate_gene_list 28(FT_short).
	DrawTopSNPTest2DMapForOneRM.py -i 2316 -l 28 -x /tmp/rm_2316_l1_C2.png -g  -C 2

Description:
	2008-10-23
		draw 2D map based on top snp test results (-log pvalue). row is no_of_top_snps(score cutoff). col is minimum distance to associate SNP to genes.
	
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib as mpl; mpl.use("Agg")
from matplotlib import rcParams
rcParams['font.size'] = 6
rcParams['legend.fontsize'] = 4
#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
rcParams['axes.labelsize'] = 4
rcParams['axes.titlesize'] = 6
rcParams['xtick.labelsize'] = 2
rcParams['ytick.labelsize'] = 4
rcParams['lines.linewidth'] = 0.3	#2008-10-31 make the linewidth of boxplot smaller
#rcParams['patch.linewidth'] = 0.5
from matplotlib.backends.backend_agg import FigureCanvasAgg

import getopt, csv, math
import Numeric, cPickle
from pymodule import PassingData, importNumericArray, write_data_matrix, SNPData
import Stock_250kDB
from sets import Set
from pymodule.DrawMatrix import drawMatrix, drawLegend, drawContinousLegend, get_font, combineTwoImages, Value2Color
import pylab

num = importNumericArray()


class DrawTopSNPTest2DMapForOneRM(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("min_distance", 1, int): [20000, 'm', 1, 'minimum distance allowed from the SNP to gene. Deprecated.'],\
							("get_closest", 0, int): [0, 'g', 0, 'only get genes closest to the SNP within that distance'],\
							('min_MAF', 1, float): [0, 'n', 1, 'minimum Minor Allele Frequency.'],\
							('max_pvalue_per_gene', 0, int): [0, 'a', 0, 'take the most significant among all SNPs associated with one gene'],\
							('no_of_top_snps', 1, int): [100, 'f', 1, 'For test_result_type=2. how many number of top snps based on score or -log10(pvalue). Deprecated.'],\
							('results_id', 1, ):[None, 'i', 1, 'data from which ResultsMethod/ResultsByGene to draw'],\
							('list_type_id', 1, ):[None, 'l', 1, 'data from which list_type to look at'],\
							("test_type_id", 1, int): [15, 'y', 1, 'which type of tests. check db table analysis_method. likely be 14,15 etc.'],\
							("results_type", 1, int): [1, 'w', 1, 'data from Which result to output. 1: ResultsMethod, 2: ResultsByGene'],\
							('null_distribution_type_id', 0, int):[1, 'C', 1, 'Type of null distribution. 1=original, 2=permutation, 3=random gene list. in db table null_distribution_type'],\
							("allow_two_sample_overlapping", 1, int): [0, '', 0, 'whether to allow one SNP to be assigned to both candidate and non-candidate gene group'],\
							('font_path', 1, ):['/usr/share/fonts/truetype/freefont/FreeSerif.ttf', 'e', 1, 'path of the font used to draw labels'],\
							('font_size', 1, int):[20, 's', 1, 'size of font, which determines the size of the whole figure.'],\
							("output_fname", 0, ): [None, 'o', 1, 'Filename to store data matrix'],\
							("fig_fname", 1, ): [None, 'x', 1, 'File name for the figure'],\
							("no_of_ticks", 1, int): [15, 't', 1, 'Number of ticks on the legend'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-10-23
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def getTopSNPTestType_id_ls(self, get_closest, min_MAF, allow_two_sample_overlapping, results_type, \
							test_type_id, null_distribution_type_id):
		"""
		"""
		sys.stderr.write("Getting TopSNPTestType_id_ls ...")
		TopSNPTestType_id_ls = []
		rows = Stock_250kDB.CandidateGeneTopSNPTestRMType.query.filter_by(get_closest=get_closest).\
			filter(Stock_250kDB.CandidateGeneTopSNPTestRMType.min_MAF>=min_MAF-0.0001).\
			filter(Stock_250kDB.CandidateGeneTopSNPTestRMType.min_MAF<=min_MAF+0.0001).\
			filter_by(allow_two_sample_overlapping = allow_two_sample_overlapping).\
			filter_by(results_type=results_type).\
			filter_by(test_type_id=test_type_id).\
			filter_by(null_distribution_type_id=1)
		for row in rows:
			TopSNPTestType_id_ls.append(row.id)
		sys.stderr.write("Done.\n")
		return TopSNPTestType_id_ls
	
	def get_no_of_top_snps_info(cls, db, from_where_clause):
		"""
		2008-11-04
			#there's a chance it occurs twice due to float difference in min_score
		2008-10-23
		"""
		sys.stderr.write("Getting no_of_top_snps_info ...")
		rows = db.metadata.bind.execute("select distinct t.no_of_top_snps, t.min_score %s order by no_of_top_snps"%from_where_clause)
		id_ls = []
		id2index = {}
		label_ls = []
		no_of_separators = 0
		for row in rows:
			if row.no_of_top_snps not in id2index:	#there's a chance it occurs twice due to float difference in min_score
				id2index[row.no_of_top_snps] = len(id_ls)
				id_ls.append(row.no_of_top_snps)
				label_ls.append('%s %s'%(row.no_of_top_snps, row.min_score))
		list_info = PassingData()
		list_info.id2index = id2index
		list_info.id_ls = id_ls
		list_info.label_ls = label_ls
		sys.stderr.write("Done.\n")
		return list_info
	get_no_of_top_snps_info = classmethod(get_no_of_top_snps_info)
	
	def get_min_distance_info(cls, db, from_where_clause):
		"""
		2008-10-23
		"""
		sys.stderr.write("Getting min_distance_info ...")
		rows = db.metadata.bind.execute("select distinct t.min_distance %s order by min_distance"%from_where_clause)
		id_ls = []
		id2index = {}
		label_ls = []
		no_of_separators = 0
		for row in rows:
			id2index[row.min_distance] = len(id_ls)
			id_ls.append(row.min_distance)
			label_ls.append('%s'%(row.min_distance))
		list_info = PassingData()
		list_info.id2index = id2index
		list_info.id_ls = id_ls
		list_info.label_ls = label_ls
		sys.stderr.write("Done.\n")
		return list_info
	
	get_min_distance_info = classmethod(get_min_distance_info)
	
	def get_data_matrix(cls, db, row_info, col_info, from_where_clause, need_other_values=False,\
					null_distribution_type_id=2):
		"""
		2008-11-04
			get data_matrix_candidate_sample_size_null & data_matrix_candidate_gw_size_null also if need_other_values=True
		"""
		sys.stderr.write("Getting data matrix ...")
		data_matrix = num.zeros([len(row_info.id2index), len(col_info.id2index)], num.float)
		data_matrix[:] = -1
		max_no_of_null_data = 100
		if need_other_values:
			data_matrix_candidate_sample_size = num.zeros([len(row_info.id2index), len(col_info.id2index)], num.float)
			data_matrix_candidate_sample_size[:] = -1
			data_matrix_non_candidate_sample_size = num.zeros([len(row_info.id2index), len(col_info.id2index)], num.float)
			data_matrix_non_candidate_sample_size[:] = -1
			data_matrix_candidate_gw_size = num.zeros([len(row_info.id2index), len(col_info.id2index)], num.float)
			data_matrix_candidate_gw_size[:] = -1
			data_matrix_non_candidate_gw_size = num.zeros([len(row_info.id2index), len(col_info.id2index)], num.float)
			data_matrix_non_candidate_gw_size[:] = -1
			
			data_matrix_candidate_sample_size_null = num.zeros([len(row_info.id2index), len(col_info.id2index), max_no_of_null_data], num.float)
			data_matrix_candidate_sample_size_null[:] = -1
			data_matrix_candidate_gw_size_null = num.zeros([len(row_info.id2index), len(col_info.id2index), max_no_of_null_data], num.float)
			data_matrix_candidate_gw_size_null[:] = -1
		else:
			data_matrix_candidate_sample_size = None
			data_matrix_non_candidate_sample_size = None
			data_matrix_candidate_gw_size = None
			data_matrix_non_candidate_gw_size = None
			data_matrix_candidate_sample_size_null = None
			data_matrix_candidate_gw_size_null = None
		rows = db.metadata.bind.execute("select t.id, t.no_of_top_snps, t.min_distance, t.pvalue, t.candidate_sample_size, \
			t.non_candidate_sample_size, t.candidate_gw_size, t.non_candidate_gw_size %s"%from_where_clause)
		min_value = None
		max_value = None
		for row in rows:
			row_index = row_info.id2index[row.no_of_top_snps]
			col_index = col_info.id2index[row.min_distance]
			if row.pvalue>0:
				data_value = -math.log10(row.pvalue)
				if min_value==None:
					min_value = data_value
				elif data_value<min_value:
					min_value = data_value
				
				if max_value==None:
					max_value=data_value
				elif data_value>max_value:
					max_value =data_value
			else:
				data_value = -2	#0 pvalue
			data_matrix[row_index, col_index] = data_value
			if need_other_values:
				data_matrix_candidate_sample_size[row_index, col_index] = row.candidate_sample_size
				data_matrix_non_candidate_sample_size[row_index, col_index] = row.non_candidate_sample_size
				data_matrix_candidate_gw_size[row_index, col_index] = row.candidate_gw_size
				data_matrix_non_candidate_gw_size[row_index, col_index] = row.non_candidate_gw_size
				null_datas = db.metadata.bind.execute("select candidate_sample_size, candidate_gw_size from %s where observed_id=%s and null_distribution_type_id=%s"%\
											(Stock_250kDB.TopSNPTestRMNullData.table.name, row.id, null_distribution_type_id))
				i = 0
				for null_data in null_datas:
					data_matrix_candidate_sample_size_null[row_index, col_index, i] = null_data.candidate_sample_size
					data_matrix_candidate_gw_size_null[row_index, col_index, i] = null_data.candidate_gw_size
					i+=1
					if i>=max_no_of_null_data:	#no more than this
						break
		sys.stderr.write("Done.\n")
		return_data = PassingData()
		return_data.data_matrix = data_matrix
		return_data.data_matrix_candidate_sample_size = data_matrix_candidate_sample_size
		return_data.data_matrix_non_candidate_sample_size = data_matrix_non_candidate_sample_size
		return_data.data_matrix_candidate_gw_size = data_matrix_candidate_gw_size
		return_data.data_matrix_non_candidate_gw_size = data_matrix_non_candidate_gw_size
		return_data.min_value = min_value
		return_data.max_value = max_value
		return_data.data_matrix_candidate_sample_size_null = data_matrix_candidate_sample_size_null
		return_data.data_matrix_candidate_gw_size_null = data_matrix_candidate_gw_size_null
		return return_data
	get_data_matrix = classmethod(get_data_matrix)
	
	def plotSubCurve(cls, rdata, no_of_top_snps_info, min_distance_info, which_min_distance, no_of_rows=2, legend_ls=[], \
					patch_ls=[], pd=None):
		"""
		2008-10-28
			called by plotCurve()
		"""
		plot_kw = {'linewidth':0.5,\
				'markerfacecolor': 'w',\
				'markersize': 1}
		if len(legend_ls)==0:
			fill_in_legend = True
		else:
			fill_in_legend = False
		which_figure = 2*which_min_distance + 1
		if pd.left_1st_ax is None:
			pd.left_1st_ax = pylab.subplot(no_of_rows, 2, which_figure, frameon=False)
		else:
			pylab.subplot(no_of_rows, 2, which_figure, frameon=False, sharex=pd.left_1st_ax)
		pylab.title('min-distance=%s'%min_distance_info.label_ls[which_min_distance])
		pylab.xlabel('score cutoff')
		pylab.grid(True, alpha=0.3)
		score_cutoff_ls = []
		candidate_ratio_ls = []
		non_candidate_ratio_ls = []
		candidate_vs_non_candidate_ratio_ls = []
		pvalue_ls = []
		boxplot_data_ls = []
		score_cutoff_for_boxplot_ls = []
		emp_pvalue = []
		for i in range(rdata.data_matrix.shape[0]):
			if rdata.data_matrix[i][which_min_distance]>=0:	#not NA (missing)
				min_score = float(no_of_top_snps_info.label_ls[i].split(' ')[1])
				if pd.min_x==None:
					pd.min_x = min_score
				elif min_score<pd.min_x:
					pd.min_x = min_score
				if pd.max_x is None:
					pd.max_x = min_score
				elif min_score>pd.max_x:
					pd.max_x = min_score
				score_cutoff_ls.append(min_score)
				candidate_sample_size = rdata.data_matrix_candidate_sample_size[i][which_min_distance]
				candidate_gw_size = rdata.data_matrix_candidate_gw_size[i][which_min_distance]
				candidate_ratio = candidate_sample_size/candidate_gw_size
				candidate_ratio_ls.append(candidate_ratio)
				
				non_candidate_ratio = rdata.data_matrix_non_candidate_sample_size[i][which_min_distance]/rdata.data_matrix_non_candidate_gw_size[i][which_min_distance]
				non_candidate_ratio_ls.append(non_candidate_ratio)
				
				candidate_vs_non_candidate_ratio_ls.append(candidate_ratio/non_candidate_ratio)
				pvalue_ls.append(rdata.data_matrix[i][which_min_distance])
				
				no_of_top_snps = rdata.data_matrix_candidate_sample_size[i][which_min_distance] + rdata.data_matrix_non_candidate_sample_size[i][which_min_distance]
				no_of_total_snps = rdata.data_matrix_candidate_gw_size[i][which_min_distance] + rdata.data_matrix_non_candidate_gw_size[i][which_min_distance]
				boxplot_data = []
				no_of_hits = 0 
				for k in range(rdata.data_matrix_candidate_sample_size_null.shape[2]):
					candidate_sample_size_null = rdata.data_matrix_candidate_sample_size_null[i][which_min_distance][k]
					candidate_gw_size_null = rdata.data_matrix_candidate_gw_size_null[i][which_min_distance][k]
					if candidate_sample_size_null==-1:	#NA from this point on
						break
					
					non_candidate_sample_size_null = no_of_top_snps-candidate_sample_size_null
					non_candidate_gw_size_null = no_of_total_snps - candidate_gw_size_null
					if candidate_gw_size_null==0 or non_candidate_sample_size_null<=0 or non_candidate_gw_size_null==0:	#will cause zero denominator
						continue
					candidate_ratio_null = candidate_sample_size_null/candidate_gw_size_null
					non_candidate_ratio_null = non_candidate_sample_size_null/non_candidate_gw_size_null
					boxplot_data.append(candidate_ratio_null/non_candidate_ratio_null)
					if candidate_ratio_null>=candidate_ratio:
						no_of_hits += 1
				if len(boxplot_data)>2:
					score_cutoff_for_boxplot_ls.append(min_score)
					boxplot_data_ls.append(boxplot_data)
					if no_of_hits==0:
						emp_pvalue.append(4)
					else:
						emp_pvalue.append(-math.log(no_of_hits/float(k)))
		pylab.ylabel("candidate/non-candidate")
		if len(boxplot_data_ls)>0:
			pylab.boxplot(boxplot_data_ls, notch=1, positions=score_cutoff_for_boxplot_ls, widths=0.05)
		
		a = pylab.plot(score_cutoff_ls, candidate_vs_non_candidate_ratio_ls, 'g-', **plot_kw)
		
		if fill_in_legend:
			legend_ls.append("candidate-ratio/non-candidate-ratio")
			patch_ls.append(a[0])

		
		#draw pvalue in the same axe but with scale on the right side
		ax2 = pylab.twinx()
		ax2.set_ylabel('-log10(pvalue)')
		a = ax2.plot(score_cutoff_ls, pvalue_ls, 'r-', **plot_kw)
		if fill_in_legend:
			legend_ls.append("-log10(pvalue)")
			patch_ls.append(a[0])
		#2008-10-31 draw the empirical pvalue
		a = pylab.plot(score_cutoff_for_boxplot_ls, emp_pvalue, 'c-', **plot_kw)
		if fill_in_legend:
			legend_ls.append("empirical pvalue")
			patch_ls.append(a[0])
		
		#2008-10-28 draw the right side
		which_figure = 2*which_min_distance + 2
		if pd.right_1st_ax is None:
			pd.right_1st_ax = pylab.subplot(no_of_rows, 2, which_figure, frameon=False)
		else:
			pylab.subplot(no_of_rows, 2, which_figure, frameon=False, sharex=pd.right_1st_ax)
		pylab.grid(True, alpha=0.3)
		pylab.title('min-distance=%s'%min_distance_info.label_ls[which_min_distance])
		pylab.xlabel('score cutoff')
		pylab.ylabel("Ratio")
		a = pylab.plot(score_cutoff_ls, candidate_ratio_ls, 'co-', **plot_kw)
		if fill_in_legend:
			legend_ls.append("candidate ratio")
			patch_ls.append(a[0])
		a = pylab.plot(score_cutoff_ls, non_candidate_ratio_ls, 'bo-', **plot_kw)
		if fill_in_legend:
			legend_ls.append("non-candidate ratio")
			patch_ls.append(a[0])
		
		ax3 = pylab.twinx()
		ax3.plot(score_cutoff_ls, pvalue_ls, 'ro-', **plot_kw)
		ax3.set_ylabel('-log10(pvalue)')
		return pd
	plotSubCurve=classmethod(plotSubCurve)
		
	def plotCurve(cls, rdata, no_of_top_snps_info, min_distance_info, output_fname=None, need_svg=False):
		"""
		2008-10-29
			draw candidate ratio, non-candidate ratio, pvalues etc against score cutoff at different distances.
			
		"""
		sys.stderr.write("Plotting curves ... ")
		
		pylab.clf()
		#fig = Figure()
		fig = pylab.gcf()
		#canvas = FigureCanvasAgg(fig)
		
		#calculate the number of rows needed according to how many score_rank_data, always two-column
		pylab.subplots_adjust(left=0.08, right=0.92,bottom = 0.05, hspace = 0.3)
		no_of_rows = rdata.data_matrix.shape[1]
		legend_ls = []
		patch_ls = []
		pd = PassingData(left_1st_ax = None,
								right_1st_ax = None,
								min_x = None,
								max_x = None)
		for i in range(rdata.data_matrix.shape[1]):
			pd = cls.plotSubCurve(rdata, no_of_top_snps_info, min_distance_info, i, no_of_rows=no_of_rows, \
							legend_ls=legend_ls, patch_ls=patch_ls, pd=pd)
		
		#put a main title and single legend for all plots
		ax = pylab.axes([0.1, 0.1, 0.9, 0.9], frameon=False)
		ax.set_xticks([])
		ax.set_yticks([])
		ax.legend(patch_ls, legend_ls, loc='upper right', handlelen=0.02)
		pd.left_1st_ax.set_xlim([pd.min_x-0.2, pd.max_x+0.2])
		pd.right_1st_ax.set_xlim([pd.min_x-0.2, pd.max_x+0.2])
		"""
		#put a main title and single legend for all plots
		ax = pylab.axes([0.1, 0.1, 0.8, 0.85], frameon=False)
		ax.set_xticks([])
		ax.set_yticks([])
		title = '%s by %s'%(self.results_id, self.list_type_id)
		ax.set_title(title)
		"""
		
		if output_fname:
			pylab.savefig('%s'%output_fname, dpi=300)
			if need_svg:
				pylab.savefig('%s.svg'%output_fname, dpi=300)
		sys.stderr.write("Done.\n")
	
	plotCurve = classmethod(plotCurve)
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup()
		session = db.session
		
		TopSNPTestType_id_ls = self.getTopSNPTestType_id_ls(self.get_closest, self.min_MAF, self.allow_two_sample_overlapping, self.results_type, \
							self.test_type_id, self.null_distribution_type_id)
		if not TopSNPTestType_id_ls:
			sys.stderr.write("No TopSNPTestType matches the input requirements. Exit.\n")
			sys.exit(3)
		TopSNPTestType_id_ls_str = map(str, TopSNPTestType_id_ls)
		from_where_clause = "from %s t, %s y where t.type_id=y.id and t.results_id=%s and t.list_type_id=%s and y.id in (%s)"%\
			(Stock_250kDB.CandidateGeneTopSNPTestRM.table.name, Stock_250kDB.CandidateGeneTopSNPTestRMType.table.name,\
			self.results_id, self.list_type_id, ','.join(TopSNPTestType_id_ls_str))
		
		no_of_top_snps_info = self.get_no_of_top_snps_info(db, from_where_clause)
		min_distance_info = self.get_min_distance_info(db, from_where_clause)
		rdata = self.get_data_matrix(db, no_of_top_snps_info, min_distance_info, from_where_clause, need_other_values=True, \
									null_distribution_type_id=self.null_distribution_type_id)
		
		header = ['no_of_top_snps', ''] + min_distance_info.label_ls
		strain_acc_list = no_of_top_snps_info.label_ls
		category_list = no_of_top_snps_info.label_ls
		
		if SNPData.isDataMatrixEmpty(rdata.data_matrix):
			sys.stderr.write("Nothing fetched from database.\n")
			sys.exit(3)
		if self.output_fname:
			write_data_matrix(rdata.data_matrix, self.output_fname, header, strain_acc_list, category_list)
		
		if self.fig_fname:
			font = get_font(self.font_path, font_size=self.font_size)	#2008-08-01
			value2color_func = lambda x: Value2Color.value2HSLcolor(x, rdata.min_value, rdata.max_value)
			im_legend = drawContinousLegend(rdata.min_value, rdata.max_value, self.no_of_ticks, value2color_func, font)
			#im.save('%s_legend.png'%self.fig_fname_prefix)
			im = drawMatrix(rdata.data_matrix, value2color_func, no_of_top_snps_info.label_ls,\
						min_distance_info.label_ls, with_grid=1, font=font)
			im = combineTwoImages(im, im_legend, font=font)
			im.save(self.fig_fname)
		
		output_fname_prefix='%s_curve.png'%os.path.splitext(self.fig_fname)[0]
		self.plotCurve(rdata, no_of_top_snps_info, min_distance_info, output_fname_prefix)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DrawTopSNPTest2DMapForOneRM
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()