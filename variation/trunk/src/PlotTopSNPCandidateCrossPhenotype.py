#!/usr/bin/env python
"""

Examples:
	
	#emma on all flowering phenotype. candidate-ratio/non-candidate-ratio
	PlotTopSNPCandidateCrossPhenotype.py -m 1000 -q 7 -l 28 -A 1-7,39-61,80-82, -x /Network/Data/250k/tmp-yh/TopSNPCandidate_m1000_q7_FT_C1B1.png -C1 -B 1
	
	#emma on all disease phenotypes on list 130 (disease resistance) with output in both picture and table matrix format.
	m=5000; q=7; l=130;
	PlotTopSNPCandidateCrossPhenotype.py  -m$m -q $q -l $l -A  9-13,32-38,65-74 
	-x /Network/Data/250k/tmp-yh/TopSNPCandidateCrossPhenotype/figures/TopSNPCandidate_m$m\_q$q\_disease_phenotype_list$l\.png 
	-o /Network/Data/250k/tmp-yh/TopSNPCandidateCrossPhenotype/Disease_m$m\_q$q\_list$l\_ratio.tsv -C1 -B1
	
	#2009-10-4 wilcox (q=1) on all flowering time, call method 32, gw-looping permutation pvalue (C=2), enrichment ratio (B=1)
	m=20000; q=1; l=145; C=2; B=1; j=32; ~/script/variation/src/PlotTopSNPCandidateCrossPhenotype.py -m$m -e $q -l $l
	-A 1-7,39-59,80-82 
	-x /tmp/TopSNPCandidate_m$m\_q$q\_flower_phenotype_j$j\_list$l\_C$C\B$B\.png 
	-o /tmp/Flower_j$j\_m$m\_q$q\_list$l\_C$C\B$B\.tsv -C$C -B$B -j $j -u yh
	
Description:
	2008-11-11
		draw 3D bar chart showing candidate_sample_size, candidate-ratio or other variables along the cutoff and phenotype axises.
	
	if output_fname is specified, matrix of all data will be dumped.
	
	if analysis_method_id=6(RF), score_cutoff_take_log is set to True.
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import pylab

from DrawTopSNPTest2DMapForOneRM import DrawTopSNPTest2DMapForOneRM
from matplotlib import rcParams
rcParams['font.size'] = 4
rcParams['legend.fontsize'] = 8
#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
rcParams['axes.labelsize'] = 4
rcParams['axes.titlesize'] = 4
rcParams['xtick.labelsize'] = 4
rcParams['ytick.labelsize'] = 4
rcParams['lines.linewidth'] = 0.3	#2008-10-31 make the linewidth of boxplot smaller

import csv
import numpy
from pymodule import PassingData, getListOutOfStr, write_data_matrix
import StringIO
from MpiGeneListRankTest import MpiGeneListRankTest
import matplotlib.axes3d as axes3d
import Stock_250kDB
from sets import Set

class PlotTopSNPCandidateCrossPhenotype(DrawTopSNPTest2DMapForOneRM):
	__doc__ = __doc__
	option_default_dict = DrawTopSNPTest2DMapForOneRM.option_default_dict.copy()
	option_default_dict.update({('data_type', 1, int): [1, 'B', 1, 'which tupe of data to draw. 1:ratio, 2:pvalue, 3:candidate_sample_size, 4:non_candidate_sample_size']})
	def __init__(self,  **keywords):
		DrawTopSNPTest2DMapForOneRM.__init__(self, **keywords)
	
	def init_3d_plot(self):
		pylab.clf()
		ax = axes3d.Axes3D(pylab.figure(), azim=60,elev=50)	#pylab.gcf() doesn't work. have to use pylab.figure()
		ax.set_xlabel('cutoff')
		ax.set_ylabel('result')
		ax.set_zlabel('data')
		return ax
	
	color_ls = ['b', 'g','r', 'c', 'm', 'y']
	def plot_one_bar(self, ax, rdata, no_of_top_snps_info, min_distance_info, min_distance, result_index=0, data_type=1, output_fname=None, \
							need_svg=False, title='', commit=0, preset_xlim =None, max_no_of_jump_pts=5, score_cutoff_take_log=False):
		"""
		2008-11-13
			set alpha=0.8 (was 0.5)
			set linewidth=0.5 (was 0)
		"""
		if min_distance not in min_distance_info.id2index:
			return False
		which_min_distance = min_distance_info.id2index[min_distance]
		score_cutoff_ls = []
		candidate_ratio_ls = []
		non_candidate_ratio_ls = []
		candidate_vs_non_candidate_ratio_ls = []
		pvalue_ls = []
		boxplot_data_ls = []
		score_cutoff_for_boxplot_ls = []
		emp_pvalue = []
		
		score_cutoff_jump_pt_ls = []
		candidate_sample_size_jump_pt_ls = []
		prev_candidate_sample_size = None
		
		candidate_sample_size_ls = []
		non_candidate_sample_size_ls = []
		for i in range(rdata.data_matrix.shape[0]):
			if rdata.data_matrix[i][which_min_distance]>=0:	#not NA (missing)
				min_score = float(no_of_top_snps_info.label_ls[i].split(' ')[1])
				if score_cutoff_take_log:
					min_score = math.log10(min_score+1)
				"""
				if pd.min_x==None:
					pd.min_x = min_score
				elif min_score<pd.min_x:
					pd.min_x = min_score
				if pd.max_x is None:
					pd.max_x = min_score
				elif min_score>pd.max_x:
					pd.max_x = min_score
				"""
				candidate_sample_size = rdata.data_matrix_candidate_sample_size[i][which_min_distance]
				if prev_candidate_sample_size==None:
					score_cutoff_jump_pt_ls.append(min_score)
					candidate_sample_size_jump_pt_ls.append(candidate_sample_size)
					prev_candidate_sample_size = candidate_sample_size
				else:
					if candidate_sample_size!=prev_candidate_sample_size:
						if len(score_cutoff_jump_pt_ls)<max_no_of_jump_pts:
							score_cutoff_jump_pt_ls.append(min_score)
							candidate_sample_size_jump_pt_ls.append(candidate_sample_size)
						prev_candidate_sample_size = candidate_sample_size
						
				candidate_gw_size = rdata.data_matrix_candidate_gw_size[i][which_min_distance]
				candidate_ratio = candidate_sample_size/candidate_gw_size
				
				non_candidate_sample_size = rdata.data_matrix_non_candidate_sample_size[i][which_min_distance]
				non_candidate_ratio = non_candidate_sample_size/rdata.data_matrix_non_candidate_gw_size[i][which_min_distance]
				
				log_pvalue = rdata.data_matrix[i][which_min_distance]
				if log_pvalue==-2:	#-2 is put there when pvalue=0 and can't do log10
					log_pvalue = 5
				
				if non_candidate_ratio!=0:	#if non_candidate_ratio is 0. candidate_ratio/non_candidate_ratio is inf. cause plot error.
					candidate_ratio_ls.append(candidate_ratio)
					non_candidate_ratio_ls.append(non_candidate_ratio)
					candidate_vs_non_candidate_ratio_ls.append(candidate_ratio/non_candidate_ratio)
					score_cutoff_ls.append(min_score)
					pvalue_ls.append(log_pvalue)
					candidate_sample_size_ls.append(candidate_sample_size)
					non_candidate_sample_size_ls.append(non_candidate_sample_size)
					
		if data_type==1:
			data_ls = candidate_vs_non_candidate_ratio_ls
			
		elif data_type==2:
			data_ls = pvalue_ls
		elif data_type==3:
			data_ls = candidate_sample_size_ls
		elif data_type==4:
			data_ls = non_candidate_sample_size_ls
		
		if score_cutoff_ls:
			color_index = result_index%len(self.color_ls)
			if score_cutoff_take_log:
				bar_width = (max(score_cutoff_ls)-min(score_cutoff_ls))/float(len(score_cutoff_ls))
			else:
				bar_width = 0.1
			ax.bar(score_cutoff_ls, data_ls, z=result_index, dir='y',color=self.color_ls[color_index], width=bar_width, alpha=0.8, linewidth=0.5)
			
			#ax.set_yticks([30,20,10,0], ['30z', '20z', '10z', '0z'])	#2008-11-11 it reduces 3d into 1 line
			ax.text3D(0,result_index, 0, title, size=4)	#add title for this line of bars
			return (score_cutoff_ls, data_ls)
		else:
			return False
	
	def get_score_cutoff2index(self, data_to_output_ls):
		"""
		2008-11-11
			get a global score_cutoff2index out of all score_cutoff in data_to_output_ls
		"""
		score_cutoff_set = Set()
		for score_cutoff_ls, data_ls in data_to_output_ls:
			score_cutoff_set.update(Set(score_cutoff_ls))
		
		score_cutoff_ls = list(score_cutoff_set)
		score_cutoff_ls.sort()	#in ascending
		score_cutoff_ls.reverse()	#in descending
		score_cutoff2index = {}
		for i in range(len(score_cutoff_ls)):
			score_cutoff2index[score_cutoff_ls[i]] = i	#descending order. bigger score_cutoffs, lower index
		return score_cutoff_ls, score_cutoff2index
	
	def output_data(self, data_to_output_label_ls, data_to_output_ls, min_distance, output_fname):
		"""
		2008-11-11
			data_to_output_ls is a list of (score_cutoff_ls, data_ls). each score_cutoff_ls might be a bit different from each other.
			1. get score_cutoff2index out of all score_cutoffs in descending order
			2. each row is same score_cutoff. column is data_ls of one result from analysis_method on phenotype.
			3. first column is score cutoffs.
			4. 2nd column is min_distance. 3rd and so forth columns are data.
		"""
		sys.stderr.write("Outputting data matrix ...")
		score_cutoff_ls, score_cutoff2index = self.get_score_cutoff2index(data_to_output_ls)
		header = ['score_cutoff', 'min_distance'] + data_to_output_label_ls
		no_of_cols = len(data_to_output_label_ls)
		data_matrix = numpy.zeros([len(score_cutoff2index), no_of_cols], numpy.float)
		data_matrix[:] = -1
		for j in range(no_of_cols):
			sub_score_cutoff_ls, data_ls = data_to_output_ls[j]
			for i in range(len(sub_score_cutoff_ls)):
				score_cutoff = sub_score_cutoff_ls[i]
				data = data_ls[i]
				row_index = score_cutoff2index[score_cutoff]
				data_matrix[row_index][j] = data
		
		category_list = [min_distance]*len(score_cutoff2index)
		write_data_matrix(data_matrix, output_fname, header, score_cutoff_ls, category_list)
		sys.stderr.write("Done.\n")
		
	def run(self):
		"""
		2008-11-11
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup()
		session = db.session
		
		param_obj = PassingData(call_method_id=self.call_method_id, \
								analysis_method_id=getattr(self, 'analysis_method_id', None),\
								analysis_method_id_ls=getattr(self, 'analysis_method_id_ls', None),\
								phenotype_method_id_ls=getattr(self, 'phenotype_method_id_ls', None),\
								list_type_id_ls=self.list_type_id_ls, \
								results_type=self.results_type)
		params_ls = MpiGeneListRankTest.generate_params(param_obj)
		
		ax = self.init_3d_plot()
		result_index = 0
		data_to_output_label_ls = []
		data_to_output_ls = []
		for results_id, list_type_id in params_ls:
			rm = Stock_250kDB.ResultsMethod.get(results_id)
			list_type = Stock_250kDB.GeneListType.get(list_type_id)
			title = '%s on %s_%s (%s)'%\
				(rm.analysis_method.short_name, rm.phenotype_method.id, rm.phenotype_method.short_name, results_id)
			#a short output label
			output_label = '%s_%s (%s)'%\
				(rm.phenotype_method.id, rm.phenotype_method.short_name, results_id)
			phenotype_label = '%s_%s'%\
				(rm.phenotype_method.id, rm.phenotype_method.short_name)
			TopSNPTestType_id_ls = self.getTopSNPTestType_id_ls(self.get_closest, self.min_MAF, self.allow_two_sample_overlapping, self.results_type, \
								self.test_type_id, self.null_distribution_type_id)
			if self.commit:
				rows = Stock_250kDB.CandidateVsNonRatioPlot.query.filter_by(type_id=TopSNPTestType_id_ls[0]).\
					filter_by(results_id=results_id).filter_by(list_type_id=list_type_id)
				if rows.count()>0:
					row = rows.first()
					sys.stderr.write('%s already in db (%s of them) with first id=%s.\n'%(title, rows.count(), row.id))
					continue
			
			if not TopSNPTestType_id_ls:
				sys.stderr.write("No TopSNPTestType matches the input requirements. Exit.\n")
				sys.exit(3)
			TopSNPTestType_id_ls_str = map(str, TopSNPTestType_id_ls)
			from_where_clause = "from %s t, %s y where t.type_id=y.id and t.results_id=%s and t.list_type_id=%s and y.id in (%s)"%\
				(Stock_250kDB.CandidateGeneTopSNPTestRM.table.name, Stock_250kDB.CandidateGeneTopSNPTestRMType.table.name,\
				results_id, list_type_id, ','.join(TopSNPTestType_id_ls_str))
			
			no_of_top_snps_info = self.get_no_of_top_snps_info(db, from_where_clause)
			min_distance_info = self.get_min_distance_info(db, from_where_clause)
			rdata = self.get_data_matrix(db, no_of_top_snps_info, min_distance_info, from_where_clause, need_other_values=True, \
										null_distribution_type_id=self.null_distribution_type_id)
			if rm.analysis_method_id=='6':	#For random forest, take log and determine bar length according to score_cutoff_ls
				score_cutoff_take_log = True
			else:
				score_cutoff_take_log = False
			return_code = self.plot_one_bar(ax, rdata, no_of_top_snps_info, min_distance_info, self.min_distance, result_index=result_index, data_type=self.data_type, \
							output_fname=None, \
							need_svg=False, title=phenotype_label, commit=0, preset_xlim =None, score_cutoff_take_log=score_cutoff_take_log)
			if return_code:
				data_to_output_label_ls.append(output_label)
				data_to_output_ls.append(return_code)
				result_index += 1
		if self.fig_fname:
			pylab.savefig(self.fig_fname, dpi=300)
			pylab.savefig('%s.svg'%os.path.splitext(self.fig_fname)[0], dpi=300)
		#pylab.show()
		if self.output_fname and data_to_output_ls:
			self.output_data(data_to_output_label_ls, data_to_output_ls, self.min_distance, self.output_fname)
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PlotTopSNPCandidateCrossPhenotype
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()