#!/usr/bin/env python
"""
Examples:
	./src/PlotQCCall.py -i /mnt/nfs/NPUTE_data/QC_avg_output

Description:
	program to plot MpiQCCall.py's avg.csv output.
"""

import sys, os, math
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))


class PlotQCCall(object):
	keyword_candidate_ls = ['strain or snp', 'after_imputation', 'min_call_probability', 'max_call_mismatch_rate', 'max_call_NA_rate', 'max_snp_mismatch_rate',
						'max_snp_NA_rate', 'npute_window_size']
	candidate_plot_var_name_ls = keyword_candidate_ls[1:]+ ['avg_mismatch_rate', 'avg_NA_rate', \
								'no_of_total_accessions_filtered', 'no_of_total_snps_filtered', 'no_of_total_snps_removed', 'no_of_monomorphic_snps_removed']
	__doc__ = __doc__  + 'candidates for keywords: ' + repr(keyword_candidate_ls) + '\nCandidate plot variables: ' + repr(candidate_plot_var_name_ls)
	
	option_default_dict = {('input_dir', 1, ): ['/mnt/nfs/NPUTE_data/QC_avg_output', 'i', 1, ],\
							('output_fname', 0, ): [None, 'o', 1, 'not given, the program would construct an output filename based on arguments.', ],\
							('var_to_plot_indices', 1, ):['127', 't', 1, 'index corresponding to candidate plot variable list. -h to see the list'],\
							('strain or snp', 1,): ['strain', 's', 1,],\
							('after_imputation', 0, ): [1, 'f'],\
							('min_call_probability', 0, ): [None, 'y', 1, 'minimum probability for a call to be non-NA if there is a 3rd column for probability.', ],\
							('max_call_mismatch_rate', 0, ): [None, 'x', 1, 'maximum mismatch rate of an array call_info entry. used to exclude bad arrays.'],\
							('max_call_NA_rate', 0, ): [0.4, 'a', 1, 'maximum NA rate of an array call_info entry. used to exclude bad arrays.'],\
							('max_snp_mismatch_rate', 0, ): [0.2, 'w', 1, 'maximum snp error rate, used to exclude bad SNPs', ],\
							('max_snp_NA_rate', 0, ): [0.4, 'v', 1, 'maximum snp NA rate, used to exclude SNPs with too many NAs', ],\
							('npute_window_size', 0, ): [50, 'n', 1, 'NPUTE window size',],\
							('plot_type', 1, int): [1, 'p', 1, '1=plot3D, 2=contour',],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-05-13
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		plot_func_dict = {1:PlotQCCall._plot3D,\
						2:PlotQCCall._plotContour}
		self._plot = plot_func_dict[self.plot_type]
		
		
	#2008-05-13 read in MpiQCCall.py output
	def get_qccall_results(self, input_dir):
		import os,sys,csv
		from variation.src.MpiQCCall import MpiQCCall
		from pymodule import PassingData
		var_name_ls = ['strain or snp', 'after_imputation'] + MpiQCCall.common_var_name_ls
		avg_var_name_pair_ls, partial_header_avg = MpiQCCall.generate_avg_variable_names(MpiQCCall.avg_var_name_ls)
		var_name_ls += partial_header_avg
		
		files = os.listdir(input_dir)
		passingdata_ls = []
		no_of_objects = len(files)
		for i in range(no_of_objects):
			sys.stderr.write("\t%d/%d: from %s ... \n"%(i+1, no_of_objects, files[i]))
			filename = os.path.join(input_dir, files[i])
			reader = csv.reader(open(filename))
			try:
				reader.next()
			except:
				import traceback
				sys.stderr.write('\terror in reading this file. ignored.\n')
				traceback.print_exc()
	  			print sys.exc_info()
	  			del reader
				continue
			for row in reader:
				passingdata = PassingData()
				for i in range(len(var_name_ls)):
					var_name = var_name_ls[i]
					if i!=0:
						value = float(row[i])
					else:	#the first column is strain or snp, no float conversion
						value = row[i]
					setattr(passingdata, var_name, value)	#
				#two new variables record no of accessions/snps lost
				passingdata.no_of_total_accessions_filtered = passingdata.no_of_accessions_filtered_by_mismatch + passingdata.no_of_accessions_filtered_by_na
				passingdata.no_of_total_snps_filtered = passingdata.no_of_snps_filtered_by_mismatch +\
					passingdata.no_of_snps_filtered_by_na
				passingdata.no_of_total_snps_removed = passingdata.no_of_total_snps_filtered +\
					passingdata.no_of_monomorphic_snps_removed
				
				passingdata_ls.append(passingdata)
			del reader
		return passingdata_ls
	
	def filter_passingdata(cls, passingdata):
		"""
		if any avg/std = -1, it's not enough data, exclude it
		"""
		import os,sys,csv
		from variation.src.MpiQCCall import MpiQCCall
		avg_var_name_pair_ls, partial_header_avg = MpiQCCall.generate_avg_variable_names(MpiQCCall.avg_var_name_ls)
		good_data = 1
		for avg_var_name_pair in avg_var_name_pair_ls:
			ls_var_name, avg_var_name, std_var_name = avg_var_name_pair
			if getattr(passingdata, avg_var_name)==-1:
				good_data = 0
				break
			if getattr(passingdata, std_var_name)==-1:
				good_data = 0
				break
		return good_data
	filter_passingdata = classmethod(filter_passingdata)
	
	def plot_variables(self, passingdata_ls, plot_var_name_ls, output_fname, **keywords):
		"""
		2008-05-13
			keywords is to restrict the data
		"""
		plot_var_value_ls = []
		for i in range(len(plot_var_name_ls)):
			plot_var_value_ls.append([])
		
		for passingdata in passingdata_ls:
			good_data = self.filter_passingdata(passingdata)
			for key, value in keywords.iteritems():
				if getattr(passingdata, key)!=value:
					good_data=0
					break
			if good_data==1:
				for i in range(len(plot_var_name_ls)):
					var_name = plot_var_name_ls[i]
					plot_var_value = getattr(passingdata, var_name)
					plot_var_value_ls[i].append(plot_var_value)
		self._plot(plot_var_value_ls, plot_var_name_ls, output_fname)
		return plot_var_value_ls
	
	#plot_variables  = classmethod(plot_variables)
	
	def _plotContour(cls, plot_var_value_ls, plot_var_name_ls, output_fname_prefix):
		import pylab
		import matplotlib.axes3d as p3
		pylab.clf()
		pylab.contour(plot_var_value_ls[0], plot_var_value_ls[1], plot_var_value_ls[2])
		pylab.xlabel(plot_var_name_ls[0])
		pylab.ylabel(plot_var_name_ls[1])
		pylab.title(plot_var_name_ls[2])
		#pylab.hot()
		pylab.colorbar()
		#fig.add_axes(ax)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=100)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=100)
		pylab.show()
	_plotContour = classmethod(_plotContour)
	
	
	def _plot3D(cls, plot_var_value_ls, plot_var_name_ls, output_fname_prefix):
		import pylab
		import matplotlib.axes3d as p3
		pylab.clf()
		fig=pylab.figure()
		ax = p3.Axes3D(fig)
		# plot3D requires a 1D array for x, y, and z
		# ravel() converts the 100x100 array into a 1x10000 array
		#ax.plot3D(plot_var_value_ls[0], plot_var_value_ls[1], plot_var_value_ls[2])
		ax.scatter3D(plot_var_value_ls[0], plot_var_value_ls[1], plot_var_value_ls[2])
		ax.set_xlabel(plot_var_name_ls[0])
		ax.set_ylabel(plot_var_name_ls[1])
		ax.set_zlabel(plot_var_name_ls[2])
		#fig.add_axes(ax)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=100)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=100)
		pylab.show()
	_plot3D = classmethod(_plot3D)
	
	def run(self):
		if self.debug:	#serial debug
			import pdb
			pdb.set_trace()
		passingdata_ls = self.get_qccall_results(self.input_dir)
		
		#6 free variables, control 4 and show 2 + final strain/snp_mismatch or no_of_total_accessions_removed or no_of_total_snps_removed

		keywords = {}
		output_fname_ls = []
		for keyword_candidate in self.keyword_candidate_ls:
			value = getattr(self, keyword_candidate)
			if value !='' and value is not None:
				keywords[keyword_candidate] = value
				output_fname_ls.append('%s_%s'%(keyword_candidate, value))
		
		plot_var_name_ls = []
		for var_index in self.var_to_plot_indices[:3]:
			var_index = int(var_index)
			plot_var_name_ls.append(self.candidate_plot_var_name_ls[var_index])
			output_fname_ls.append(self.candidate_plot_var_name_ls[var_index])
		#plot_var_name_ls = ['min_call_probability', 'max_call_mismatch_rate', 'avg_mismatch_rate']
		
		if not self.output_fname:
			self.output_fname = '_'.join(output_fname_ls)
		
		self.plot_variables(passingdata_ls, plot_var_name_ls, self.output_fname, **keywords)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PlotQCCall
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	
