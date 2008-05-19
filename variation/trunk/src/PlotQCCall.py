#!/usr/bin/env python
"""
Examples:
	#default x axs='min_call_probability', y-axis='max_call_mismatch_rate', z-axis='avg_mismatch_rate'
	./src/PlotQCCall.py -i /mnt/nfs/NPUTE_data/QC_avg_output

	#replace the z axis with no_of_total_accessions_filtered
	./src/PlotQCCall.py -i /mnt/nfs/NPUTE_data/QC_avg_output -t 1_2_9
	
	#replace the z axis with no_of_total_snps_filtered
	./src/PlotQCCall.py -i /mnt/nfs/NPUTE_data/QC_avg_output -t 1_2_14
	
	#avg_mismatch_rate (7) ~ min_call_probability (1) * npute_window_size (6)
	./src/PlotQCCall.py -i /mnt/nfs/NPUTE_data/QC_avg_output -x 0.15 -w 0.25 -n -1 -t 1_6_7
	
	#avg_mismatch_rate (7) ~ max_call_mismatch_rate (2) * npute_window_size (6)	#min_call_probability not fixed
	./src/PlotQCCall.py -i /mnt/nfs/NPUTE_data/QC_avg_output -x -1 -w 0.25 -n -1 -t 2_6_7
	
	#avg_mismatch_rate (7) ~ max_call_mismatch_rate (2) * max_call_NA_rate (3)
	./src/PlotQCCall.py -i /mnt/nfs/NPUTE_data/QC_avg_output -x -1 -w 0.25 -a -1 -n 50 -t 2_3_7
	
	#avg_mismatch_rate (7) ~ no_of_total_accessions_filtered (9) * no_of_total_snps_removed (13) [after imputation]
	./src/PlotQCCall.py -i /mnt/nfs/NPUTE_data/QC_avg_output -v -1 -w -1 -a -1 -n 50 -t 9_13_7 
	
	#avg_mismatch_rate (7) ~ no_of_total_accessions_filtered (9) * no_of_total_snps_removed (13) [before imputation]
	./src/PlotQCCall.py -i /mnt/nfs/NPUTE_data/QC_avg_output -v -1 -w -1 -a -1 -n 50 -t 9_13_7 -f 0
	
Description:
	program to plot MpiQCCall.py's avg.csv output.
	Basically pick 3 variable to draw 3D plots while using other variables to restrict data points (subsetting data).
	
	doc/QC_parameter/QC_parameter_tuning.tex has more example commandlines and corresponding figures.
	
	-1 in arguments invalidates that argument. not used in subsetting data
"""

import sys, os, math
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))


class PlotQCCall(object):
	keyword_candidate_ls = ['strain or snp', 'after_imputation', 'min_call_probability', 'max_call_mismatch_rate', 'max_call_NA_rate', 'max_snp_mismatch_rate',
						'max_snp_NA_rate', 'npute_window_size']
	candidate_plot_var_name_ls = keyword_candidate_ls[1:]+ ['avg_mismatch_rate', 'avg_NA_rate', \
								'no_of_total_accessions_filtered', 'no_of_accessions_filtered_by_mismatch', 'no_of_accessions_filtered_by_na', \
								'no_of_total_snps_filtered', 'no_of_total_snps_removed', 'no_of_snps_filtered_by_mismatch', 'no_of_snps_filtered_by_na',\
								'no_of_monomorphic_snps_removed']
	candidate_plot_var_name_index_tuple_ls = [(candidate_plot_var_name_ls[i], i) for i in range(len(candidate_plot_var_name_ls))]
	__doc__ = __doc__  + 'candidates for keywords: ' + repr(keyword_candidate_ls) + '\nCandidate plot variables and indices: ' + repr(candidate_plot_var_name_index_tuple_ls)
	
	option_default_dict = {('input_dir', 1, ): ['/mnt/nfs/NPUTE_data/QC_avg_output', 'i', 1, ],\
							('output_fname', 0, ): [None, 'o', 1, 'not given, the program would construct an output filename based on arguments.', ],\
							('output_for_qhull_fname', 0, ): [None, 'q', 1, 'qhull output'],\
							('var_to_plot_indices', 1, ):['1_2_7', 't', 1, '3 indices separated by _, corresponding to candidate plot variable list. -h to see the list'],\
							('strain or snp', 1, ): ['strain', 's', 1,],\
							('after_imputation', 0, int): [1, 'f'],\
							('min_call_probability', 0, float): [-1, 'y', 1, 'minimum probability for a call to be non-NA if there is a 3rd column for probability.', ],\
							('max_call_mismatch_rate', 0, float): [-1, 'x', 1, 'maximum mismatch rate of an array call_info entry. used to exclude bad arrays.'],\
							('max_call_NA_rate', 0, float,): [0.4, 'a', 1, 'maximum NA rate of an array call_info entry. used to exclude bad arrays.'],\
							('max_snp_mismatch_rate', 0, float,): [0.2, 'w', 1, 'maximum snp error rate, used to exclude bad SNPs', ],\
							('max_snp_NA_rate', 0, float,): [0.4, 'v', 1, 'maximum snp NA rate, used to exclude SNPs with too many NAs', ],\
							('npute_window_size', 0, int): [50, 'n', 1, 'NPUTE window size',],\
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
		"""
		var_name_ls = ['strain or snp', 'after_imputation'] + MpiQCCall.common_var_name_ls
		avg_var_name_pair_ls, partial_header_avg = MpiQCCall.generate_avg_variable_names(MpiQCCall.avg_var_name_ls)
		var_name_ls += partial_header_avg
		"""
		files = os.listdir(input_dir)
		passingdata_ls = []
		no_of_objects = len(files)
		var_name_ls = []
		for i in range(no_of_objects):
			sys.stderr.write("\t%d/%d: from %s ... \n"%(i+1, no_of_objects, files[i]))
			filename = os.path.join(input_dir, files[i])
			reader = csv.reader(open(filename))
			try:
				row = reader.next()
				if len(var_name_ls)==0:
					var_name_ls = row
			except:
				if self.debug:
					import traceback
					traceback.print_exc()
	  				sys.stderr.write('%s\n'%sys.exc_info())
	  			sys.stderr.write('\terror in reading this file. ignored.\n')
	  			del reader
				continue
			for row in reader:
				passingdata = PassingData()
				for i in range(len(var_name_ls)):
					var_name = var_name_ls[i]
					if var_name!='strain or snp':
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
		return passingdata_ls, var_name_ls
	
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
			if hasattr(passingdata, avg_var_name) and getattr(passingdata, avg_var_name)==-1:
				good_data = 0
				break
			if hasattr(passingdata, std_var_name) and getattr(passingdata, std_var_name)==-1:
				good_data = 0
				break
		return good_data
	filter_passingdata = classmethod(filter_passingdata)
	
	def output_passingdata_in_qhull_format(self, passingdata_ls, candidate_plot_var_name_ls, plot_var_name_ls, output_for_qhull_fname):
		"""
		2008-05-16
			to prepare the data (3D in three variables chosen in plot_var_name_ls) in qhull format
		"""
		sys.stderr.write("Outputting in qhull format to %s ..." %output_for_qhull_fname)
		import csv
		writer = csv.writer(open(output_for_qhull_fname, 'w'), delimiter=' ')
		writer.writerow([len(plot_var_name_ls)])
		writer.writerow([len(passingdata_ls)])
		for passingdata in passingdata_ls:
			one_row = []
			for var_name in plot_var_name_ls:
				one_row.append(getattr(passingdata, var_name))
			writer.writerow(one_row)
		
		writer.writerow(['#'] + candidate_plot_var_name_ls)
		
		for passingdata in passingdata_ls:
			one_row = ['#']
			for var_name in candidate_plot_var_name_ls:
				one_row.append(getattr(passingdata, var_name))
			writer.writerow(one_row)
		del writer
		sys.stderr.write(" Done.\n")
		
	def plot_variables(self, passingdata_ls, plot_var_name_ls, output_fname, **keywords):
		"""
		2008-05-13
			keywords is to restrict the data
		"""
		plot_var_value_ls = []
		for i in range(len(plot_var_name_ls)):
			plot_var_value_ls.append([])
		self.good_passingdata_ls = []
		
		for passingdata in passingdata_ls:
			good_data = self.filter_passingdata(passingdata)
			for key, value in keywords.iteritems():
				if hasattr(passingdata, key) and getattr(passingdata, key)!=value:
					good_data=0
					break
			if good_data==1:
				for i in range(len(plot_var_name_ls)):
					var_name = plot_var_name_ls[i]
					plot_var_value = getattr(passingdata, var_name)
					plot_var_value_ls[i].append(plot_var_value)
				self.good_passingdata_ls.append(passingdata)
		self._plot(plot_var_value_ls, plot_var_name_ls, output_fname)
		return self.good_passingdata_ls
	
	#plot_variables  = classmethod(plot_variables)
	
	def _plotContour(cls, plot_var_value_ls, plot_var_name_ls, output_fname_prefix):
		"""
		2008-05-18
			this function is not runnable yet. pylab.contour()'s 3 arguments are prepared wrongly.
		"""
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
	
	def on_click_showing_passingdata(cls, event):
		"""
		2008-05-15
			inspired by examples/pick_event_demo.py from matplotlib source code
			
			useless as Axes3D doesn't support pick event, check plone doc, log/python/plot
		"""
		good_passingdata_ls = getattr(cls, good_passingdata_ls)
		if good_passingdata_ls:
			ind = event.ind
        	print "picked: ",ind
        	for plot_var_name in candidate_plot_var_name_ls:
        		plot_var_value = getattr(good_passingdata_ls[ind], plot_var_name)
        		print plot_var_name, plot_var_value
	
	on_click_showing_passingdata = classmethod(on_click_showing_passingdata)
	
	def _plot3D(cls, plot_var_value_ls, plot_var_name_ls, output_fname_prefix, picker_function=on_click_showing_passingdata):
		import pylab
		import matplotlib.axes3d as p3
		pylab.clf()
		fig=pylab.figure()
		fig.canvas.mpl_connect('pick_event', picker_function)
		ax = p3.Axes3D(fig)
		# plot3D requires a 1D array for x, y, and z
		# ravel() converts the 100x100 array into a 1x10000 array
		#ax.plot3D(plot_var_value_ls[0], plot_var_value_ls[1], plot_var_value_ls[2])
		ax.scatter3D(plot_var_value_ls[0], plot_var_value_ls[1], plot_var_value_ls[2], picker=5)
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
		passingdata_ls, total_var_name_ls = self.get_qccall_results(self.input_dir)
		
		#6 free variables, control 4 and show 2 + final strain/snp_mismatch or no_of_total_accessions_removed or no_of_total_snps_removed

		keywords = {}
		output_fname_ls = []
		for keyword_candidate in self.keyword_candidate_ls:
			if not hasattr(self, keyword_candidate):
				continue
			value = getattr(self, keyword_candidate)
			chosen=1
			if isinstance(value, float) or isinstance(value, int):
				if value<0:
					chosen = 0
			elif value =='' or value is None:
				chosen = 0
			if chosen:
				keywords[keyword_candidate] = value
				
				output_fname_ls.append('%s=%s'%(keyword_candidate.replace('_', ' '), value))
		
		plot_var_name_ls = []
		var_indices = self.var_to_plot_indices.split('_')
		for var_index in var_indices[:3]:
			var_index = int(var_index)
			plot_var_name_ls.append(self.candidate_plot_var_name_ls[var_index])
			output_fname_ls.append(self.candidate_plot_var_name_ls[var_index].replace('_', ' '))
		#plot_var_name_ls = ['min_call_probability', 'max_call_mismatch_rate', 'avg_mismatch_rate']
		
		if not self.output_fname:
			self.output_fname = ', '.join(output_fname_ls)
		
		good_passingdata_ls = self.plot_variables(passingdata_ls, plot_var_name_ls, self.output_fname, **keywords)
		if self.output_for_qhull_fname:
			self.output_passingdata_in_qhull_format(good_passingdata_ls, total_var_name_ls, plot_var_name_ls, self.output_for_qhull_fname)

def recoverQhullOutput(qhull_input_fname, qhull_output_fname, output_fname):
	"""
	2008-05-16
		pass the qhull input/output and pick the hull points and output corresponding complete stat data to output_fname
	"""
	import os, sys, csv
	from sets import Set
	inf1 = open(qhull_output_fname)
	point_index_set = Set()
	for line in inf1:
		point_index = int(line[:-1])
		point_index_set.add(point_index)
	del inf1
	
	reader = csv.reader(open(qhull_input_fname), delimiter=' ')	#outputted by output_passingdata_in_qhull_format()
	writer = csv.writer(open(output_fname, 'w'))
	point_index = 0
	for row in reader:
		if row[0]=='#':	#outputted by output_passingdata_in_qhull_format
			if point_index == 0:
				writer.writerow(row[1:])	#keep the header
			elif point_index-1 in point_index_set:	#0 is the header, 1 is the actual info
				writer.writerow(row[1:])
			point_index += 1
	del reader, writer

"""
qhull_input_fname = '/tmp/qhull.in'
qhull_output_fname = '/tmp/qhull.out'
output_fname = '/tmp/qhull.out.passingdata'
recoverQhullOutput(qhull_input_fname, qhull_output_fname, output_fname)
"""

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PlotQCCall
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	
