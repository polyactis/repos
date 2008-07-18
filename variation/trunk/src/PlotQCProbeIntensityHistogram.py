#!/usr/bin/env python

"""
Examples:
	#function 1
	PlotQCProbeIntensityHistogram.py -i 396 -o ~/Desktop/array_intensity_histogram/396_QC.png -m 50000

	#function 2
	PlotQCProbeIntensityHistogram.py -i 498 -o ~/Desktop/array_intensity_histogram/498.png  -p ~/script/affy/250k_test/250kprobe_subset.txt


Description:
	1. check intensity histogram of QC probes given an array id
	2. check intensity histogram of probes given a probe subset file
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, stat, getopt
import traceback, gc, subprocess


class PlotQCProbeIntensityHistogram:
	__doc__ = __doc__
	option_default_dict = {('array_id', 1, int): [None, 'i'],\
						('cdf_fname', 1, ):[os.path.expanduser('~/script/variation/genotyping/250ksnp/data/atSNPtilx520433_rev2/Full/atSNPtil_geno/LibFiles/atSNPtil_geno.cdf')],\
							('figure_ofname', 1, ): ['', 'o', 1, 'Figure output filename', ],\
							('input_dir', 1, ): ['/Network/Data/250k/db/raw_data/', 'n', 1, 'Directory where cel files are stored.' ],\
							('no_of_bins', 1, int): [100, 'f', 1, 'Number of bins in the histogram plot'],\
							('max_intensity', 1, int,):[25000, 'm', 1, 'Maximum intensity to put into plot' ],\
							('probe_subset_fname', 0, ):['', 'p', 1, 'probe subset filename, if given, will drawn histogram of those probes rather than QC probes. ~/script/affy/250k_test/250kprobe_subset.txt' ],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  **keywords):
		"""
		2008-07-05
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		self.getIntensityFuncDict = {1: self.getQCIntensityLs,
									2: self.getIntensityLs}
		if os.path.isfile(self.probe_subset_fname):
			self.func_type = 2
		else:
			self.func_type = 1
		
	def getIntensityLs(self, intensity_output_fname, max_intensity):
		"""
		2007-07-18
		"""
		sys.stderr.write("Getting intensity from %s ... "%(intensity_output_fname))
		reader = csv.reader(open(intensity_output_fname), delimiter='\t')
		reader.next()
		intensity_ls = []
		for row in reader:
			for i in range(1,len(row)):
				intensity = float(row[i])
				if intensity<=max_intensity:
					intensity_ls.append(math.log10(intensity))
		sys.stderr.write("Done.\n")
		return intensity_ls
	
	def getQCIntensityLs(self, qc_intensity_fname, max_intensity):
		"""
		2007-07-18
		"""
		reader = csv.reader(open(qc_intensity_fname), delimiter='\t')
		reader.next()
		intensity_ls = []
		for row in reader:
			intensity = float(row[3])
			if intensity<=max_intensity:
				intensity_ls.append(math.log10(intensity))
		return intensity_ls
	
	def _plot(self, intensity_ls, array_id, figure_ofname, no_of_bins=25, max_intensity=25000):
		"""
		2008-07-05
		"""
		import pylab
		pylab.clf()
		pylab.title("Array: %s"%array_id)
		pylab.hist(intensity_ls, no_of_bins)
		pylab.savefig(figure_ofname, dpi=100)
		#pylab.show()

	def run(self):
		array_cel_fname = os.path.join(self.input_dir, '%s_raw_data.cel'%self.array_id)
		
		intensity_output_fname = '/tmp/%s_%s.tsv'%(self.array_id, self.func_type)
		
		if os.path.isfile(self.probe_subset_fname):
			command = os.path.expanduser('~/script/affy/sdk/calvin_files/ReadatSNPtilgeno')
			cmd_p = subprocess.Popen([command, '-i', array_cel_fname, '-d', self.cdf_fname, '-e', repr(self.array_id), '-p', self.probe_subset_fname, '-o', intensity_output_fname ], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		else:
			OutputQCIntensity_command = os.path.expanduser('~/script/affy/sdk/calvin_files/OutputQCIntensity')
			cmd_p = subprocess.Popen([OutputQCIntensity_command, '-i', array_cel_fname, '-d', self.cdf_fname, '-o', intensity_output_fname ], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		cmd_p.wait()
		cmd_out = cmd_p.stdout.read()
		print cmd_out
		cmd_err = cmd_p.stderr.read()
		sys.stderr.write(cmd_err)
		
		intensity_ls = self.getIntensityFuncDict[self.func_type](intensity_output_fname, self.max_intensity)
		
		self._plot(intensity_ls, self.array_id, self.figure_ofname, self.no_of_bins, self.max_intensity)
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PlotQCProbeIntensityHistogram
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()