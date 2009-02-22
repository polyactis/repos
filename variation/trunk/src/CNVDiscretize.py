#!/usr/bin/env python
"""

Examples:
	CNVDiscretize.py -i /Network/Data/250k/tmp-yh/CNV/call_method_17_CNV_array_intensity_norm_chr4_GADA_out_amp.tsv  -o /Network/Data/250k/tmp-yh/CNV/call_method_17_CNV_array_intensity_norm_chr4_GADA_out_discrete_m.3125.tsv 
	
	CNVDiscretize.py -i /Network/Data/250k/tmp-yh/CNV/call_method_17_CNV_array_intensity_norm_chr$i\_GADA_M10_discrete_m.3125.tsv -o /Network/Data/250k/tmp-yh/CNV/discrete_m.3125_chromosomal_plot_M10/call_method_17_chr$i  -y2
	
Description:
	Program to 
	1. discretize CNV amplitude file outputted by GADA into 1 (normal & duplication) or -1 (deletion).
	2. make segmental plots in which each plot contains consecutive 1000 probes. y-axis is probe value summarizing across all accessions.
		input could be output of RunGADA.py or from step 1.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import getopt, numpy
from variation.src.common import nt2number, number2nt
from pymodule import ProcessOptions, SNPData, PassingData
from DB_250k2Array import DB_250k2Array
import pylab

class CNVDiscretize(object):
	__doc__ = __doc__
	"""
	2009-2-12
	"""
	option_default_dict = {("input_fname", 1, ): [None, 'i', 1, 'StrainXProbe intensity matrix file (amplitude output by RunGADA.py)'],\
							('output_fname',1,): ['', 'o', 1, 'output filename'],\
							('max_del_intensity', 1, float):[-0.03125, 'm',1, 'maximum intensity to call deletion'],\
							('run_type',1, int): [1, 'y', 1, 'running type. 1: Discretize the CNV amplitude into 1 or -1 according to max_del_intensity. 2: Plot chromosomal-wide probe value aggregated across all samples'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2009-2-12
		"""
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self):
		cnvIntensityData = SNPData(input_fname=self.input_fname, turn_into_array=1, ignore_2nd_column=1, matrix_data_type=float)
		probe_pos_ls = []
		avg_intensity_ls = []
		
		if self.run_type == 1:
			newDataMatrix = numpy.ones(cnvIntensityData.data_matrix.shape, numpy.int)
		
		for j in range(cnvIntensityData.data_matrix.shape[1]):
			probe_id = cnvIntensityData.col_id_ls[j]
			probe_id = probe_id.split('_')
			probe_id = map(int, probe_id)
			probe_pos_ls.append(probe_id[1])
			avg_intensity_ls.append(numpy.sum(cnvIntensityData.data_matrix[:,j]))
			if self.run_type==1:
				for i in range(cnvIntensityData.data_matrix.shape[0]):
					if cnvIntensityData.data_matrix[i][j]<=self.max_del_intensity:
						newDataMatrix[i][j] = -1
		
		if self.run_type==1:
			newData = SNPData(row_id_ls=cnvIntensityData.row_id_ls, col_id_ls=cnvIntensityData.col_id_ls, data_matrix=newDataMatrix)
			newData.tofile(self.output_fname)
		elif self.run_type==2:
			block_size = 1000
			no_of_probes = len(probe_pos_ls)
			no_of_blocks = no_of_probes/block_size
			for i in range(no_of_blocks):
				if i*block_size>no_of_probes:
					break
				start_index = i*block_size
				end_index = min((i+1)*block_size, no_of_probes)
				fname = '%s_%s_%s.png'%(self.output_fname, probe_pos_ls[start_index], probe_pos_ls[end_index])
				pylab.clf()
				pylab.plot(probe_pos_ls[start_index:end_index], avg_intensity_ls[start_index:end_index], '.', markersize=4, alpha=0.4)
				pylab.xlabel('chromosome position')
				pylab.ylabel('sum intensity')
				pylab.savefig(fname, dpi=300)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CNVDiscretize
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()