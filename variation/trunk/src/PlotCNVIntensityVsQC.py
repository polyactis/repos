#!/usr/bin/env python
"""

Examples:
	
	PlotCNVIntensityVsQC.py -i /Network/Data/250k/tmp-yh/CNV/call_method_17_CNV_array_intensity_norm_chr4_GADA_out_amp.tsv -q /Network/Data/250k/tmp-yh/CNV/2010_CNV_probe_qc_deletion.tsv -o /Network/Data/250k/tmp-yh/CNV/no_of_deletion_vs_intensity/
	
Description:
	Program to generate QC data from a SNP data for CNV data. 3 StrainXProbe files generated:
		mismatch matrix file. containing the no_of_mismatches for each probe. -2 is NA.
		insertion matrix file. containing the no_of_insertions for each probe. -2 is NA.
		deletion matrix file. containing the no_of_deletions for each probe. -2 is NA.
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
from CNVNormalize import CNVNormalize
import Stock_250kDB

class PlotCNVIntensityVsQC(object):
	__doc__ = __doc__
	"""
	2009-2-12
	"""
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("input_fname", 1, ): [None, 'i', 1, 'StrainXProbe intensity matrix file'],\
							('qc_fname', 1, ): [None, 'q', 1, 'StrainXProbe QC no of mismatch/deletion/insertion(es) matrix'],\
							('output_dir',1,): ['', 'o', 1, 'output directory to hold no_of_mismatches vs intensity scatter plots'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2009-2-12
		"""
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def getBeforeGADAIntensityData(self, input_fname):
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		
		data_matrix, probe_id_ls, chr_pos_ls, header = CNVNormalize.get_input(input_fname)
		
		col_id_ls = []
		for chr_pos in chr_pos_ls:
			col_id_ls.append('%s_%s'%(chr_pos[0], chr_pos[1]))
		
		ecotype_id_ls = []
		for array_id in header[1:-2]:
			array = Stock_250kDB.ArrayInfo.get(int(array_id))
			if array:
				ecotype_id = array.maternal_ecotype_id
				
			else:
				ecotype_id = -1
			ecotype_id_ls.append('%s'%ecotype_id)
		cnvIntensityData = SNPData(row_id_ls=ecotype_id_ls, col_id_ls=col_id_ls, data_matrix=data_matrix.transpose())
		return cnvIntensityData
	
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		cnvIntensityData = self.getBeforeGADAIntensityData(self.input_fname)
		#cnvIntensityData = SNPData(input_fname=self.input_fname, turn_into_array=1, ignore_2nd_column=1, matrix_data_type=float)
		
		qcData = SNPData(input_fname=self.qc_fname, turn_into_array=1, ignore_2nd_column=1)
		if not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
		
		for probe_id in qcData.col_id_ls:
			if probe_id in cnvIntensityData.col_id2col_index:
				cnv_col_index = cnvIntensityData.col_id2col_index[probe_id]
				qc_col_index = qcData.col_id2col_index[probe_id]
				count_ls = []
				intensity_ls = []
				for i in range(len(qcData.row_id_ls)):
					row_id = qcData.row_id_ls[i]
					if qcData.data_matrix[i][qc_col_index]>=0 and row_id in cnvIntensityData.row_id2row_index:
						cnv_row_index = cnvIntensityData.row_id2row_index[row_id]
						count = qcData.data_matrix[i][qc_col_index]
						count_ls.append(count)
						intensity_ls.append(cnvIntensityData.data_matrix[cnv_row_index][cnv_col_index])
				count_set = set(count_ls)
				if len(count_set)>0 and count_set!=set([0]):
					pylab.clf()
					ax = pylab.axes([0.1, 0.1, 0.8, 0.8], frameon=False)
					ax.grid(True, alpha=0.3)
					pylab.plot(count_ls, intensity_ls, '.', markersize=5, alpha=0.4)
					pylab.xlabel('count')
					pylab.ylabel('CNV probe intensity')
					pylab.ylim([-1,1])
					xlim = list(ax.get_xlim())
					xlim[0] -= 1
					xlim[1] += 1
					ax.set_xlim(xlim)
					pylab.title(probe_id)
					pylab.savefig(os.path.join(self.output_dir, '%s.png'%probe_id), dpi=300)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PlotCNVIntensityVsQC
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()