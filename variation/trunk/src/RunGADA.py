#!/usr/bin/env python
"""

Examples:
	RunGADA.py -i -o
	RunGADA.py -i ./call_method_17_CNV_array_intensity_chr4_line_no_888148_1107622_norm.tsv -o call_method_17_CNV_array_intensity_chr4_line_no_888148_1107622_norm_w1_34.out -w1,34 -u yh
	
Description:
	2008-12-05 program to run GADA (Pique-Regi2008) to infer states/amplitude of CNV probes.
		it takes output from CNVNormalize.py as input and outputs 2 files,
		file 1 containing StrainXSNP matrix with -1(Loss),0,+1(Gain) as states.
		file 2 (with '_amp' embedded in file 1's name) containing StrainXSNP matrix with actual amplitude values.
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv, numpy, traceback, subprocess
from pymodule import figureOutDelimiter, getListOutOfStr
from CNVNormalize import CNVNormalize
from PlotGroupOfSNPs import PlotGroupOfSNPs
import StringIO
import Stock_250kDB
from sets import Set

class RunGADA(CNVNormalize):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname', 1, ): ['', 'i', 1, 'CNV intensity matrix, probe X arrays. 1st column is probe id. 2nd last col is chr. last col is pos.', ],\
							('output_fname', 1, ): ['', 'o', 1, 'self-explanatory', ],\
							('which_array_id_ls', 0, ): [None, 'w', 1, 'list of array ids indicating which arrays(s) format: 0,1-3', ],\
							('tmp_input_fname', 1, ): ['/tmp/GADA_input', '', 1, 'temporary file to store GADA input'],\
							('tmp_output_fname', 1, ): ['/tmp/GADA_output', '', 1, 'temporary file to store GADA output', ],\
							('GADA_path', 1, ): [os.path.expanduser('~/script/variation/bin/GADA/GADA'), '', 1, 'path to the GADA program'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		"""
		2008-12-05
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
		self.which_array_id_ls = getListOutOfStr(self.which_array_id_ls, data_type=int)
	def prepareGADAinput(self, data_matrix, which_column, tmp_input_fname):
		"""
		"""
		if self.report:
			sys.stderr.write("Preparing column %s input for GADA ..."%which_column)
		
		of = open(tmp_input_fname, 'w')
		no_of_rows, no_of_cols = data_matrix.shape
		
		for i in range(no_of_rows):
			of.write('%s\n'%data_matrix[i][which_column])
		del of
		if self.report:
			sys.stderr.write("Done.\n")
	
	def _GADA(self, GADA_path, tmp_input_fname, tmp_output_fname):
		"""
		"""
		if self.report:
			sys.stderr.write("Running GADA ...")
		commandline = '%s -a 0.8 -T 2 -M 2 -s -0.4 -b 0.0 -c'%GADA_path
		inf = open(tmp_input_fname)
		command_handler = subprocess.Popen(commandline, shell=True, stdin=inf, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		stdout_content, stderr_content = command_handler.communicate()
		del inf
		if stderr_content:
			sys.stderr.write('stderr of %s: %s \n'%(commandline, stderr_content))
		if self.report:
			sys.stderr.write("Done.\n")
		of = open(tmp_output_fname, 'w')
		of.write(stdout_content)
		del of
		return stdout_content
	
	
	def findOutWhichColumn(cls, col_id_ls, which_array_id_set):
		"""
		2008-12-05
		"""
		col_index_to_return_ls = []
		for i in range(len(col_id_ls)):
			col_id = int(col_id_ls[i])
			if col_id in which_array_id_set:
				col_index_to_return_ls.append(i)
		return col_index_to_return_ls
	findOutWhichColumn = classmethod(findOutWhichColumn)
	
	def output_header(self, writer, chr_pos_ls):
		"""
		2008-12-05
		"""
		sys.stderr.write("Outputting header ...")
		header = ['ecotype_id', 'array_id', ]
		for chr_pos in chr_pos_ls:
			chr_pos = '_'.join(chr_pos)
			header.append(chr_pos)
		writer.writerow(header)
		sys.stderr.write("Done.\n")
	state2number = {'G': 1, 'N':0, 'L':-1}
	
	def output_GADA_output(self, writer, writer_amp, GADA_output, array_id, ecotype_id):
		"""
		2008-12-05
		"""
		sys.stderr.write("Outputting GADA output ...")
		GADA_output = StringIO.StringIO(GADA_output)
		counter = 0
		data_row = [ecotype_id, array_id]
		data_row_amp = [ecotype_id, array_id]
		for line in GADA_output:
			if line[0]=='#':
				continue
			else:
				counter += 1
			if counter>1:	#skip the first line of real data. it's header "Start   Stop    Lenght  Ampl    State"
				row = line[:-1].split('\t')
				start = int(row[0])
				stop = int(row[1])
				state = row[-1]
				amp = float(row[3])
				if amp>0.5:
					state = 'G'
				elif amp<-0.5:
					state='L'
				else:
					state = 'N'
				state_number = self.state2number[state]
				for i in range(start-1, stop):
					data_row.append(state_number)
					data_row_amp.append(amp)
		writer.writerow(data_row)
		writer_amp.writerow(data_row_amp)
		sys.stderr.write("Done.\n")
	
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		
		data_matrix, probe_id_ls, chr_pos_ls, header = self.get_input(self.input_fname)
		array_id_ls = header[1:-2]
		array_id_ls = map(int, array_id_ls)
		if self.which_array_id_ls:
			col_index_ls = self.findOutWhichColumn(array_id_ls, Set(self.which_array_id_ls))
		else:
			col_index_ls = range(data_matrix.shape[1])
		
		output_amp_fname = '%s_amp.tsv'%os.path.splitext(self.output_fname)[0]
		writer_amp  = csv.writer(open(output_amp_fname, 'w'), delimiter='\t')
		writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
		self.output_header(writer, chr_pos_ls)
		self.output_header(writer_amp, chr_pos_ls)
		
		for col_index in col_index_ls:
			array_id = array_id_ls[col_index]
			array = Stock_250kDB.ArrayInfo.get(array_id)
			if array:
				ecotype_id = array.maternal_ecotype_id
			else:
				continue
			self.prepareGADAinput(data_matrix, col_index, self.tmp_input_fname)
			GADA_output = self._GADA(self.GADA_path, self.tmp_input_fname, self.tmp_output_fname)
			self.output_GADA_output(writer, writer_amp, GADA_output, array_id, ecotype_id)
		del writer, writer_amp

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = RunGADA
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()